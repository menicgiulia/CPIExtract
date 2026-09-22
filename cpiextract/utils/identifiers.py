'''Retrieves protein or compond identifiers from servers.'''

import pubchempy as pcp
import pandas as pd
import time
import re

from ..servers.BiomartServer import BiomartServer
from ..servers.MyGeneServer import MyGeneServer
from ..servers.PubChemServer import PubChemServer


def _pubchem_call_with_retry(func, *args, max_retries: int = 3, backoff_base: float = 1.0,
                              verbose: bool = False, **kwargs):
    """
    Calls a PubChem-hitting function (pubchempy under the hood) with retry-with-backoff.
    PubChem's PUG-REST service is known to intermittently return errors/timeouts under
    load - without this, a single transient failure surfaces identically to "this
    compound doesn't exist" (both currently fall through to the same broad except:
    below), which is misleading and wastes the person's time chasing a bad input that
    was never actually the problem. Retries max_retries times total, waiting
    backoff_base * 2^attempt seconds between attempts (1s, 2s, 4s by default) before
    giving up and letting the final exception propagate to the caller's own handling.
    """
    last_exception = None
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            last_exception = e
            if attempt < max_retries - 1:
                wait = backoff_base * (2 ** attempt)
                if verbose:
                    print(f"  PubChem call failed (attempt {attempt + 1}/{max_retries}): "
                          f"{e} - retrying in {wait:.1f}s...")
                time.sleep(wait)
    raise last_exception


# Uniprot format checks:
#   Format 1 (reviewed Swiss-Prot): [OPQ][0-9][A-Z0-9]{3}[0-9]  e.g. P11473
#   Format 2 (TrEMBL):              [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}  e.g. A0A000
#   Optional isoform suffix:        -\d+  e.g. P11473-2
UNIPROT_REGEX = re.compile(
    r'^([OPQ][0-9][A-Z0-9]{3}[0-9]'
    r'|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(-\d+)?$'
)

def protein_identifiers(input_id: int | str, gene_server=None, server_select='mygene') -> pd.DataFrame:
    """
    Retrieves protein identifiers from Biomart.

    Sends a search query to Biomart using the input_id parameter as input.

    Parameters
    ----------
    input_id : int/string
        Identifier of the protein(s) to search on the server. It can be one of the following types: \\
        Entrezgene_id (int), hgnc_id, ensembl_peptide_id, ensembl_gene_id, uniprotswissprot, chembl or hgnc_symbol

    Returns
    -------
    DataFrame
        Dataframe of the protein(s) containing the following values: \\
        uniprotswissprot, entrezgene_id, chembl, ensembl_peptide_id, ensembl_gene_id, hgnc_symbol, gene_type, description
    """

    if gene_server is not None:
        pass
    elif server_select == 'mygene':
        gene_server = MyGeneServer()
    elif server_select == 'biomart':
        gene_server = BiomartServer()
    else:
        raise ValueError(f"server_select must be 'mygene' or 'biomart', got '{server_select}'")

    try:
        # Check if input_id is an integer
        if isinstance(input_id, int):
            inputtype = 'entrezgene_id'
        # Find input type based on input size and/or words
        elif input_id.find('HGNC:') == 0:
            inputtype='hgnc_id'
        elif input_id.find('ENSP') == 0 and len(input_id) == 15:
            inputtype='ensembl_peptide_id'
        elif input_id.find('ENSG') == 0 and len(input_id) == 15:
            inputtype='ensembl_gene_id'
        elif UNIPROT_REGEX.match(str(input_id)):
            inputtype='uniprotswissprot'
        elif input_id.find('CHEMBL') == 0:
            inputtype='chembl'
        else:
            inputtype='hgnc_symbol'
    except:
        raise TypeError("Input needs to be entrezgene_id, hgnc_id, ensembl_peptide_id, ensembl_gene_id, uniprotswissprot or chembl id. \
                        If error persists, then likely input identifier does not exist on ensembl.")
    
    ensembl=gene_server
    attributes = ['uniprotswissprot','entrezgene_id','chembl','ensembl_peptide_id','ensembl_gene_id']
    columns = ['uniprot','entrez','chembl','ensembl_peptide_id','ensembl_gene_id']
    input_protein = ensembl.search(inputtype, input_id, attributes, columns)

    # If query returns no result and input type is 'uniprotswissprot' try to search again using 'hgnc_symbol' 
    if len(input_protein)==0 and inputtype=='uniprotswissprot':
        inputtype = 'hgnc_symbol'
        
        input_protein = ensembl.search(inputtype, input_id, attributes, columns)

    # Add HGNC symbol, gene type, description, HGNC ID to output DataFrame
    attributes = ['hgnc_symbol', 'gene_biotype', 'description', 'hgnc_id']
    columns = ['hgnc_symbol', 'gene_type', 'description', 'hgnc_id']
    hgnc_protein = ensembl.search(inputtype, input_id, attributes, columns)
    input_protein = pd.concat([input_protein, hgnc_protein], axis=1) 

    # HGNC data is gene-level (one value) — fill down to all peptide rows
    hgnc_cols = ['hgnc_symbol', 'gene_type', 'description', 'hgnc_id']
    for col in hgnc_cols:
        if col in input_protein.columns:
            input_protein[col] = input_protein[col].ffill()
    
    return input_protein


def compound_identifiers(input_id: int | str | dict[str, str | int], verbose: bool = False,
                          max_retries: int = 3, retry_backoff: float = 1.0) -> pd.DataFrame:
    """
    Retrieves compound identifiers and synonyms from PubChem.

    Sends a search query to PubChem using the input_id parameter as input.

    Parameters
    ----------
    input_id : int/dictionary/string
        Identifier of the compound(s) to search on the API. It can be one of the following types: \\
        CID (int), InChI, inchikey or smiles
    verbose : bool, default False
        If True, prints a message when the input resolves to a salt (or other
        multi-covalent-component structure) and gets replaced with its PubChem parent
        compound, and also prints each PubChem retry attempt as it happens. Desalting
        failures and PubChem call failures on the final retry attempt are always
        printed regardless of this flag, since they indicate something going wrong
        rather than a routine diagnostic.
    max_retries : int, default 3
        Number of attempts for each individual PubChem call before giving up and
        letting the failure propagate. PubChem's API intermittently errors/times out
        under load; this retries with exponential backoff rather than failing
        immediately on the first transient issue.
    retry_backoff : float, default 1.0
        Base delay in seconds between retry attempts (doubles each attempt: 1s, 2s,
        4s by default).

    Returns
    -------
    DataFrame
        Dataframe of the compound(s) containing the following values: \\
        cid, synonyms, inchi, inchikey, smiles, connectivity_smiles, iupac_name, molecular_formula, molecular_weight
    """

    if isinstance(input_id, dict):
        input_compound = pd.DataFrame.from_dict(input_id)
    else:
        try:
            # Tracks the input and which resolved value to compare it against
            original_identifier = None
            identifier_type = None

            if isinstance(input_id, int):
                cids = _pubchem_call_with_retry(pcp.get_cids, input_id, namespace='cid', domain='compound',
                                                 cids_type='parent', as_dataframe=False,
                                                 max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                c = _pubchem_call_with_retry(pcp.Compound.from_cid, cids[0],
                                              max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                original_identifier, identifier_type = input_id, 'cid'
            
            else:
                # Find input is an inchi
                if input_id.find('InChI=') == 0: 
                    cids = _pubchem_call_with_retry(pcp.get_cids, input_id, namespace='inchi', searchtype=None,
                                                     domain='compound', cids_type='parent', as_dataframe=False,
                                                     max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    c = _pubchem_call_with_retry(pcp.Compound.from_cid, cids[0],
                                                  max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    original_identifier, identifier_type = input_id, 'inchi'
                # Find input is a inchikey
                elif len(input_id) == 27 and input_id.find('-') == 14: 
                    cids = _pubchem_call_with_retry(pcp.get_cids, input_id, namespace='inchikey', searchtype=None,
                                                     domain='compound', cids_type='parent', as_dataframe=False,
                                                     max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    c = _pubchem_call_with_retry(pcp.Compound.from_cid, cids[0],
                                                  max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    original_identifier, identifier_type = input_id, 'inchikey'
                # If not any of the above, assume input is smiles
                else: 
                    cids = _pubchem_call_with_retry(pcp.get_cids, input_id, namespace='smiles', searchtype=None,
                                                     domain='compound', cids_type='parent', as_dataframe=False,
                                                     max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    c = _pubchem_call_with_retry(pcp.Compound.from_cid, cids[0],
                                                  max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
                    original_identifier, identifier_type = input_id, 'smiles'

            # CID case: comparing directly against c.cid
            if identifier_type == 'cid' and c.cid != original_identifier and verbose:
                print(f"Desalting triggered: input CID {original_identifier} resolved to a "
                      f"different compound (CID {c.cid}) - likely a salt or other "
                      f"multi-component structure.")

            # API PubChem Limit. conservative at 0.5 while 0.25 is the hard minimum 
            time.sleep(0.4)

            # Replace c.to_series() with PubChemServer to unify property keys
            pcs = PubChemServer()
            cid_str = str(c.cid)
            api_columns = pcs.get_columns(list(pcs.properties.values()))
            data = _pubchem_call_with_retry(pcs.get_compounds, cid_str, api_columns, namespace='cid',
                                             max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
            data = data.rename(columns=pcs.properties)

            # InChI/InChIKey/SMILES case: comparing the input identifier against record
            if (identifier_type is not None and identifier_type != 'cid'
                    and identifier_type in data.columns):
                resolved_identifier = data[identifier_type].iloc[0]
                if resolved_identifier != original_identifier and verbose:
                    print(f"Desalting triggered: input {identifier_type} "
                          f"'{original_identifier}' resolved to a different compound "
                          f"(CID {c.cid}, {identifier_type}='{resolved_identifier}') "
                          f"- likely a salt or other multi-component structure.")

            # Fetch synonyms
            synonyms_df = _pubchem_call_with_retry(pcs.get_synonyms, cid_str,
                                                    max_retries=max_retries, backoff_base=retry_backoff, verbose=verbose)
            synonyms = synonyms_df['Synonym'].iloc[0] if len(synonyms_df) > 0 else []
            #print(synonyms)
            # Add IUPAC identifiers into the synonyms list for search in other databases (e.g. ChEMBL)
            synonyms.extend([
                data['inchi'].iloc[0],
                data['inchikey'].iloc[0],
                data['smiles'].iloc[0],
                data['connectivity_smiles'].iloc[0],
                data['iupac_name'].iloc[0],
                ])
            data['synonyms'] = [synonyms]
            data['cid'] = c.cid
            data['input_id']=input_id

            # Convert dictionary to Dataframe
            input_compound = data
            cols = ['input_id'] + [c for c in input_compound.columns if c != 'input_id']
            input_compound = input_compound[cols]

        except Exception as e:
            raise TypeError(
                f"Failed to resolve compound identifier via PubChem: {e}\n"
                f"If this is a server/connection error (e.g. 'ServerBusyError', HTTP 503, timeout), "
                f"PubChem's own service is likely the cause. With service errors, CPIExtract's automatic "
                f"retries were exhausted, might need to wait for PubChem server stability. "
                f"Otherwise, confirm the input is a valid CID, InChI, InChIKey, or SMILES string."
            ) from e
    
    if 'inchikey' in input_compound.columns:
        input_compound = input_compound.copy()
        # Handle single row case
        if len(input_compound) == 1:
            inchikey_val = input_compound['inchikey'].iloc[0]
            input_compound['inchikey_fb'] = pcs.get_inchikey_first_block(inchikey_val)
        # Handle multiple rows case
        else:
            input_compound['inchikey_fb'] = input_compound['inchikey'].apply(
                lambda x: pcs.get_inchikey_first_block(x) if pd.notna(x) else None)
    else:
        # If no inchikey column exists, set to None
        input_compound['inchikey_fb'] = None

    return input_compound