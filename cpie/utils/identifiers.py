from ..servers.BiomartServer import BiomartServer
import pubchempy as pcp
import pandas as pd

def protein_identifiers(input_id):
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
        elif len(input_id) == 6:
            inputtype='uniprotswissprot'
        elif input_id.find('CHEMBL') == 0:
            inputtype='chembl'
        else:
            inputtype='hgnc_symbol'
    except:
        raise TypeError("Input needs to be entrezgene_id, hgnc_id, ensembl_peptide_id, ensembl_gene_id, uniprotswissprot or chembl id. \
                        If error persists, then likely input identifier does not exist on ensembl.")
    
    ensembl = BiomartServer()
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

    return input_protein


def compound_identifiers(input_id):
    """
    Retrieves compound identifiers and synonyms from PubChem.

    Sends a search query to PubChem using the input_id parameter as input.

    Parameters
    ----------
    input_id : int/dictionary/string
        Identifier of the compound(s) to search on the API. It can be one of the following types: \\
        CID (int), InChI, inchikey or smiles

    Returns
    -------
    DataFrame
        Dataframe of the compound(s) containing the following values: \\
        cid, synonyms, inchi, inchikey, isomeric_smiles, canonical_smiles, iupac_name, molecular_formula, molecular_weight
    """

    if isinstance(input_id, dict):
        input_compound = pd.DataFrame.from_dict(input_id)
    else:
        try:
            if isinstance(input_id, int):
                # Retrieve info from PubChem using CID input
                c = pcp.Compound.from_cid(input_id) 
            
            else:
                # Find input is an inchi
                if input_id.find('InChI=') == 0: 
                    cids = pcp.get_cids(input_id, namespace='inchi', searchtype=None, as_dataframe=False)
                    c = pcp.Compound.from_cid(cids[0])
                # Find input is a inchikey
                elif len(input_id) == 27 and input_id.find('-') == 14: 
                    cids = pcp.get_cids(input_id, namespace='inchikey', searchtype=None, as_dataframe=False) 
                    c = pcp.Compound.from_cid(cids[0])
                # If not any of the above, assume input is smiles
                else: 
                    cids = pcp.get_cids(input_id, namespace='smiles', searchtype=None, as_dataframe=False) 
                    c = pcp.Compound.from_cid(cids[0])

            # Select info to keep from pubchem API into series
            data = c.to_series(properties=['isomeric_smiles','canonical_smiles','inchi','inchikey','synonyms','iupac_name','molecular_formula','molecular_weight']) 
        
            # Add IUPAC identifiers into the synonyms list for search in other databases (e.g. ChEMBL)
            data['synonyms'].append(data['inchi'])
            data['synonyms'].append(data['inchikey'])
            data['synonyms'].append(data['isomeric_smiles'])
            data['synonyms'].append(data['canonical_smiles'])
            data['synonyms'].append(data['iupac_name'])
            if isinstance(input_id, int):
                data['cid'] = input_id
            else:
                data['cid'] = cids[0]
            
            # Convert dictionary to Dataframe
            input_compound = data.to_frame().transpose()

        # Error if input is not the proper format or not found in PubChem
        except: 
            raise TypeError("Input needs to be CID, InChI, InChIKey, or SMILES. If error persists, then likely input identifier does not exist on PubChem.")
    
    return input_compound