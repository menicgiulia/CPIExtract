import time
from chembl_webresource_client.new_client import new_client
import pandas as pd

class ChEMBLServer(object):
    
    @staticmethod
    def identify_chembl_ids(compound: pd.DataFrame, pubchem_ids=True):
        chembl_total = []
        if pubchem_ids:
            # Search for all the chembl ids for the compound
            for _, row in compound.iterrows():
                syns = row['synonyms']
                if syns:
                    chembl_total.extend([syn for syn in syns if (syn and 'CHEMBL' in syn)])
        else:
            # Collect a list of ChEMBL ids from the PubChem names
            chemblEx=[] 
            chemblSyn=[] 
            # Search for all the chembl ids of synonyms for the compound
            for _, row in compound.iterrows():
                syns = row['synonyms']
                for syn in syns:
                    try:
                        molecule = new_client.molecule # This uses the chembl webresource client to find chembl ids
                        molS = molecule.filter(molecule_synonyms__molecule_synonym__iexact=syn)\
                                                .only(['molecule_chembl_id'])
                        molE = molecule.filter(pref_name__iexact=syn).only(['molecule_chembl_id'])
                        # Add all Exact hits
                        if len(molE) > 0:
                            chemblEx.append(molE[0]['molecule_chembl_id']) 

                        # Add all Synonym hits
                        if len(molS) > 0:
                            chemblSyn.append(molS[0]['molecule_chembl_id'])
                    except:
                        # API not working, wait 5 seconds
                        time.sleep(5)
            chembl_total = chemblEx + chemblSyn
            
        # If search is with ChEMBL or PubChem search was not successful 
        if not pubchem_ids or (chembl_total == [] and pubchem_ids):  
            chemblKey=[] 
            chemblSmi=[]  
            # Find the inchikey mapping to ChEMBL id
            inchikey = compound['inchikey'][0]
            if pd.notna(inchikey):
                try:
                    molecule = new_client.molecule
                    mols = molecule.filter(molecule_structures__standard_inchi_key=inchikey).only(['molecule_chembl_id'])
                    # Add all inchikey hits
                    if len(mols) > 0:
                        chemblKey.append(mols[0]['molecule_chembl_id'])
                except:
                    time.sleep(5)

            # Find the SMILES mapping to ChEMBL id
            smiles = compound['isomeric_smiles'][0]
            if pd.notna(smiles):
                try:
                    molecule = new_client.molecule
                    # Add all smiles hits
                    mols=molecule.filter(molecule_structures__canonical_smiles__connectivity=smiles)\
                                            .only(['molecule_chembl_id'])
                    if len(mols) > 0:
                        chemblSmi.append(mols[0]['molecule_chembl_id'])
                except:
                    time.sleep(5)

            # Merge all hits 
            chembl_total.extend(chemblKey + chemblSmi)

        # Obtain a single list of unique ids
        chembl_total = list(set(chembl_total))
        if len(chembl_total) > 0:
            return chembl_total
        else:
            return [None]

    @staticmethod
    def molecule_search(ids, inchi, chembl_id):

        attempts = 0
        max_attempts = 10
        while True:
            try:
                molecule = new_client.molecule
                molecule.set_format('json')
                record_via_client = molecule.filter(molecule_chembl_id__in=ids).only(['molecule_chembl_id','molecule_structures'])
                for record in record_via_client:
                    if record['molecule_structures']:
                        inchi.append(record['molecule_structures']['standard_inchi'])#
                        chembl_id.append(record['molecule_chembl_id'])
                break
            except: 
                # Avoid the chembl API error, and give up after 10 attempts
                time.sleep(5)
                attempts += 1
                if attempts >= max_attempts:
                    print(f"attempts ({max_attempts}) reached. Chembl API does not work.")
                    break