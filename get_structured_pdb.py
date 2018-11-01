import sys, requests, json

# this script should provide a function for getting the 'best' structure
# from the pdb (for a provided uniprotid). Let's go with the single best for now.


def get_structured_pdb(uniprotid):
    '''This function takes a uniprot ID and fetches the 'best'
    structure available from the PDB'''

    # for this we will use the PDBe REST API.
    # This request is based on SIFTS
    # for more information on that see: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html

    requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/' + uniprotid

    r = requests.get(requestURL, headers={"Accept": "application/json"})
    # we check that this is an entry at pdb
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # read the text as a json object
    # this will return a dictionary, with the uniprotid as the only key
    # the value will be a list where each pdb structure is an element
    structure_info = json.loads(r.text)
    structure_list = structure_info[uniprotid]

    # for now let's just take the first structure in the list
    # PDBe has it ordered so that the 'best' is first
    # each element in the structure list is a dictionary.
    best_structure = structure_list[0]['pdb_id']
    chain_id = structure_list[0]['chain_id']
    # coverage seems like a key parameter, let's record that too
    coverage = structure_list[0]['coverage']

    # donwload the structure
    # for now i will just download the entire structures
    # isolating chains, and cleaning it up for calculations, are going to have to wait.
    structureURL = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(best_structure)
    r = requests.get(structureURL)
    # this way the model will be saved as uniprotAC_template_provider.pdb
    with open('experimental_structures/{}.pdb'.format(best_structure), 'w') as f:
        f.write(r.text)

    return best_structure, chain_id, coverage


if __name__ == '__main__':
    print(get_structured_pdb(sys.argv[1]))
