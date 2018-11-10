#! python3
import requests
imp
# this file defines the structure class.
# it will be used by the albumin


class structure:

    def __init__(self, UniprotAC):
        # the structure will belong to a uniprot AC
        # it might as well know which one
        self.uniprotac = UniprotAC

    # the most important method for the structureclass, is looking through the pdb for structures
    def get_pdb(self, uniprotac):
        '''this method looks in the pdb for the best experimental structure
        and downloads it'''

        # for this we will use the PDBe REST API.
        # This request is based on SIFTS
        # for more information on that see: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html

        requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/' + uniprotac

        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at pdb
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # read the text as a json object
        # this will return a dictionary, with the uniprotid as the only key
        # the value will be a list where each pdb structure is an element
        structure_info = json.loads(r.text)
        structure_list = structure_info[uniprotac]

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

        # we can do coverage too in the form of an alignment to self.sequence

    # it will have a path, as the file will be downloaded
    self.path = ''
    # it will have a sequence
    self.sequence = 0
    # it will cover a certain amount of the sequence
    self.coverage = 0

    # it will contain an alignment to the 'master' fasta seq

    # method for cleaning
    # method for mutfiling
