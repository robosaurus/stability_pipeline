#! python3
import json
import requests

class albumin:

    def __init__(self, uniprotAC):
        self.uac = uniprotAC

    def get_sequenced(self):
        '''this method gets downloads the fasta file of the uniprot AC,
            and parses the file to store the sequence as self.uniprot_seq'''
        # the following information should be gotten from a lookup in uniprot

        # use the good old restful API
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/" + self.uac
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # the response_body a json object, let's parse it as a python dict.
        self.response_dict = json.loads(r.text)

        # there is a bunch of usefule information here. I will put some of it in variables
        # for easy acces
        # this self.sequence will be a master sequence of sorts.
        # this is the numbering we will use for the database and it will be the sequence
        # we will align the structure sequnces to
        self.sequence = self.response_dict['sequence']['sequence']  # get the sequence from uniprot
        self.length = self.response_dict['sequence']['length']  # get from uniprot

        # amelie specifically asked us to look out for isoforms.
        # so here they are, though actually reacting on it not implemented yet:
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/{}/isoforms".format(self.uac)
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        # the response_body a json object, let's parse it as a python dict.
        self.isoforms_list = json.loads(r.text)


    def get_pdb(self):
        '''this method looks in the pdb for the best experimental structure
        and downloads it'''

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
 
                # we can do coverage too in the form of an alignment to self.sequence

    def get_swismodel(self):
        # same as pdb, but for homology models.
        pass
    
    def make_mutfiles(self, structure):
        # outputs the mutfiles for a given structures
        pass
