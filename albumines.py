#! python3
import json
import requests
import sys


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


    def get_swismodel(self):
        # same as pdb, but for homology models.
        pass
