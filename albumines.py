#! python3
import json
import requests
import sys
import os


class albumin:

    def __init__(self, uniprotAC):
        self.uniprotAC = uniprotAC

    def get_sequenced(self):
        '''this method gets downloads the fasta file of the uniprot AC,
            and parses the file to store the sequence as self.uniprot_seq'''
        # the following information should be gotten from a lookup in uniprot

        # use the good old restful API
        requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/" + self.uniprotAC
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
        #requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/{}/isoforms".format(self.uac)
        #r = requests.get(requestURL, headers={"Accept": "application/json"})
        #if not r.ok:
        #    r.raise_for_status()
        #    sys.exit()
        ## the response_body a json object, let's parse it as a python dict.
        #self.isoforms_list = json.loads(r.text)

    def get_variants(self):
        '''This method gets the variation data from uniprot and stores it under self.exac_variants and self.clinvar_variants.
        For now it only records ClinVar and ExAC entries'''

        # let's put the variants in these dictionaries.
        self.exac_variants = {}
        self.clinvar_variants = {}

        # documentation of the proteins API: https://www.ebi.ac.uk/proteins/api/doc/#!/variation/getVariation

        requestURL = "https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=100&accession=" + self.uniprotAC
        r = requests.get(requestURL, headers={"Accept": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        responseBody = r.text

        # the response is a json file, with the variation information
        # read it as a python dict
        variation_dictionary_list = json.loads(responseBody)
        # they return a 1 length list for some reason
        variation_dictionary = variation_dictionary_list[0]

        # the actual variants are stored under variation_dictionary['features'], as a list, of dicts.
        # with the following keys:
        # dict_keys(['type', 'description', 'alternativeSequence', 'begin', 'end', 'xrefs', 'wildType', 'somaticStatus', 'cytogeneticBand', 'consequenceType', 'genomicLocation', 'association', 'clinicalSignificances', 'sourceType'])

        for i in range(0, len(variation_dictionary['features'])):
            ctype = variation_dictionary['features'][i]['consequenceType']
            wtAA = variation_dictionary['features'][i]['wildType']
            # alt_seq is the just the alt_aa in case of sAA subs.
            alt_seq = variation_dictionary['features'][i]['alternativeSequence']
            # the 'begin' and 'end' value is always the same. So let's go with 'begin'.
            pos = variation_dictionary['features'][i]['begin']

            # for now we only care about exac and clinvar variants
            # they will have a key called 'xrefs'
            if 'xrefs' in variation_dictionary['features'][i]:
                # one variant can have several xrefs
                for element in variation_dictionary['features'][i]['xrefs']:
                    # check if the variant is referenced in ExAC or ClinVar
                    # and put them in the dicts.
                    # I would have like to have the substiturion as a key, and the ids as values.
                    # but since several variants can have the same sub, we have to do it the other way around.
                    if element['name'] == 'ClinVar':
                        self.clinvar_variants[element['id']] = [wtAA + pos + alt_seq, ctype]
                        # TO DO: add clinical significance, maybe
                    elif element['name'] == 'ExAC':
                        self.exac_variants[element['id']] = [wtAA + pos + alt_seq, ctype]
                        # To DO: add frequuencies
                        # this will most likely be a seperate lookup with the eXac API

        # finally dump the variant information in some json files
        # so they are accesible for the score parser.
        variant_dict = {}
        variant_dict['exac_variants'] = self.exac_variants
        variant_dict['clinvar_variants'] = self.clinvar_variants

        self.path_to_accession_folder = 'uniprot_accessions/{}'.format(self.uniprotAC)
        # check if it is a folder, otherwise make it
        if not os.path.isdir(self.path_to_accession_folder):
            os.mkdir(self.path_to_accession_folder)
        # and dump the jsons there
        with open('{}/variants.json'.format(self.path_to_accession_folder), 'w') as variant_file:
            json.dump(variant_dict, variant_file)

    def get_swismodel(self):
        # same as pdb, but for homology models.
        pass
