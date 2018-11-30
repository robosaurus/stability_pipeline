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
                        self.clinvar_variants[element['id']] = [wtAA + pos + alt_seq, ctype, variation_dictionary['features'][i]['clinicalSignificances']]
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

    def get_gene_name(self):
        '''This method gets the name of the gene for a uniprot accession and stores it under self.gene_name'''

        # first we need to get a gene name from uniprot.
        requestURL = 'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={}'.format(self.uniprotAC)
        r = requests.get(requestURL, headers={"Accept": "application/json"})

        if not r.ok:
            r.raise_for_status()
            print('entry {} not found at uniprot'.format(self.uniprotAC))
            return

        uniprot_entry_dict = json.loads(r.text)
        try:
            self.gene_name = uniprot_entry_dict[0]['gene'][0]['name']['value']
        except:
            print('gene_name_not_found')
            return
        return(self.gene_name)

    def get_clinvar_variants(self):
        '''this method gets the clinvar variant for a uniprot accession,
        by looking up the gene name in the clinvar database.
        Only variants recorded as giving a the single aa sub is recorded.'''

        # we are gonna need this dictionary:

        aminocodes = {
            "ALA": "A",
            "CYS": "C",
            "ASP": "D",
            "GLU": "E",
            "PHE": "F",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LYS": "K",
            "LEU": "L",
            "MET": "M",
            "ASN": "N",
            "PRO": "P",
            "GLN": "Q",
            "ARG": "R",
            "SER": "S",
            "THR": "T",
            "VAL": "V",
            "TRP": "W",
            "TYR": "Y"
        }


        request_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={}[gene]&retmax=10000&retmode=json'.format(self.gene_name)
        r = requests.get(request_url, headers={"accept": "application/json"})
        variant_dict = json.loads(r.text)

        number_of_entries = variant_dict["esearchresult"]["count"]
        print(number_of_entries)
        entries_list = variant_dict["esearchresult"]["idlist"]

        # put all the relevant entries in this dictionary:
        self.clinvar_sAAsubs = {}
        for entry in entries_list:
            request_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={}&retmode=json'.format(entry)
            r = requests.get(request_url)
            variant_info = json.loads(r.text)

            variant_title = variant_info['result'][entry]['title']
            variant_type = variant_info['result'][entry]['variation_set'][0]['variant_type']
            # note there is other information here, if you are interested. such as review status, and last time it was reviewed
            clinical_significance = variant_info['result'][entry]['clinical_significance']['description']
            review_status = variant_info['result'][entry]['clinical_significance']['review_status']
            trait_set = variant_info['result'][entry]['trait_set'][0]['trait_name']
            print(clinical_significance, trait_set, review_status)
            # let's tease out just the SNVs
            number_of_snvs = 0
            if variant_type == 'single nucleotide variant':
                number_of_snvs += 1
                # and from these let's look at the ones recorded as giving an aa sub.
                if '(p.' in variant_title:
                    # determine the substitution:
                    aa_change = variant_title.split()[-1]
                    from_aa = aa_change[3:6]
                    to_aa = aa_change[-4:-1]
                    # sometimes there is no change:
                    if '=' in to_aa:
                        to_aa = from_aa
                    # and in single letter code, using the dict above
                    from_aa_single = aminocodes[from_aa.upper()]
                    to_aa_single = aminocodes[to_aa.lower()]

                    position = variant_title.split()


                    print(from_aa, to_aa)
                    self.clinvar_sAAsubs[entry] = {
                        'title': variant_title,
                        'clinical_significance': clinical_significance,
                        'review_status': review_status,
                        'trait_set': trait_set
                    }



    def pdb_map(self):
        '''this method should make something like a coverage map, of the uniprot ac,
        in terms of pdb structures, and what residues they cover'''

        # first get all the pdb structures associated with the accession.
        # and we will use the sifts API
        requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/' + self.uniprotAC

        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at pdb
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        all_pdbs = json.loads(r.text)
