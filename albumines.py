#! python3
import json
import requests
import sys
import os
import re


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
        Only variants recorded as giving a the single aa sub is recorded.
        Note that this method is quite slow to complete, as it relies on
        making thousands of http requests through the clinvar API'''

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
            "TYR": "Y",
            "TER": "stop"
        }


        request_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={}[gene]&retmax=10000&retmode=json'.format(self.gene_name)
        r = requests.get(request_url, headers={"accept": "application/json"})
        variant_dict = json.loads(r.text)

        number_of_entries = variant_dict["esearchresult"]["count"]
        print(number_of_entries)
        entries_list = variant_dict["esearchresult"]["idlist"]

        entry_numbo = 0
        # put all the relevant entries in this dictionary:
        self.clinvar_sAAsubs = {}
        for entry in entries_list:
            entry_numbo = entry_numbo + 1
            print('processing entry {}'.format(entry_numbo))
            request_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={}&retmode=json'.format(entry)
            r = requests.get(request_url)
            try:
                # sometimes this fails for some reason
                variant_info = json.loads(r.text)
            except json.decoder.JSONDecodeError:
                print('couldn\'t get info on variant{}'.format(entry))
                continue

            variant_title = variant_info['result'][entry]['title']
            variant_type = variant_info['result'][entry]['variation_set'][0]['variant_type']
            # note there is other information here, if you are interested. such as review status, and last time it was reviewed
            clinical_significance = variant_info['result'][entry]['clinical_significance']['description']
            review_status = variant_info['result'][entry]['clinical_significance']['review_status']
            trait_set = variant_info['result'][entry]['trait_set'][0]['trait_name']
            # let's tease out just the SNVs
            number_of_snvs = 0
            if variant_type == 'single nucleotide variant':
                number_of_snvs = number_of_snvs + 1
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
                    to_aa_single = aminocodes[to_aa.upper()]

                    # we use regex to find the variant position
                    # match one or more digits in the aa change
                    digits = re.compile(r'\d+')
                    position = digits.search(aa_change).group()

                    # now we are ready to record the mutation
                    # if it gives an aa_sub and not a stop codon
                    if from_aa_single!=to_aa_single and to_aa_single!="stop":
                        mut_spec = from_aa_single+position+to_aa_single

                        self.clinvar_sAAsubs[entry] = {
                            'title': variant_title,
                            'clinical_significance': clinical_significance,
                            'review_status': review_status,
                            'trait_set': trait_set,
                            'mutation': mut_spec
                        }

        # dump the variation dict in as a json in the uniprotac folder
        self.path_to_accession_folder = 'uniprot_accessions/{}'.format(self.uniprotAC)
        # check if it is a folder, otherwise make it
        if not os.path.isdir(self.path_to_accession_folder):
            os.mkdir(self.path_to_accession_folder)
        # and dump the jsons there
        with open('{}/clinvar_sAA_variants.json'.format(self.path_to_accession_folder), 'w') as variant_file:
            json.dump(self.clinvar_sAAsubs, variant_file)

        # and let's also make a file with some overall stats
        with open('{}/clinvar_stats.txt'.format(self.path_to_accession_folder), 'w') as stat_file:
            stat_file.write('clinvar has {} entries for {}, of these {} are snvs, and of these {} give a single aa sub \n (max entries parsed is 10.000)'.format(number_of_entries, self.gene_name,  number_of_snvs, len(self.clinvar_sAAsubs)))

        return(self.clinvar_sAAsubs)

    def get_exac_variants(self):
        '''This method gets the exac variants for a uniprot AC.
        it only works if the gene name has been found'''

        # we will use the exac restful API
        # documented here: http://exac.hms.harvard.edu/

        # the single_aa_sub entries will be stored here:
        self.exac_sAAsubs = {}

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
            "TYR": "Y",
            "TER": "stop"
        }

        request_url = 'http://exac.hms.harvard.edu//rest/awesome?query={}&service=variants_in_gene'.format(self.gene_name)
        print(request_url)
        # we must allow redirects, since we are looking up gene names
        r = requests.get(request_url, headers={"accept": "application/json"}, allow_redirects=True)
        print(r.url)
        variant_dict = json.loads(r.text)
        number_of_variants = len(variant_dict)
        # let's record the number of missense
        number_of_missense = 0
        # let's parse the entries in turn.
        # thankfully all the information we need is contained in this dict.
        # so we avoid requesting info for every single variant (I am looking at you Clinvar)
        for entry in variant_dict:
            if entry["category"] == "missense_variant":
                number_of_missense += 1
                # stops and indels are included here. So to get just the single AA subs
                if entry["major_consequence"] == "missense_variant":
                    # the protein level consequence is stored here
                    HGVSp = entry["HGVSp"]
                    # like for clinvar we parse the entry:
                    # the single_aa_sub entries will be stored here:
                    aa_change = HGVSp.split('.')[-1]
                    from_aa = aa_change[0:3]
                    to_aa = aa_change[-3:]

                    # and in single letter code, using the dict above
                    from_aa_single = aminocodes[from_aa.upper()]
                    to_aa_single = aminocodes[to_aa.upper()]

                    # we use regex to find the variant position
                    # match one or more digits in the aa change
                    digits = re.compile(r'\d+')
                    position = digits.search(aa_change).group()
                    mut_spec = from_aa_single + position + to_aa_single
                    print(mut_spec)

                    # make an entry in the exac variant dict
                    self.exac_sAAsubs[entry["variant_id"]] = {
                        'allele_count': entry["allele_count"],
                        'allele_freq': entry["allele_freq"],
                        'HGVSp': entry["HGVSp"],
                        'mut': mut_spec
                    }
        # Let's dump it as a json in the proper folder again
        self.path_to_accession_folder = 'uniprot_accessions/{}'.format(self.uniprotAC)
        # check if it is a folder, otherwise make it
        if not os.path.isdir(self.path_to_accession_folder):
            os.mkdir(self.path_to_accession_folder)
        # and dump the jsons there
        with open('{}/exac_sAA_variants.json'.format(self.path_to_accession_folder), 'w') as variant_file:
            json.dump(self.exac_sAAsubs, variant_file)

        # and let's also make a file with some overall stats
        with open('{}/exac_stats.txt'.format(self.path_to_accession_folder), 'w') as stat_file:
            stat_file.write('exac has {} entries for {}, of these {} are missense variants, and {} give a single aa sub\n'.format(number_of_variants, self.gene_name, str(number_of_missense), len(self.exac_sAAsubs)))

        return(self.exac_sAAsubs)

    def pdb_map(self):
        '''this method should make something like a coverage map, of the uniprot ac,
        in terms of pdb structures, and what residues they cover.
        It should also look in the swissprot repo'''

        # first get all the pdb structures associated with the accession.
        # and we will use the sifts API

        requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/' + self.uniprotAC
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at pdb
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        all_pdbs = json.loads(r.text)
        print(len(all_pdbs))
        print(len(all_pdbs))
