#! python3
import json
import requests
import sys
import os
import re
import numpy as np
import pandas as pd
from operator import itemgetter
from itertools import count, groupby
import rosetta_paths


class albumin:

    def __init__(self, uniprotAC, output_path=rosetta_paths.default_output_path):
        self.uniprotAC = uniprotAC

        # sometimes the is a / and sometimes there is not
        if output_path[-1] == '/':
            self.out_path = output_path[:-1]
        else:
            self.out_path = output_path

        self.out_path = output_path

        # let's make sure the right directories exist
        if not os.path.isdir('{}/{}/'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}/'.format(self.out_path, self.uniprotAC))

        if not os.path.isdir('{}/{}/experimental_structures'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}/experimental_structures'.format(self.out_path, self.uniprotAC))

        if not os.path.isdir('{}/{}/cleaned_structures'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}/cleaned_structures'.format(self.out_path, self.uniprotAC))

        if not os.path.isdir('{}/{}/homology_models'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}/homology_models'.format(self.out_path, self.uniprotAC))

        if not os.path.isdir('{}/{}/predictions'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}/predictions'.format(self.out_path, self.uniprotAC))

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
            return('ERROR')
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
            print('processing entry clinvar entry {}'.format(entry_numbo))
            request_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={}&retmode=json'.format(entry)
            # sometimes the request fails, and we don't want to crash everything
            try:
                r = requests.get(request_url)
            except:
                print('failed to get variant{}'.format(entry))
                continue

            try:
                # and also sometimes this fails for some reason
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
        self.path_to_accession_folder = '{}/{}'.format(self.out_path, self.uniprotAC)
        # check if it is a folder, otherwise make it
        if not os.path.isdir(self.path_to_accession_folder):
            os.mkdir(self.path_to_accession_folder)
        # and dump the json there
        self.path_to_clinvar_variants = '{}/clinvar_sAA_variants.json'.format(self.path_to_accession_folder)
        with open(self.path_to_clinvar_variants, 'w') as variant_file:
            json.dump(self.clinvar_sAAsubs, variant_file)

        # and let's also make a file with some overall stats
        with open('{}/clinvar_stats.txt'.format(self.path_to_accession_folder), 'w') as stat_file:
            stat_file.write('clinvar has {} entries for {}, of these {} are snvs, and of these {} give a single aa sub \n (max entries parsed is 10.000)'.format(number_of_entries, self.gene_name,  number_of_snvs, len(self.clinvar_sAAsubs)))

        return(self.path_to_clinvar_variants)

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
                number_of_missense = number_of_missense + 1
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
        self.path_to_accession_folder = '{}/{}'.format(self.out_path, self.uniprotAC)
        # check if it is a folder, otherwise make it
        if not os.path.isdir(self.path_to_accession_folder):
            os.mkdir(self.path_to_accession_folder)
        # and dump the jsons there
        self.path_to_exac_variants = '{}/exac_sAA_variants.json'.format(self.path_to_accession_folder)
        with open(self.path_to_exac_variants, 'w') as variant_file:
            json.dump(self.exac_sAAsubs, variant_file)

        # and let's also make a file with some overall stats
        with open('{}/exac_stats.txt'.format(self.path_to_accession_folder), 'w') as stat_file:
            stat_file.write('exac has {} entries for {}, of these {} are missense variants, and {} give a single aa sub\n'.format(number_of_variants, self.gene_name, str(number_of_missense), len(self.exac_sAAsubs)))

        return(self.path_to_exac_variants)

    def pdb_map(self):
        '''this method looks for the experimental structure SIFTS mappings for a uniprot accessions,
        It writes all the maps to a file, and returns a list of recomended structures.
        (chosen by length, and for maximum coverage).
        note that the coverage map is only an approximation, and is based on the mapping info
        in sifts, not actual alignments. YMMV'''

        # first get all the pdb structures associated with the accession.
        # and we will use the sifts API

        requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/' + self.uniprotAC
        try:
            r = requests.get(requestURL, headers={"Accept": "application/json"})
            # we check that this is an entry at pdb
        except requests.exceptions.HTTPError:
            print('somethings wrong with the http request')
            print('no pdb structures for {}'.format(self.uniprotAC))

        all_pdbs = json.loads(r.text)
        number_of_structures = len(all_pdbs[self.uniprotAC]['PDB'])

        print('has {} mappings'.format(number_of_structures))
        # lets make a list of dictionaries, for each mapping
        mapping_list = []
        # this variable is just used to keep track of what entry is being processed
        mapping_number = 0
        # this dictionary records the structure_id + chain_id as keys, and the length as values
        structure_map_len_dict = {}
        # this dictionary records
        for mapping in all_pdbs[self.uniprotAC]['PDB']:
            pdb_id = mapping
            # all_pdbs[self.uniprotAC]['PDB'][mapping] is actually a list of dicts, so let's loop through the entries
            for element in all_pdbs[self.uniprotAC]['PDB'][mapping]:
                mapping_info = element
                # remember where you come from
                mapping_info['pdb_id'] = pdb_id
                chain_id = mapping_info['chain_id']
                structure_id = '{}_{}'.format(pdb_id, chain_id)
                # record the length of the mapping
                uniprot_start = mapping_info['unp_start']
                uniprot_end = mapping_info['unp_end']
                map_length = int(uniprot_end) - int(uniprot_start)
                pdb_map_len = mapping_info['end']['residue_number'] - mapping_info['start']['residue_number'] + 1
                mapping_info['map_length'] = map_length
                # record the mapping in the length in the structure_map_len_dict
                structure_map_len_dict[structure_id] = map_length
                mapping_info['pdb_map_len'] = pdb_map_len
                # check if this is the longest map
                # and add it to the list of mappings:
                mapping_list.append(mapping_info)
                # increment the mapping number
                mapping_number = mapping_number + 1

        # now choosing which mappings to go by, and what structures to use is not trivial.
        # ideally we want complete coverage of the uniprot sequence. But only once.

        # Let's try Inti's pandas solution
        # make a dataframe with every uniprot residue index

        df = pd.Series(np.arange(1, len(self.sequence)+1)).to_frame('uniprot')
        for structure in mapping_list:
            # determine the structure id
            structure_id = '{}_{}'.format(structure['pdb_id'], structure['chain_id'])
            # make a new dataframe with the structure id name
            # and populate it with nan values
            df[structure_id] = np.nan
            # take the partial of the original dataframe, that corresponds to the mapped residues
            # and make a copy
            # !!! now this is not completely halal
            # but we map the entire length of the pdb to the sequuntial numbers in the unp numbering scheme
            # because unp segments are often longer than pdb, segments, because of gaps in the alignment.
            # i cannot take that into account here,
            # and it also means that the structural coverage dataframe is not completely accurate.
            # but for structure selection it will have to do.
            partial = df.query('uniprot >= {} and uniprot <= {}'.format(structure['unp_start'], structure['unp_start'] + structure['pdb_map_len'] - 1)).copy()
            # in this partial dataframe, add a colomn for our structure
            # and assign the pdb numbers to those
            partial[structure_id] = np.arange(structure['start']['residue_number'], structure['end']['residue_number']+1)
            # and then update the dataframe
            df.update(partial)

        # now we have a dataframe with all the mappings.
        # based on this we recommend some structures,
        # that together gives the highest coverage
        total_coverage_df = df
        # and let's write the dataframe to a file, in case someone wants to have a look:
        # check if there is a folder first
        if not os.path.isdir('{}/{}'.format(self.out_path, self.uniprotAC)):
            os.mkdir('{}/{}'.format(self.out_path, self.uniprotAC))

        total_coverage_df.to_csv(path_or_buf='{}/{}/experimental_coverage_map.csv'.format(self.out_path, self.uniprotAC))

        # Let's start with the longest mapping, and then, looping through the rest, add
        # add structures if they cover residues in the uniprot sequence not previously covered
        # keep the ordering in this list. it's a list of doubles (structure_id, length)
        length_structure_order = []
        for entry in sorted(structure_map_len_dict.items(), key=itemgetter(1), reverse=True):
            length_structure_order.append(entry)

        # now as Inti proposed, we loop through each residue in the total_coverage_df
        # find the longest map that covers this residue,
        # check if it is part of the recommended structures list, and add it if it is not.
        recomended_structures = []
        for index, row in total_coverage_df.iterrows():
            # loop through the structures that have this residue present
            list_of_structures_covering_this_res = []
            # get the colouns that have a value:
            # This is obviously not the pandas way,
            # but i do not know pandas
            for cell in range(1, len(row)):
                # check if it has a value
                if not np.isnan(row[cell]):
                    # if it does, add it to the list of structures that cover this residue
                    list_of_structures_covering_this_res.append(total_coverage_df.dtypes.index[cell])
            # now find the longest map in this list:
            # we will do that by looping through the list_of_longest_maps,
            # and finding the first match:
            while True:
                for structure in length_structure_order:
                    if structure[0] in list_of_structures_covering_this_res:
                        # check if it is already in the list of recomended structures else add it
                        if structure[0] not in recomended_structures:
                            recomended_structures.append(structure[0])
                        # and break the while loop
                        break
                # otherwise there is no coverage, and we break anyway
                break

        # save the list of recomended structures as self.recomended_experimental_structures
        self.recomended_experimental_structures = recomended_structures
        self.total_coverage_df = total_coverage_df
        return(recomended_structures)

    def swiss_model_map(self):
        '''This method looks for homology models in the swissmodel repository,
        It is designed to be invoked after the pdb_map method, and looks for homology models,
        That cover areas with no experimental structure.'''

        # first go through the experimental coverage map and determine what residues are without experimental structure coverage
        without_structure = []

        for index, row in self.total_coverage_df.iterrows():
            # loop through the structures that have this residue present
            # get the colouns that have a value:
            # This is obviously not the pandas way,
            # but i do not know pandas
            # this variable watches for structure
            no_structure = True
            for cell in range(1, len(row)):
                # check if it has a value
                if not np.isnan(row[cell]):
                    # then this residue is covered
                    no_structure = False
            # if no structure has been found, we add the residue to
            # the list of residues without structure
            if no_structure is True:
                without_structure.append(index+1)

        # now convert this list into a a series of ranges, with some itertools black magic
        G = (list(x) for _, x in groupby(without_structure, lambda x, c=count(): next(c)-x))
        ranges_without_coverage_string = (",".join("-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G))
        # now there might no be anythn here. so only make this list if there is
        if len(ranges_without_coverage_string.split(',')) > 0:
            list_of_ranges_without_coverage = ranges_without_coverage_string.split(',')
        else:
            list_of_ranges_without_coverage = 0

        print('list_of_ranges_without_coverage')
        print(list_of_ranges_without_coverage)

        # let's keep the accepted models in this list
        list_of_homology_models = []

        # now loop through the ranges, and ask for structures from swissmodel
        # demanding a model for every range seems a little excessive.
        # some residues are just missing from the protein altogether, and the uniprot sequence is just inaccurate
        # let's go with ranges longer than 10 residues
        for stretch in list_of_ranges_without_coverage:
            if '-' in stretch:
                # then this is a proper stretch
                print(stretch)
                start, end = stretch.split('-')
                # no tiny stretches
                if int(end)-int(start) > 10:
                    # ask for a swissmodel covering this range
                    # and sort it by sequence id
                    requestURL = 'https://swissmodel.expasy.org/repository/uniprot/' + self.uniprotAC + '.json?range='+stretch+'?sort=seqid&provider=swissmodel'
                    r = requests.get(requestURL, headers={"Accept": "application/json"})
                    # we check that this is an entry at swismodel. Otherwise this will raise a 404 error.
                    print(requestURL)
                    if not r.ok:
                        r.raise_for_status()
                        continue

                    model_info = json.loads(r.text)
                    pdb_file_url = model_info['result']['structures'][-1]['coordinates']
                    # if we don't already have this in the list, append it
                    if pdb_file_url not in list_of_homology_models:
                        list_of_homology_models.append(pdb_file_url)

        print('here is the list of model urls')
        print(list_of_homology_models)

    def best_swiss_model(self):

        requestURL = 'https://swissmodel.expasy.org/repository/uniprot/' + self.uniprotAC + '.pdb?provider=swissmodel&sort=seqsim'
        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at swismodel. Otherwise this will raise a 404 error.
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # read the text as a json object
        # this will return a dictionary with 2 keys: 'query' and 'result'
        print(r.text)
