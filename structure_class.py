#! python3
import requests
import json
import sys
import subprocess
import os
# and import the file paths to rosetta and utility scripts
import rosetta_paths
import pdb_to_fasta_seq


# this file defines the structure class.
# it will be used by the albumin


class structure:

    def __init__(self, UniprotAC='no_ac', output_path='output/'):
        # the structure will belong to a uniprot AC
        # it might as well know which one
        self.uniprotac = UniprotAC

        # sometimes i put a backslash and sometimes they don't
        if output_path[-1] == '/':
            self.out_path = output_path[:-1]
        else:
            self.out_path = output_path

        # let's make the right directories
        if not os.path.isdir('{}/{}/experimental_structures'.format(self.out_path, self.uniprotac)):
            os.mkdir('{}/{}/experimental_structures'.format(self.out_path, self.uniprotac))

        if not os.path.isdir('{}/{}/cleaned_structures'.format(self.out_path, self.uniprotac)):
            os.mkdir('{}/{}/cleaned_structures'.format(self.out_path, self.uniprotac))

        if not os.path.isdir('{}/{}/homology_models'.format(self.out_path, self.uniprotac)):
            os.mkdir('{}/{}/homology_models'.format(self.out_path, self.uniprotac))

        print(rosetta_paths.path_to_rosetta)

        # the most important method for the structureclass, is getting structures from the pdbe
    def get_pdb(self, structureID_chainID):
        structure_id, chain_id = structureID_chainID.split('_')

        structureURL = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(structure_id)
        r = requests.get(structureURL)
        # this way the model will be saved as uniprotAC_pdbname.pdb
        # and this naming convention will be stored as self.sys_name
        self.sys_name = '{}_{}'.format(self.uniprotac, structure_id)
        path_to_pdbfile = self.out_path + '/experimental_structures/{}.pdb'.format(self.sys_name)
        with open(path_to_pdbfile, 'w') as pdb_file:
            pdb_file.write(r.text)

        self.chain_id = chain_id
        self.path = path_to_pdbfile

    def get_best_pdb(self):
        '''this method looks in the pdb for the best experimental structure
        and downloads it, 'best' as defined by the pdb themselves'''

        # for this we will use the PDBe REST API.
        # This request is based on SIFTS
        # for more information on that see: https://www.ebi.ac.uk/pdbe/api/doc/sifts.html

        requestURL = 'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/' + self.uniprotac

        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at pdb
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        # read the text as a json object
        # this will return a dictionary, with the uniprotid as the only key
        # the value will be a list where each pdb structure is an element
        structure_info = json.loads(r.text)
        structure_list = structure_info[self.uniprotac]

        # for now let's just take the first structure in the list
        # PDBe has it ordered so that the 'best' is first
        # each element in the structure list is a dictionary.
        self.best_structure = structure_list[0]['pdb_id']
        self.pdb_id = structure_list[0]['pdb_id']
        self.chain_id = structure_list[0]['chain_id']
        # coverage seems like a key parameter, let's record that too
        self.coverage = structure_list[0]['coverage']

        # donwload the structure
        # for now just download the entire structures
        structureURL = 'http://www.rcsb.org/pdb/files/{}.pdb'.format(self.best_structure)
        r = requests.get(structureURL)
        # this way the model will be saved as uniprotAC_pdbname.pdb
        # and this naming convention will be stored as self.sys_name
        self.sys_name = '{}_{}'.format(self.uniprotac, self.best_structure)
        path_to_pdbfile = '{}/experimental_structures/{}.pdb'.format(self.out_path, self.sys_name)
        with open(path_to_pdbfile, 'w') as pdb_file:
            pdb_file.write(r.text)

        self.path = path_to_pdbfile

    # We will not be interested in everything in that pdb, only certain chains
    # and we need to clean up the pdb before running Rosetta on it.
    # luckily rosetta provides just such a tool, for chain isolation and clean up.
    def clean_up_and_isolate(self, path_to_pdb, chains):
        '''this method takes a pdb, and cleans it up for rosetta, and isolates the
        relevant chains. It uses the clean_pdb.py tool bundled with the release
        version of Rosetta'''

        # let's just call the script from the shell.
        # it may seem a little clumsy to call a python script
        # from a shell command spawned from another python script.
        # but this way i do not have to modify their script :)
        # also, they use python2

        # first we find the clean_pdb.py script
        # the path is imported from rosetta_paths.py module
        path_to_clean_pdb = rosetta_paths.path_to_clean_pdb

        shell_command = 'python2 {} {} {}'.format(path_to_clean_pdb, path_to_pdb, chains)
        print('here is some output from the clean_pdb.py script')
        subprocess.call(shell_command, cwd='{}/cleaned_structures/'.format(self.out_path), shell=True)
        print('end of output from clean_pdb.py')

        self.path_to_cleaned_pdb = '{}/{}/cleaned_structures/{}_{}.pdb'.format(self.out_path, self.uniprotac, self.sys_name, chains)
        # this script also produces a fasta, that we can convinently use to set the sequence
        # of the structure
        path_to_cleaned_fasta = '{}/{}cleaned_structures/{}_{}.fasta'.format(self.out_path, self.uniprotac, self.sys_name, chains)
        fasta_file = open(path_to_cleaned_fasta, 'r')
        fasta_lines = fasta_file.readlines()
        fasta_file.close()
        # let's store the fasta sequence in self.fasta_seq
        self.fasta_seq = ''
        # the first line is always the '>' line
        for line in fasta_lines[1:]:
            self.fasta_seq = self.fasta_seq + line.strip()

        # now here should be some code for aligning the structure sequence to the uniprot sequence
        # determine the coverage (is this enough structure?)
        # and determine a numbering
        return(self.path_to_cleaned_pdb)

    def get_swiss_model(self, get_model_number=0):
        '''this method gets the best available homology model,
        from the swissmodel repository.
        Swissmodel provides a list of models, the get_model_number variable
        determines which model index to fetch (counting from 0). The default is to get the first (the 'best')'''

        # we will rquest models from the swissmodel repository, through their restful API
        # documented here: https://swissmodel.expasy.org/docs/repository_help
        # or checkout the more fun version of the docs here:
        # https://swissmodel.expasy.org/docs/smr_openapi

        requestURL = 'https://swissmodel.expasy.org/repository/uniprot/{}.json?provider=swissmodel'.format(self.uniprotac)

        r = requests.get(requestURL, headers={"Accept": "application/json"})
        # we check that this is an entry at swissmodel
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        structure_info = json.loads(r.text)
        # all the structures are available as a list here:
        list_of_models = structure_info['result']['structures']
        print('swissmodel has {} models in the repo for uniprot accession {}'.format(len(list_of_models), self.uniprotac))
        print('getting the one at index {}'.format(get_model_number))
        model_info = list_of_models[get_model_number]
        # the model_info is a dictionary with the following keys:
        # dict_keys(['similarity', 'gmqe', 'oligo-state', 'crc64', 'coverage', 'alignment', 'md5', 'from', 'qmean_norm', 'coordinates', 'to', 'identity', 'template', 'provider', 'qmean', 'method'])

        # first let's record some information on the structure
        # firstly the sequence id:
        self.seq_id = model_info['identity']
        # coverage seems like a key parameter
        self.coverage = model_info['coverage']
        # and of course the template that the model is based on
        self.template = model_info['template']
        # they even provide an alignment. nifty :)
        self.alignment = model_info['alignment']
        # let's get some model evaluation metrics too. The qmean, and the normalised qmean.
        # I assume this is qmean4
        # for information on qmean evaluation see:
       # Benkert, P., Tosatto, S.C.E. and Schomburg, D. (2008). "QMEAN: A comprehensive scoring function for model quality assessment." Proteins: Structure, Function, and Bioinformatics, 71(1):261-277.
        self.qmean = model_info['qmean']
        self.qmean_norm = model_info['qmean_norm']

        # we assign the structure a systematic name, based on uniprotAC, template used and the sequence ID,
        self.sys_name = '{}_{}_sm{}'.format(self.uniprotac, self.template, self.seq_id)
        print(self.sys_name)

        path_to_pdbfile = 'homology_models/{}.pdb'.format(self.sys_name)

        model_url = model_info['coordinates']
        r = requests.get(model_url)
        with open(path_to_pdbfile, 'w') as f:
            f.write(r.text)

        self.path = path_to_pdbfile

        # Note that the swissmodel repo is constantly updated, and models can dissapear from their servers.

    def make_mutfiles(self):
        '''this function makes Rosetta mutfiles, specifying all the possible AA substitutions
        for the structure. The .clean_up_and_isolate() method should be run first.
        The input is a fasta sequence, and the output is a folder in rosetta_runs/<structure_name>/mutfiles/'''

        # it will run in it's own folder. that we can make a rosetta_runs/self.sys_name folder
        # this is where we will put the files
        self.path_to_run_folder = 'rosetta_runs/{}'.format(self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(self.path_to_run_folder):
            os.mkdir(self.path_to_run_folder)

        path_to_mutfiles = '{}/mutfiles/'.format(self.path_to_run_folder)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(path_to_mutfiles):
            os.mkdir(path_to_mutfiles)

        # write the files
        # each residue will be a single file, with the 20 possible AA subs (including itself)
        # rosetta counts from 1 and not 0 i think
        # also watch out for the extra character in self.fasta_seq.
        # writing mutfiles for:
        print(self.fasta_seq)
        for residue_number in range(1, len(self.fasta_seq)):
            # pad the mutfile numbers with 0'es, to make it more traceable, since the sbatch makes an ls call.
            mutfile = open(path_to_mutfiles+'mutfile{:0>5}'.format(str(residue_number)), 'w')
            mutfile.write('total 20\n')
            # and then a line for each type of AA
            for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
                mutfile.write('1\n')
                mutfile.write(self.fasta_seq[residue_number-1] + ' ' + str(residue_number) + ' ' + AAtype + '\n')
            mutfile.close()

    def rosetta_sbatch_relax(self, structure_path='defaults to self.path_to_cleaned'):
        '''this function writes an sbatch script, specifying the relaxation of the sctructure.
        it returns the path to the sbatch script.
        The flags for the relaxation are taken from rosetta_parameters/relax_flagfile'''

        structure_path = self.path_to_cleaned_pdb

        self.path_to_run_folder = '{}/uniprot_accessions/{}/rosetta_runs/{}'.format(self.out_path, self.uniprotac, self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(self.path_to_run_folder):
            os.makedirs(self.path_to_run_folder)

        path_to_sbatch = '{}/rosetta_relax.sbatch'.format(self.path_to_run_folder)
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=pre_relax_rosetta
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching rosetta relax
        {}bin/relax.linuxgccrelease -s {} -relax:script {}/cart2.script @{}/relax_flagfile
    '''.format(rosetta_paths.path_to_rosetta, structure_path, rosetta_paths.path_to_parameters, rosetta_paths.path_to_parameters))
        sbatch.close()
        print(path_to_sbatch)
        return path_to_sbatch

    def rosetta_relax(self):
        '''this function sruns a rosetta relaxation of the structure,
        it follows the recomended protocol for a pre-relaxation prior to
        a cartesian_ddg run'''

        # it will run in it's own folder. that we can make a rosetta_runs/self.sys_name folder
        self.path_to_run_folder = '{}/{}/rosetta_runs/{}'.format(self.out_path, self.uniprotac, self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(self.path_to_run_folder):
            os.mkdir(self.path_to_run_folder)

        # we will submit the sbatch from this folder.
        # so the paths should be relative to that
        self.path_to_cleaned_pdb = '../../cleaned_structures/{}_{}.pdb'.format(self.sys_name, self.chain_id)

        # and then we run the relax app
        # from the appropriate rosetta run folder
        # the relax protocol is taken from the rosetta docs: https://www.rosettacommons.org/docs/latest/cartesian-ddG
        # and it is dependent on the cart2.script file in the rosetta_parameters folder
        # remember to change nstructs to 20, when you are done testing
        shell_command = 'srun --mem-per-cpu=5000M {}/bin/relax.linuxgccrelease -s {} -use_input_sc \
-constrain_relax_to_start_coords \
-ignore_unrecognized_res \
-nstruct 1 \
-relax:coord_constrain_sidechains  \
-relax:cartesian \
-score:weights ref2015_cart \
-relax:min_type lbfgs_armijo_nonmonotone \
-out:suffix _bn15_calibrated \
-relax:script ../../rosetta_parameters/cart2.script'.format(self.path_to_rosetta, self.path_to_cleaned_pdb)
        print('calling to the shell:{}'.format(shell_command))
        subprocess.call(shell_command, shell=True,  cwd=self.path_to_run_folder)

    def write_rosetta_cartesian_sbatch(self):
        '''this function writes an sbatch file.
        it should only be called after clean, make_mutfiles, and rosetta_relax'''

        # you should change to a more reliable way of specifying the structure. Maybe a selection among the 20 in rosetta_relax
        # the memory spec is based on some brief testing. 200 M is not enough 1000 is. (based on a ~800 residue protein)
        # turns out it is not enough at all!!!!
        # change nstructs to 3, when you are done testing.
        path_to_sbatch = '{}/rosetta_cartesian_saturation_mutagenesis.sbatch'.format(self.path_to_run_folder)
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=Rosetta_cartesian_ddg
#SBATCH --array=0-{}
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab
LST=(`ls mutfiles/mutfile*`)
OFFSET=0
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta
{}/bin/cartesian_ddg.linuxgccrelease -database {} -s {} -fa_max_dis 9.0 -ddg::dump_pdbs true -ddg:iterations 1 -ddg:mut_file ${{LST[$INDEX]}} -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID -score:weights beta_nov16_cart -ddg:mut_only -ddg:bbnbrs 1 -beta_cart -ddg:mut_only
    '''.format(len(self.fasta_seq), self.path_to_rosetta, self.path_to_rosetta[0:-7]+'/database/', '*_bn15_calibrated*.pdb'))
        sbatch.close()
        print(path_to_sbatch)
        return path_to_sbatch

    def write_parse_ddg_sbatch(self):
        '''this function writes an sbatch script, with a call to parse_rosetta_ddgs.
        This is the easiest way to submit a job with a dependency'''

        score_sbatch_path = '{}/parse_ddgs.sbatch'.format(self.path_to_run_folder)
        score_sbatch = open(score_sbatch_path, 'w')
        score_sbatch.write('''#!/bin/sh
#SBATCH --job-name=collect_rosetta_ddgs
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --partition=sbinlab

# This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs
# it will output a file in the prediction_files/ folder.
python3 parse_rosetta_ddgs.py {} {} {} {} {}
'''.format(self.path_to_mapping_dict, self.sys_name, self.chain_id, self.fasta_seq, self.uniprotac))
        score_sbatch.close()

        return score_sbatch_path

    def muscle_align_to_uniprot(self):
        '''this method uses the muscle application to make an alignment between the structure sequence
        and the uniprot sequence'''

        # first determine the sequence of the structure
        self.fasta_seq = pdb_to_fasta_seq.pdb_to_fasta_seq(self.path_to_cleaned_pdb)
        # first make a fasta that contains the structure sequence and the uniprot sequence
        # first check that the folder exists
        if not os.path.isdir('{}/{}/{}'.format(self.out_path, self.uniprotac, self.sys_name)):
            os.mkdir('{}/{}/{}'.format(self.out_path, self.uniprotac, self.sys_name))
        # then write the fasta file, that will be used as input for muscle
        with open('{}/{}/{}/fasta_file.fasta'.format(self.out_path, self.uniprotac, self.sys_name), 'w') as fasta_file:
            fasta_file.write('did this work?')

    def align_pdb_to_uniprot(self):
        '''This method gets the alignment of the pdb residues to the uniprot sequence.
        The alignment is fetched from SIFTS.
        A dictionary is stored under self.pdb_to_uniprot, that has the chain IDs as keys and a dictionary as value.
        This dictionary has the pdb numbers as keys, and the uniprot numbering as value.'''

        # in order to maintain the uniprot numbering scheme we need an alignment of the pdb residues
        # to the uniprot sequence. We get it from SIFTS, the same system that helped us find pdbs in
        # the first place.
        #        http://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/:accession
        # sifts also has an xml based api. That has some of the same information.
        # but they seem in the process of phasing that out. So let's go with the new one.

        requestURL = 'http://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/{}'.format(self.pdb_id)
        r = requests.get(requestURL)
        # we check that this is an entry at in the database
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        mapping_dict = json.loads(r.text)
        # the uniprot mappings will be here:
        # mapping_dict[self.pdb_id]['UniProt']
        # it is a dictionary with each uniprot ID present in the structure, as a key
        # let's put he ones we are currently working on here:
        uniprot_segments = mapping_dict[self.pdb_id]['UniProt'][self.uniprotac]['mappings']
        # this is a list of segments from the structure that are present in the uniprot.
        # each element is a dictionary with the following keys:
        # dict_keys(['entity_id', 'end', 'chain_id', 'start', 'unp_end', 'unp_start', 'struct_asym_id'])
        # lets make a dictionary, with the pdb numbering as keys, and the uniprot numbering as values
        # because of chains, and the possibility of the same number in different chains,
        # it will be a dictionary of chain dictionaries.
        self.pdb_to_uniprot = {}
        for mapping in uniprot_segments:
            print('new mapping')
            # make a chain dict in the pdb_to_uniprot dict:
            self.pdb_to_uniprot[mapping['chain_id']] = {}
            print(mapping)
            # for residue number in pdb mapping:
            for residue in range(int(mapping['start']['residue_number']), int(mapping['end']['residue_number'])+1):
                # for this residue the unprot numbering is:
                uniprot_number = int(mapping['unp_start']) + residue - 1
                # assign the number in the chain dictionary to the corresponding uniprot ID
                self.pdb_to_uniprot[mapping['chain_id']][residue] = uniprot_number

        # finally we write the mapping to a file, so the ddg_parser can find it.
        self.path_to_mapping_dict = '{}/pdb_to_uniprot.txt'.format(self.path_to_run_folder)
        with open(self.path_to_mapping_dict, 'w') as mapping_file:
            json.dump(self.pdb_to_uniprot, mapping_file)

        def gently_clean_pdb():
            # here is a method for cleaning the pdb, that is less drastic than clean and isolate.
            # and we need to keep track of the numbering.
            pass
