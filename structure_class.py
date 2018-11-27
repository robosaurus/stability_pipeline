#! python3
import requests
import json
import sys
import subprocess
import os
import xml.etree.ElementTree as ET
# and some functions for parsing ddg_predictions
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg

# this file defines the structure class.
# it will be used by the albumin


class structure:

    def __init__(self, UniprotAC):
        # the structure will belong to a uniprot AC
        # it might as well know which one
        self.uniprotac = UniprotAC

    # also this class should know the path to rosetta:
    path_to_rosetta = '/groups/sbinlab/software/Rosetta_2018_Oct_d557f8/source'

    # the most important method for the structureclass, is looking through the pdb for structures
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
        path_to_pdbfile = 'experimental_structures/{}.pdb'.format(self.sys_name)
        with open(path_to_pdbfile, 'w') as f:
            f.write(r.text)

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

        shell_command = 'python2 ../clean_pdb.py ../{} {}'.format(path_to_pdb, chains)
        print('here is some output from the clean_pdb.py script')
        subprocess.call(shell_command, cwd='cleaned_structures/', shell=True)
        print('end of output from clean_pdb.py')

        # this script also produces a fasta, that we can convinently use to set the sequence
        # of the structure
        path_to_cleaned_fasta = 'cleaned_structures/{}_{}.fasta'.format(self.sys_name, chains)
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

    def rosetta_relax(self):
        '''this function sruns a rosetta relaxation of the structure,
        it follows the recomended protocol for a pre-relaxation prior to
        a cartesian_ddg run'''

        # it will run in it's own folder. that we can make a rosetta_runs/self.sys_name folder
        self.path_to_run_folder = 'rosetta_runs/{}'.format(self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(self.path_to_run_folder):
            os.mkdir(self.path_to_run_folder)

        # so the paths should be relative to that
        self.path_to_cleaned_pdb = '../../cleaned_structures/{}_{}.pdb'.format(self.sys_name, self.chain_id)

        # and then we run the relax app
        # from the appropriate rosetta run folder
        # the relax protocol is taken from the rosetta docs: https://www.rosettacommons.org/docs/latest/cartesian-ddG
        # and it is dependent on the cart2.script file in the rosetta_parameters folder
        # remember to change nstructs to 20, when you are done testing
        shell_command = 'srun {}/bin/relax.linuxgccrelease -s {} -use_input_sc \
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
        '''this function writes an sbatch file, and submits it to slurm.
        it should only be called after clean, make_mutfiles, and rosetta_relax'''

        # you should change to a more reliable way of specifying the structure. Maybe a selection among the 20 in rosetta_relax
        # the memory spec is based on some brief testing. 200 M is not enough 1000 is. (based on a ~800 residue protein)
        # change nstructs to 3, when you are done testing.
        path_to_sbatch = '{}/rosetta_cartesian_saturation_mutagenesis.sbatch'.format(self.path_to_run_folder)
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=Rosetta_cartesian_ddg
#SBATCH --array=0-{}%256
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem 1000
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

    def parse_rosetta_ddgs(self, path_to_mapping, sys_name, chain_id, fasta_seq, exac_variants='', clinvar_variants=''):
        '''This function parses the result of a Rosetta ddg submission,
        It returns a dictionary with the variants as keys, and the ddgs as values.
        It only works if the sbatch job has finished.'''
        # you need to make parse rosetta ddgs into a standalone script.
        path_to_run_folder = 'rosetta_runs/{}'.format(sys_name)

        # first lets cat all the *.ddg files, into a single text
        rosetta_summary_file = '{}_{}.rosetta_cartesian.dgs'.format(sys_name, chain_id)
        # Leave out the wildtype, just in case someone forgets the -ddg:muts_only flag
        shell_command = 'cat *.ddg | grep -v WT > {}'.format(rosetta_summary_file)
        print('calling to the shell:')
        print(shell_command)
        subprocess.call(shell_command, cwd=path_to_run_folder, shell=True)

        # and then we read the pdb_to_uniprot mapping file
        with open(path_to_mapping, 'r') as mapping_file:
            pdb_to_uniprot = json.load(mapping_file)

        print(pdb_to_uniprot)

        # and then we parse the file, with the functions imported from parse_cartesian_functions.py
        rosetta_cartesian_ddgs_dict = ddgs_from_dg(rosetta_cartesian_read('{}/{}'.format(path_to_run_folder, rosetta_summary_file), fasta_seq))
        # Now we just need to print it nicely into a file
        # there has got to be a more elegant way to do this...
        # ACDEFGHIKLMNPQRSTVWY
        # first open a file to write to
        scorefile = open('prediction_files/{}_{}_ddg.txt'.format(sys_name, chain_id), 'w')
        # write the header
        scorefile.write('#Rosetta cartesian_ddg stability predictions for {}\n'.format(sys_name))
        scorefile.write('#sequence is {}\n'.format(fasta_seq))
        scorefile.write('UAC_pos\t A \t C \t D \t E \t F \t G \t H \t I \t K \t L \t M \t N \t P \t Q \t R \t S \t T \t V \t W \t Y \n')
        scorefile_line = '{}' + '\t {:.3}'*20 + '\n'
        for i in range(1, len(self.fasta_seq.strip()) + 1):
            print(i)
            print(pdb_to_uniprot[chain_id][str(i)])
            try:
                scorefile.write(scorefile_line.format(pdb_to_uniprot[chain_id][str(i)], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'A'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'C'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'D'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'E'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'F'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'G'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'H'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'I'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'K'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'L'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'M'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'N'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'P'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'Q'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'R'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'S'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'T'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'V'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'W'], rosetta_cartesian_ddgs_dict[fasta_seq[i-1]+str(i)+'Y']))
            except(KeyError):
                print('missing DATA! ', i)

        if exac_variants != '':
            scorefile.write('Exac variants for uniprot Accesion {}:\n'.format(sys_name.split()[0]))
            for key in exac_variants:
                scorefile.write(key)
                scorefile.write(': ')
                for element in exac_variants[key]:
                    scorefile.write(element)
                    scorefile.write(' ')
                scorefile.write('\n')
        if clinvar_variants != '':
            scorefile.write('clinvar variants for uniprot Accesion {}:\n'.format(sys_name.split()[0]))
            for key in clinvar_variants:
                scorefile.write(key)
                scorefile.write(': ')
                for element in clinvar_variants[key]:
                    scorefile.write(element)
                    scorefile.write(' ')
                scorefile.write('\n')

        scorefile.close()

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
python3 parse_rosetta_ddgs.py {} {} {} {}
'''.format(self.path_to_mapping_dict, self.sys_name, self.chain_id, self.fasta_seq))
        score_sbatch.close()

        return score_sbatch_path

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
