#! python3
import requests
import json
import sys
import subprocess
import os

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
    def get_pdb(self):
        '''this method looks in the pdb for the best experimental structure
        and downloads it'''

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
        self.chain_id = structure_list[0]['chain_id']
        # coverage seems like a key parameter, let's record that too
        self.coverage = structure_list[0]['coverage']

        # donwload the structure
        # for now i will just download the entire structures
        # isolating chains, and cleaning it up for calculations, are going to have to wait.
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
        # it may seem a little clumsy to call a python cript
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
            self.fasta_seq = self.fasta_seq + line

        # now here should be some code for aligning the structure sequence to the uniprot sequence
        # determine the coverage (is this enough structure?)
        # and determine a numbering

    def make_mutfiles(self):
        '''this function makes Rosetta mutfiles, specifying all the possible AA substitutions
        for the structure. The .clean_up_and_isolate() method should be run first.
        The input is a fasta sequence, and the output is a folder in mutfiles/'''

        # this is where we will put the files
        path_to_mutfiles = 'mutfiles/{}_mutfiles/'.format(self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(path_to_mutfiles):
            os.mkdir(path_to_mutfiles)

        # write the files
        # each residue will be a single file, with the 20 possible AA subs (including itself)
        # rosetta counts from 1 and not 0 i think
        for residue_number in range(1, len(self.fasta_seq)+1):
            mutfile = open(path_to_mutfiles+'mutfile{}'.format(str(residue_number)), 'w')
            mutfile.write('total 20\n')
            # and then a line for each type of AA
            for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
                mutfile.write('1\n')
                mutfile.write(self.fasta_seq[residue_number-1] + ' ' + str(residue_number) + ' ' + AAtype + '\n')
            mutfile.close()

    # now we are about ready to submit some rosetta jobs
    # i am not sure how to go about that.
    # slurm maybe?
    # or should it be run locally?
    # let us start with a local boy
    # there is the relaxation of the structure first.

    def rosetta_relax(self):
        '''this function should srun a rosetta relaxation of the structure'''

        # it will run in it's own folder. that we can make a rosetta_runs/self.sys_name folder
        self.path_to_run_folder = 'rosetta_runs/{}'.format(self.sys_name)
        # check if the folder exists, otherwise make it
        if not os.path.isdir(self.path_to_run_folder):
            os.mkdir(self.path_to_run_folder)

        # so the paths should be relative to that
        self.path_to_cleaned_pdb = '../../cleaned_structures/{}_{}.pdb'.format(self.sys_name, self.chain_id)

        # first we put the pdbs to be relaxed in a list called lst
        # because that is how the relax application wants it
        lst_file = open(self.path_to_run_folder + '/lst', 'w')
        lst_file.write(self.path_to_cleaned_pdb + '\n')
        lst_file.close()

        # and then we run the relax app
        # from the appropriate rosetta run folder
        shell_command = 'srun {}/bin/relax.linuxgccrelease -ex1 -ex2 -flip_HNQ -no_optH false -use_input_sc -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -out:suffix _bn16_calibrated -beta -score:weights beta_nov16_cart.wts -ddg::legacy false -optimize_proline -in:file:s {}'.format(self.path_to_rosetta, self.path_to_cleaned_pdb)
        print('calling to the shell:{}'.format(shell_command))
        subprocess.call(shell_command, shell=True,  cwd=self.path_to_run_folder)
