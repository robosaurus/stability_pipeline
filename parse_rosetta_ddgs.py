import json
import sys
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg


# this file is a script for parseing the results of a rosetta ddg, satu7ration mut run.
# it needs to be launched by itself, so we can submit it through slurm, and tell it to wait.


def parse_rosetta_ddgs(path_to_mapping, sys_name, chain_id, fasta_seq, exac_variants='', clinvar_variants=''):
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
    for i in range(1, len(fasta_seq.strip()) + 1):
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


# and since we need to call this from the shell
if __name__ == '__main__':
    parse_rosetta_ddgs(path_to_mapping=sys.argv[1], sys_name=sys.argv[2], chain_id=sys.argv[3], fasta_seq=sys.argv[4], exac_variants='', clinvar_variants='')
