import json
import sys
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg


# this file is a script for parsing the results of a rosetta ddg, saturation mut run.
# it needs to be launched by itself, so we can submit it through slurm, and tell it to wait.


def parse_rosetta_ddgs(sys_name, chain_id, fasta_seq, uniprotAC, out_path, path_to_uniprot_index_list):
    '''This function parses the result of a Rosetta ddg submission,
        It returns a dictionary with the variants as keys, and the ddgs as values.
        It only works if the sbatch job has finished.'''
    # you need to make parse rosetta ddgs into a standalone script.
    path_to_run_folder = '{}/{}/rosetta_runs/{}'.format(out_path, uniprotAC, sys_name)

    # first lets cat all the *.ddg files, into a single text
    rosetta_summary_file = '{}_{}.rosetta_cartesian.dgs'.format(sys_name, chain_id)
    # Leave out the wildtype, just in case someone forgets the -ddg:muts_only flag
    shell_command = 'cat *.ddg | grep -v WT > {}'.format(rosetta_summary_file)
    print('calling to the shell:')
    print(shell_command)
    subprocess.call(shell_command, cwd=path_to_run_folder, shell=True)

    # drop this, and use the list instead
    # and then we read the pdb_to_uniprot mapping file
    # with open(path_to_mapping, 'r') as mapping_file:
    #     pdb_to_uniprot = json.load(mapping_file)

    # read the mapping list:
    # this list translates a structure residue index to the corresponding uniprot index
    uniprot_index_list = []
    with open(path_to_uniprot_index_list) as list_file:
        map_as_string = list_file.readlines()
        for number in map_as_string.split():
            uniprot_index_list.append(int(number.strip()))

    print(uniprot_index_list)

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
    for i in range(0, len(fasta_seq.strip())):
        print(i)
        print(uniprot_index_list[str(i)])
        try:
            scorefile.write(scorefile_line.format(uniprot_index_list[i], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'A'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'C'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'D'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'E'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'F'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'G'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'H'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'I'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'K'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'L'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'M'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'N'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'P'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'Q'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'R'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'S'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'T'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'V'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'W'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i)+'Y']))
        except(KeyError):
            print('missing DATA! on residue', i+1)

    # now parse the variant information and put it at the end.
    path_to_exac_variants = '{}/{}/exac_sAA_variants.json'.format(out_path, uniprotAC)
    with open(path_to_exac_variants, 'r') as variant_file:
        variant_dict = json.load(variant_file)

    scorefile.write('Exac variants for uniprot Accession {}:\n'.format(uniprotAC))
    for key in variant_dict:
        scorefile.write(key)
        scorefile.write(': ')
        for element in variant_dict[key]:
            scorefile.write(element)
            scorefile.write(' ')
        scorefile.write('\n')

    # and the clinvar variants
    path_to_clinvar_variants = '{}/{}/clinvar_sAA_variants.json'.format(out_path, uniprotAC)
    with open(path_to_clinvar_variants, 'r') as variant_file:
        variant_dict = json.load(variant_file)
    scorefile.write('clinvar variants for uniprot Accession {}:\n'.format(uniprotAC))
    for key in variant_dict:
        scorefile.write(key)
        scorefile.write(': ')
        for element in variant_dict[key]:
            scorefile.write(element)
            scorefile.write(' ')
        scorefile.write('\n')

    scorefile.close()


# and since we need to call this from the shell
if __name__ == '__main__':
    parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[2], fasta_seq=sys.argv[3], uniprotAC=sys.argv[4], path_to_uniprot_index_list=sys.argv[5])
