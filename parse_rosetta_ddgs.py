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
    path_to_run_folder = '{}/{}/rosetta_runs/{}_{}'.format(out_path, uniprotAC, sys_name, chain_id)
    print('the path to run folder is')
    print(path_to_run_folder)

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
        map_as_string = list_file.readlines()[0]
        for number in map_as_string.split(','):
            uniprot_index_list.append(int(number.strip()))

    print(uniprot_index_list)

    # and then we parse the file, with the functions imported from parse_cartesian_functions.py
    rosetta_cartesian_ddgs_dict = ddgs_from_dg(rosetta_cartesian_read('{}/{}'.format(path_to_run_folder, rosetta_summary_file), fasta_seq))
    # Now we just need to print it nicely into a file
    # there has got to be a more elegant way to do this...
    # ACDEFGHIKLMNPQRSTVWY
    # first open a file to write to
    scorefile = open('{}/{}/predictions/{}_{}_ddg.txt'.format(out_path, uniprotAC, sys_name, chain_id), 'w')
    # write the header
    # print some stuff for debugging
    for key in rosetta_cartesian_ddgs_dict:
        print(key, rosetta_cartesian_ddgs_dict[key])

    # ok, let's make a uniprot numbering version of the rosetta ddgs dictionary
    uniprot_numbering_ddgs_dict = {}
    for key in rosetta_cartesian_ddgs_dict:
        position = int(key[1:-1])
        uniprot_position = uniprot_index_list[position-1]
        uniprot_numbering_ddgs_dict[key[0]+str(uniprot_position)+key[-1]] = rosetta_cartesian_ddgs_dict[key]

    scorefile.write('#Rosetta cartesian_ddg stability predictions for {}\n'.format(sys_name))
    scorefile.write('#sequence is {}\n'.format(fasta_seq))
    scorefile.write('UAC_pos\t A \t C \t D \t E \t F \t G \t H \t I \t K \t L \t M \t N \t P \t Q \t R \t S \t T \t V \t W \t Y \n')
    scorefile_line = '{}' + '\t {:.3}'*20 + '\n'
    for i in range(0, len(fasta_seq.strip())):
        print(i)
        print(uniprot_index_list[i])
        try:
            scorefile.write(scorefile_line.format(uniprot_index_list[i], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'A'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'C'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'D'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'E'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'F'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'G'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'H'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'I'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'K'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'L'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'M'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'N'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'P'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'Q'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'R'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'S'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'T'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'V'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'W'], rosetta_cartesian_ddgs_dict[fasta_seq[i]+str(i+1)+'Y']))
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
            print(element)
            print(variant_dict[key][element])
            scorefile.write(str(element)+':')
            scorefile.write(' ')
            scorefile.write(str(variant_dict[key][element]))
            scorefile.write(' ')
        if variant_dict[key]['mut'] in uniprot_numbering_ddgs_dict:
            mutation = variant_dict[key]['mut']
            scorefile.write('predicted_ddg: {}'.format(uniprot_numbering_ddgs_dict[mutation]))
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
            print(element)
            print(variant_dict[key][element])
            scorefile.write(str(element)+':')
            scorefile.write(' ')
            scorefile.write(str(variant_dict[key][element]))
            scorefile.write(' ')
        if variant_dict[key]['mutation'] in uniprot_numbering_ddgs_dict:
            mutation = variant_dict[key]['mutation']
            scorefile.write('predicted_ddg: {}'.format(uniprot_numbering_ddgs_dict[mutation]))
        scorefile.write('\n')

    scorefile.close()


# and since we need to call this from the shell
if __name__ == '__main__':
    parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[2], fasta_seq=sys.argv[3], uniprotAC=sys.argv[4], path_to_uniprot_index_list=sys.argv[5], out_path=sys.argv[6])
