import sys
import os

def make_rosetta_mutfiles(fasta_file):
    '''This function creates rosetta mutfiles, specifying each possible
    amino acid substitution for each position in the sequence.
    Each residue will have 1 mutfile with 20 mutation specs
     lets make a single mutfile for each residue
    The input is a FASTA file. And the output will be a folder with
    a similar name'''

    # First read the fasta file
    fasta = open(fasta_file, 'r')
    fasta_text = fasta.readlines()
    fasta.close()

    fasta_seq = ''
    for line in fasta_text[1:]:
        fasta_seq += line.strip()

    print(len(fasta_seq))
    name = fasta_file.split('.')[0]
    # this is where we will put the files
    path_to_mutfiles = 'mutfiles/{}_mutfiles/'.format(name)
    # check if the folder exists, otherwise make it
    if not os.path.isdir(path_to_mutfiles):
        os.mkdir(path_to_mutfiles)

    for residue_number in range(1, len(fasta_seq)+1):
        mutfile = open(path_to_mutfiles+'mutfile{}'.format(str(residue_number)), 'w')
        mutfile.write('total 20\n')
        # and then a line for each type of AA
        for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
            mutfile.write('1\n')
            mutfile.write(fasta_seq[residue_number-1] + ' ' + str(residue_number) + ' ' + AAtype + '\n')
        mutfile.close()


if __name__ == '__main__':
    make_rosetta_mutfiles(sys.argv[1])
