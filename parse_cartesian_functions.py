import numpy as np
import scipy

def rosetta_cartesian_read(pathtofile, protein_seq='KLHKEPATLIKAIDGDTVKLMYKGQPMTFRLLLVDTPETKHPKKGVEKYGPEASAFTKKMVENAKKIEVEFDKGQRTDKYGRGLAYIYADGKMVNEALVRQGLAKVAYVYKPNNTHEQHLRKSEAQAKKEKLNIWSEDNADSGQ'):
    score_file = open(pathtofile, "r")
    score_data = score_file.readlines()
    score_file.close()

    # this is just a dictionary we need to convert the three letter AAs from
    # the rosetta cartesian output into the oneletter code, from everywhere
    # else

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

    cartesian_scores = {}

    for line in score_data:
        score_fields = line.split()
        description = score_fields[2]
        three_letter_code = description[-4:-1]
        one_letter = aminocodes[three_letter_code]
        res_number = description[4:-4]
        dg = float(score_fields[3])
        key = protein_seq[int(res_number)-1]+res_number+one_letter
        if key in cartesian_scores:
            cartesian_scores[key].append(dg)
        else:
            cartesian_scores[protein_seq[int(res_number)-1]+res_number
                             + one_letter] = [dg]

    return cartesian_scores


def low_3_avg_dg_from_dg(dictionary_of_dGs):
    # this function takes a dictionary of dGs from the rosetta silent read
    # function, and turns into a dictionary with the values being the average
    # of the three lowest scores
    avg_low_three = {}
    for key in dictionary_of_dGs:
        list_of_values = []
        for item in dictionary_of_dGs[key]:
            list_of_values.append(float(item))
        list_of_values.sort()
        avg_low_three[key] = np.mean(list_of_values[0:3])
    return avg_low_three


def ddgs_from_dg(dictionary_of_dGs):
    # This function takes a dictionary of dGs, and returns a dictionary of ddGs
    # (the dictionary needs to have the "wt" mutations, ofcourse)
    # we will use the mean of the three lowest dg scores, and subtract the
    # mean of the three lowest dgs of the wts

    # this dictionary will hold all the dGs for the wt mutations. The value to
    # subtract from the corresponding actual mutations

    wt_dGs = {}
    for entry in dictionary_of_dGs:
        if entry[0] == entry[-1]:
            residue_number = entry[1:-1]
            for item in dictionary_of_dGs[entry]:
                if residue_number in wt_dGs:
                    wt_dGs[residue_number].append(float(item))
                else:
                    wt_dGs[residue_number] = [float(item)]

    # now we have all the values, let's reduce this to a single number.
    # Let's go with the average of the 50 its. A previous version used
    # the lowest three. But i don't roll like that anymore
    ddgs = {}
    dgs_as_floats = {}
    for mutation in dictionary_of_dGs:
        dgs_as_floats[mutation] = []
        for value in dictionary_of_dGs[mutation]:
            dgs_as_floats[mutation].append(float(value))

    for mutation in dictionary_of_dGs:
        residue_number = mutation[1:-1]
        ddgs[mutation] = np.mean(dgs_as_floats[mutation])-np.mean(wt_dGs[residue_number])

    return ddgs


def foldx_reader(path_to_1stn_fxout_avg_mat, AAsequence):
    aa_sequence = AAsequence
    foldx_matrix = open(path_to_1stn_fxout_avg_mat, "r")
    foldx_data = foldx_matrix.readlines()
    foldx_matrix.close()

    foldx_ddg_dict = {}

    for line in foldx_data:
        line_fields = line.split()
        if line_fields[0] != "AA" or line_fields[0][0] != '#':
            mutateto = line_fields[0]
            residue_number = 1
            for item in line_fields[1:]:
                mutatefrom = aa_sequence[residue_number-1]
                foldx_ddg_dict[mutatefrom + str(residue_number)+mutateto] = float(item)
                residue_number += 1

    return foldx_ddg_dict


def foldx_reader_wrappers_delight(path_to_foldx_ddgs, AAseq):
    '''This is a function for reading the output of saturation mutagenesis,
    performed by foldX using the wrappers delight wrapper.
    It takes an .ddgs summary file, and a string containing the
    amino acid sequence of the protein, and returns a dictionary
    with the mutations as keys, and te ddg values as values'''

    foldx_file = open(path_to_foldx_ddgs, 'r')
    foldx_data = foldx_file.readlines()
    foldx_file.close()

    # this is the order of the fields in the wrappers delight
    # ddg summary  file
    field_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                   'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    ddg_dictionary = {}

    # you should really fix this. And parse the first few lines too!

    for line in foldx_data[2:]:
        line_fields = line.split()
        # the residue index is given in the first field
        residue_index = int(line_fields[0])
        # the wild type residue is found by looking at the AAseq that has been input
        mutatefrom = AAseq[residue_index-1]

        for field_index, value in enumerate(line_fields[1:]):
            mutateto = field_order[field_index]
            ddg_dictionary[mutatefrom + str(residue_index) + mutateto] = float(value)

    return ddg_dictionary
