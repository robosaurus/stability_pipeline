from albumines import albumin
from structure_class import structure
import sys

# this is a version for of the pipeline invokation function
# that only writes the sbatch scripts, but does not submit them


def predict_stability_no_submit(uniprot_accesion, out_path):
        '''This function runs a uniprot accesion number through the stability pipeline.
        It will ask for variants associated with the accession in clinvar and exac,
        It will make an experimental coverage map, and select suitable experimental structures.
        for each structure it will:
        write sbatch script for rosetta relaxation of the structure,
        write sbatch script detailing each possible single aa substitution in the structure,
        write an sbatch script that collects and parses the ddg values.
        The sbatch files will NOT be submitted to slurm. If you want to also submit the jobs,
        you can use the sister function in predict_stability.py'''

        # first make an instance of the albumin class
        albumine_instance = albumin(uniprot_accesion, output_path=out_path)
        # get the sequence and length
        albumine_instance.get_sequenced()
        # get the exac_variants and clinvar variants
        # but first the gene name
        albumine_instance.get_gene_name()
        albumine_instance.get_exac_variants()
        # getting the clinvar variants take a lot of time, as it relies on making thousands of http requests
        albumine_instance.get_clinvar_variants()
        # now look for experimental structures through sifts, and find a list of recomended structures,
        # for maximum coverage:
        print(albumine_instance.uniprotAC, albumine_instance.gene_name)
        structure_list = albumine_instance.pdb_map()
        print('experimental structures', structure_list)

        # and the loop through the structures in the structure list, and prepare for calculations
        for element in structure_list:
                # make an instance of the structure class
                structure_instance = structure(uniprot_accesion, output_path=albumine_instance.out_path)
                # get the pdb
                structure_instance.get_pdb(element)
                # and clean up and isolate
                # this will isolate a single chain the the pdb
                structure_instance.clean_up_and_isolate(structure_instance.path, structure_instance.chain_id)
                # now use the muscle application to make a proper alignment between the
                # the uniprot sequence and the cleaned structure sequence
                # this will allow us to translate the pdb numbering to uniprot numbering
                structure_instance.muscle_align_to_uniprot(albumine_instance.sequence)
                # now make some mutfiles.
                # for now this simply specifies all possible single AA mutations for the
                # given structure and chain.
                # it would be cool to change it so it only mutates the residues that map to the uniprot.
                structure_instance.make_mutfiles()
                # and the write the sbatch file specifying a relaxation:
                path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax()
                # and the sbatch for the actual ddg calculations:
                # in order for the clinvar variants, and the exac variants to be included,
                # we need to add them as variables to the structure instance
                path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_sbatch()
                # and finally the parse ddgs sbatch
                path_to_parse_ddgs_sbatch = structure_instance.write_parse_ddg_sbatch()

                # If there are no experimental structures, instead try get a homology model
                list_of_homology_models = []
                if structure_list == []:
                        homology_model = structure(albumine_instance.uniprotAC, output_path=albumine_instance.out_path)
                        try:
                                homology_model.get_swiss_model()
                                list_of_homology_models.append(homology_model)
                        except:
                                print('no models available')

                                # if this worked, prepare it for a rosetta run.
                                # this means doing all the stuff, we just did for the experimental structures

                        for model in list_of_homology_models:
                                # clean up and isolate
                                # this will isolate a single chain the the pdb
                                model.clean_up_and_isolate(model.path, 'ignorechain')
                                # since swissmodels don't have chain specs, we set self.chain_id to an empty string
                                model.chain_id = ''
                                # now use the muscle application to make a proper alignment between the
                                # the uniprot sequence and the cleaned structure sequence
                                # this will allow us to translate the pdb numbering to uniprot numbering
                                model.muscle_align_to_uniprot(albumine_instance.sequence)
                                # now make some mutfiles.
                                # for now this simply specifies all possible single AA mutations for the
                                # given structure and chain.
                                # it would be cool to change it so it only mutates the residues that map to the uniprot.
                                model.make_mutfiles()
                                # and the write the sbatch file specifying a relaxation:
                                path_to_relax_sbatch = model.rosetta_sbatch_relax()
                                # and the sbatch for the actual ddg calculations:
                                # in order for the clinvar variants, and the exac variants to be included,
                                # we need to add them to the structure instance
                                path_to_ddg_calc_sbatch = model.write_rosetta_cartesian_sbatch()
                                # and finally the parse ddgs sbatch
                                path_to_parse_ddgs_sbatch = model.write_parse_ddg_sbatch()


if __name__ == '__main__':
        uniprot_accesion = sys.argv[1]
        out_path = sys.argv[2]
        predict_stability_no_submit(uniprot_accesion, out_path)
