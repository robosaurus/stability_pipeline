from albumines import albumin
from structure_class import structure
import subprocess
import sys



def predict_stability_for_ac(uniprot_accesion, out_path):
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
                print(path_to_parse_ddgs_sbatch)
                # and submit the sbatch files:
                # first we relax:
                relax_call = subprocess.Popen('sbatch rosetta_relax.sbatch', stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)
                # and get the slurm id information
                relax_process_id_info = relax_call.communicate()
                # the actual id is here:
                relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
                # submit the ddg calculations, with the relaxation as a dependency
                cart_ddg_call = subprocess.Popen('sbatch --dependency=afterany:{} rosetta_cartesian_saturation_mutagenesis.sbatch'.format(relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)
                # again we get the slurm id
                cart_ddg_process_id_info = cart_ddg_call.communicate()
                cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]
                # aaaand submit the final piece, the parsin of the results:
                parse_results_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_ddgs.sbatch'.format(cart_ddg_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)

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
                                print(path_to_parse_ddgs_sbatch)
                                # and submit the sbatch files:
                                # first we relax:
                                relax_call = subprocess.Popen('sbatch rosetta_relax.sbatch', stdout=subprocess.PIPE, shell=True, cwd=model.path_to_run_folder)
                                # and get the slurm id information
                                relax_process_id_info = relax_call.communicate()
                                # the actual id is here:
                                relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
                                # submit the ddg calculations, with the relaxation as a dependency
                                cart_ddg_call = subprocess.Popen('sbatch --dependency=afterany:{} rosetta_cartesian_saturation_mutagenesis.sbatch'.format(relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=model.path_to_run_folder)
                                # again we get the slurm id
                                cart_ddg_process_id_info = cart_ddg_call.communicate()
                                cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]
                                # aaaand submit the final piece, the parsing of the results:
                                parse_results_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_ddgs.sbatch'.format(cart_ddg_process_id), stdout=subprocess.PIPE, shell=True, cwd=model.path_to_run_folder)


if __name__ == '__main__':
        uniprot_accesion = sys.argv[1]
        out_path = sys.argv[2]
        predict_stability_for_ac(sys.argv[1], sys.argv[2])
