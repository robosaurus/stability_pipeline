from albumines import albumin
from structure_class import structure
import subprocess

for uniprotAC in ['P43246']: #, 'P40692', 'P35557', 'P51587', 'P60484', 'Q8NEA6']:
	# make a albumin instance
	alb_ins = albumin(uniprotAC)
	# get the sequence and length
	alb_ins.get_sequenced()
	# get the exac_variants and clinvar variants
	# but first the gene name
	alb_ins.get_gene_name()
	alb_ins.get_exac_variants()
	#alb_ins.get_clinvar_variants()
	# now look for experimental structures through sifts, and find a list of recomended structures,
	# for maximum coverage:
	print(alb_ins.uniprotAC, alb_ins.gene_name)
	structure_list = alb_ins.pdb_map()
	print(structure_list)
	# and the loop through the structures in the structure list, and prepare for calculations
	for element in structure_list:
		# make an instance of the structure class
		struc_ins = structure(uniprotAC, output_path=alb_ins.out_path)
		# get the pdb
		struc_ins.get_pdb(element)
		# and clean up and isolate
		# this will isolate a single chain the the pdb
		struc_ins.clean_up_and_isolate(struc_ins.path, struc_ins.chain_id)
		# now use the muscle application to make a proper alignment between the 
		# the uniprot sequence and the cleaned structure sequence
		# this will allow us to translate the pdb numbering to uniprot numbering
		struc_ins.muscle_align_to_uniprot(alb_ins.sequence)
		# now make some mutfiles.
		# for now this simply specifies all possible single AA mutations for the
		# given structure and chain.
		# it would be cool to change it so it only mutates the residues that map to the uniprot.
		struc_ins.make_mutfiles()
		# and the write the sbatch file specifying a relaxation:
		path_to_relax_sbatch = struc_ins.rosetta_sbatch_relax()
		# and the sbatch for the actual ddg calculations:
		# in order for the clinvar variants, and the exac variants to be included,
		# we need to add them to the structure instance
		path_to_ddg_calc_sbatch = struc_ins.write_rosetta_cartesian_sbatch()
		# and finally the parse ddgs sbatch
		path_to_parse_ddgs_sbatch = struc_ins.write_parse_ddg_sbatch()
		print(path_to_parse_ddgs_sbatch)
		## and submit the sbatch files:
		# first we relax:
		relax_call = subprocess.Popen('sbatch rosetta_relax.sbatch', stdout=subprocess.PIPE, shell=True, cwd=struc_ins.path_to_run_folder)
		# and get the slurm id information
		relax_process_id_info = relax_call.communicate()
		# the actual id is here:
		relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
		# submit the ddg calculations, with the relaxation as a dependency
		cart_ddg_call = subprocess.Popen('sbatch --dependency=afterany:{} rosetta_cartesian_saturation_mutagenesis.sbatch'.format(relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=struc_ins.path_to_run_folder)
		# again we get the slurm id
		cart_ddg_process_id_info = cart_ddg_call.communicate()
		cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]
		# aaaand submit the final piece, the parsin of the results:
		parse_results_call= subprocess.Popen('sbatch --dependency=afterany:{} parse_ddgs.sbatch'.format(cart_ddg_process_id), stdout=subprocess.PIPE, shell=True, cwd=struc_ins.path_to_run_folder)
