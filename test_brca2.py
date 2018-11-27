from albumines import albumin
from structure_class import structure
import subprocess

# let's test the pipe one brca2
# make a albumin instance
P51587 = albumin('P51587')
# get the sequence and length
P51587.get_sequenced()
# get the exac_variants and clinvar variants
P51587.get_variants()
# make a structure instance to hold the 'best' structure from the pdb
best_pdb = structure('P51587')
# get the structure from the pdb
best_pdb.get_best_pdb()
# clean it up, and get the sequence
best_pdb.clean_up_and_isolate(best_pdb.path, best_pdb.chain_id)
# make mutfiles
best_pdb.make_mutfiles()
# relax the structure
#best_pdb.rosetta_relax()
# align the pdb to uniprot, and write the mapping_dict
best_pdb.align_pdb_to_uniprot()
# write the sbatch files.
calc_sbatch_path = best_pdb.write_rosetta_cartesian_sbatch()
print(calc_sbatch_path)
score_sbatch_path = best_pdb.write_parse_ddg_sbatch()
# and submit the sbatch files:
sbatch_call = subprocess.Popen('sbatch rosetta_cartesian_saturation_mutagenesis.sbatch', stdout=subprocess.PIPE, shell=True, cwd=best_pdb.path_to_run_folder)
# this will give us standard out from the sbatch submission
sbatch_process_ID_info = sbatch_call.communicate()
# this is where the process Id is
calc_sbatch_process_ID = str(sbatch_process_ID_info[0]).split()[3][0:-3]
print('the sbatch process id is {}'.format(calc_sbatch_process_ID))
sbatch_score_command = 'sbatch --dependency=afterany:{} {}'.format(calc_sbatch_process_ID, score_sbatch_path)
subprocess.Popen(sbatch_score_command, shell=True)
