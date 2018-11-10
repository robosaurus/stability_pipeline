from albumines import albumin
from structure_class import structure

# make a albumin instance
P43246 = albumin('p43246')
# get the sequence and length
P43246.get_sequenced()

# print(P43246.length)
# print(P43246.sequence)

# make a structure instance to hold the 'best' structure from the pdb
best_pdb = structure('P43246')
# get the structure from the pdb
best_pdb.get_pdb()
# clean it up, and get the sequence
best_pdb.clean_up_and_isolate(best_pdb.path, best_pdb.chain_id)
# make mutfiles
best_pdb.make_mutfiles()

best_pdb.rosetta_relax()

# best_pdb.write_sbatch_rosetta_cartesian (makes an sbatch file, and submits it)

# best_pdb.parse_ddgs (this will store all the ddgs, under self. And also write a file self.sys_name.ddgs)

# i think the next part will be running the actual calculations on the
# DEiC cluster
# this should take a structure instance as input
# relax
# wait for relaxation
# write an sbatch file
# submit sbatch file
# wait for sbatch job.
# parse the results,  make a file similar to ones for wrappers delight
# these files will later be parsed and logged in the database
