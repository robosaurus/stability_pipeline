from albumines import albumin
from structure_class import structure
# let's try with this small protein
# P62945
# ok, another small protein
# O00631
# make a albumin instance
O00631 = albumin('O00631')
# get the sequence and length
O00631.get_sequenced()

print(O00631.length)
print(O00631.sequence)

# make a structure instance to hold the 'best' structure from the pdb
best_pdb = structure('O00631')
# get the structure from the pdb
best_pdb.get_pdb()
# clean it up, and get the sequence
best_pdb.clean_up_and_isolate(best_pdb.path, best_pdb.chain_id)
#make mutfiles
best_pdb.make_mutfiles()

#best_pdb.rosetta_relax()
#best_pdb.submit_rosetta_cartesian()

# best_pdb.write_sbatch_rosetta_cartesian (makes an sbatch file, and submits it)
#best_pdb.parse_rosetta_ddgs()

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
