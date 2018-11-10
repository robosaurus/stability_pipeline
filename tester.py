from albumines import albumin
from structure_class import structure

# make a albumin isntance
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
