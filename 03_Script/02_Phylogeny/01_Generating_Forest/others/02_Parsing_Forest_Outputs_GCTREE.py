print("starting Python")

import os
import re
import pickle
import pandas as pd
import copy
import numpy as np


os.environ['QT_QPA_PLATFORM']='offscreen'

# Checking whether inputs are in folder 

path = os.getcwd()
print(path)
#files = [f for f in os.listdir('.') if os.path.isfile(f)]
#for f in files:
#	print(f)
#if re.search('Light', path):
#	max_index_exploration  = 20000
#else:
#	max_index_exploration = 2
max_index_exploration = 700000 #71	
	
	
Clone_Size = int(re.search(r'_n(\d+)_Id', path).group(1))

# Choosing the right mutability file 

if re.search('MOUSE|Mouse|mouse|mus|Mus|MUS|Murine|murine', path):
	if re.search('Light', path):
		selected_model = "MKS5F"
	else:
		selected_model = "MH_RS5A"
else:
	if re.search('Light', path):
		selected_model = "HLS5F"
	else:
		selected_model = "HHS5F"
#elif re.search('Light', path):
#	selected_model = "HLS5F"
#else:
#	selected_model = "HHS5F"

mutability_file = 	"../../../../../../01_Reference/Substitution_Models/Mutability__" + selected_model + ".csv"
substitution_file = 	"../../../../../../01_Reference/Substitution_Models/Substitution__" + selected_model + ".csv"
print(mutability_file)

#TO REMOVE AFTER
#mutability_file = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/01_Reference/Substitution_Models/Mutability__S5F.csv"
#substitution_file = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/01_Reference/Substitution_Models/Substitution__S5F.csv"

#mutability_file = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/01_Reference/Substitution_Models/Mutability__HHS5F.csv"
#substitution_file = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/01_Reference/Substitution_Models/Substitution__HHS5F.csv"


# Retrieving galton parameters 
galton_params = open('default.forest_summary2.log').readline().rstrip()
print(galton_params)
galton_params = eval(galton_params) 
# https://stackoverflow.com/questions/8494514/converting-string-to-tuple
print( "GALSTON = ", galton_params)




# definining mutability parsimony functions 
import gctree.mutation_model as mm
mut_model = mm.MutationModel(mutability_file=mutability_file, substitution_file=substitution_file)
mutability_distance = mm._mutability_distance(mut_model) #, splits=splits

def mutability_parsimony(ctree):
    return sum(mutability_distance(n.up.sequence, n.sequence) for n in ctree.tree.iter_descendants())




# definining mutability parsimony functions 

import gctree.isotyping as iso
with open("isotypemap.txt", "r") as fh:
	isotypemap = dict(map(lambda x: x.strip(), line.split(",")) for line in fh)


def isotype_parsimony(ctree):
    return sum(iso.isotype_distance(list(n.up.isotype.keys())[0], list(n.isotype.keys())[0]) for n in ctree.tree.iter_descendants())

def hamming_distance(seq1: str, seq2: str) -> int:
    r"""Hamming distance between two sequences of equal length.

    Args:
        seq1: sequence 1
        seq2: sequence 2
    """
    return sum(x != y for x, y in zip(seq1, seq2))


#####################################
# Parsing the forest and ascribing scores 
#####################################


#with open('default.inference.parsimony_forest.p', 'rb') as fh:
#    forest = pickle.load(fh)

print("looking for error")
#forest.add_isotypes(isotypemap=None, isotypemap_file="isotypemap.txt", idmap=None, idmap_file="idmap.txt", isotype_names=["IgM", "IgD", "IgG3", "IgG1", "IgG2","IgG2A","IgG2B", "IgE", "IgA"])
#print("not in isotype")

#summary_table = pd.DataFrame(columns = ['Tree_Index', 'Mut_Pars', 'GcTree_ll']) 
summary_table = []
all_trees = []
topo_list = []


 ##################################
# TEMPORARY: for GCTree benchmark
##################################

with open( f"tree_16.p",  "rb") as f:
	tree = pickle.load(f) 
#with open( f"tree_16.p",  "rb") as f:
#	tree = pickle.load(f) 
print(tree)
# Sometimes, nodes are attached to parents that are more mutated, even if tree is parsimonious. These trees are discarded
Ancestral_Seq =  tree.tree.search_nodes(name="Ancestor")[0]
#print("Ancestral seq")
#print(Ancestral_Seq.sequence)
aberrant_tree = False
all_aberrations = 0
print("STARTING LOOP")
for node in tree.tree.traverse():
	print(node.name)
	#if node.name == "Ancestor":
	#	continue
	
	#print("TREE TRAVERSE----------------")
	curent_nMut = hamming_distance(Ancestral_Seq.sequence, node.sequence ) 
	##parent_nMut = hamming_distance(Ancestral_Seq.sequence, node.up.sequence ) 
	
	dist_to_NCA = curent_nMut
	node.add_feature("distance_to_NCA", dist_to_NCA  )
	
	if curent_nMut == 1 and node.up.name != "Ancestor" :
			aberrant_tree = True
			#break
			
	#print(curent_nMut -parent_nMut)
	##if curent_nMut -parent_nMut  <0 :
		#aberrant_tree = True
		##all_aberrations = all_aberrations + 1
		#print("ABERRATION")
		#print(node.up.name)
		##new_parent = node.up.up
		#node.detach()
		#new_parent.add_child(node)
		#print(node.up.name)
		#if all_aberrations > max_aberration:
		#	aberrant_tree = True
		#	break
		#else:
		#	max_aberration = all_aberrations
	

	
gctree_score = tree.ll( galton_params[0] , galton_params[1] )[0]
idx = 16.2
summary_table.append([idx, gctree_score, mutability_parsimony(tree) , isotype_parsimony(tree), 0 , all_aberrations ] )

#tree.render(f"tree_{idx}.svg", idlabel=True) #to comment to save time? 
tree.newick(f"tree_{idx}.nk")
with open(f"tree_{idx}.p", "wb") as f:
	pickle.dump(tree, f)


summary_table = pd.DataFrame(summary_table, columns=('Tree_Index', 'GcTree_ll', 'Mut_Pars', 'Iso_Pars','Pars_Penalty', 'all_aberrations' )) 
summary_table.to_csv('treestat_custom.txt', header=True, index=False, sep='\t', mode='a')



colormap = tree.feature_colormap("distance_to_NCA", cmap='turbo')
tree.render(f"tree_{idx}.svg", colormap=colormap,  idlabel=True)

#####################################
#Adding extra trees if few topologies available 
#####################################
#summary_table.columns = ['Tree_Index', 'GcTree_ll', 'Mut_Pars', 'Iso_Pars']


