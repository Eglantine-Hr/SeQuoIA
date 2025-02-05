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

max_index_exploration = 2500000 	
	
	
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



# Retrieving galton parameters 
galton_params = open('default.forest_summary2.log').readline().rstrip()
print(galton_params)
galton_params = eval(galton_params) 

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


with open('default.inference.parsimony_forest.p', 'rb') as fh:
    forest = pickle.load(fh)

print("looking for error")
forest.add_isotypes(isotypemap=None, isotypemap_file="isotypemap.txt", idmap=None, idmap_file="idmap.txt", isotype_names=["IgM", "IgD", "IgG3", "IgG1", "IgG2","IgG2A","IgG2B", "IgE", "IgA"])
print("not in isotype")

#summary_table = pd.DataFrame(columns = ['Tree_Index', 'Mut_Pars', 'GcTree_ll']) 
summary_table = []
all_trees = []
topo_list = []

#to add new tree in the forest, the trees must have minimal complexity for regraft 
max_real_internal_nodes = 0 

table_length = 0

max_aberration = 2

for idx, tree in enumerate(forest, start = 0 ): 
    print(idx)
    #if len(all_trees) > max_index_exploration:
    if idx > max_index_exploration or table_length > 1000 :
    	print("stopping tree space browsing")
    	break
    #if idx < 1:
    #	continue 
    
    # Sometimes, nodes are attached to parents that are more mutated, even if tree is parsimonious. These trees are discarded
    Ancestral_Seq =  tree.tree.search_nodes(name="Ancestor")[0]
    #print("Ancestral seq")
    #print(Ancestral_Seq.sequence)
    aberrant_tree = False
    all_aberrations = 0
    for node in tree.tree.traverse():
    	if node.name == "Ancestor":
    		continue
    	#print("TREE TRAVERSE----------------")
    	curent_nMut = hamming_distance(Ancestral_Seq.sequence, node.sequence ) 
    	parent_nMut = hamming_distance(Ancestral_Seq.sequence, node.up.sequence ) 
    	
    	if curent_nMut == 1 and node.up.name != "Ancestor" :
    			#aberrant_tree = True
    			all_aberrations = all_aberrations + 1
    			new_parent = node.up.up
    			node.detach()
    			new_parent.add_child(node)
    			if all_aberrations > max_aberration:
    				aberrant_tree = True
    				break
    			
    	#print(curent_nMut -parent_nMut)
    	if curent_nMut -parent_nMut  <0 :
    		#aberrant_tree = True
    		all_aberrations = all_aberrations + 1
    		#print("ABERRATION")
    		#print(node.up.name)
    		new_parent = node.up.up
    		node.detach()
    		new_parent.add_child(node)
    		#print(node.up.name)
    		if all_aberrations > max_aberration:
    			aberrant_tree = True
    			break
    		#else:
    		#	max_aberration = all_aberrations
    	
    #
    if aberrant_tree: 
    	continue 
    print("Tree not aberrant") 
    if all_aberrations < max_aberration :
    	 max_aberration = all_aberrations
    
    print(all_aberrations)
    table_length = table_length +1
    gctree_score = tree.ll( galton_params[0] , galton_params[1] )[0]
    print(idx, gctree_score , mutability_parsimony(tree), isotype_parsimony(tree) )
    
    #summary_table.append([idx, gctree_score, mutability_parsimony(tree) , isotype_parsimony(tree), 0 ] )
    if gctree_score in topo_list and os.path.getsize("default.inference.parsimony_forest.p") > 20000: 
    	#more than 30 kb (30000) => all trees within a topology will not be processed 
    	print("large clone: all topologies will not be considered")
    	continue
    	
    already_explored = False
    #if len(all_trees) ==0:
    #	all_trees.append(tree)		
    all_trees.append(tree)
    
    #if os.path.getsize("default.inference.parsimony_forest.p") < 30000:
    #			already_explored = False
    #else:			
    #	for already_added in all_trees:
    #			TreeDist = tree.tree.compare(already_added.tree, unrooted=True)['rf']
    #			if TreeDist ==0 and len(already_added.tree) == len(already_added.tree) and len(all_trees) > 500:
    #				already_explored = True
    #				break
    
    #if len(summary_table) > 5000:
    #	for already_added in all_trees:
    #		TreeDist = tree.tree.compare(already_added.tree, unrooted=True)['rf']
    #		if TreeDist ==0 and len(already_added.tree) == len(already_added.tree):
    #			already_explored = True
    #			break
    	
    
    #if already_explored == False :
    if already_explored == False and isotype_parsimony(tree) ==0 :
                    	all_trees.append(tree)
                    	summary_table.append([idx, gctree_score, mutability_parsimony(tree) , isotype_parsimony(tree), 0 , all_aberrations ] )
                    	topo_list+= [gctree_score]
                    	tree.render(f"tree_{idx}.svg", idlabel=True) #to comment to save time? 
                    	tree.newick(f"tree_{idx}.nk")
                    	with open(f"tree_{idx}.p", "wb") as f:
                   	 	pickle.dump(tree, f)
    
    


#####################################
#Adding extra trees if few topologies available 
#####################################




summary_table = pd.DataFrame(summary_table, columns=('Tree_Index', 'GcTree_ll', 'Mut_Pars', 'Iso_Pars','Pars_Penalty', 'all_aberrations' )) 
summary_table.to_csv('treestat_custom_original.txt', header=True, index=False, sep='\t', mode='a')



new_index = idx  #last index 
all_topologies = np.unique( np.array(summary_table)[:,1] )
print(all_topologies)
print("TOPOLOGY !!!!!")
print(len(topo_list) )



max_iter = 1
cur_iter = 0 


if len(topo_list) < 50000:
    new_trees = []
    for current_topo in all_topologies:
        #limiting tree space exploration 
         
        if cur_iter > max_iter:
        	break
        cur_iter = cur_iter + 1
        
        new_index = new_index + 1 
        # Getting index of tree to prune and regraft -----------------------------------------
        #print(current_topo)
        #base_tree = summary_table.loc[summary_table['GcTree_ll'] ==  current_topo  ]
        #print(base_tree)
        #base_tree = base_tree.loc[base_tree['Iso_Pars'] ==  base_tree['Iso_Pars'].min()  ]
        #print(base_tree)
        #base_tree = base_tree.loc[base_tree['Mut_Pars'] ==  base_tree['Mut_Pars'].min()  ]
        #current_index = int(float(base_tree.iloc[0]['Tree_Index']))
       
        #print(current_index)
        #base_tree = all_trees[ current_index ]
        
        print("NEW")
        #base_tree = summary_table.loc[summary_table['GcTree_ll'] ==  summary_table['GcTree_ll'].max()  ] #all_aberrations
        base_tree = summary_table.loc[summary_table['all_aberrations'] ==  summary_table['all_aberrations'].min()  ]
        base_tree = base_tree.loc[summary_table['Iso_Pars'] ==  base_tree['Iso_Pars'].min()  ]
        base_tree = base_tree.loc[summary_table['GcTree_ll'] ==  base_tree['GcTree_ll'].max()  ]
        print(base_tree)
        current_index = int(base_tree['Tree_Index'].iloc[0]) #base_tree._get_value(0, 'Tree_Index') 
        print("current index")
        print(current_index)
        
        
        #Modyfying tree -----------------------------------------
        #current_index = 0 
        ###with open( f"tree_{current_index}.p",  "rb") as f:
        with open( f"tree_{current_index}.p",  "rb") as f:
            current_tree = pickle.load(f)
        
        #with open( f"default.inference.1.p",  "rb") as f:
        #    current_tree = pickle.load(f)
        
        #for leaf in current_tree: #terminal leaves, not internal nodes need to add .tree to get ete3 format from collapsed tree
        
        leaf_num = 0 #to keep track and modify gctree ll 
        
        for leaf in current_tree.tree.traverse():
            leaf_num = leaf_num + 1
            print("LEAF NAME")
            nameOfLeaf = leaf.name
            print(nameOfLeaf)
            if leaf.is_root() or not nameOfLeaf.startswith("seq") : #or leaf.up.is_root() # or not nameOfLeaf.startswith("seq"):
                continue
            leaf_parent = leaf.up
            print("parent") 
            print(leaf_parent.name)
            #old condition. If replacing grandparent by parent, we can unresolve nodes related to ancestor 
            #if leaf_parent.is_root():
            #    continue
            # here we discard trees with not real intermediate sequences => useless 
            #if not leaf_parent.name.startswith("seq") or not leaf_parent.name.startswith("seq"):
            #    continue
            
            
            # do not move if the sequence is virtual   
            #if not nameOfLeaf.startswith("seq"): #re.match(r"[Ss]eq", nameOfLeaf) == None 
            #    print(str(nameOfLeaf) + " is not real" )
            #    continue
            #if not nameOfLeaf.startswith("seq"): #re.match(r"[Ss]eq", nameOfLeaf) == None 
            #    print(str(nameOfLeaf) + " is not real" )
            #    continue
        
            
            # SPR
            ####with open( f"tree_{current_index}.p",  "rb") as f:
            #with open( f"tree_{current_index}.p",  "rb") as f:
            #with open( f"default.inference.1.p",  "rb") as f:
            with open( f"tree_{current_index}.p",  "rb") as f:
                current_tree = pickle.load(f) #not elegant, but need to reopen the file otherwise trees stored in leaves are overwritten
            trees_for_SPR = [ current_tree ] #current_tree.tree.copy()
            
            
            leafToMove = trees_for_SPR[0].tree.search_nodes(name=nameOfLeaf)[0]
            
            
            
           
            
            # Get parent an grand parent names -----------------------------------------
            
            #alternative version: children still in place while current node is moved => unresoving node 
            leafToMove_children = leafToMove.children
            children_names = [leafToMove_children[i].name for i in range(len(leafToMove_children)) ]
            print(children_names)
            
            
            #NEW # if very abundant node, moving node without children is not considered (unlikely to be clonal burst, should be there for some time) 
            if leafToMove.abundance < 5: 
                #trees_for_SPR = trees_for_SPR[0] 
                if leafToMove_children != [] :
                	print("another tree will be aded")
                	#
                	#
                	#
                	with open( f"tree_{current_index}.p",  "rb") as f:
                    		new_tree_2 = pickle.load(f) #not elegant, but need to reopen the file otherwise trees stored in leaves are overwritten
                    		trees_for_SPR.append(new_tree_2)
            
            
            
            
            #for new_tree in  trees_for_SPR:
            
            
            
            
            
            print("length of treee for SPR")
            print(len(trees_for_SPR))
            for k in range(len(trees_for_SPR)):
                print("k is equal to " + str(k))
                new_tree = trees_for_SPR[k]
                leafToMove = new_tree.tree.search_nodes(name=leaf.name)[0]
                leafToMove_children = leafToMove.children #actualize with right tree
                
                leaf_parent = leafToMove.up
                leaf_parent_name = leafToMove.up.name
                
                
                
                
                if leafToMove.up.is_root():
                    leaf_grandparent = leafToMove.up
                else:
                    leaf_grandparent = leafToMove.up.up
                #leaf_grandparent = leafToMove.up.up
                ##leaf_grandparent = new_tree.tree.search_nodes(name="Ancestor")[0] 
                # do not attach to virtual node
                #if  re.match(r"[Ss]eq|[Aa]ncestor", leaf_grandparent.name) == None :
                #    leaf_grandparent = leaf_grandparent.up
                #    print(str(leaf.name) + "was moved up")
                
                
                #if new_tree == trees_for_SPR[0]:
                if k == 0:
                    #leaving children in place
                    leafToMove.detach()
                    
                    # do not attach to virtual node
                    if  re.match(r"[Ss]eq|[Aa]ncestor", leaf_grandparent.name) == None :
                        leaf_grandparent = leaf_grandparent.up
                        print(str(leaf_grandparent.name) + "was moved up")
                    
                    
                    # do not attach to virtual node
                    if  re.match(r"[Ss]eq|[Aa]ncestor", leaf_parent.name) == None :
                    	leaf_parent = leaf_parent.up
               
                    
                    grandparent_name = leaf_grandparent.name    
                    leaf_grandparent.add_child(leafToMove)
                    MovedLeaf = new_tree.tree.search_nodes(name=leaf.name)[0]
                    if MovedLeaf.abundance < 1: 
                        for childn in children_names:
                            childn = new_tree.tree.search_nodes(name=childn)[0]
                            childn.detach()
                            leaf_parent.add_child(childn)
                    else:
                        for childn in children_names:
                            childn = new_tree.tree.search_nodes(name=childn)[0]
                            childn.detach()
                            MovedLeaf.add_child(childn)
                else:
                    print(nameOfLeaf)
                    print(new_index)
                    

                    print("DELETE STRATEGY") #unresolving the node so that parent and child become sisters without affecting progeny
                    leafToMove_copy = leafToMove.copy()
                    print(leafToMove_copy.get_ascii(show_internal=True))
                    print(leafToMove_copy.children[0].name)
                    while len(leafToMove_copy.children) > 0:
                    	leafToMove_copy.children[0].detach()
                    	#leaf_parent.add_child(leafToMove_copy.children[0])
                    print(leafToMove_copy.get_ascii(show_internal=True))
                    leafToMove.delete()
                    leaf_parent.add_child(leafToMove_copy)
                    print(new_tree.tree.get_ascii(show_internal=True))
                    
                    #moving current node children as well 
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                 
                #removing unecessary intermediaries -----------------------------------------
                #leaf_parent = new_tree.tree.search_nodes(name=leaf_parent_name)[0]
                #if len(leaf_parent.children) < 2 and re.match(r"[Ss]eq|[Aa]ncestor", leaf_parent.name) == None :
                #    leaf_parent.delete()
                    
                #http://etetoolkit.org/docs/latest/reference/reference_tree.html
                # Adding to New tree list if does not correspond to tree that was already added-----------------------------------------
                already_explored = False
                for already_added in new_trees:
                    #TreeDist = new_tree.tree.compare(already_added.tree, unrooted=True)['rf']
                    TreeDist = already_added.tree.compare(new_tree.tree, unrooted=True)['rf']
                    #print(already_added.tree)
                    #print(new_tree.tree)
                    print("distance trees =  ")
                    print(TreeDist)
                    if TreeDist ==0 and len(new_tree.tree) == len(already_added.tree) and grandparent_name != "Ancestor" and k == 0:
                        already_explored = True
                        print("already explored")
                        break
                # 
                if already_explored == False:
                    
                    new_trees.append( new_tree )
                    #print("adding this tree")
                    #print(new_tree.tree)
                    #print(new_trees)
                    new_index = new_index + 1 
                    new_tree.newick(f"tree_{new_index}.nk")
                    new_tree.render(f"tree_{new_index}.svg", idlabel=True)
                    
                    with open(f"tree_{new_index}.p", "wb") as f:
                        pickle.dump(tree, f)
                    
                    #compute statistics for this new tree
                    #summary_table.append([new_index, current_tree.ll( galton_params[0] , galton_params[1] )[0], mutability_parsimony(current_tree) , isotype_parsimony(current_tree) ])
                    GcTree_ll = new_tree.ll( galton_params[0] , galton_params[1] )[0] - leaf_num*10**(-7) - k**10**(-6)  #to avoid having the same value as original trees
                   
                    print(new_index)
                    print(GcTree_ll)
                    print( mutability_parsimony(new_tree)  )
                    print(new_tree.tree.get_ascii(show_internal=True))
                    print(isotype_parsimony(new_tree) )
                    new_row = {'Tree_Index':new_index, 'GcTree_ll':GcTree_ll, 'Mut_Pars': mutability_parsimony(new_tree), 'Iso_Pars':isotype_parsimony(new_tree),'Pars_Penalty':-100, 'all_aberrations':0 }
                    
                    summary_table = summary_table.append(new_row, ignore_index=True)
                    print('bug2')
                    #df = pd.concat([df, pd.DataFrame([insert_row])])
                    print(summary_table)
            
            # New case: attach to root if not solvable
            #    if k == 0 and re.match(r"[Aa]ncestor", leaf_parent_name ) != None and re.match(r"[Aa]ncestor", grandparent_name ) != None :
            #        leafToMove = new_tree.tree.search_nodes(name=nameOfLeaf)[0]
            #        leafToMove.detach()
            #        tree_root = leafToMove = new_tree.tree.search_nodes(name="Ancestor")[0]
            #        tree_root.add_child(leafToMove)
            #        
            #        new_index = new_index + 1 
            #        new_tree.newick(f"tree_{new_index}.nk")
            #        new_tree.render(f"tree_{new_index}.svg", idlabel=True)
            #        
            #        with open(f"tree_{new_index}.p", "wb") as f:
            #            pickle.dump(tree, f)
                
            





#####################################
#Saving outputs 
#####################################

    
#summary_table = pd.DataFrame(summary_table, columns=('Tree_Index', 'GcTree_ll', 'Mut_Pars', 'Iso_Pars' )) 
summary_table.to_csv('treestat_custom.txt', header=True, index=False, sep='\t', mode='a')


print( "pars scores done") 



#####################################
#def hamming_distance(seq1: str, seq2: str):
#    return sum(x != y for x, y in zip(seq1, seq2))    

#def regular_parsimony(ctree):
#    return sum(hamming_distance(n.up.sequence, n.sequence) for n in ctree.tree.iter_descendants())
    

