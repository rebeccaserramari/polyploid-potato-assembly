import sys
import numpy as np
import copy
import glob
import os
import pandas as pd

from collections import defaultdict


def main():
	#assign each node ID to a component 
	path = "/home/rebecca/work/hifi-potato/altusHifi_270121/altusGFA/spud_altus_hifiasm.r_utg.noseq.components/"
	node_to_comp = {}
	for filename in glob.glob(os.path.join(path, '*.gfa')):
		with open(os.path.join(os.getcwd(), filename)) as f:
			component = filename.split(path)[1].split('.gfa')[0].strip('/')
			for i,line in enumerate(f):
				if line[0] == "S":
					node = line.strip().split('\t')[1]
					node_to_comp[node] = component
	
	print("number of nodes assigned to a component: ", len(node_to_comp.keys()))
	
	#read in the graph file to get a list of all nodes
	#read in the dosage file to extract nodes with dosage 1
	
	dosagefile = "/home/rebecca/work/hifi-potato/coverage/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth-coverage-cns.csv"
	gfafile = "/home/rebecca/work/hifi-potato/altusHifi_270121/altusGFA/spud_altus_hifiasm.r_utg.noseq.gfa"
	#gfafile = sys.argv[1]
	
	allnodes = []
	with open(gfafile) as gfa:
		for i,l in enumerate(gfa):
			if l[0] == 'S':
				allnodes.append(l.strip().split('\t')[1])
				
	haplotigs, diplotigs, triplotigs, tetraplotigs, replotigs = [],[],[],[], []
	nodes = []
	with open(dosagefile) as dosage:
		for i,line in enumerate(dosage):
			if i==0:
				continue
			else:
				node = line.strip().split(',')[0]
				nodes.append(node)
				
				if int(line.strip().split(',')[2]) == 1 or int(line.strip().split(',')[2]) == 0 :
					haplotigs.append(node)
				if int(line.strip().split(',')[2]) == 2:
					diplotigs.append(node)
				if int(line.strip().split(',')[2]) == 3:
					triplotigs.append(node)
				if int(line.strip().split(',')[2]) == 4:
					tetraplotigs.append(node)
				if int(line.strip().split(',')[2]) > 4:
					replotigs.append(node)
	print("allnodes: ", len(allnodes))	
	#nodes that have not been assigned to a component should be singletons
	singletons = [node for node in allnodes if not node in node_to_comp.keys()]
	print("number of single nodes: ", len(singletons))
	assert(len(singletons)+len(node_to_comp) == len(allnodes))
	for node in singletons:
		node_to_comp[node] = 'single'+node
	print(node_to_comp['utg020110l'])				
	print("nodes from count file: ", len(nodes))				
	print("haplotigs: ", len(haplotigs))
	print("diplotigs: ", len(diplotigs))
	print("triplotigs: ", len(triplotigs))
	print("tetraplotigs: ", len(tetraplotigs))
	print("replotigs: ", len(replotigs))
	
	#add nodes with dosage 0 if those are not in dosage file
	haplotigs.extend([node for node in allnodes if not node in nodes])
	print("haplotigs: ", len(haplotigs))
	
	return();
	#compute correlation
	countfile = "/home/rebecca/work/hifi-potato/allnodes_merged_shortreadcounts_k71_bcalm_filtered_0.1_notangle.tsv"
	
	df = pd.read_csv(countfile, sep='\t')
	
	df = df.dropna(axis='columns')
	#df.iloc[:,2:] = df.iloc[:,2:].apply(lambda x: [y if y <= 1 else 1 for y in x])
	
	norm = df.iloc[:,2:]
	norm_t = norm.T
	
	print("computing correlation")
	corrtable = norm_t.corr(method='spearman')
	print("correlation matrix computed")
	
	#build seed clusters like before, but with r = 0.9 --> how many clusters?
	corr_values = defaultdict(list)
	
	for i in range(len(corrtable.columns)-1):
		corr_values[i] = list(corrtable.iloc[:,i])
	
	index_to_cluster = defaultdict(list)
	index_to_nodes = defaultdict(list)
	
	possible_starts = [i for i in corr_values.keys() if float(max(corrtable.iloc[:,i]))>0.5]
	if len(possible_starts) > 0:
		print("first possible start: ", possible_starts[0])
	#print("possible starts: ", possible_starts)
	print("number of possible start positions: ", len(possible_starts))
	print("number of nodes with no high-conf corr: ", len(corr_values)-len(possible_starts))
	clustersets = []
	for startindex in possible_starts:
		testcluster = []
		maxlist = [i for i in range(len(corr_values[startindex])) if corr_values[startindex][i]>0.5]
		firstmax_id = corr_values[startindex].index(max(corr_values[startindex]))
		firstmax = max(corr_values[startindex])
		testcluster.append(startindex)
		testcluster.extend(maxlist)
		for cur in maxlist:
			col = corr_values[cur]
			testcluster.extend([i for i in range(len(col)) if col[i]>0.5])
			testcluster = list(set(testcluster))
		
		index_to_cluster[startindex] =  testcluster
		clustersets.append(set(testcluster))
	print("number of clustersets: ", len(clustersets))
	grouped_sets = {}
	for cur_set in sorted(clustersets, key=len,reverse=True):
		is_subset = False
		for group, sets in grouped_sets.items():
			if not cur_set.isdisjoint(set().union(*sets)):
				sets.append(cur_set)
				is_subset = True
		if not is_subset:
			grouped_sets[tuple(cur_set)] = [cur_set]
	
	print("number of grouped sets: ", len(grouped_sets.keys()))
	union_sets = {}
	all_ids = []
	for group, sets in grouped_sets.items():
		#print("group: ", group, " set: ", set().union(*sets), " length: ", len(set().union(*sets)))
		union_sets[group] = tuple(set().union(*sets))
		all_ids.extend(list(set().union(*sets)))
	
	print("number of union sets: ", len(union_sets.keys()))
	print("number of all ids: ", len(all_ids))
		
		
	#TODO in sets are currently only indices, no nodes!
#	index = 0
#	nodeclusters = {}
#	for k, first in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
#		nodes = set([haplo_df.iloc[s,:]['node'] for s in first])	
#		nodeclusters[index] = nodes 
#		index += 1	
		
		
	#transform sets of nodes into sets of component IDs/singleton IDs
	comp_sets = {}
	for group, indices in union_sets.items():
		comps = []	
		for index in indices:
			node = df.loc[[index],'node'].values[0]
			comps.append(node_to_comp[node])
		comp_sets[group] = list(set(comps))
	print("len component clusters: ", len(comp_sets))

	samesets = []
	longer_sets= []
	merged_compcluster = {}
	visited = set()
	for testgroup, testcomps in comp_sets.items():
	    if testgroup not in visited:
	        merged_compcluster[testgroup] = []
	        same = []
	        same_single = []
	        cs = extract_components(testcomps)
	        same.extend(cs)
	        ss = extract_singletons(testcomps)
	        same_single.extend(ss)
	        for group, comps in comp_sets.items():       
	            if group != testgroup:
	                new_cs = extract_components(comps)
	                new_ss = extract_singletons(comps)
	                if len(new_cs) > 1:
	                    longer_sets.append(group)               
	                if len(new_cs) > 0:                    
	                    intersect_frac = float(len(set(new_cs).intersection(set(cs)))/len(new_cs))
	                    if intersect_frac >0.5:
	                        same.extend(new_cs)                        
	                        visited.add(group)
	                else:
	                    intersect_frac = float(len(set(new_ss).intersection(set(ss)))/len(new_ss))
	                    if intersect_frac >0.5:
	                        same_single.extend(new_ss)                        
	                        visited.add(group)
	                   # samesets.append(new_cs)
	        merged_compcluster[testgroup] = list(set(same))+list(set(same_single))
	  
	print("number of sets with only singletons: ", len([c for key,c in comp_sets.items() if len(extract_components(c))==0]))
	print("number of sets with more than one component: ", len(longer_sets))      
	print("length of merged compluster: ", len(merged_compcluster))
	print([(len(i),len(j)) for i,j in merged_compcluster.items()])
	print("clusters: ")
	for i,c in merged_compcluster.items():
	    if len(c) > 0:
	        print(sorted([int(co.split('component')[1]) for co in c if 'component' in co]), "singletons: ", [co.split('single')[1] for co in c if 'single' in co])
	allcomponents = [i for key,c in merged_compcluster.items() for i in c]
	print("all components: ", len(allcomponents))
	print("all components without doubles: ", len(set(allcomponents)))
	#compute pairwise intersections:
	for key,c in merged_compcluster.items():
	    for key2, c2 in merged_compcluster.items():
	        if key2 != key:
	            intersect = len(set(c).intersection(set(c2)))
	            if intersect != 0 and len(c2)!= 136 and len(c)!= 136:
	                print("intersect, len c, len c2: ", intersect, len(c), len(c2))
	singletongroups = [i for i in merged_compcluster.keys() if len(i) ==1]
	print("length of groups with only one item: ", len(singletongroups))
	print([comp_sets[i] for i in singletongroups])
	#print(len(same))
	#print(samesets) 
	
	#compute list of nodesets for each merged group
	clusters = {}
	for group,compgroups in merged_compcluster.items():
	    nodes = []
	    for g in compgroups:
	        nodes.extend(comp_to_nodes[g])
	    clusters[group] = list(set(nodes))
	
	print("number of clusters: ", len(clusters))    
	print([len(c) for key,c in clusters.items()]) 
	print([len(key) for key,c in clusters.items() if len(c)==5807])

	
if __name__ == '__main__':
	main()	

	

