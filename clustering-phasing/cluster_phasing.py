import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import copy
from collections import defaultdict


def write_clusters_to_file(clusters, outfile):
	with open(outfile, 'w') as outf:
		for key,nodes in clusters.items():
			outf.write(str(key)+'\t'+','.join(nodes)+'\n')

def print_cluster(clusters):
	for key, el in clusters.items():
		print('cluster ',key,', length ',len(el), ': ')
		print(', '.join(el))
		
def node_to_cluster(elems,nodeclusters):
	ass = []
	for el in elems:
		assigned = False
		for key, nodes in nodeclusters.items():
			if el in nodes:
				ass.append(key)
				assigned = True
		if not assigned:
			ass.append('-')
	return(ass)

def add_higher_dosage_nodes(clusters, new_clusters, nodes, dosage_df, haplo_df, haplo_norm_t):
	
	new_df = dosage_df.loc[dosage_df['node'].isin(nodes)]
	new_df = new_df.dropna(axis='columns')
	haplo_df = haplo_df.dropna(axis='columns')
	new_norm = new_df.iloc[:,2:]
	new_norm_t = new_norm.T
	node_to_correlated, node_to_negative = {}, {}
	for node in nodes:
		
		if len(new_df.loc[new_df['node']==node].index.values) == 0:
			continue
		n_index = new_df.loc[new_df['node']==node].index.values[0]
		node_to_correlated[node] = []
		node_to_negative[node] = []
		maxcorr, mincorr = 0,0
		maxclus, minclus = -1, -1
		for i,clusnodes in clusters.items():
			correlated, negative = 0,0
			for clusnode in clusnodes:
				c_index = haplo_df.loc[haplo_df['node']==clusnode].index.values[0]
				corr = new_norm_t.loc[:,n_index].corr(haplo_norm_t.loc[:, c_index], method='spearman')
				#if corr > 0.4:
				#	correlated += 1
				if (corr > 0.1):
					correlated += 1
				if corr > maxcorr:
					maxcorr = corr
					maxclus = i
				#if corr < -0.3:
				#	negative += 1
				if (corr < -0.1):
					negative += 1
				if corr < mincorr:
					mincorr = corr
					minclus = i
			node_to_correlated[node].append((i, correlated, round(float(correlated/len(clusnodes))*100, 2)))
			node_to_negative[node].append((i, negative, round(float(negative/len(clusnodes))*100,2)))
	#	print("node: ", node, " maxcorr, maxclus: ", maxcorr, maxclus, " mincorr, minclus: ", mincorr, minclus)	
	print("node_to_correlated: ", node_to_correlated)
	print("node_to_negative: ", node_to_negative)
	cluster_to_nodes = {}
	for key,el in node_to_correlated.items():
		common = []
		different = []
		for i in range(len(el)):
			entry = el[i]
			if entry[2] > 50:
				common.append(entry[0])
			neg_entry = node_to_negative[key][i]
			if neg_entry[2] > 50:
				different.append(neg_entry[0])
		#assert that no clusters have both too many negative and positive entries (should not happen)		
		assert(len(set(common).intersection(set(different))) == 0)
		print("node, common, different: ", key, common, different)
		if not (tuple(common) in cluster_to_nodes.keys() ):
			cluster_to_nodes[tuple(common)] = []
		cluster_to_nodes[tuple(common)].append(key)
		
		for com in common:
			new_clusters[com].add(key)
	print("cluster to nodes: ", cluster_to_nodes)
	return(new_clusters)
	
def compute_clustering(haplo_df, haplo_norm_t):
	print("computing correlation")
	haplo_corrtable = haplo_norm_t.corr(method='spearman')
	for i in haplo_corrtable.index:
	    haplo_corrtable.loc[i, i] = 0.0
	print("correlation matrix computed")
	
	#haplo_corrtable = haplo_corrtable.dropna(axis=0,how='all').dropna(axis='columns',how='all')
	print("number of columns: ", len(haplo_corrtable.columns))
	print("number of rows: ", len(haplo_corrtable))
		
	corr_values = defaultdict(list)
	
	for i in range(len(haplo_corrtable.columns)-1):
		#corr_values[i] = list(haplo_corrtable.iloc[i+1:,i])
		corr_values[i] = list(haplo_corrtable.iloc[:,i])
	
	index_to_cluster = defaultdict(list)
	index_to_nodes = defaultdict(list)
	
	#startindex = 2
	#if not whole column is nan:
	#possible_starts = [i for i in corr_values.keys() if len([el for el in corr_values[i] if not str(el)=='nan'])>0 and max([el for el in corr_values[i] if not str(el)=='nan'])>0.5]
	#possible_starts = [i for i in corr_values.keys() if float(max(corr_values[i]))>0.5]
	possible_starts = [i for i in corr_values.keys() if float(max(haplo_corrtable.iloc[:,i]))>0.5]
	if len(possible_starts) > 0:
		print("first possible start: ", possible_starts[0])
	print("possible starts: ", possible_starts)
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
		#for cur in maxlist:
		#	col = corr_values[cur]
		#	testcluster.extend([i for i in range(len(col)) if col[i]>0.5])
		#	testcluster = list(set(testcluster))
		
		index_to_cluster[startindex] =  testcluster
		clustersets.append(set(testcluster))
	print("number of clustersets: ", len(clustersets))
	print("size of clustersets: ", [len(i) for i in clustersets])
	grouped_sets = {}
	for cur_set in sorted(clustersets, key=len,reverse=True):
		is_subset = False
		for group, sets in grouped_sets.items():
			#not only common elements, but number of matching elements must be half the size of the set
			#if len(common) >= 0.5*len(set().union(*sets))     
			common = cur_set.intersection(set().union(*sets))
			if len(common) >= (0.9*len(cur_set)):            
			#if not cur_set.isdisjoint(set().union(*sets)):
				sets.append(cur_set)
				is_subset = True
		if not is_subset:
			grouped_sets[tuple(cur_set)] = [cur_set]
	
	print("number of grouped sets: ", len(grouped_sets.keys()))
	union_sets = {}
	all_ids = []
	for group, sets in grouped_sets.items():
		print("group: ", group, " set: ", set().union(*sets), " length: ", len(set().union(*sets)))
		union_sets[group] = tuple(set().union(*sets))
		all_ids.extend(list(set().union(*sets)))
	
	print("number of union sets: ", len(union_sets.keys()))
	print("number of all ids: ", len(all_ids))
    
	#TODO test
#	union_sets = {}
#	for cur in clustersets:
#		union_sets[tuple(cur)] = tuple(cur)
        
	#define the number of haplotype clusters (homologous linkage groups)
	#define clusters where a pair of nodes has high negative correlation
	countfirst = -1
	negative_corr = {}
	positive_corr = {}
	for k, first in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
	#for k, first in sorted(clustersets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
		countfirst += 1
		countsecond = -1
		for j, second in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
		#for j, second in sorted(clustersets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
			countsecond += 1
			negative, positive =find_correlation(first, second, haplo_corrtable) 
			if negative:
				print("hom linkage first: ", countfirst, len(first), " second: ", countsecond, len(second))
				if countfirst not in negative_corr.keys():			
					negative_corr[countfirst] = [countsecond]
				else:
					negative_corr[countfirst].append(countsecond)
			#else:
			#	negative_corr[countfirst] = []
			elif positive:
				print("positive: ", countfirst, len(first), " second: ", countsecond, len(second))	
				if countfirst not in positive_corr.keys():			
					positive_corr[countfirst] = [countsecond]
				else:
					positive_corr[countfirst].append(countsecond)
			if can_join(first, second, haplo_corrtable, haplo_df):
				print("can be joined: ", countfirst, len(first), " second: ", countsecond, len(second))
	#assign the small groupsets to the largest clusters
	merged, visited = merge_clusters(union_sets,negative_corr, positive_corr)
	#merged, visited = merge_clusters(clustersets,negative_corr, positive_corr)
	print("merged: ", merged)		
	print("used in clustering: ", visited)
	
	index = 0
	nodeclusters = {}
	for k, first in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
		nodes = set([haplo_df.iloc[s,:]['node'] for s in first])	
		nodeclusters[index] = nodes 
		index += 1
	
	merged_nodeclusters = {}
	for ind, nodes in nodeclusters.items():
		if ind in merged.keys():
			for j in merged[ind]:
				if ind in merged_nodeclusters.keys():
					merged_nodeclusters[ind] = merged_nodeclusters[ind].union(nodeclusters[j])
				else:
					merged_nodeclusters[ind] =nodeclusters[j]
		else:
			if ind not in visited:
				merged_nodeclusters[ind] = nodeclusters[ind]
			
			
	print("merged nodeclusters: ",  [(i, len(merged_nodeclusters[i])) for i in merged_nodeclusters.keys()])			
	
	largest_merged_nodeclusters = {key:val for (key, val) in sorted(merged_nodeclusters.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True)[:4] }
	
			
		
	remaining = [i for i in corr_values.keys() if not i in all_ids]
	print("unmapped nodes: ", len(remaining))
	unmappable = []
	unmappable_ind = []
	for rem in remaining:
		col = [i for i in corr_values[rem] if not str(i)=='nan']
		highest = haplo_corrtable.iloc[:,rem].nlargest(3).values
		highest_ind = haplo_corrtable.iloc[:,rem].nlargest(3).index.values #highest ind are the real indices (loc)
		highest_ind_nodes = haplo_df.loc[highest_ind,'node'].values
		assigned = False
		for key, nodes in largest_merged_nodeclusters.items():	
			if (set(highest_ind_nodes).issubset(nodes)):
			#if (any(elem != '-' for elem in clus_assigns_nodes) and all(elem == clus_assigns_nodes[0] for elem in clus_assigns_nodes)):		
				largest_merged_nodeclusters[key].add(haplo_df.iloc[rem,:]['node'])
				assigned = True
		if not assigned:
			clus_assigns_nodes = node_to_cluster(highest_ind_nodes, largest_merged_nodeclusters)
			
			cleaned = [i for i in clus_assigns_nodes if i!='-']
			if len(cleaned) >0 and all(elem==cleaned[0] for elem in cleaned):
				largest_merged_nodeclusters[cleaned[0]].add(haplo_df.iloc[rem,:]['node'])
			else:
				print("rem: ", rem, " node:", haplo_df.iloc[rem,:]['node'], "highest: ", [i for i in zip(highest, highest_ind, highest_ind_nodes,clus_assigns_nodes)])		
				unmappable.append(haplo_df.iloc[rem,:]['node'])
				unmappable_ind.append(rem)
	
	print("unmappable: ", len(unmappable))	
	print("unmappable: ", unmappable)
	print("unmappable indices: ", unmappable_ind)	
		
	print("second iteration: ")
	next_unmappable = []
	next_unmappable_ind = []
	for rem in unmappable_ind:
		col = [i for i in corr_values[rem] if not str(i)=='nan']
		highest = haplo_corrtable.iloc[:,rem].nlargest(3).values
		highest_ind = haplo_corrtable.iloc[:,rem].nlargest(3).index.values #highest ind are the real indices (loc)
		highest_ind_nodes = haplo_df.loc[highest_ind,'node'].values
		assigned = False
		for key, nodes in largest_merged_nodeclusters.items():	
			if (set(highest_ind_nodes).issubset(nodes)):
				largest_merged_nodeclusters[key].add(haplo_df.iloc[rem,:]['node'])
				assigned = True
		if not assigned:
			clus_assigns_nodes = node_to_cluster(highest_ind_nodes, largest_merged_nodeclusters)
			cleaned = [i for i in clus_assigns_nodes if i!='-']
			if len(cleaned) >0 and all(elem==cleaned[0] for elem in cleaned):
				largest_merged_nodeclusters[cleaned[0]].add(haplo_df.iloc[rem,:]['node'])
			else:
				print("rem: ", rem, " node:", haplo_df.iloc[rem,:]['node'], "highest: ", [i for i in zip(highest, highest_ind, highest_ind_nodes,clus_assigns_nodes)])
				next_unmappable.append(haplo_df.iloc[rem,:]['node'])
				next_unmappable_ind.append(rem)
	
	print("remaining unmappable: ", len(next_unmappable_ind))
	print("remaining nodes: ", next_unmappable)
	
	print("clusters: ")
	print_cluster(largest_merged_nodeclusters)  	
	print("merged: ", merged)	
	return(largest_merged_nodeclusters)				

def find_correlation(first, second, haplo_corrtable):
	negative = False
	positive = False	
	positivecount = 0
	for f in first:
		for s in second:
			corr = haplo_corrtable.iloc[f,s]
			corr_diag = haplo_corrtable.iloc[s,f]
			assert(corr == corr_diag)
			if (corr < -0.3):
				negative = True
			if (corr > 0.3):
				positivecount += 1
	if (positivecount >= 1):
		positive = True			
	return((negative, positive))

def can_join(first, second, haplo_corrtable, haplo_df):
	unmappable_ind = []
	can_join = True
	second_nodes = [haplo_df.iloc[s,:]['node'] for s in second]
	for f in first:
		highest = haplo_corrtable.iloc[:,f].nlargest(3).values
		highest_ind = haplo_corrtable.iloc[:,f].nlargest(3).index.values #highest ind are the real indices (loc)
		highest_ind_nodes = haplo_df.loc[highest_ind,'node'].values
								
		if (highest_ind_nodes.any() not in second_nodes):
			can_join = False
	return(can_join)	

def contains_contradicting(pos_corr_list,neg_corr):
	contradicting = False	
	for el in pos_corr_list:
		for j in pos_corr_list:
			
			if el in neg_corr.keys() and j in neg_corr[el]:
				contradicting = True

	return(contradicting)
		
def merge_clusters(union_sets, neg_corr, pos_corr):
	visited = set()
	merged = {}
	for x in pos_corr.keys():
		if (x not in visited):	
			for el in pos_corr[x]:
				if x in neg_corr.keys():
					if not any(j in neg_corr[x] for j in pos_corr[el]) and not contains_contradicting(pos_corr[x],neg_corr):
						print("x: ", x, "el: ", el, pos_corr[el], neg_corr[x])
						visited.add(el)
						if (x not in merged.keys()):
							merged[x] = {el}
						else:
							merged[x].add(el)	
	#	visited.add(x)
	return(merged, visited)

#indices_to_merge = [[3, 16, 22, 23, 59, 63, 85, 90, 102, 118, 175, 177, 206]]
#singlenodes =  singletons = [['utg003443l', 'utg005132l', 'utg005927l', 'utg000428l', 'utg006912l', 'utg002559l', 'utg000509l', 'utg011993l', 'utg010078l', 'utg000416l', 'utg006226l', 'utg016459l', 'utg001938l', 'utg012503l', 'utg010155l', 'utg007053l', 'utg007114l']]
#rg = ['ch09']

#indices_to_merge = [[6, 9, 13, 43, 88, 115, 117, 154, 159, 172, 173]]
#singlenodes = [['utg014961l', 'utg014782l', 'utg014154l', 'utg016845l', 'utg003691l', 'utg016896l', 'utg012891l', 'utg011195l', 'utg006552l', 'utg001880l', 'utg003558l', 'utg011990l', 'utg010789l', 'utg012732l', 'utg019238l', 'utg019189l', 'utg016599l', 'utg019428l', 'utg018234l', 'utg000444l', 'utg000019l', 'utg011708l', 'utg007707l', 'utg001908l', 'utg014822l', 'utg013608l', 'utg001290l', 'utg005633l', 'utg017116l', 'utg018035l', 'utg007528l', 'utg011207l', 'utg009803l', 'utg012089l', 'utg008367l', 'utg007249l']]
#rg = ['ch03']

#indices_to_merge = [[2, 4, 7, 17, 26, 38, 39, 45, 49, 66, 83, 87, 94, 95, 97, 98, 101, 104, 106, 107, 113, 119, 123, 127, 131, 134, 135, 137, 138, 141, 142, 147, 150, 151, 152, 153, 157, 163, 165, 189, 192, 198, 201]]
#singlenodes = [['utg015925l', 'utg001014l', 'utg009123l', 'utg006440l', 'utg000434l', 'utg001774l', 'utg005817l', 'utg001152l', 'utg012255l', 'utg016897l', 'utg008840l']]
#rg = ['ch04']

#indices_to_merge = [[15, 30, 37, 47, 100, 183]]
#singlenodes = [['utg013155l', 'utg002339l', 'utg017901l']]
#rg = ['ch08']

#indices_to_merge = [[11, 29, 40, 50, 52, 61, 105, 110, 133]]
#singlenodes = [['utg017634l', 'utg007181l', 'utg014816l', 'utg009512l', 'utg001630l', 'utg000843l', 'utg017173l', 'utg015327l', 'utg013350l', 'utg000081l', 'utg002110l', 'utg000824l', 'utg003113l']]
#rg = ['ch05']
 
#indices_to_merge = [[18, 56]]
#singlenodes = [['utg012616l', 'utg009557l', 'utg017031l', 'utg002019l', 'utg002269l']]
#rg = ['ch06']    
 
#indices_to_merge = [[5, 20, 31, 67, 72, 74, 80, 149, 162, 191]]
#singlenodes = [['utg006878l', 'utg015358l', 'utg000282l', 'utg015928l', 'utg001518l', 'utg009301l', 'utg017192l', 'utg007146l', 'utg009443l', 'utg017375l', 'utg015872l', 'utg000431l', 'utg004009l', 'utg001960l', 'utg001502l', 'utg008371l', 'utg013903l', 'utg013270l', 'utg004012l', 'utg000840l']]
#rg = ['ch07']

#indices_to_merge = [[19, 32, 33, 34, 36, 42, 44, 54, 64, 73, 77, 78, 84, 91, 93, 96, 108, 114, 132, 140, 160, 167, 169, 185, 190, 194, 200]]
#singlenodes = [['utg002937l', 'utg002474l', 'utg010747l', 'utg017657l', 'utg009236l', 'utg006399l', 'utg007565l', 'utg009931l', 'utg000120l', 'utg004470l', 'utg002233l', 'utg000439l', 'utg015404l', 'utg004502l', 'utg008745l', 'utg006261l', 'utg012359l', 'utg012495l', 'utg009852l', 'utg006362l', 'utg006270l', 'utg004916l', 'utg000455l', 'utg003285l', 'utg008652l', 'utg001527l', 'utg001257l', 'utg003010l', 'utg004331l', 'utg008055l', 'utg010169l', 'utg011792l', 'utg005249l', 'utg000060l', 'utg012604l', 'utg012109l', 'utg002177l', 'utg002494l', 'utg009022l', 'utg007914l']]
#rg = ['ch01']

indices_to_merge = [[8, 21, 57, 60, 68, 103, 111]]
singlenodes = ['utg018065l', 'utg019760l', 'utg019310l', 'utg017238l', 'utg010961l', 'utg009403l', 'utg001713l', 'utg002523l']
rg = ['ch10']

#indices_to_merge = [[1, 12, 24, 35, 41, 53, 62, 65, 70, 79, 82, 89]]
#singlenodes = [['utg009919l', 'utg000588l', 'utg009242l', 'utg002803l', 'utg003108l', 'utg002442l', 'utg001239l', 'utg001037l', 'utg007560l', 'utg003416l', 'utg003107l', 'utg010245l', 'utg001863l', 'utg003950l']]
#rg = ['ch11']

#indices_to_merge = [[8, 14, 27, 51, 55, 81, 92, 99, 126, 161, 186]]
#singlenodes = [['utg001854l', 'utg007538l', 'utg002239l', 'utg013231l', 'utg006475l', 'utg005734l', 'utg015144l', 'utg018286l', 'utg013781l', 'utg007561l', 'utg005822l', 'utg003505l', 'utg002487l', 'utg017515l', 'utg016159l', 'utg016006l', 'utg000696l', 'utg017599l', 'utg016319l', 'utg017743l', 'utg000103l', 'utg008703l', 'utg018397l', 'utg014366l', 'utg001929l', 'utg008186l']]
#rg = ['ch12']

#indices_to_merge = [[10, 31, 48, 58, 69, 71, 76, 109, 120, 139, 164, 207]]
#singlenodes = [['utg015800l', 'utg011606l', 'utg017381l', 'utg012407l', 'utg008651l', 'utg000141l', 'utg012543l', 'utg010108l', 'utg010554l', 'utg004900l', 'utg018167l', 'utg018562l', 'utg016239l', 'utg018825l', 'utg000494l', 'utg015570l', 'utg014989l', 'utg015397l', 'utg008781l', 'utg009710l', 'utg017474l', 'utg006400l', 'utg011567l', 'utg015979l', 'utg013853l', 'utg003928l', 'utg007282l', 'utg009925l', 'utg003995l', 'utg002067l', 'utg000539l', 'utg001759l', 'utg008002l', 'utg007469l', 'utg007476l', 'utg002343l', 'utg010048l', 'utg014831l', 'utg013149l', 'utg002232l', 'utg011683l', 'utg000376l', 'utg002017l', 'utg018995l', 'utg014091l', 'utg005937l', 'utg000873l', 'utg011791l']]
#rg = ['ch02']

#clusterpath = "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/k71/components_clusters_bcalm/"
gfapath = "/home/rebecca/work/hifi-potato/altusHifi_270121/altusGFA/spud_altus_hifiasm.r_utg.components/"
wholedosagefile =  "/home/rebecca/work/hifi-potato/allnodes_merged_shortreadcounts_k71_bcalm_filtered_0.1_notangle.tsv"
outpath =  "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_kmercounts_030222"
cutoff = 0.1
if outpath[-1] != '/':
	outpath += '/'

comp_to_phased_lengths = {}
comp_to_lengths = {}
component_lengths = {}
for rg_index in range(len(indices_to_merge)):
	readgroup = rg[rg_index]	
	allnodes = []
	comp_to_phased_lengths[readgroup] = 0
	comp_to_lengths[readgroup] = 0
	component_lengths[readgroup] = []
	for i in indices_to_merge[rg_index]:
		gfafile = gfapath+'component'+str(i)+'.gfa'
		with open(gfafile) as gfa:
			for line in gfa:
				if line[0] == 'S':
					node = line.strip().split('\t')[1]
					allnodes.append(node)
		print('i, len allnodes: ', i, len(allnodes))
	print('number of all nodes: ', len(allnodes))
	print('number of set allnodes: ', len(set(allnodes)))
	
	if len(singlenodes) > 0:
		singletons = singlenodes[rg_index]
	else:
		singletons = []
	print('number of singletons: ', len(singletons))
	print('number of set singletons: ', len(set(singletons)))
	allnodes += singletons
	print('number of set allnodes after adding singletons: ', len(set(allnodes)))

	nodefile = outpath+"allnodes_RG"+str(readgroup)+".txt"
	with open(nodefile, 'w') as nodef:
		for n in allnodes:
			nodef.write(n+',')
				
	
	outfile = outpath + 'shortreadcounts_k71_bcalm_filtered_0.1_RG'+str(readgroup)+'.tsv'
	nodes_in_countfile = 0
	with open(outfile, 'w') as outf:
		with open(wholedosagefile) as dosage:
				for i,l in enumerate(dosage):
					parts = l.strip().split('\t')
					if i == 0:
						outf.write(l)                        
					if l[0] == 'u':
						nodename = parts[0]
						if nodename in allnodes:
							nodes_in_countfile += 1
							outf.write(l)
print("nodes found in count file: ", nodes_in_countfile)  

dosage_file = outfile
print(dosage_file)
df = pd.read_csv(dosage_file, sep='\t')

df = df.dropna(axis='columns')

norm = df.iloc[:,2:]
norm_t = norm.T

print("computing correlation")
corrtable = norm_t.corr(method='spearman')
print("correlation matrix computed")

dosagefile = "/home/rebecca/work/hifi-potato/coverage/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth-coverage-cns.csv"
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

                
print("nodes: ", len(nodes))				
print("haplotigs: ", len(haplotigs))
print("diplotigs: ", len(diplotigs))
print("triplotigs: ", len(triplotigs))
print("tetraplotigs: ", len(tetraplotigs))
print("replotigs: ", len(replotigs))

print("all nodes: ", len(allnodes))
dosage1_nodes = [node for node in allnodes if node in haplotigs]
print("dosage 1: ", len(dosage1_nodes))
#add the nodes that are not in the depth file ( = coverage 0)
dosage1_nodes.extend([node for node in allnodes if not node in nodes])
print("dosage 1 after: ", len(dosage1_nodes))

haplotigdf = df.loc[df['node'].isin(dosage1_nodes)]
haplotignorm = haplotigdf.iloc[:,2:]
haplotignorm_t = haplotignorm.T

print("computing correlation")
corrtable = haplotignorm_t.corr(method='spearman')
print("correlation matrix computed")

result_clusters = compute_clustering(haplotigdf, haplotignorm_t)

print("clusters computed")
print("adding higher dosages")

new_clusters = copy.deepcopy(result_clusters)
new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, diplotigs, df, haplotigdf, haplotignorm_t)

new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, triplotigs, df, haplotigdf, haplotignorm_t)

new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, tetraplotigs, df, haplotigdf, haplotignorm_t)

print_cluster(new_clusters)

clusterlist = []
for key,el in new_clusters.items():
    clusterlist.append(el)
    
print("union cluster 1 and 2: ", len(clusterlist[0].intersection(clusterlist[1]))) 
print("union cluster 1 and 2: ", ','.join(clusterlist[0].intersection(clusterlist[1])))
print("union cluster 1 and 3: ", len(clusterlist[0].intersection(clusterlist[2]))) 
print("union cluster 1 and 4: ", len(clusterlist[0].intersection(clusterlist[3])))
print("union cluster 1 and 4: ", ','.join(clusterlist[0].intersection(clusterlist[3])))
print("haplotigs in 1 and 4: ", len([i for i in clusterlist[0].intersection(clusterlist[3]) if i in haplotigs]))
print("haplotigs in 1 and 4: ", ','.join([i for i in clusterlist[0].intersection(clusterlist[3]) if i in haplotigs]))
print("union cluster 2 and 3: ", len(clusterlist[1].intersection(clusterlist[2]))) 
print("union cluster 2 and 4: ", len(clusterlist[1].intersection(clusterlist[3]))) 
print("union cluster 3 and 4: ", len(clusterlist[2].intersection(clusterlist[3]))) 

cluster_outfile = '/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_clusters_030222/clusters_'+rg[0]+'.tsv'
print(cluster_outfile)
write_clusters_to_file(new_clusters, cluster_outfile)

