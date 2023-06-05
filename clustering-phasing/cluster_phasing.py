import sys
import numpy as np
import pandas as pd
import copy
from collections import defaultdict


def write_clusters_to_file(clusters, outfile):
	i = 0
	with open(outfile, 'w') as outf:
		for key,nodes in clusters.items():
			outf.write(str(i)+'\t'+','.join(nodes)+'\n')
			i += 1

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
		
		index_to_cluster[startindex] =  testcluster
		clustersets.append(set(testcluster))
	print("number of clustersets: ", len(clustersets))
	print("size of clustersets: ", [len(i) for i in clustersets])
	grouped_sets = {}
	for cur_set in sorted(clustersets, key=len,reverse=True):
		is_subset = False
		for group, sets in grouped_sets.items():
			#not only common elements, but number of matching elements must be half the size of the set
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
    
        
	#define the number of haplotype clusters (homologous linkage groups)
	#define clusters where a pair of nodes has high negative correlation
	countfirst = -1
	negative_corr = {}
	positive_corr = {}
	for k, first in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
		countfirst += 1
		countsecond = -1
		for j, second in sorted(union_sets.items(), key=lambda kv: (len(kv[1]), kv[0]), reverse=True):
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


#indices_to_merge = [[8, 21, 57, 60, 68, 103, 111]]
#singlenodes = [['utg018065l', 'utg019760l', 'utg019310l', 'utg017238l', 'utg010961l', 'utg009403l', 'utg001713l', 'utg002523l']]
#rg = ['ch10']


clusterfile = sys.argv[1]
gfapath = sys.argv[2]
kmercountfile =  sys.argv[3]
outpath =  sys.argv[4]
dosagefile = sys.argv[5]
if outpath[-1] != '/':
	outpath += '/'

#read in from file
with open(clusterfile) as cfile:
	for i,l in enumerate(cfile):
		assert(l[0] == '[')
		#files should contain only one line		
		assert(i < 1)
		indices_to_merge = eval(l.split('singletons:')[0])
		singlenodes = eval(l.strip().split('singletons:')[1])
readgroup = clusterfile.split('_')[-1].split('.')[0]
cutoff = 0.1

stats_outfile = outpath+"stats_"+readgroup+'.tsv'

comp_to_phased_lengths = {}
comp_to_lengths = {}
component_lengths = {}

allnodes = []
comp_to_phased_lengths[readgroup] = 0
comp_to_lengths[readgroup] = 0
component_lengths[readgroup] = []

for i in indices_to_merge:
	if gfapath[-1] == '/':
		gfafile = gfapath+'component'+str(i)+'.gfa'
	else:
		gfapath+= '/'
		gfafile = gfapath+'component'+str(i)+'.gfa'
	with open(gfafile) as gfa:
		for line in gfa:
			if line[0] == 'S':
				node = line.strip().split('\t')[1]
				allnodes.append(node)


if len(singlenodes) > 0:
	singletons = singlenodes
else:
	singletons = []

allnodes.extend(singletons)

nodefile = outpath+"allnodes_"+str(readgroup)+".txt"
with open(nodefile, 'w') as nodef:
	for n in allnodes:
		nodef.write(n+',')

outfile = outpath + 'shortreadcounts_k71_bcalm_filtered_0.1_'+str(readgroup)+'.tsv'
nodes_in_countfile = []
with open(outfile, 'w') as outf:
	with open(kmercountfile) as kmercounts:
			for i,l in enumerate(kmercounts):
				parts = l.strip().split('\t')
				if i == 0:
					outf.write(l)                        
				if l[0] == 'u':
					nodename = parts[0]
					if nodename in allnodes:
						nodes_in_countfile.append(nodename)
						outf.write(l)

dosage_file = outfile
print(dosage_file)
df = pd.read_csv(dosage_file, sep='\t')

df = df.dropna(axis='columns')

norm = df.iloc[:,2:]
norm_t = norm.T

print("computing correlation")
corrtable = norm_t.corr(method='spearman')
print("correlation matrix computed")


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


dosage1_nodes = [node for node in nodes_in_countfile if node in haplotigs]
dosage2_nodes = [node for node in nodes_in_countfile if node in diplotigs]
dosage3_nodes = [node for node in nodes_in_countfile if node in triplotigs]
dosage4_nodes = [node for node in nodes_in_countfile if node in tetraplotigs]

haplotigdf = df.loc[df['node'].isin(dosage1_nodes)]
haplotignorm = haplotigdf.iloc[:,2:]
haplotignorm_t = haplotignorm.T

print("computing correlation")
corrtable = haplotignorm_t.corr(method='spearman')
print("correlation matrix computed")

#write stats
with open(stats_outfile,'w') as stats:
    stats.write("total length\t"+str(sum(node_to_lengths[n] for n in allnodes)))
    stats.write('\n')
    stats.write("phase-informative length\t"+str(sum(node_to_lengths[n] for n in nodes_in_countfile)))
    stats.write('\n')
    stats.write("phase-dosage-informative length\t"+str(sum(node_to_lengths[n] for n in nodes_in_countfile if n in nodes)))
    stats.write('\n')
    stats.write("dosage1 length\t"+str(sum(node_to_lengths[n] for n in dosage1_nodes)))
    stats.write('\n')
    stats.write("dosage2 length\t"+str(sum(node_to_lengths[n] for n in dosage2_nodes)))
    stats.write('\n')
    stats.write("dosage3 length\t"+str(sum(node_to_lengths[n] for n in dosage3_nodes)))
    stats.write('\n')
    stats.write("dosage4 length\t"+str(sum(node_to_lengths[n] for n in dosage4_nodes)))
    stats.write('\n')


result_clusters = compute_clustering(haplotigdf, haplotignorm_t)

print("clusters computed")
print("adding higher dosages")

new_clusters = copy.deepcopy(result_clusters)
new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, diplotigs, df, haplotigdf, haplotignorm_t)

new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, triplotigs, df, haplotigdf, haplotignorm_t)

new_clusters = add_higher_dosage_nodes(result_clusters, new_clusters, tetraplotigs, df, haplotigdf, haplotignorm_t)

print_cluster(new_clusters)
    
#write stats into file before last step of adding unphased nodes to the set
with open(stats_outfile, 'a') as stats:
    i=0
    for n,clus in new_clusters.items():
        length = sum([lengths[node] for node in clus])
        stats.write(str(i)+"_before"+"\t"+str(length/1000000)+"Mb")
        stats.write('\n')
        i+=1
"""
examine which nodes are still unphased:
1. check which nodes of allnodes are not in one of the clusters
2. for each of these nodes, compute correlation to the four clusters
"""
#allnodes = nodes_in_countfile
nodes_nocounts = [n for n in allnodes if not n in nodes_in_countfile]
dosage1_nodes = [node for node in allnodes if node in haplotigs]

print("nodes in clusters: ", sum([len(val) for key,val in new_clusters.items()]))
l = [val for key,val in new_clusters.items()]
nodes_inclusters = list(set.union(*l))

dosage1_remaining = [n for n in dosage1_nodes if n not in nodes_inclusters and n in nodes_in_countfile]
dosage1_remaining_old = [n for n in dosage1_nodes if not n in nodes_inclusters]

# ^idea: nodes without k-mer counts cannot be recovered

#create dosage1_remaining_tocorrs
df = pd.read_csv(dosage_file, sep='\t')
df = df.dropna(axis='columns')
norm = df.iloc[:,2:]
norm_t = norm.T
print("computing correlation")
corrtable = norm_t.corr(method='spearman')
print("correlation matrix computed")

dosage1_remaining_to_corrs = {}

for f in corrtable.columns:
    f_node = df.loc[[f],'node'].values[0]
    if f_node in dosage1_remaining:
        if f_node not in dosage1_remaining_to_corrs:
            dosage1_remaining_to_corrs[f_node] = []
            values = corrtable.loc[:,f].values
            nodes = [df.loc[[i],'node'].values[0] for i in corrtable.loc[:,f].index.values]
            assert(len(values) == len(nodes))
            for i in range(len(nodes)):
                pair = (nodes[i], values[i])
                dosage1_remaining_to_corrs[f_node].append(pair)

#check correlation of dosage1_remaining to the clusters

rem_to_group = {}
for node in dosage1_remaining:
    rem_to_group[node] = []
    for i, clus in new_clusters.items():
        corr = compute_corr([node], clus, dosage1_remaining_to_corrs) 
        rem_to_group[node].append((i, corr))
    el = sorted(rem_to_group[node], key=lambda t:t[1], reverse=True)[:4]   
    print("node: ", node, node_to_lengths[node], "el: ", el, "diplotig" if node in diplotigs else "triplotig" if node in triplotigs else "tetraplotig")
    
    for (n,val) in el[:1]:
        if val > 1:
            new_clusters[n].add(node)

print_cluster(new_clusters)

with open(stats_outfile, 'a') as stats:
    i=0
    for n,clus in new_clusters.items():
        length = sum([lengths[node] for node in clus])
        stats.write(str(i)+"_after"+"\t"+str(length/1000000)+"Mb")
        stats.write('\n')
        i+=1           

dosage2_remaining = [n for n in dosage2_nodes if not n in nodes_inclusters and n in nodes_in_countfile]
print("length of dosage2 nodes remaining: ", len(dosage2_remaining))
print([(n, str(lengths[n]/1000000)) for n in dosage2_remaining])
print(sum(lengths[n]/1000000 for n in dosage2_remaining))
print(sum(lengths[n]/1000000 for n in dosage1_remaining))

with open(stats_outfile, 'a') as stats:
    stats.write("dos1_remaining length\t"+str(sum(lengths[n]/1000000 for n in dosage1_remaining)))
    stats.write('\n')
    stats.write("dos2_remaining length\t"+str(sum(lengths[n]/1000000 for n in dosage2_remaining)))
    stats.write('\n')
        

cluster_outfile = outpath+"clusters_"+readgroup+'_phased.tsv'

print(cluster_outfile)
write_clusters_to_file(new_clusters, cluster_outfile)

