import sys

#hifi-potato/allnodes_merged_shortreadcounts_k71_bcalm_filtered_0.1.tsv
countfile = sys.argv[1]
gfafile = sys.argv[2]

node_to_lengths = {}
with open(gfafile) as gfa:
	for l in gfa:
		if (l[0] == 'S'):
			parts = l.strip().split('\t')
			node = parts[1]
			length = len(parts[2])
			node_to_lengths[node] = length

lengths = []
nodes_withkmers = []
#count first and second column of countfile
with open(countfile) as kmerfile:
	for i,l in enumerate(kmerfile):
		if i>0:
			node = l.strip().split('\t')[0]			
			kmers = int(l.strip().split('\t')[1])
			if kmers > 0:
				nodes_withkmers.append(node)
				lengths.append(node_to_lengths[node])	
print("{} nodes with kmers".format(len(nodes_withkmers)))

covered_len = [node_to_lengths[node] for node in nodes_withkmers]
total_len = sum(node_to_lengths.values())
print("total sequence: ", total_len)
print("sequence covered by nodes with kmers: ", sum(covered_len))
print("fraction: ",float(sum(covered_len)/total_len))
			
			