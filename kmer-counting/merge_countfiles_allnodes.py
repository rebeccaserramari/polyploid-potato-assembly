import os
import sys
import glob
from collections import defaultdict

outfile = sys.argv[1]
#oldfile = sys.argv[2]
path = '/home/rebecca/mounts/cluster/project/serramar/potato-hifi/altusHifi/altusGFA_270121/wholegenome_altus_k71_unique_counting/count_files/'

#allnodes_HZK_11-0343_S192_shortreadcounts_k71_bcalm.txt

#nodes_to_counts_tot = defaultdict(list)
#nodes_to_counts_bin = defaultdict(list)
#samples = []
#with open(oldfile) as old:
#	for i,line in enumerate(old):
#		if i==0:
#			samples= line.strip().split('\t')[2:]
samples = sys.argv[2]
print("number of samples: ", len(samples))


for sample in samples:
	filename = path+"allnodes_"+sample+"_shortreadcounts_k71_bcalm.txt"
	with open(filename) as f:
		sample = filename.split(path)[1].split('allnodes_')[1].split('_shortreadcounts')[0]
		print(sample)
		for i,line in enumerate(f):			
			node_name = line.strip().split('\t')[0]
			kmer_number = line.strip().split('\t')[1]
			kmer_count_bin = line.strip().split('\t')[2]
			kmer_count_tot = line.strip().split('\t')[3]
			key = node_name+"\t"+kmer_number
			nodes_to_counts_tot[key].append(kmer_count_tot)
			nodes_to_counts_bin[key].append(kmer_count_bin)

print("number of samples: ", len(samples) )
   		
with open(outfile, 'w') as outf:
	outf.write('node'+'\t')
	outf.write('kmer'+'\t')
	for sample in samples:
		outf.write(sample)
		outf.write('\t')
	outf.write('\n')
	for el in nodes_to_counts_tot:
		outf.write(el)
		outf.write('\t')
		for sample in nodes_to_counts_tot[el]:
			outf.write(sample)
			outf.write('\t')
		outf.write('\n')   		

with open(outfile.split('.tsv')[0]+'_binary.tsv', 'w') as outf:
	outf.write('node'+'\t')
	outf.write('kmer'+'\t')
	for sample in samples:
		outf.write(sample)
		outf.write('\t')
	outf.write('\n')
	for el in nodes_to_counts_bin:
		outf.write(el)
		outf.write('\t')
		for sample in nodes_to_counts_bin[el]:
			outf.write(sample)
			outf.write('\t')
		outf.write('\n')   		
