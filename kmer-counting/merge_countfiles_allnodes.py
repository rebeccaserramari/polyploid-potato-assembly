import sys

outfile = sys.argv[1]
samples = sys.argv[2]
path = sys.argv[3]
k = sys.argv[4]

for sample in samples:
	filename = path+"allnodes_"+sample+"_shortreadcounts_k"+k+".txt"
	with open(filename) as f:
		for i,line in enumerate(f):			
			node_name = line.strip().split('\t')[0]
			kmer_number = line.strip().split('\t')[1]
			kmer_count_bin = line.strip().split('\t')[2]
			kmer_count_tot = line.strip().split('\t')[3]
			key = node_name+"\t"+kmer_number
			nodes_to_counts_tot[key].append(kmer_count_tot)
			nodes_to_counts_bin[key].append(kmer_count_bin)
   		
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
