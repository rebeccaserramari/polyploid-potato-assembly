import sys

clusterfile = sys.argv[1]
outfile = sys.argv[2]

node_to_haplo = {}
with open(clusterfile) as cfile:
	for i,l in enumerate(cfile):
		for node in l.strip().split("\t")[1].split(','):
			if node not in node_to_haplo.keys():
				node_to_haplo[node] = []
			node_to_haplo[node].append(l.strip().split("\t")[0])
			
with open(outfile, 'w') as outf:
	outf.write('node,haplotypes\n')
	for node,haplolist in node_to_haplo.items():
		outf.write(node+','+';'.join(haplolist)+'\n')
			
			