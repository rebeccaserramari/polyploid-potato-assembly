import sys

chrom = sys.argv[1]
clusterfile = sys.argv[2]
nodefile = sys.argv[3]
outfile = sys.argv[4] 
#clusterfile = "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_clusters/clusters_"+chrom+"_new.tsv"
#nodefile = "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_kmercounts_030222/allnodes_RG"+chrom+".txt"
#outfile = "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_clusters/clusters_"+chrom+"_new_withunphased.tsv"

phasednodes = []
with open(clusterfile) as cf:
	for i,l in enumerate(cf):
		cnodes = l.strip().split('\t')[1].split(',')
		phasednodes.extend(cnodes)
print("{} phased nodes".format(len(phasednodes)))
print("{} phased nodes".format(len(set(phasednodes))))

allnodes = []
with open(nodefile) as nf:
	for l in nf:
		allnodes = l.strip().split(',')
allnodes = [n for n in allnodes if len(n)>0]		
print("{} nodes".format(len(allnodes)))				
print("{} nodes".format(len(set(allnodes))))				

unphased = [n for n in allnodes if not n in phasednodes]
print("{} unphased nodes".format(len(unphased)))
print("{} unphased nodes".format(len(set(unphased))))

with open(outfile,'w') as outf:
	with open(clusterfile) as cf:
		for l in cf:
			outf.write(l)
	outf.write("unphased\t"+','.join(unphased))
			
