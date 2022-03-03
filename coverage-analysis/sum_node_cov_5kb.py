import sys
from collections import defaultdict

covfile = sys.argv[1]
outfile = sys.argv[2]
nodefile = sys.argv[3]

covdict = defaultdict(list)

relevant_nodes = set()
with open(nodefile) as nodef:
	for i,line in enumerate(nodef):
		relevant_nodes.add(line.strip())

count = 0
seen_nodes = set()
with open(covfile) as f:
	for i,line in enumerate(f):
		parts = line.strip().split('\t')
		node = parts[0]
		if (node in relevant_nodes):
			if (node not in seen_nodes):
				count += 1
				print(count, ': ', node)
			seen_nodes.add(node)
			covdict[node].append(sum([int(parts[i]) for i in range(2,len(parts))]))

avg_covdict = defaultdict(list)
total_avg = dict()
for node in covdict.keys():
	if (len(covdict[node]) > 5000):
		i = 5000
		curr = 0
		while i < len(covdict[node]):
			frag = covdict[node][curr:i]
			avg = float(sum(frag))/len(frag)
			avg_covdict[node].append(avg)
			curr = i
			i += 5000
		lastfrag = covdict[node][curr:]
		lastavg = sum(lastfrag)/len(lastfrag)
		avg_covdict[node].append(lastavg)
	else:
		avg = sum(covdict[node])/len(covdict[node])
		avg_covdict[node].append(avg)
	total_avg[node] = sum(covdict[node])/len(covdict[node])

with open(outfile, 'w') as outf:
	for node in avg_covdict.keys():
		outf.write(node+'\t'+str(total_avg[node])+'\t'+str(len(covdict[node]))+'\t'+','.join([str(i) for i in avg_covdict[node]])+'\n')
