import sys

infile = sys.argv[1]
outfile = sys.argv[2]
cutoff = float(sys.argv[3])
if (len(sys.argv) > 4):
	targetnodelist = sys.argv[4].strip().split(',')
else:
	targetnodelist = []

sep, ending = "",""
if '.tsv' in infile:
	sep = '\t'
	ending = '.tsv'
if '.csv' in infile:
	sep = ','
	ending = '.csv'
	

with open(outfile+'_'+str(cutoff)+ending, 'w') as outf:
	with open(infile) as f:
		for i,l in enumerate(f):
			counts = []
			#print("line: ", i)
			if i == 0:
				outf.write(l)
				continue
			parts = l.strip().split(sep)
			node = parts[0]
			kmers = int(parts[1])
			if len(targetnodelist) > 0:			
			
				if node not in targetnodelist:
					outf.write(l)
				else:
					for c in parts[2:]:
						if int(c) < cutoff *kmers:
							counts.append(0)
						else:
							counts.append(c)
					outf.write(node+sep+str(kmers)+sep+sep.join([str(i) for i in counts])+'\n')		
			else:
				for c in parts[2:]:
					if int(c) < cutoff *kmers:
						counts.append(0)
					else:
						counts.append(c)
				outf.write(node+sep+str(kmers)+sep+sep.join([str(i) for i in counts])+'\n')	