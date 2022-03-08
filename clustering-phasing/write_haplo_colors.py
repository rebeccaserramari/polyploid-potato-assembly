import sys

infile = sys.argv[1]
outfile = infile.split('.csv')[0]+"_corrected.csv"


hap_to_col = {0:"red",1:"green",2:"blue",3:"orange"}

haplos = []
with open(infile) as inf:
	for i,l in enumerate(inf):
		if i>0:
			haplos.extend(l.strip().split(',')[1].split(';'))				
haplos = sorted(list(set(haplos)))

assert(len(haplos) == 4)

#transform the existing haplos into 0,1,2,3
with open(outfile, 'w') as outf:
	with open(infile) as inf:
		for i,l in enumerate(inf):
			if i == 0:
				outf.write(l)
			if i>0:
				haplolist = l.strip().split(',')[1].split(';')
				haplolist = [str(haplos.index(i)) for i in haplolist]				
				outf.write(l.strip().split(',')[0]+','+';'.join(haplolist))
				outf.write('\n')
		
for haplo in range(4):
	colfile = outfile.split('_new_nodes_corrected.csv')[0]+"_colors_h"+str(haplo)+".csv"
	with open(colfile, 'w') as colf:
		colf.write("node,color\n")
		with open(outfile) as inf:
			for i,l in enumerate(inf):
				if i>0:
					if str(haplo) in l.strip().split(',')[1].split(';'):				
						colf.write(l.strip().split(',')[0]+','+hap_to_col[haplo])
						colf.write('\n')

