import sys

infile = sys.argv[1]
outpath = sys.argv[2]


hap_to_col = {0:"red",1:"green",2:"blue",3:"orange"}

haplos = []
with open(infile) as inf:
	for i,l in enumerate(inf):
		if i>0:
			haplos.extend(l.strip().split(',')[1].split(';'))				
haplos = sorted(list(set(haplos)))

assert(len(haplos) == 4)
		
for haplo in range(4):
	with open(outpath+"_colors_h"+str(haplo)+".csv", 'w') as outf:
		outf.write("node,color\n")
		with open(infile) as inf:
			for i,l in enumerate(inf):
				if i>0:
					if str(haplo) in l.strip().split(',')[1].split(';'):				
						outf.write(l.strip().split(',')[0]+','+hap_to_col[haplo])
						outf.write('\n')

