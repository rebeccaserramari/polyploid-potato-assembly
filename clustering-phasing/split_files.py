import sys

clusterfile = sys.argv[1]
outfiles = sys.argv[2:]
print(outfiles)
#read in from file
with open(clusterfile) as cfile:
		for i,l in enumerate(cfile):
			assert(l[0] == '[')
			components = l.strip().split('singletons:')[0].strip('][').split(', ')
			singletons = l.strip().split('singletons:')[1].strip('][').split(', ')
			if i < len(outfiles):
				
				with open(outfiles[i], 'w') as rgfile:
					rgfile.write(l)
			

