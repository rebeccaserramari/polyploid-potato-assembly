import sys
from statistics import mean

def calculate_N50(list_of_lengths):
	"""Calculate N50 for a sequence of numbers.
 
	Args:
		list_of_lengths (list): List of numbers.
 
	Returns:
		float: N50 value.
 
	"""
	if len(list_of_lengths)==0:
		print("list is empty. Cannot compute N50.")
		return
	else:
		tmp = []
		for tmp_number in set(list_of_lengths):
					tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
		tmp.sort()
	 
		if (len(tmp) % 2) == 0:
			median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
		else:
			median = tmp[int(len(tmp) / 2)]
	 
		return median

def calculate_N50_alternative(contiglist, totallength):
	
	#sort contiglist in descending order
	#add all lengths from sorted list
	#when sum exceeds totallength/2, output the last added contig
	sorted_contigs = sorted(contiglist, reverse=True)
	print("longest contig: ", str(sorted_contigs[0]/1000000)+" Mb")
	lengthsum = 0
	for i,c in enumerate(sorted_contigs):
		lengthsum += c
		if lengthsum > totallength/2:
			return(c)
	return(0)	
		
def print_stats(somelist):
	print('length: ', len(somelist))
	if len(somelist) > 0:
		print('total: ', sum(somelist))
		print('mean: ', mean(somelist))
		print('max: ', max(somelist))
		print('min: ', min(somelist))
		#print('N50: ', calculate_N50(somelist))

if __name__ == '__main__':
	gfafile = sys.argv[1]
	
	lengths = []
	with open(gfafile) as gfa:
		for i,l in enumerate(gfa):			
			parts = l.strip().split('\t')
			if l[0] == 'S':
				nodename = parts[1]
				nodelength = len(parts[2])
				lengths.append(nodelength)
		print_stats(lengths)		