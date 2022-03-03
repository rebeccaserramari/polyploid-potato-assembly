import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

from statistics import mean
from matplotlib.colors import LogNorm
from collections import defaultdict


def estimate_cn(co,m):
	res = 0	
	if co >= 0 and co <= m+(0.5*m):
		res = 1
	elif co > m+(0.5*m) and co <= 2.5*m:	
		res = 2
	elif co > 2.5*m and co <= 3.5*m:
		res = 3
	elif co > 3.5*m and co <= 4.5*m:
		res = 4
	else:
		res = 5
	return(res)
		

covfile = sys.argv[1]
outfile = sys.argv[2]

covs = []
fragcovs = []
lengths = []
uniquelengths = []
fraglengths = []
node_to_cns = defaultdict(list)

node_to_length = dict()
double_interval_covs = []

node_to_avg_cn = {}
node_to_double = {}
nodecount = 0
with open(covfile) as covf:
	for i,line in enumerate(covf):
		nodecount +=1
		node = line.strip().split('\t')[0]
		avgcov = float(line.strip().split('\t')[1])
		length = int(line.strip().split('\t')[2])
		uniquelength=0
		if len(line.strip().split('\t'))> 4:
			uniquelength = int(line.strip().split('\t')[3])
			avg_cov_list = [float(i) for i in line.strip().split('\t')[4].split(',')]
		else:		
			avg_cov_list = [float(i) for i in line.strip().split('\t')[3].split(',')]
		covs.append(avgcov)
		lengths.append(length)
		uniquelengths.append(uniquelength)
		curr = 0
		double_intervals = []
		
		#collect average coverages over two interval lengths 
		while curr < len(avg_cov_list)-1:
			avg = float(avg_cov_list[curr]+avg_cov_list[curr+1])/2
			double_intervals.append(avg)
			curr += 2
		if (len(avg_cov_list)%2 == 1):
			double_intervals.append(avg_cov_list[-1])
		node_to_double[node] = double_intervals
		double_interval_covs.extend(double_intervals)

		#coverages in all intervals 
		fragcovs.extend(avg_cov_list)
		fraglengths.extend([length for i in range(len(avg_cov_list))])

assert(len([i for i in covs if i<0]) ==0)


m = mean(fragcovs)
avg_m = mean(covs)
double_interval_avg_m = mean(double_interval_covs)
all_cns = []
node_to_double_cns = {}

#estimate dosages (copy numbers) for each double interval length
for node, doub in node_to_double.items():
	cns = []
	for co in doub:
		cns.append(estimate_cn(co,double_interval_avg_m))
	node_to_double_cns[node] = cns		
	
node_to_avg_cov = {}		
with open(covfile) as covf:
	for i,line in enumerate(covf):
		nodename = line.strip().split('\t')[0]
		avgcov = float(line.strip().split('\t')[1])
		length = int(line.strip().split('\t')[2])
		node_to_length[nodename] = length
		#estimate dosage based on previously computed average
		node_to_avg_cn[nodename] = estimate_cn(avgcov, avg_m)
		node_to_avg_cov[nodename] = avgcov
		uniquelength=0
		if len(line.strip().split('\t')) > 4:
			uniquelength = int(line.strip().split('\t')[3])
			avg_cov_list = [float(i) for i in line.strip().split('\t')[4].split(',')]
		else:
			avg_cov_list = [float(i) for i in line.strip().split('\t')[3].split(',')]
		cns = []		
		cos = []
		for co in avg_cov_list:
			copynum = estimate_cn(co, m)
			cns.append(copynum)
			cos.append(round(co,2))
		node_to_cns[nodename] = cns
		all_cns.extend(cns)
		
with open(outfile,'w') as outf:
	outf.write('node,avg-coverage,avg-dosage,windowed-dosages\n')	
	for node in node_to_double.keys():
		outf.write(node+','+str(round(node_to_avg_cov[node],2))+','+str(node_to_avg_cn[node])+','+';'.join([str(i) for i in node_to_double_cns[node]]))
		outf.write('\n')

