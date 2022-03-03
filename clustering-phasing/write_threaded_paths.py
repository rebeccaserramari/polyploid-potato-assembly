#!/usr/bin/python
import sys
from statistics import mean
from compute_graphstats import calculate_N50, print_stats, calculate_N50_alternative 
from write_phasedpaths import compute_pathlength, write_sequences, basepair_length
from write_paths import find_BFS_SP

redfile = sys.argv[3]
red_nodes = []
color = ""
with open(redfile) as redf:
	for i,line in enumerate(redf):
		if i>0:
			parts = line.strip().split(',')
			#assert(parts[1] == 'red')
			red_nodes.append(parts[0])
			color = parts[1]
print("number of already phased nodes: ", len(red_nodes))			

def fill_nodes(graph, right_children, left_children):
	new_phased_nodes = []	
	for node in graph:
		#if node is phased and has only one neighbour to either direction, neighbour can also be phased
		if node in red_nodes:
			if len(right_children[node]) == 1:
				new_phased_nodes.append(right_children[node][0][0])
			if len(left_children[node]) == 1:
				new_phased_nodes.append(left_children[node][0][0])	
	return(new_phased_nodes)

def find_simple_blocks(graph, right_children, left_children, dist, nodeseqs, red_nodes):		
	colored = []	
	for node in graph:
		#if node x is phased and has only two neighbours to either direction, both unphased, and both have the same phased, unique neighbour neighbour y,
		#choose the shortest of both optional paths from x to y and fill the intermediate length with Ns
		if node in red_nodes:
			for children_type in [right_children, left_children]:
				if len(children_type[node]) == 2:
					#right_c = [n[0] for n in right_children[node]]
					right_c = children_type[node]				
					#if right_c both not in red_nodes
					first = right_c[0]
					second = right_c[1]
					if children_type==right_children:
						direc='-'
					else:
						direc='+'
					if first[0] not in red_nodes and second[0] not in red_nodes and first[0] not in colored and second[0] not in colored:
						#check whether both of the right_c nodes have the same neighbour
						#assert that the left neighbour of them is 'node'
						
						if first[1] == '-':
							children_first = left_children[first[0]]
							assert((node,direc) in right_children[first[0]])
						else:
							
							children_first = right_children[first[0]]
							#print("node, children_first, first[1]: ", node, children_first, first[1])						
							assert((node,direc) in left_children[first[0]])
						if second[1] == '-':
							children_second = left_children[second[0]]
							assert((node,direc) in right_children[second[0]])
						else:
							children_second = right_children[second[0]]
							assert((node,direc) in left_children[second[0]])
						
						if len(children_first)==1 and len(children_second)==1:
							#only in case the neighbour to the other side is their only neighbour, and phased again
													
							y = children_first[0][0]
							#print("node, right_c, children_first, children_second, y: ",node, right_c, children_first, children_second, y)
							if (y == children_second[0][0]):
								if y in red_nodes:
									#check sequences of right_c[0] and [1], take shortest, but also overlaps?
									#check unique sequences of right_c[0] and right_c[1]
									unique_0 = dist[(node,right_c[0][0])]+dist[(right_c[0][0],y)] - len(nodeseqs[node])
									unique_1 = dist[(node,right_c[1][0])]+dist[(right_c[1][0],y)] - len(nodeseqs[node])
									if unique_0 < unique_1:
										shortest = right_c[0][0]
									else:
										shortest = right_c[1][0]
									colored.append(shortest)	
					
	return(colored)


def BFS_SP(graph, start, goal,right_children, left_children):
	visited = []
	queue = [[start]]
     
	#desired node is reached
	if start == goal:
		#print("Same Node")
		return
     
#	while queue:
#		path = queue.pop(0)
#		node = path[-1]
#         
#		if node not in visited:
#			print(node)
#			neighbours = graph[node]
#			phased_neighbours = False 
#			if len(neighbours) >0:
#				for neighbour in neighbours:
#					if neighbour in phased_nodes:
#						new_path = list(path)
#						new_path.append(neighbour)
#						queue.append(new_path)
#						phased_neighbours = True
#					else:
#						if neighbour == neighbours[-1]:
#							new_path = list(path)
#							new_path.append(neighbour)
#							queue.append(new_path)
#					
#					if neighbour == goal:
#				#		print("Shortest path = ", *new_path)
#						return(new_path)
#			if not phased_neighbours:
#				print("no neighbour is part of haplotype: ", node)
#			visited.append(node)

	path = queue.pop(0)
	node = path[-1]
	direction = sys.argv[4]		
	while True:
		#print(visited)
		if node not in visited:	
			visited.append(node)
		#	print(node)
			if direction == '-':
				children = left_children[node]
			elif direction == '+':
				children = right_children[node]
			
		#	print("node: ", node)
		#	print("dir: ", direction)
		#	print("children: ", children)					
				
			#children = graph[node]
			phased_found = False
			for child, childdir in children:
		#		print("child: ", child)
				if child in phased_nodes:
					path.append(child)
					phased_found = True
					visited.append(node)
					node = child
					direction = childdir
		#			print("child was phased: ", child)
				if child == goal:
					return(path)
			
			if not phased_found:
		#		print("round 2: ",node)
				has_phasedchild = False
				for child, childdir in children:
					if childdir == '-':
						grandchildren = left_children[child]
					elif childdir == '+':
						grandchildren = right_children[child]
		#			print(grandchildren)
					#grandchildren = graph[child]
					for gc, gcdir in grandchildren:
						#if gc not in visited and gc in phased_nodes:
						if gc in phased_nodes:						
							has_phasedchild = True
							path.append(child)
							path.append(gc)
							visited.append(node)
							visited.append(child)
		#					print("visited: ", visited)
							node = gc
							direction = gcdir
							break
						if gc == goal:
							return(path)
				if not has_phasedchild:
					print("need to make a break here")					
					return(path)
					
	return

def reverse_complement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	return ''.join([complement[base] for base in seq[::-1]])


def get_sequence(node, sign, dist, nodeseqs):
	#node in forward direction	
	if sign == '+':
		sequence = nodeseqs[node][:dist]
	else:
		assert(sign == '-')
		sequence = reverse_complement(nodeseqs[node])[:dist]
	return(sequence)

def basepair_length(gfapath, nodeseqs):
	print("in basepair length")
	distances = {}
	sequences = {}
	signs = {}
	dist = 0
	with open(gfapath) as gfa:
		for i,line in enumerate(gfa):
		#	if (i%10000 == 0):
		#		print("line {}".format(i))
			parts = line.strip().split('\t')
			if parts[0] == 'L':
				if (i%1000==0):
					print("L line {}".format(i))
				if len(parts) > 6:
					dist = int(parts[6].split(':')[-1])
				else:
					dist = len(nodeseqs[parts[1]]) - int(parts[5].split('M')[0])
							
				distances[(parts[1],parts[3])] = dist
				signs[(parts[1], parts[3])] = parts[2]
		#		sequences[(parts[1],parts[3])] = get_sequence(parts[1], parts[2],dist, nodeseqs)
				#compute overlap between both 
	return(distances, sequences, signs)
	
def write_adjacency(gfafile):
	node_to_list = {}
	node_to_left = {}
	node_to_right = {}
	with open(gfafile) as gfa:
		for i,line in enumerate(gfa):
			parts = line.strip().split('\t')
			if parts[0] == 'S':
			#	if parts[1] != "utg012726l":
				if parts[1] != "utg019940l":
					node_to_list[parts[1]] = []
					node_to_left[parts[1]] = []
					node_to_right[parts[1]] = []
				
			if parts[0] == 'L':
				left = parts[1]
				right = parts[3]
			#	if left != "utg012726l" and right != "utg012726l":
				if left != "utg019940l" and right != "utg019940l":
					node_to_list[left].append(right)
					direction = parts[2]
					if direction == '+':
						node_to_right[left].append((right,parts[4]))
					elif direction == '-':
						node_to_left[left].append((right,parts[4]))
	return(node_to_list, node_to_right, node_to_left)	
	
def compute_pathlength(nodes, distances):
	dist = 0	
	for i in range(len(nodes)-1):
		if (nodes[i],nodes[i+1]) not in distances.keys():
			return(-1)
		curdist = distances[(nodes[i],nodes[i+1])]
	#	print("len: ",nodes[i], nodes[i+1], curdist)
		dist += curdist	
	return(dist)

def compute_pathsequence(nodes, sequences):
	seq = ""
	for i in range(len(nodes)-1):
		seq += sequences[(nodes[i],nodes[i+1])]	
	return(seq)


def write_sequences(gfafile):
	nodeseqs = {}
	with open(gfafile) as f:
		for i,line in enumerate(f):
		#	if (i%10000==0):
		#		print("line {}".format(i))
			if line[0] == 'S':
				name = line.strip().split('\t')[1]
				seq = line.strip().split('\t')[2]
				nodeseqs[name] = seq
	return(nodeseqs)



def write_phased(all_phased):
	haplofile = "/home/rebecca/work/hifi-potato/wholegenome_kmercounts/readgroup_clusters/clusters_ch03_colors_h0_blocks_comp9_13.csv"
	with open(haplofile, 'w') as outf:
		outf.write("node,color\n")		
		for node in all_phased:
			outf.write(node+',red\n')
	return(haplofile)

def find_path(node, red_right_children, red_left_children):
	#find path of red nodes
	path = [node]
	branching_nodes = []
	print("in find path, node: ", node)
	assert((len(red_right_children[node])==0) or (len(red_left_children[node])==0))
	if len(red_right_children[node])==0:
		direction = '-'
	else:
		direction = '+'
	while True:		
	#	print("node: ", node, " red_right_children[node]: ", red_right_children[node], " red_left_children[node]: ", red_left_children[node] )
		if direction == '-':
			children = red_left_children[node]
		else:
			assert(direction == '+')
			children = red_right_children[node]
		if len(children)>1:
			branching_nodes.append(node)	
			#TODO: test... What should actually happen when two branches are phased?
			print("path: ", path)			
			return(path,branching_nodes)
		if len(children) == 0:
			assert(node in end_nodes)
			return(path,branching_nodes)
		for child,childdir in children:
			path.append(child)
			node=child
			direction = childdir

def path_exists(graph, start, goal, grey_nodes):
	visited = []
	queue = [[start]]
     
	#desired node is reached
	if start == goal:
		#print("Same Node")
		return
     
	while queue:
		path = queue.pop(0)
		node = path[-1]
         
		if node not in visited:
				neighbours = graph[node]
             
				for neighbour in neighbours:
					if neighbour in grey_nodes:
						new_path = list(path)
						new_path.append(neighbour)
						queue.append(new_path)
                 
					if neighbour == goal:
				#		print("Shortest path = ", *new_path)
						new_path = list(path)
						new_path.append(goal)
						if len(new_path) == 2:
							assert(new_path[0] == start and new_path[1] == goal)
							return(False)
						else:
							return(new_path)
				visited.append(node)
	return(False)	

def phased_path_exists(graph, start, goal, grey_nodes):
	visited = []
	queue = [[start]]
     
	#desired node is reached
	if start == goal:
		#print("Same Node")
		return
     
	while queue:
		path = queue.pop(0)
		node = path[-1]
         
		if node not in visited:
			neighbours = graph[node]
			phased_neighbours = False 
			for neighbour in neighbours:
				if neighbour in grey_nodes:
					new_path = list(path)
					new_path.append(neighbour)
					queue.append(new_path)
					phased_neighbours = True
				
					
				if neighbour == goal:
			#		print("Shortest path = ", *new_path)
					new_path = list(path)
					new_path.append(goal)
					return(new_path)
			if not phased_neighbours:
				new_path = list(path)
				new_path.append(neighbours[0])
				queue.append(new_path)		
			visited.append(node)
	return(False)	

def compute_paths(node_to_extendedpaths, node, node_to_path, node_to_end):
	path = node_to_extendedpaths[node][0]
	start = path[0]
	end = path[-1]
#	print("node, start, end: ", node, start, end)
	#only look at grey parts to not double the start and end nodes
	if len(path) > 2:
		path = path[1:-1]
	#get path that belongs to start node
	if start in node_to_path.keys():
		startpath = node_to_path[start][::-1]
#		print("if startpath: ", startpath)
	else:
		startpath = node_to_path[node_to_end[start]]
#		print("else startpath: ", startpath)
	node_to_startpath[node] = startpath	
	#get path that belongs to end node
	if end in node_to_path.keys():
		endpath = node_to_path[end]
#		print("if endpath: ", endpath)
	else:
		endpath = node_to_path[node_to_end[end]][::-1]
#		print("else endpath: ", endpath)
	node_to_endpath[node] = endpath		
	#concatenate
	wholepath = startpath+path+endpath
	return(wholepath)

def merge_pairs(p1, p2,  p1_path, p2_path):
	#p1 is the new c_pair, p2 is ip from the list p	
	
	#for the case that both boundaries are already part of the path, the new path is fully integrated within the old one, no need to merge
	if p1[0] in p2_path and p1[1] in p2_path:
		print("equal in merge_pairs")
		return(p2, p2_path)
	if p1[0] == p2[1]:
		first = p2[0]
		second = p1[1]
		new_path = p2_path + p1_path[1:]
	elif p1[1] == p2[0]:
		first = p1[0]
		second = p2[1]
		new_path = p1_path + p2_path[1:]
	elif p1[0] == p2[0]:
		first = p1[1]
		second = p2[1]
		new_path = p1_path[::-1] + p2_path[1:]
	elif p1[1] == p2[1]:
		first = p1[0]
		second = p2[0]
		new_path = p1_path + p2_path[::-1][1:]
	return((first, second), new_path)	
			
def merge_pairs_multiple(pair,merge_with, pair_to_path, node_to_extendedpaths):
	p1 = merge_with[0]
	p2 = merge_with[1]
	p1_path = pair_to_path[p1]
	p2_path = pair_to_path[p2]
	mid_path = node_to_extendedpaths[pair[0]][0]
	print("merge multiple: ", pair, mid_path)
	if pair[0] in p1:
		if pair[0] == p1[0]:
			start = p1[1]
			p1_path = p1_path[::-1]
		#	mid_path = mid_path[::-1]
		else:
			assert(pair[0] == p1[1])
			start = p1[0]
			p1_path = p1_path
			mid_path = mid_path
		assert(pair[1] in p2)
		if pair[1] == p2[0]:
			end = p2[1]
			p2_path = p2_path
		else:
			end = p2[0]	
			p2_path = p2_path[::-1]
	elif pair[0] in p2:
		if pair[0] == p2[0]:
			start = p2[1]
			p2_path = p2_path
			mid_path = mid_path[::-1]
		else:
			assert(pair[0] == p2[1])
			start = p2[0]
			p2_path = p2_path[::-1]
			mid_path = mid_path[::-1]
		assert(pair[1] in p1)
		if pair[1] == p1[0]:
			end = p1[1]
			p1_path = p1_path[::-1]
		else:
			end = p1[0]		
			p1_path = p1_path
	merged_path = p1_path+ mid_path[1:-1]+p2_path
#	if start == merged_path[-1] and end == merged_path[0]:
#		merged_path = merged_path[::-1]		
	return((start, end), merged_path)

if __name__=='__main__':
	gfafile = sys.argv[1]
	outfile = sys.argv[2]
	allnodes_file = sys.argv[4]
	write_outputlengths, write_sequence = False, False
	if sys.argv[5] == 'write_output':
		write_outputlengths = True
	if sys.argv[6] == 'write_sequence':
		write_sequence = True
	haplotype = sys.argv[3].split('.')[0].split('_')[-1]
	chromosome = sys.argv[3].split('.')[0].split('_')[-3]
	allnodes = []	
	testnodes_all = []
	with open(allnodes_file) as nf:
		for l in nf:
			if l.strip().split('\t')[0] == 'unphased':
				#allnodes = l.strip().split(',')
				allnodes = l.strip().split('\t')[1].split(',')
			testnodes_all.extend(l.strip().split('\t')[1].split(','))	
	allnodes = [n for n in allnodes if len(n)>0]		
	
	print("create graph")
	graph, right_children, left_children = write_adjacency(gfafile)
	print("size of graph: ", len(graph))
	print("graph created")
	nodeseqs = write_sequences(gfafile)
	print("sequences written")
	full_chrom_length = sum(len(nodeseqs[n]) for n in testnodes_all)
	print("full length of chrom: ", full_chrom_length)
#	pl = []
#	with open("chr10_pathlengths.txt") as pfile:
#		for i,l in enumerate(pfile):
#			pl.append(int(l.strip()))
##	full_chrom_n50 = calculate_N50_alternative(pl, full_chrom_length)
#	full_chrom_n50 = calculate_N50_alternative(pl, 195328000)
#	print("full chrom n50: ", full_chrom_n50)		
	
	unphased_length = sum(len(nodeseqs[n]) for n in allnodes)
	print("unphased length: ", unphased_length)
	print("number of red nodes: ", len(red_nodes))
	print("number of red nodes set: ", len(set(red_nodes)))
	red_len = sum(len(nodeseqs[n]) for n in red_nodes)
	print("phased length: ", red_len)
	total_length = red_len+ 0.25*unphased_length
#	total_length = red_len+ 0.5*unphased_length
#	total_length = red_len + unphased_length
	print("total length: ", total_length)	
	
	bpdists, sequences, signs = basepair_length(gfafile, nodeseqs)
	print("after base pair length")
#	all_red_lengths = [len(nodeseqs[i]) for i in red_nodes]	
#	print_stats(all_red_lengths)
#	print("N50 all red: ", calculate_N50(all_red_lengths))	
		
	red_nodes = [i for i in red_nodes if i in graph]

	#compute initial stats
	red_lengths = [len(nodeseqs[i]) for i in red_nodes]
	print("sum(red_lengths): ", sum(red_lengths))	
	print_stats(red_lengths)
	print("N50 graph red: ", calculate_N50(red_lengths))
	n50_1 = calculate_N50_alternative(red_lengths, total_length)	
	print("N50 graph red, alternative: ", n50_1)
	
	#identify red and grey nodes
	grey_nodes = [i for i in graph if not i in red_nodes]
	print("{} red, {} grey nodes and {} nodes in total".format(len(red_nodes),len(grey_nodes),len(graph)))	

	#elongate blocks where neighbours are unambiguous	
	print("fill nodes")	
	new_phased_nodes = fill_nodes(graph, right_children, left_children)
#	new_phased_nodes.remove("utg001436l")
	original_red = red_nodes
	#update red and grey nodes
	red_nodes = red_nodes + new_phased_nodes	
#	red_nodes.remove("utg002130l")
#	red_nodes.append("utg003415l")
	grey_nodes = [i for i in graph if not i in red_nodes]
	#update colors in file
	updated_colorfile = sys.argv[3].split('.')[0]+'_blocks.csv'	
	with open(updated_colorfile,'w') as upd_c:
		upd_c.write("node,color\n")
		for node in red_nodes:
			upd_c.write(node+','+color+'\n')
	print("{} red and {} grey nodes after extension".format(len(red_nodes), len(grey_nodes)))		
	print("length of red nodes after extension: ", sum(len(nodeseqs[n]) for n in set(red_nodes)))
	#fill simple unphased bubbles in between
	bluecolored = find_simple_blocks(graph, right_children, left_children, bpdists, nodeseqs, red_nodes)			
	print("{} unphased bubbles inbetween filled".format(len(bluecolored)))
	
	#update colors in file		
	updated_colorfile_bubbles = sys.argv[3].split('.')[0]+'_blocks_bubbles.csv'	
	with open(updated_colorfile_bubbles,'w') as upd_cb:
		upd_cb.write("node,color\n")
		for node in red_nodes:
			upd_cb.write(node+','+color+'\n')
		for node in bluecolored:
			upd_cb.write(node+','+"blue\n")
		
	red_blue_nodes = red_nodes + bluecolored
	grey_nodes = [i for i in graph if not i in red_blue_nodes]
	print("{} red/blue and {} grey nodes".format(len(red_blue_nodes),len(grey_nodes)))
	print("red/blue nodes: ", len(set(red_blue_nodes)))
	
	
	#separate into grey and red graph
	red_graph, red_right_children, red_left_children = {}, {}, {}
	for node in graph:
		if node in red_blue_nodes:
			red_graph[node] = [ i for i in graph[node] if i in red_blue_nodes]
			red_right_children[node] = [i for i in right_children[node] if i[0] in red_blue_nodes]
			red_left_children[node] = [i for i in left_children[node] if i[0] in red_blue_nodes]

	grey_graph, grey_right_children, grey_left_children = {}, {}, {}
	for node in graph:
		if node not in red_blue_nodes:
			grey_graph[node] = [ i for i in graph[node] if i in grey_nodes]
			grey_right_children[node] = [ i for i in right_children[node] if i[0] in grey_nodes]
			grey_left_children[node] = [ i for i in left_children[node] if i[0] in grey_nodes]

	#identify end nodes of red blocks
	end_nodes = []	
	for node in red_graph:
	#	rc_red = [i[0] for i in red_right_children[node] if i[0] in red_nodes]	
	#	lc_red = [i[0] for i in red_left_children[node] if i[0] in red_nodes]
		rc_red = red_right_children[node]	
		lc_red = red_left_children[node]	
		if len(rc_red)==0 or len(lc_red)==0:
			end_nodes.append(node)
	print("{} end nodes".format(len(end_nodes)))			
	
	#for red nodes: start path finding from there and stop when extension not possible --> should end in an end node
	node_to_path = {}
	node_to_end = {}
	pathlengths = []
	allnodes_to_path = {}
	node_to_pathlength = {}
	for i,node in enumerate(end_nodes):
		#print(node)
		if node not in node_to_end.keys():
			
			path, branching_nodes = find_path(node, red_right_children, red_left_children)	
			pathlength = compute_pathlength(path, bpdists)	
			if pathlength == 0:
				assert(len(path)==1)
				pathlength = len(nodeseqs[node])
				#pathlengths.append(pathlength)
				#print("EMPTY: ", path)	
			pathlengths.append(pathlength)
			node_to_path[node] = path
			
		#	node_to_path[path[-1]] = path[::-1]
			for n in path:
				allnodes_to_path[n] = path
			node_to_pathlength[node] = pathlength
			node_to_end[node] = path[-1]
			
			if path[-1] not in node_to_end:	
				node_to_end[path[-1]] = node
		#print("i ", i, "node: ", node, "path: ", path, "branching: ", branching_nodes)

	print("{} problematic paths".format(len([i for i in pathlengths if i==-1])))
	pathlengths = [i for i in pathlengths if not i==-1]
	print("found paths for each end node. ")
	print("{} paths".format(len(node_to_path)))
	nodes_in_paths = set([i for n in node_to_path.keys() for i in node_to_path[n] ])
	print("number of nodes part of paths: ", sum(len(node_to_path[n]) for n in node_to_path.keys()))
	print("number of nodes part of paths: ", len(nodes_in_paths))
	print("length of nodes part of paths: ", sum(len(nodeseqs[n]) for n in nodes_in_paths))
	print("length of blue nodes: ", sum(len(nodeseqs[n]) for n in bluecolored))
	print("{} pathlengths".format(len(pathlengths)))
	print("sum pathlengths: ", sum(pathlengths))
	remaining_nodes = [i for i in red_blue_nodes if not i in nodes_in_paths]
	print("remaining nodes: ", len(remaining_nodes))
	remaining_length = sum(len(nodeseqs[i]) for i in remaining_nodes)
	print("remaining length: ", remaining_length)
	print("remaining + phased: ", remaining_length+sum(pathlengths))
	print("path stats: ")
	print_stats(pathlengths)
	print("N50: ", calculate_N50(pathlengths))
	n50_2 = calculate_N50_alternative(pathlengths, total_length)	
	print("N50 alternative: ", n50_2)	
	
	
	#TODO: look at BOTH ends of a red block here
	#try connecting each end node to another end node to connect the blocks
	node_to_extendedpaths = {}
	unique_paths, zero_paths, multi_paths = [], [], []
#	for i,node1 in enumerate(node_to_path.keys()):
	used_ends = []
	for i,node1 in enumerate(end_nodes):
		if node1 not in used_ends:
		#	print("node ",i,": ",node1)	
			node_to_extendedpaths[node1] = []
	#		for node in node_to_path:
			print("node1, node_to_end[node1]: ", node1, node_to_end[node1])
			for node in end_nodes:
				if node != node_to_end[node1]:
					p = path_exists(graph,node1,node, grey_nodes)
					if p:
						used_ends.append(p[-1])
						node_to_extendedpaths[node1].append(p)
		#	print("number of paths: ", len(node_to_extendedpaths[node1]))	
			if len(node_to_extendedpaths[node1])==1:
				unique_paths.append(node1)
			if len(node_to_extendedpaths[node1])==0:
				zero_paths.append(node1)
			if len(node_to_extendedpaths[node1])>5:
				multi_paths.append(node1)		

	print("number of node_to_extendedpaths: ", len(node_to_extendedpaths))
	for key,val in node_to_extendedpaths.items():
		if len(val) > 0:
			print(key,": ",len(val), ": ", val)
		
	print("number of nodes with one unique path: ", len(unique_paths))
	print(unique_paths)
	print("first 5 unique paths: ", [node_to_extendedpaths[n] for n in unique_paths[:5]])
	print("number of nodes without a path: ", len(zero_paths))
	zero_and_single = [n for n in zero_paths if len(right_children[n])==0 and len(left_children[n])==0]
	zero_notsingle = [n for n in zero_paths if len(right_children[n])!=0 and len(left_children[n])!=0]
	print("number of nodes without a path because they are singletons: ", len(zero_and_single))
	print("nodes without path, but no singletons: ", zero_notsingle)
	print("number of nodes with > 5 paths: ", len(multi_paths))
	print(multi_paths)
	
#	#add an exception for the case that no paths can be connected, i.e. unique_paths is empty
#	if len(unique_paths) == 0:
#		for i,node in enumerate(end_nodes):
#			if node in node_to_path:
#				node_to_extendedpaths[node] = node_to_path[node]
#			else:
#				assert(node_to_end[node] in node_to_path)
#				node_to_extendedpaths[node] = node_to_path[node]	
#		unique_paths = end_nodes	
	
	#compute stats for the grey intermediate paths
	grey_pathlengths = []	
	for node in unique_paths:
		path = node_to_extendedpaths[node][0]
		if len(path) > 2:
			path = path[1:-1]
		greypathlength = compute_pathlength(path, bpdists)	
		if greypathlength == 0:
			print("EMPTY: ", path)	
		grey_pathlengths.append(greypathlength)
		
	grey_pathlengths = [i for i in grey_pathlengths if not i==-1]		
	print("{} grey path lengths".format(len(grey_pathlengths)))
	print_stats(grey_pathlengths)
#	print("N50 :", calculate_N50(grey_pathlengths))
	
	#connect paths with one unique grey path:
	whole_lengths = []
	whole_paths = []
	node_to_startpath = {}
	node_to_endpath = {}
	endpairs_to_paths = {}
	for node in unique_paths:
		path = node_to_extendedpaths[node][0]
		start = path[0]
		end = path[-1]
		print("node, start, end: ", node, start, end)
		#only look at grey parts to not double the start and end nodes
		if len(path) > 2:
			path = path[1:-1]
		#get path that belongs to start node
		if start in node_to_path.keys():
			startpath = node_to_path[start][::-1]
			print("if startpath: ", startpath)
		else:
			startpath = node_to_path[node_to_end[start]]
			print("else startpath: ", startpath)
		node_to_startpath[node] = startpath	
		#get path that belongs to end node
		if end in node_to_path.keys():
			endpath = node_to_path[end]
			print("if endpath: ", endpath)
		else:
			endpath = node_to_path[node_to_end[end]][::-1]
			print("else endpath: ", endpath)
		node_to_endpath[node] = endpath		
		#concatenate
		wholepath = startpath+path+endpath
		print("wholepath: ", wholepath)
		whole_paths.append(wholepath)
		endpairs_to_paths[(start, end)] = wholepath
		endpairs_to_paths[(end, start)] = wholepath
		#get length of concatenated path
		wholepathlength = compute_pathlength(wholepath, bpdists)
		print("wholepathlength: ", wholepathlength)	
		whole_lengths.append(wholepathlength)

	whole_lengths = [i for i in whole_lengths if not i==-1]
	print("{} whole lengths".format(len(whole_lengths)))
	print(zero_paths)		
	print_stats(whole_lengths)
	print("N50: ", calculate_N50(whole_lengths))
	
	#check which nodes are present in the long paths
	nodes_in_wholepaths = set([i for n in whole_paths for i in n ])
	print("nodes part of full paths: ", len(nodes_in_wholepaths))
	
	#for those nodes which are not present in the paths, get their path from the set of unconnected paths
	nodes_outside_paths = set([n for n in red_blue_nodes if n not in nodes_in_wholepaths])
	print("nodes outside full paths: ", len(nodes_outside_paths))
	print("nodes outside full paths: ", nodes_outside_paths)
	path_set = []
	path_set_lengths = set()
	special_ns = []
	for n in nodes_outside_paths:
		if n in allnodes_to_path.keys():
			n_path = allnodes_to_path[n]
			path_set_lengths.add(node_to_pathlength[n_path[0]])
			if n_path not in path_set:
				path_set.append(n_path)
			print("node, path: ",n, n_path)
		else:
			special_ns.append(n)
	
	print("nodes that were not part of a path: ", special_ns)		
	print("number of path_set_lengths: ", len(path_set_lengths))		
	print("number of path_set: ", len(path_set))		
	finalpaths = whole_paths + path_set
	print("number of paths that are not singletons: ", len([p for p in finalpaths if len(p) > 1]))	
	print("number of paths that are singletons: ", len([p for p in finalpaths if len(p) == 1]))	
	finallengths = whole_lengths + list(path_set_lengths)		
	print("sum of finallengths: ", sum(finallengths))
	print_stats(finallengths)
	print("N50: ", calculate_N50(finallengths))
	n50_3 = calculate_N50_alternative(finallengths, total_length)
	print("N50 alternative: ", n50_3)	
	print("N50s: ", n50_1, n50_2, n50_3)
#	print("N50 factors: ", float(n50_2/n50_1), float(n50_3/n50_2), float(n50_3/n50_1))

	print("length of node_to_startpath.keys(): ", len(node_to_startpath.keys()))
	#for whole_paths: check which endpath of one node is equal to startpath of other node
	#start with one startpath that is not equal to any endpath 
	for node, sp in node_to_startpath.items():
		for node2, ep in node_to_endpath.items():		
			if sp == ep or sp == node_to_startpath[node2][::-1]:
				print("equal: ", node, node2)
		
	#merge the long paths
	#pairs of connected end nodes
	c_pairs = []	
	for node in unique_paths:
		path = node_to_extendedpaths[node][0]
		start = path[0]
		end = path[-1]	
		#avoiding to add pairs (a,b) and (b,a)		
		if (end,start) not in c_pairs:	
			c_pairs.append((start, end))
	# c_pairs = [(x,y), (z,a),(w,v)]
	#pairs of path beginning and end nodes
	p_pairs = []	
	pair_to_path = {}
	for n,p in node_to_path.items():
		p_pairs.append((p[0],p[-1]))
		pair_to_path[(p[0],p[-1])] = p
	#p_pairs = [(w,x),(y,z),(a,b)]	
	p = []
	cur_p = []
	resulting_merged_paths = []
	for pair in c_pairs:
		print("c_pair: ", pair)
		pair_path = compute_paths(node_to_extendedpaths, pair[0], node_to_path, node_to_end)
		#extend pair = (x,y) by end(x) or end(y)
		#pair = (node_to_end[pair[0]], node_to_end[pair[1]])
		#TODO: test if contains pair OR extended pair ...
		
		#is there a pair in p which includes pair?	
		included = False
		merge_with = []
		for ip in p:
			
			if any(x == y for x, y in zip(pair, ip)) or any(x == y for x, y in zip(pair[::-1], ip)):
				#pair is included
				#...merge
				included = True
				merge_with.append(ip)
						
			else:
				cur_p.append(ip)
		print("len(merge_with): ", len(merge_with))			
		if len(merge_with) == 1:
			ip = merge_with[0]		
			old_pair = pair
			pair, mergedpath = merge_pairs(pair,ip, node_to_extendedpaths[pair[0]][0], pair_to_path[ip])
			print("pair, old_pair: ", pair, old_pair)
			if pair == ip:
				print("equal pair: ", pair)
				
			else:			
				if pair[0] in old_pair:
					
					if pair[0] in node_to_path:
						print("pair 0: ", node_to_path[pair[0]])
						mergedpath = node_to_path[pair[0]][::-1] + mergedpath[1:]
					else:
						print("pair 0 not: ", node_to_path[node_to_end[pair[0]]])
						assert(node_to_end[pair[0]] in node_to_path)
						mergedpath = node_to_path[node_to_end[pair[0]]] + mergedpath[1:]
					pair = (node_to_end[pair[0]],pair[1])
				else:
					assert(pair[1] in old_pair)
				#	pair = (pair[0], node_to_end[pair[1]])
					if pair[1] in node_to_path:
						print("pair 1: ", node_to_path[pair[1]])
						mergedpath += node_to_path[pair[1]][1:]
					else:
						print("pair 1 not: ", node_to_path[node_to_end[pair[1]]])
						assert(node_to_end[pair[1]] in node_to_path)
						mergedpath += node_to_path[node_to_end[pair[1]]][::-1][1:]		
					pair = (pair[0], node_to_end[pair[1]])		
		
		if len(merge_with) == 2:
			pair, mergedpath = merge_pairs_multiple(pair,merge_with, pair_to_path, node_to_extendedpaths)
				
								
		#if cannot be merged to any existing path, append to list
		#if it was merged, 'pair' should contain the already merged pair and be added too	
		#TODO: extend only here?
		if not included:
			mergedpath = endpairs_to_paths[pair[0],pair[1]]
			pair = (node_to_end[pair[0]], node_to_end[pair[1]])	
				
		print("pair after: ", pair)
		#ensure that paths are written from start to end (in the right order of the tuple)
		if pair[0] == mergedpath[-1] and pair[-1] == mergedpath[0]:
			mergedpath = mergedpath[::-1]
		pair_to_path[pair] = mergedpath
		print(mergedpath)
		#TODO: how to store only the resulting paths?
		
		cur_p.append(pair)
		p= cur_p
		cur_p = []	
		print("p: ", p)	
		print("len p: ", len(p))
	pathlengths = []
	paths = []
	for i,pair in enumerate(p):
		print(pair, pair_to_path[pair][0], pair_to_path[pair][-1])
	#	assert(pair[0] == pair_to_path[pair][0] and pair[1] == pair_to_path[pair][-1])
		resulting_merged_paths.append(pair_to_path[pair])
		
		
		
		path = find_BFS_SP(graph, pair[0], pair[1])
		pathlength = compute_pathlength(path, bpdists)	
		if pathlength == 0:
			assert(len(path)==1)
			pathlength = len(nodeseqs[path[0]])		
		pathlengths.append(pathlength)	
		paths.append(path)		
#	print(pair_to_path[('utg004280l', 'utg007923l')])	
	testpath = phased_path_exists(graph, 'utg004280l', 'utg006698l', red_blue_nodes)
	print(testpath)		
	print(pathlengths)
	print(sum(pathlengths))
	if pathlengths:
		print("mean: ", mean(pathlengths))
	print([(k,v) for k,v in pair_to_path.items()][:2])
	
	#check which nodes are present in the final paths
#	nodes_in_finalpaths = set([i for n in paths for i in n ])
	nodes_in_finalpaths = set([i for n in resulting_merged_paths for i in n ])
	print("nodes part of final paths: ", len(nodes_in_finalpaths))
	
	#for those nodes which are not present in the paths, get their path from the set of unconnected paths
	nodes_outside_finalpaths = set([n for n in red_blue_nodes if n not in nodes_in_finalpaths])
	print("nodes outside final paths: ", len(nodes_outside_finalpaths))
	print("nodes outside final paths: ", nodes_outside_finalpaths)
	path_set = []
	path_set_lengths = set()
	special_ns = []
	for n in nodes_outside_finalpaths:
		if n in allnodes_to_path.keys():
			n_path = allnodes_to_path[n]
			path_set_lengths.add(node_to_pathlength[n_path[0]])
			if n_path not in path_set:
				path_set.append(n_path)
			#print("node, path: ",n, n_path)
		else:
			special_ns.append(n)
	
	print("paths not in final paths: ", len(path_set_lengths))
	print("paths not in final paths: ", path_set_lengths)
	print("paths not in final paths: ", list(path_set)[:5])
	final_finallengths = pathlengths + list(path_set_lengths)
	print("mean of all: ", mean(final_finallengths))
	print("no. of all: ", len(final_finallengths))
	print(final_finallengths)
	n50_4 = calculate_N50_alternative(final_finallengths, total_length)
	print("N50 alternative, final final: ", n50_4)
	
	finallist = resulting_merged_paths
	print("number of final paths before: ", len(finallist))
	print("final paths: ", [(p[0],p[-1]) for p in finallist])	
	
	for n,p in node_to_path.items():
		if n not in nodes_in_finalpaths:
			finallist.append(p)
	
	print("number of final paths after: ", len(finallist))
	print("final paths: ", [(p[0],p[-1]) for p in finallist])	
#	print("final paths: ", finallist)	
	
	#compute the sequences of the paths
	seqs = []
	seqs_filled = []
	for p in finallist:
		pathseq = ""
		pathseq_filled = ""
		if len(p) > 1:		
			for i in range(len(p)-1):
				start = p[i]
				end = p[i+1]
				if (start, end) not in bpdists:
					print("error in: ", p)
				dist = bpdists[(start, end)]
				sign = signs[(start, end)]
				seq = get_sequence(start, sign, dist, nodeseqs)
				if start in red_nodes:
					pathseq+= seq
					pathseq_filled += seq
				else:
					pathseq += ('N'*len(seq)) 
					if len(right_children[start])==2 and len(left_children[start])==2:
						pathseq_filled += seq
					else:
						pathseq_filled += ('N'*len(seq))
		else:
			node = p[0]
			seq = nodeseqs[node]
			pathseq = seq
			pathseq_filled = seq
		seqs.append(pathseq)		
		seqs_filled.append(pathseq_filled)	
	
	print("length of resulting seqs: ", [len(i) for i in seqs])
	print("fraction of 'N's included: ", [round(len([a for a in i if a=='N'])/len(i),4) for i in seqs])	
	print("summed length: ", str(sum(len(i) for i in seqs)/1000000)+ " Mb")
	print("fraction of 'N's of total length: ", round(sum(len([a for a in i if a=='N']) for i in seqs)/sum(len(i) for i in seqs), 4))	

	print("fraction of 'N's included, filled seqs: ", [round(len([a for a in i if a=='N'])/len(i),4) for i in seqs_filled])	
	print("fraction of 'N's of total length, filled seqs: ", round(sum(len([a for a in i if a=='N']) for i in seqs_filled)/sum(len(i) for i in seqs_filled), 4))
	
	n50_5 = calculate_N50_alternative([len(i) for i in seqs], total_length)
	print("N50 alternative: ", str(n50_5/1000000)+" Mb")
	
#	testlist = [12316690, 1599670, 27893808, 5317165, 82892, 29149, 32913, 38784, 36788, 39407, 30882, 31231, 30664, 30254, 38646, 35008, 35213]		
#	testn50 = calculate_N50_alternative(testlist, total_length)
#	print("N50 test: ", str(testn50/1000000)+" Mb")
	
	with open("phased_nodes_colors_" + chromosome+ "_" + haplotype +".csv",'w') as of:
		of.write("node,color\n")
		for n in nodes_in_finalpaths:
			if n in [i for tup in p for i in tup]:
				of.write(n+",#DE3163\n")
			else:
				of.write(n+",#2ecc71\n")	
		for n in nodes_outside_finalpaths:
			of.write(n+",#DFFF00\n")
	
	if write_outputlengths:			
		with open(chromosome+"_pathlengths.txt","a") as plfile:
			for p in [len(i) for i in seqs]:
				plfile.write(str(p)+'\n')
			
	if write_sequence:
		with open(chromosome+"_v0.1.fasta","a") as fasta:
			for i,p in enumerate(seqs):
				fasta.write(">"+"Altus_"+chromosome+"_"+ haplotype+"_"+"haplotig"+str(i)+'\n')				
				fasta.write(p+'\n')
		with open(chromosome+"_v0.2.fasta","a") as fasta:
			for i,p in enumerate(seqs_filled):
				fasta.write(">"+"Altus_"+chromosome+"_"+ haplotype+"_"+"haplotig"+str(i)+'\n')				
				fasta.write(p+'\n')								
#	
#	print("write phased nodes to file")	
#	write_phased(all_phased)	
#	
#	print("collect sequences")
#	nodeseqs = write_sequences(gfafile)
#	bpdists, sequences = basepair_length(gfafile, nodeseqs)
#
	
#	print("compute shortest path")
#	shortest = BFS_SP(graph, start, end, right_children, left_children)
#	print(shortest)
#	print("compute path length")
#	bp = compute_pathlength(shortest, bpdists)
#	print("length: ", bp)
#	print("compute path sequence")								
#	seq = compute_pathsequence(shortest, sequences)
#	
#	#add last node to path
#	seq += nodeseqs[end]
#	
#	with open(outfile,'a') as out:
#		out.write('>'+start+'_to_'+end+'\n')
#		out.write(seq+'\n')
		


	
	
				
	