#!/usr/bin/python
import sys


def find_BFS_SP(graph, start, goal):
	visited = []
	queue = [[start]]
     
	#desired node is reached
	if start == goal:
		print("Same Node")
		return(queue.pop(0))
     
	while queue:
		path = queue.pop(0)
		node = path[-1]
         
		if node not in visited:
				neighbours = graph[node]
             
				for neighbour in neighbours:
					new_path = list(path)
					new_path.append(neighbour)
					queue.append(new_path)
                 
					if neighbour == goal:
				#		print("Shortest path = ", *new_path)
						return(new_path)
				visited.append(node)
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
	distances = {}
	sequences = {}
	dist = 0
	with open(gfapath) as gfa:
		for line in gfa:
			parts = line.strip().split('\t')
			if parts[0] == 'L':
				if len(parts) > 6:
					dist = int(parts[6].split(':')[-1])
				else:
					dist = len(nodeseqs[parts[1]]) - int(parts[5].split('M')[0])
				distances[(parts[1],parts[3])] = dist
				sequences[(parts[1],parts[3])] = get_sequence(parts[1], parts[2],dist, nodeseqs)
				#compute overlap between both 
	return(distances, sequences)
	
def write_adjacency(gfafile):
	node_to_list = {}
	with open(gfafile) as gfa:
		for i,line in enumerate(gfa):
			parts = line.strip().split('\t')
			if parts[0] == 'S':
				node_to_list[parts[1]] = []
			if parts[0] == 'L':
				left = parts[1]
				right = parts[3]
				node_to_list[left].append(right)
	return(node_to_list)	
	
def compute_pathlength(nodes, distances):
	dist = 0	
	for i in range(len(nodes)-1):
		dist += distances[(nodes[i],nodes[i+1])]	
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
			if line[0] == 'S':
				name = line.strip().split('\t')[1]
				seq = line.strip().split('\t')[2]
				nodeseqs[name] = seq
	return(nodeseqs)


if __name__=='__main__':
	gfafile = sys.argv[1]
	outfile = sys.argv[2]
	#start and end node of path
	start = sys.argv[3]
	end = sys.argv[4]
	
	print("create graph")
	graph = write_adjacency(gfafile)
	print("collect sequences")
	nodeseqs = write_sequences(gfafile)
	bpdists, sequences = basepair_length(gfafile, nodeseqs)

	print("compute shortest path")
	shortest = find_BFS_SP(graph, start,end)
	print(shortest)
	print("compute path length")
	bp = compute_pathlength(shortest, bpdists)
	print("compute path sequence")								
	seq = compute_pathsequence(shortest, sequences)
	
	#add last node to path
	seq += nodeseqs[end]
	
	with open(outfile,'w') as out:
		out.write('> path'+'\n')
		out.write(seq+'\n')
		
			
	