import numpy
import pysam

def process_groups(groups):
	if not groups:
		return None
	pg = []
	for group in groups.split(","):
		ids = [int(x) for x in group.split(":")]
		if len(ids) == 2:
			pg.append(range(ids[0], ids[1] + 1))
		else:
			pg.append(ids)
	return pg

def split_interval(interval):
	chrom,coords = interval.split(":")
	start,end = [int(x) for x in coords.replace(",","").split("-")]
	return chrom,start,end

def bam2numreads(bamfile):
	result = pysam.idxstats(bamfile)
	return numpy.sum([int(row.strip().split("\t")[2]) for row in result])
	
def _treesort(order, nodeorder, nodecounts, tree):
	# From the Pycluster library, Michiel de Hoon
		# Find the order of the nodes consistent with the hierarchical clustering
	# tree, taking into account the preferred order of nodes.
	nNodes = len(tree)
	nElements = nNodes + 1
	neworder = numpy.zeros(nElements)
	clusterids = numpy.arange(nElements)
	for i in range(nNodes):
		i1 = tree[i].left
		i2 = tree[i].right
		if i1 < 0:
			order1 = nodeorder[-i1-1]
			count1 = nodecounts[-i1-1]
		else:
			order1 = order[i1]
			count1 = 1
		if i2 < 0:
			order2 = nodeorder[-i2-1]
			count2 = nodecounts[-i2-1]
		else:
			order2 = order[i2]
			count2 = 1
		# If order1 and order2 are equal, their order is determined
		# by the order in which they were clustered
		if i1 < i2:
			if order1 < order2:
				increase = count1
			else:
				increase = count2
			for j in range(nElements):
				clusterid = clusterids[j]
				if clusterid == i1 and order1 >= order2:
					neworder[j] += increase
				if clusterid == i2 and order1 < order2:
					neworder[j] += increase
				if clusterid == i1 or clusterid == i2:
					clusterids[j] = -i-1
		else:
			if order1 <= order2:
				increase = count1
			else:
				increase = count2
			for j in range(nElements):
				clusterid = clusterids[j]
				if clusterid == i1 and order1 > order2:
					neworder[j] += increase
				if clusterid == i2 and order1 <= order2:
					neworder[j] += increase
				if clusterid == i1 or clusterid == i2:
					clusterids[j] = -i-1
	return numpy.argsort(neworder)

def sort_tree(tree, order):
	# Adapted from the Pycluster library, Michiel de Hoon
	nnodes = len(tree)
	nodeindex = 0
	nodecounts = numpy.zeros(nnodes, int)
	nodeorder = numpy.zeros(nnodes)
	nodedist = numpy.array([node.distance for node in tree])
	for nodeindex in range(nnodes):
		min1 = tree[nodeindex].left
		min2 = tree[nodeindex].right
		if min1 < 0:
			index1 = -min1-1
			order1 = nodeorder[index1]
			counts1 = nodecounts[index1]
			nodedist[nodeindex] = max(nodedist[nodeindex],nodedist[index1])
		else:
			order1 = order[min1]
			counts1 = 1
		if min2 < 0:
			index2 = -min2-1
			order2 = nodeorder[index2]
			counts2 = nodecounts[index2]
			nodedist[nodeindex] = max(nodedist[nodeindex],nodedist[index2])
		else:
			order2 = order[min2]
			counts2 = 1
		counts = counts1 + counts2
		nodecounts[nodeindex] = counts
		nodeorder[nodeindex] = (counts1*order1+counts2*order2) / counts
	# Now set up order based on the tree structure
	index = _treesort(order, nodeorder, nodecounts, tree)
	return index


