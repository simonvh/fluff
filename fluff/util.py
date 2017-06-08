import re
import sys

import numpy
import pysam
from scipy.stats import scoreatpercentile, chisquare


def split_ranges(r):
    if not r:
        return r
    p = re.compile('[:-]')
    n = []
    for part in r.split(","):
        nums = [int(x) for x in p.split(part)]
        for i in range(nums[0], nums[-1] + 1):
            n.append(i)
    return n


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
    chrom, coords = interval.split(":")
    start, end = [int(x) for x in coords.replace(",", "").split("-")]
    return chrom, start, end


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
            order1 = nodeorder[-i1 - 1]
            count1 = nodecounts[-i1 - 1]
        else:
            order1 = order[i1]
            count1 = 1
        if i2 < 0:
            order2 = nodeorder[-i2 - 1]
            count2 = nodecounts[-i2 - 1]
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
                    clusterids[j] = -i - 1
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
                    clusterids[j] = -i - 1
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
            index1 = -min1 - 1
            order1 = nodeorder[index1]
            counts1 = nodecounts[index1]
            nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index1])
        else:
            order1 = order[min1]
            counts1 = 1
        if min2 < 0:
            index2 = -min2 - 1
            order2 = nodeorder[index2]
            counts2 = nodecounts[index2]
            nodedist[nodeindex] = max(nodedist[nodeindex], nodedist[index2])
        else:
            order2 = order[min2]
            counts2 = 1
        counts = counts1 + counts2
        nodecounts[nodeindex] = counts
        nodeorder[nodeindex] = (counts1 * order1 + counts2 * order2) / counts
    # Now set up order based on the tree structure
    index = _treesort(order, nodeorder, nodecounts, tree)
    return index


def normalize_data(data, percentile=75):
    norm_data = {}
    for track, ar in data.items():
        flat = ar.flatten()
        s = scoreatpercentile(flat[~numpy.isnan(flat)], percentile)
        if s == 0:
            sys.stderr.write(
                "Error normalizing track {0} as score at percentile {1} is 0, normalizing to maximum value instead\n".format(
                    track, percentile))
            x = ar / max(flat)
        else:
            x = ar / s
            # x[x <= 0.5] = 0
            x[x >= 1.0] = 1
        norm_data[track] = x
    return norm_data


def get_absolute_scale(scale, data, per_track=False):
    try:
        scale = float(scale)
        return scale
    except:
        if type(scale) == type("") and scale.endswith("%"):
            rel_scale = float(scale[:-1])

            if per_track:
                print "Hoe"
                s = [scoreatpercentile(d, rel_scale) for d in data]
                print s
                return s
            else:
                d = numpy.array(data).flatten()
                s = scoreatpercentile(d, rel_scale)
                # Set the scale to the minimum non-zero value, otherwise
                # the plot will show nothing
                if s == 0:
                    try:
                        s = min(d[d > 0])
                    except:
                        s = 1.0
                return s


def mirror_clusters(data, labels, cutoff=0.01):
    """
    Merge mirrored profiles based on a chi2 test of the mean profiles 
    Only if the profile is mirrored over all data tracks
    Returns the labels of the two matched mirrored tracks, if there is at least one match with a p-value
    greater than the cutoff.
    If not, return (None, None)
    """
    n = len(set(labels))
    if n == 1:
        return (None, None)
    mirror = dict([(i, {}) for i in range(n)])
    for track in data.keys():
        profiles = []
        for i in range(n):
            profiles.append(numpy.mean(data[track][labels == i], 0) + 1e-10)
        for i in range(n - 1):
            for j in range(i + 1, n):
                p = chisquare(profiles[i], profiles[j][::-1])[1]
                mirror[i].setdefault(j, []).append(p)
    result = []
    for i in mirror.keys():
        for j in mirror[i].keys():
            result.append([(i, j), mirror[i][j]])
    for (i, j), ps in sorted(result, cmp=lambda a, b: cmp(numpy.mean(a[1]), numpy.mean(b[1])))[::-1]:
        # print (i,j), ps, numpy.array(ps), cutoff
        if (numpy.array(ps) >= cutoff).all():
            return (i, j)
    return (None, None)
