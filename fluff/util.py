import re
import sys

import numpy as np
import pysam
from scipy.stats import scoreatpercentile, chisquare
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import scale

# note: scoreatpercentile  become obsolete in the future.
# For Numpy 1.9 and higher, numpy.percentile provides all the functionality
# that scoreatpercentile provides. And it’s significantly faster.

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
            pg.append(list(range(ids[0], ids[1] + 1)))
        else:
            pg.append(ids)
    return pg


def split_interval(interval):
    chrom, coords = interval.split(":")
    start, end = [int(x) for x in coords.replace(",", "").split("-")]
    return chrom, start, end


def bam2numreads(bamfile):
    result = pysam.idxstats(bamfile)
    return np.sum([int(row.strip().split("\t")[2]) for row in result])


def _treesort(order, nodeorder, nodecounts, tree):
    # From the Pycluster library, Michiel de Hoon
    # Find the order of the nodes consistent with the hierarchical clustering
    # tree, taking into account the preferred order of nodes.
    nNodes = len(tree)
    nElements = nNodes + 1
    neworder = np.zeros(nElements)
    clusterids = np.arange(nElements)
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
    return np.argsort(neworder)

def normalize_data(data, percentile=75):
    norm_data = {}
    for track, ar in list(data.items()):
        flat = ar.flatten()
        s = scoreatpercentile(flat[~np.isnan(flat)], percentile)
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
                print("Hoe")
                s = [scoreatpercentile(d, rel_scale) for d in data]
                print(s)
                return s
            else:
                d = np.array(data).flatten()
                s = scoreatpercentile(d, rel_scale)
                # Set the scale to the minimum non-zero value, otherwise
                # the plot will show nothing
                if s == 0:
                    try:
                        s = min(d[d > 0])
                    except:
                        s = 1.0
                return s

def mycmp(a,b):
    """Wrap function of cmp in py2"""
    if a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0

def mirror_clusters(data, labels, cutoff=0.01):
    """
    Merge mirrored profiles based on a chi2 test of the mean profiles 
    Only if the profile is mirrored over all data tracks
    Returns the labels of the two matched mirrored tracks, if there is at least one match with a p-value
    greater than the cutoff.
    If not, return (None, None)
    """
    from functools import cmp_to_key

    n = len(set(labels))
    if n == 1:
        return (None, None)
    mirror = dict([(i, {}) for i in range(n)])
    for track in list(data.keys()):
        profiles = []
        for i in range(n):
            profiles.append(np.mean(data[track][labels == i], 0) + 1e-10)
        for i in range(n - 1):
            for j in range(i + 1, n):
                p = chisquare(profiles[i], profiles[j][::-1])[1]
                mirror[i].setdefault(j, []).append(p)
    result = []
    for i in list(mirror.keys()):
        for j in list(mirror[i].keys()):
            result.append([(i, j), mirror[i][j]])
    ### fixed for python 3 only
    key = cmp_to_key(lambda a, b:mycmp(np.mean(a[1]), np.mean(b[1])))
    for (i, j), ps in sorted(result, key=key)[::-1]:
        # print (i,j), ps, numpy.array(ps), cutoff
        if (np.array(ps) >= cutoff).all():
            return (i, j)
    return (None, None)


def cluster_profile(cluster_data, cluster_type="k", numclusters=3, dist="euclidean", random_state=None):
    """Cluster profiles for heatmap

    Takes a matrix and clusters either with kmeans or hierarchical clustering.
    Distance can be either euclidean or pearson. 

    Parameters
    ----------
    cluster_data :  array_like
        Data to cluster.

    cluster_type : str, optional
        Either 'k' for kmeans, 'h' for hierarchical or 'n' for no clustering.
        If cluster_type equals None, data is also not clustered.

    numclusters : int, optional
        Number of clusters.

    dist : str, optional
        Distance metric, either 'euclidean' or 'pearson'.

    Returns
    -------

    ind : array
        Indices of sorted input.

    labels : array 
        Cluster labels.
    """
    if dist not in ["euclidean", "pearson"]:
        raise ValueError("distance can be either 'euclidean' or 'pearson'")
    # Clustering
    if dist == "pearson":
        cluster_data = np.apply_along_axis(scale, 1, cluster_data)
    
    if cluster_type == "k":
        print("K-means clustering")
        ## K-means clustering
       
        k = KMeans(n_clusters=numclusters, random_state=random_state)
        labels = k.fit(cluster_data).labels_
        ind = labels.argsort()

        # Hierarchical clustering
    elif cluster_type == "h":
        print("Hierarchical clustering")
        a = AgglomerativeClustering(
                n_clusters=numclusters, 
                linkage="complete"
                )
        a.fit(cluster_data)
        labels = a.labels_
        c = a.n_leaves_
        t = {x:[x] for x in range(a.n_leaves_)}
        for x in a.children_:
            t[c] = t[x[0]] + t[x[1]]
            c += 1  
        ind = t[c - 1]
    # No clustering
    elif cluster_type in ["n", None]:
        ind = np.arange(len(cluster_data))
        labels = np.zeros(len(cluster_data))
    else:
        raise ValueError("Invalid value for cluster_type")

    return ind, labels


