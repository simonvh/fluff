#!/usr/bin/python
__author__ = 'george'

import multiprocessing
import os
import sys

import pysam
### External imports ###
import Pycluster
from numpy import array, hstack, arange, zeros

### My imports ###
from fluff.util import ( normalize_data, get_absolute_scale, split_ranges, 
        mirror_clusters, sort_tree )
from fluff.fluffio import load_heatmap_data, check_data
from fluff.color import parse_colors
from fluff.plot import heatmap_plot
import fluff.config as cfg

def heatmap(args):
    datafiles = args.datafiles
    for x in args.datafiles:
        if not os.path.isfile(x):
            print "ERROR: Data file '{0}' does not exist".format(x)
            sys.exit(1)
    for x in args.datafiles:
        if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
            print "Data file '{0}' does not have an index file. Creating an index file for {0}.".format(x)
            pysam.index(x)

    # Options Parser
    featurefile = args.featurefile
    datafiles = [x.strip() for x in args.datafiles]
    tracks = [os.path.basename(x) for x in datafiles]
    titles = [os.path.splitext(x)[0] for x in tracks]
    colors = parse_colors(args.colors)
    bgcolors = parse_colors(args.bgcolors)
    outfile = args.outfile
    extend_up = args.extend
    extend_down = args.extend
    fragmentsize = args.fragmentsize
    cluster_type = args.clustering[0].lower()
    merge_mirrored = args.merge_mirrored
    bins = (extend_up + extend_down) / args.binsize
    rmdup = args.rmdup
    rpkm = args.rpkm
    rmrepeats = args.rmrepeats
    ncpus = args.cpus
    distancefunction = args.distancefunction[0].lower()
    dynam = args.graphdynamics
    fontsize = args.textfontsize

    # Check for mutually exclusive parameters
    if dynam:
        if merge_mirrored:
            print "ERROR: -m and -g option CANNOT be used together"
            sys.exit(1)
        if distancefunction == 'e':
            print 'Dynamics can only be identified using Pearson correlation as metric.'
            print 'Assigning metric to Pearson correlation'
            distancefunction = 'p'

    # Warning about too much files
    if (len(tracks) > 4):
        print "Warning: Running fluff with too many files might make you system use enormous amount of memory!"
    
    # Method of clustering
    if (args.pick != None):
        pick = [i - 1 for i in split_ranges(args.pick)]
        if not all(i <= len(tracks) - 1 for i in pick):
            sys.stderr.write("You picked a non-existent file for clustering.\n")
            sys.exit(1)
    else:
        pick = range(len(datafiles))


    if not cluster_type in ["k", "h", "n"]:
        sys.stderr.write("Unknown clustering type!\n")
        sys.exit(1)
    # Number of clusters
    if cluster_type == "k" and not args.numclusters >= 2:
        sys.stderr.write("Please provide number of clusters!\n")
        sys.exit(1)
    # Distance function
    if not distancefunction in ["e", "p"]:
        sys.stderr.write("Unknown distance function!\n")
        sys.exit(1)
    else:
        if distancefunction == "e":
            METRIC = cfg.DEFAULT_METRIC
            print "Euclidean distance method"
        else:
            METRIC = "c"
            print "Pearson distance method"
    ## Get scale for each track
    tscale = [1.0 for track in datafiles]

    # Function to load heatmap data
    def load_data(featurefile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam,
                  guard=None):
        if guard is None:
            guard = []
        # Calculate the profile data
        data = {}
        regions = []
        print "Loading data"
        try:
            # Load data in parallel
            pool = multiprocessing.Pool(processes=ncpus)
            jobs = []
            for datafile in datafiles:
                jobs.append(pool.apply_async(load_heatmap_data, args=(
                featurefile, datafile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats,
                fragmentsize, dynam, guard)))
            for job in jobs:
                track, regions, profile, guard = job.get()
                data[track] = profile
        except Exception as e:
            sys.stderr.write("Error loading data in parallel, trying serial\n")
            sys.stderr.write("Error: {}\n".format(e))
            for datafile in datafiles:
                track, regions, profile, guard = load_heatmap_data(featurefile, datafile, amount_bins, extend_dyn_up,
                                                                   extend_dyn_down, rmdup, rpkm, rmrepeats,
                                                                   fragmentsize, dynam, guard)
                data[track] = profile
        return data, regions, guard

    # -g : Option to try and get dynamics
    # Extend features 1kb up/down stream
    # Cluster them in one bin
    # Cluster them in one bin
    guard = []
    amount_bins = bins
    extend_dyn_up = extend_up
    extend_dyn_down = extend_down
    if dynam:
        # load the data once to get the features which extend below 0
        guard = check_data(featurefile, extend_dyn_up, extend_dyn_down)
        extend_dyn_up = 1000
        extend_dyn_down = 1000
        amount_bins = 1

    # Load data for clustering
    data, regions, guard = load_data(featurefile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm,
                                         rmrepeats,
                                         fragmentsize, dynam, guard)
    
    # Normalize
    norm_data = normalize_data(data, cfg.DEFAULT_PERCENTILE)

    clus = hstack([norm_data[t] for i, t in enumerate(tracks) if (not pick or i in pick)])

    # Clustering
    if cluster_type == "k":
        print "K-means clustering"
        ## K-means clustering
        # PyCluster
        labels, _, nfound = Pycluster.kcluster(clus, args.numclusters, dist=METRIC)
        if not dynam and merge_mirrored:
            (i, j) = mirror_clusters(data, labels)
            while j:
                for track in data.keys():
                    data[track][labels == j] = [row[::-1] for row in data[track][labels == j]]
                for k in range(len(regions)):
                    if labels[k] == j:
                        (chrom, start, end, gene, strand) = regions[k]
                        if strand == "+":
                            strand = "-"
                        else:
                            strand = "+"
                        regions[k] = (chrom, start, end, gene, strand)
                n = len(set(labels))
                labels[labels == j] = i
                for k in range(j + 1, n):
                    labels[labels == k] = k - 1
                (i, j) = mirror_clusters(data, labels)

        ind = labels.argsort()

        # Hierarchical clustering
    elif cluster_type == "h":
        print "Hierarchical clustering"
        tree = Pycluster.treecluster(clus, method="m", dist=METRIC)
        labels = tree.cut(args.numclusters)
        ind = sort_tree(tree, arange(len(regions)))
    else:
        ind = arange(len(regions))
        labels = zeros(len(regions))

    # Load data for visualization if -g option was used
    if dynam:
        data, regions, guard = load_data(featurefile, bins, extend_up, extend_down, rmdup, rpkm, rmrepeats,
                                         fragmentsize, dynam, guard)

    f = open("{0}_clusters.bed".format(outfile), "w")
    for (chrom, start, end, gene, strand), cluster in zip(array(regions, dtype="object")[ind], array(labels)[ind]):
        if not gene:
            f.write("{0}\t{1}\t{2}\t.\t{3}\t{4}\n".format(chrom, start, end, cluster + 1, strand))
        else:
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, gene, cluster + 1, strand))
    f.close()
    # Save read counts
    readcounts = {}
    for i, track in enumerate(tracks):
        readcounts[track] = {}
        readcounts[track]['bins'] = []
        for idx, row in enumerate(data[track]):
            bins = ''
            for b in row:
                if not bins:
                    bins = '{0}'.format(b)
                else:
                    bins = '{0};{1}'.format(bins, b)
            readcounts[track]['bins'].append(bins)
    
    input_fileBins = open('{0}_readCounts.txt'.format(outfile), 'w')
    input_fileBins.write('Regions\t')
    for i, track in enumerate(titles):
        input_fileBins.write('{0}\t'.format(track))
    input_fileBins.write('\n')
    for i, track in enumerate(tracks):
        for idx in ind:
            input_fileBins.write('{0}:{1}-{2}\t'.format(regions[idx][0], regions[idx][1], regions[idx][2]))
            for i, track in enumerate(tracks):
                input_fileBins.write('{0}\t'.format(readcounts[track]['bins'][idx]))
            input_fileBins.write('\n')
        break
    input_fileBins.close()
 
    if not cluster_type == "k":
        labels = None

    scale = get_absolute_scale(args.scale, [data[track] for track in tracks])
    heatmap_plot(data, ind[::-1], outfile, tracks, titles, colors, bgcolors, scale, tscale, labels, fontsize)
