# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License
# 
import os
import sys
import tempfile

import numpy as np
import pybedtools

from fluff.track import Track

def is_equal_feature(feature, vals):
    if not vals:
        return False
    if feature.chrom != vals[0]:
        return False
    if feature.start != int(vals[1]):
        return False
    if feature.end != int(vals[2]):
        return False
    return True

def _convert_value(v):
    """
    Returns 0 if v is not specified, otherwise try int
    """
    if v:
        try:
            v = int(v)
            return v
        except ValueError:
            return v
    return 0

def load_bed_clusters(bedfile):
    """
    Reads a BED file, using the fourth column as cluster number
    Arguments: bedfile - a 4-column BED file
    Returns: a hash with cluster numbers as key, and a list of genomic locations as value
    """
    cluster_data = {}
    track = pybedtools.BedTool(bedfile)
    for f in track:
        cluster_data.setdefault(_convert_value(f.score), []).append("{0}:{1}-{2}".format(f.chrom, f.start, f.end))
    return cluster_data

def load_cluster_data(clust_file, datafiles, bins, rpkm, rmdup, rmrepeats, fragmentsize=None):
    data = {}
    for datafile in datafiles:
        result = []
        track = Track.load(datafile,
                rmdup=rmdup,
                rmrepeats=rmrepeats,
                fragmentsize=fragmentsize)
        result = track.binned_stats(clust_file,
                                  bins,
                                  split=True,
                                  rpkm=rpkm
                                  )
        data[os.path.basename(datafile)] = dict(
                [["{0}:{1}-{2}".format(vals[0], vals[1], vals[2]), [float(x) for x in vals[3:]]] for vals in result])
    return data

def load_read_counts(readCounts):
    data = {}
    indexes = {}
    titles = []
    for line in open(readCounts):
        if line.startswith('Regions'):
            idx = 0
            for datafile in line.split('\t')[1:]:
                if datafile.strip():
                    titles.append(datafile.strip())
                    data[datafile.strip()] = {}
                    indexes[idx] = datafile.strip()
                    idx += 1
        else:
            for idx, binsline in enumerate(line.split('\t')[1:]):
                if binsline.strip():
                    data[indexes[idx]][line.split('\t')[0]] = [float(x) for x in binsline.split(';')]
    return titles, data

def get_free_track(overlap, start, end, max_end, min_gap):
    first = start - min_gap * max_end
    if first < 0:
        first = 0

    for i, track in enumerate(overlap):
        if max(track[start:end]) == 0:
            track[first:end + min_gap * max_end] += 1
            return overlap, i

    overlap.append(np.zeros(max_end, dtype="i"))
    overlap[-1][first:end + min_gap * max_end] += 1
    # overlap[-1][start- min_gap * max_end:end + min_gap * max_end] += 1
    return overlap, len(overlap) - 1


def load_annotation(interval, fname, min_gap=0.05, vis="stack"):
    genes = []
    chrom, start, end = interval
    for line in open(fname):
        if not line.startswith("#") and not line.startswith("track"):
            vals = line.strip().split("\t")
            for i in [1, 2, 6, 7]:
                if len(vals) > i:
                    vals[i] = int(vals[i])
            if vals[0] == chrom:
                if vals[1] <= end and vals[2] >= start:
                    # sys.stderr.write("Adding {0}\n".format(vals[3]))
                    genes.append(vals)
    if len(genes) == 0:
        return {}
    min_start = min([gene[1] for gene in genes])
    max_end = max([gene[2] for gene in genes])
    overlap = []
    gene_tracks = {}
    for gene in sorted(genes, key=lambda x: x[1]):
        if vis == "stack":
            overlap, i = get_free_track(overlap, gene[1] - min_start, gene[2] - min_start, max_end - min_start, min_gap)
        elif vis == "merge":
            i = 0
        else:
            sys.stderr.write("Unknown visualization")
        if gene_tracks.has_key(i):
            gene_tracks[i].append(gene)
        else:
            gene_tracks[i] = [gene]
    return gene_tracks

def load_heatmap_data(featurefile, datafile, bins=100, up=5000, down=5000, rmdup=True, rpkm=False, rmrepeats=True,fragmentsize=None, dynam=False, guard=None):
    if guard is None:
        guard = []
    
    tmp = tempfile.NamedTemporaryFile(delete=False, prefix="fluff")
    regions = []
    order = {}
    count = 0
    hashcounter = 0
    if not guard and dynam:
        filt = True
    else:
        filt = False
    for i, line in enumerate(open(featurefile)):
        if line.startswith("#") or line[:5] == "track":
            hashcounter += 1
            continue
        vals = line.strip().split("\t")
        strand = "+"
        gene = ""
        if len(vals) >= 6:
            strand = vals[5]
        if len(vals) >= 4:
            gene = vals[3]
        middle = (int(vals[2]) + int(vals[1])) / 2
        start, end = middle, middle
        if strand == "+":
            start -= up
            end += down
        else:
            start -= down
            end += up
        if filt:
            if start >= 0:
                guard.append(True)
            else:
                guard.append(False)
        if not filt and start >= 0:
            if not dynam or guard[i - hashcounter]:
                regions.append([vals[0], start, end, gene, strand])
                order["{0}:{1}-{2}".format(vals[0], start, end)] = count
                count += 1
                tmp.write("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(vals[0], start, end, gene, strand))
    tmp.flush()
    track = Track.load(datafile,
            rmdup=rmdup,
            rmrepeats=rmrepeats,
            fragmentsize=fragmentsize)

    result = track.binned_stats(tmp.name, bins, split=True, rpkm=rpkm)
    # Retrieve original order
    r_data = np.array([[float(x) for x in row[3:]] for row in result])
    return os.path.basename(datafile), regions, r_data, guard  # [r_order]


def check_data(featurefile, up=5000, down=5000):
    guard = []
    for line in open(featurefile):
        if line.startswith("#") or line[:5] == "track":
            continue
        vals = line.strip().split("\t")
        strand = "+"
        
        if len(vals) >= 6:
            strand = vals[5]
        
        middle = (int(vals[2]) + int(vals[1])) / 2
        start, end = middle, middle
        if strand == "+":
            start -= up
            end += down
        else:
            start -= down
            end += up
        if start >= 0:
            guard.append(True)
        else:
            guard.append(False)
    return guard
