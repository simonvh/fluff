# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License
# 
import os
import sys
import tempfile

import HTSeq
import numpy
import pybedtools
import pysam

import fluff


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


def get_features_by_feature(track_a, track_b):
    """ Takes two BedTool tracks as input, returns overlapping features
    """
    if track_a.file_type != "bed" or track_b.file_type != "bed":
        raise ValueError, "Need BED files"
    for f in track_a:
        field_len_a = len(f.fields)
        break
    for f in track_b:
        field_len_b = len(f.fields)
        break
    i = track_a.intersect(track_b, wao=True, stream=False)
    tmp = tempfile.NamedTemporaryFile(delete=False, prefix="fluff")
    savedfile = i.saveas(tmp.name)
    tmp.flush()
    last = None
    features = []
    for line in tmp.readlines():
        vals = line.strip().split("\t")
        if field_len_a >= 6:
            feature = pybedtools.Interval(vals[0], int(vals[1]), int(vals[2]), strand=vals[5])
        else:
            feature = pybedtools.Interval(vals[0], int(vals[1]), int(vals[2]))
        if str(feature) != str(last):
            if len(features) > 0:
                if len(features) == 1 and features[0][1:3] == ['-1', '-1']:
                    yield last, []
                else:
                    yield last, features
            features = []
            last = feature
        features.append(vals[field_len_a:])
    if len(features) == 1 and features[0][1:3] == ['-1', '-1']:
        yield feature, []
    else:
        yield feature, features
    tmp.close()


class TrackWrapper():
    def __init__(self, fname):
        if fname.endswith("bam"):
            self.track = pysam.Samfile(fname, "rb")
            self.htseq_track = HTSeq.BAM_Reader(fname)
            self.ftype = "bam"
            self.chroms = self.track.references
        elif fname.endswith("bg") or fname.endswith("wig"):
            self.htseq_track = HTSeq.WiggleReader(fname, verbose=True)
            self.ftype = "wiggle"
        elif fname.endswith("gz") and os.path.exists(fname + ".tbi"):
            self.tabix_track = pysam.Tabixfile(fname)
            self.ftype = "tabix"
        else:
            self.track = pybedtools.BedTool(fname)
            self.ftype = "bed"
            # self.tabix_track = self.track.tabix()

    def __getitem__(self, window):
        """ 
        Retrieve all reads within a given window
        Arguments: window - a list or tuple containing chromosome, start, end and strand
        Returns a list of GenomicIntervals
        """
        chrom, start, end, strand = window
        if strand == None:
            strand = "."
        intervals = []
        if self.ftype == "bam":
            window = HTSeq.GenomicInterval(chrom, start, end, strand)
            for almnt in self.htseq_track[window]:
                if almnt and almnt.iv:
                    # ailmnt.iv.length = almn.iv
                    if strand == "." or strand == almnt.iv.strand:
                        intervals.append(almnt.iv)
        elif self.ftype == "bed":
            if strand == ".":
                feature = pybedtools.BedTool("{0} {1} {2}".format(chrom, start, end), from_string=True)
                s = False
            else:
                feature = pybedtools.BedTool("{0} {1} {2} 0 0 {3}".format(chrom, start, end, strand), from_string=True)
                s = True
            for read in self.track.intersect(feature, u=True, stream=True, s=s):
                intervals.append(HTSeq.GenomicInterval(chrom, read.start, read.end, str(read.strand)))
        return intervals

    def count(self, rmdup=False, rmrepeats=False):
        if self.ftype == "bam":
            if (not rmdup and not rmrepeats):
                return self.track.mapped
            c = 0
            for read in self.track:
                if (not rmdup or not read.flag & 0x0400):
                    # if (not rmrepeats or ('X0', 1) in read.tags or not 'X0' in [x[0] for x in read.tags]):
                    if (not rmrepeats) or read.mapq > 0:
                        c += 1
            return c
        if self.ftype == "bed":
            if rmrepeats:
                sys.stderr.write("Warning: rmrepeats has no result on BED files!")
            if rmdup:
                sys.stderr.write("Warning: rmdup has no result on BED files! (yet...;))")

    def read_length(self):
        if self.ftype == "bam":
            for read in self.track.fetch(until_eof=True):
                if read.alen:
                    return read.alen
        if self.ftype == "bed":
            for read in self.track:
                return read.end - read.start

    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        """ Generator 
        """
        if self.ftype == "bed":
            if rmrepeats:
                sys.stderr.write("Warning: rmrepeats has no result on BED files!")
        if self.ftype == "bam":
            for feature in track:
                min_strand = []
                plus_strand = []
                if feature.start < 0:
                    feature.start = 0
                if feature.chrom in self.chroms:
                    for read in self.track.fetch(feature.chrom, feature.start, feature.end):
                        # print read, read.mapq, read.tags
                        # if (not rmrepeats) or (('X0',1) in read.tags or not 'X0' in [x[0] for x in read.tags]):
                        if (not rmrepeats) or read.mapq > 0:
                            if read.is_reverse:
                                min_strand.append(read.pos)
                            else:
                                plus_strand.append(read.pos)

                    # Remove duplicates
                    if rmdup:
                        min_strand = sorted(set(min_strand))
                        plus_strand = sorted(set(plus_strand))
                    else:
                        min_strand = sorted(min_strand)
                        plus_strand = sorted(plus_strand)
                yield (feature, min_strand, plus_strand)
        elif self.ftype == "bed":
            for feature, features in get_features_by_feature(track, self.track):
                min_strand = []
                plus_strand = []

                for f in features:
                    if len(f) >= 6 and f[5] == "-":
                        min_strand.append(int(f[1]))
                    else:
                        plus_strand.append(int(f[1]))
                yield (feature, min_strand, plus_strand)

    def fetch_reads(self, interval, rmdup=False, rmrepeats=False):
        """ Generator 
        """
        chrom, start, end = interval

        if self.ftype == "bed":
            raise NotImplementedError

        if self.ftype == "bam":
            if chrom in self.chroms:
                for read in self.track.fetch(chrom, start, end):
                    if rmdup and (read.flag & 1024):
                        continue
                    if rmrepeats and read.mapq < 10:
                        continue
                    yield read

    def close(self):
        if self.ftype == "bam":
            self.track.close()

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
        chrom, start, end = interval
        profile = numpy.zeros(end - start, dtype="f")
        profile.fill(numpy.nan)

        if self.ftype == "bam":
            strand = {True: "-", False: "+"}

            for read in self.fetch_reads(interval, rmdup=rmdup, rmrepeats=rmrepeats):
                iv = HTSeq.GenomicInterval(chrom, read.pos, read.aend, strand[read.is_reverse])
                iv.length = fragmentsize
                region = profile[iv.start - start:iv.end - start]
                region[numpy.isnan(region)] = 0
                region += 1
        elif self.ftype == "wiggle":
            for iv, score in self.htseq_track:
                if iv.chrom == chrom:
                    if iv.start <= end and iv.end >= start:
                        if iv.start < start:
                            iv.start = start
                        if iv.end > end:
                            iv.end = end

                        profile[iv.start - start:iv.end - start] = score
        elif self.ftype == "tabix":
            for f in self.tabix_track.fetch(chrom, start, end):
                f = f.split()
                iv = HTSeq.GenomicInterval(f[0], int(f[1]) - 5, int(f[2]) + 5, ".")
                if iv.start <= end and iv.end >= start:
                    if iv.start < start:
                        iv.start = start
                    if iv.end > end:
                        iv.end = end

                    profile[iv.start - start:iv.end - start] = float(f[3])

        return profile


def _convert_value(v):
    """
    Returns 0 if v is not specified, otherwise try int
    """
    if v:
        try:
            v = int(v)
            return v
        except:
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
        result = get_binned_stats(clust_file,
                                  datafile,
                                  bins,
                                  rpkm=rpkm,
                                  rmdup=rmdup,
                                  rmrepeats=rmrepeats,
                                  fragmentsize=fragmentsize)
        result = [row.split("\t") for row in result]
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


def load_profile(interval, tracks, fragmentsize=200, rmdup=False, rmrepeats=False, reverse=False):
    profiles = []
    for track_group in tracks:
        if type(track_group) == type([]):
            profile_group = []
            for track in track_group:
                t = TrackWrapper(track)
                profile = t.get_profile(interval, fragmentsize, rmdup, rmrepeats)
                if reverse:
                    profile = profile[::-1]
                profile_group.append(profile)
        else:
            track = track_group
            t = TrackWrapper(track)
            profile_group = t.get_profile(interval, fragmentsize, rmdup, rmrepeats)
            if reverse:
                profile_group = profile_group[::-1]
        profiles.append(profile_group)
    return profiles


def get_free_track(overlap, start, end, max_end, min_gap):
    first = start - min_gap * max_end
    if first < 0:
        first = 0

    for i, track in enumerate(overlap):
        if max(track[start:end]) == 0:
            track[first:end + min_gap * max_end] += 1
            return overlap, i

    overlap.append(numpy.zeros(max_end, dtype="i"))
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


class SimpleFeature():
    def __init__(self, chrom, start, end, value, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.value = value
        self.strand = strand


class SimpleBed():
    def __init__(self, fname):
        self.f = open(fname)

    def __iter__(self):
        return self

    def next(self):
        line = self.f.readline()
        while line and (line[0] == "#" or line.startswith("track")):
            line = self.f.readline()
        if line:
            vals = line.strip().split("\t")
            start, end = int(vals[1]), int(vals[2])
            if len(vals) > 3:
                value = vals[3]
            else:
                value = 0
            if len(vals) > 5:
                if not (vals[5] is '+') or (vals[5] is '-'):
                    return SimpleFeature(vals[0], start, end, value, '+')
                else:
                    return SimpleFeature(vals[0], start, end, value, vals[5])
            elif len(vals) > 4:
                if not (vals[4] is '+') or (vals[4] is '-'):
                    return SimpleFeature(vals[0], start, end, value, '+')
                else:
                    return SimpleFeature(vals[0], start, end, value, vals[4])
            else:
                return SimpleFeature(vals[0], start, end, value, '+')
        else:
            self.f.close()
            raise StopIteration


def get_binned_stats(in_fname, data_fname, nbins, rpkm=False, rmdup=False, rmrepeats=False, fragmentsize=None,
                     split=False):
    track = TrackWrapper(data_fname)
    readlength = track.read_length()
    if not fragmentsize:
        fragmentsize = readlength
    total_reads = 1
    if rpkm:
        total_reads = track.count(rmdup, rmrepeats) / 1000000.0
    ret = []
    count = 1
    # Only use a BedTool if really necessary, as BedTools does not close open files
    # on object deletion
    if track.ftype == "bam":
        in_track = SimpleBed(in_fname)
    else:
        in_track = pybedtools.BedTool(in_fname)
    extend = fragmentsize - readlength
    for feature, min_strand, plus_strand in track.fetch_to_counts(in_track, rmdup, rmrepeats):
        binsize = (feature.end - feature.start) / float(nbins)
        row = []
        overlap = []
        min_strand = [x - (fragmentsize - readlength) for x in min_strand]
        bin_start = feature.start
        while int(bin_start + 0.5) < feature.end:
            num_reads = 0
            i = 0
            c = 0
            while i < len(min_strand) and min_strand[i] <= int(bin_start + binsize + 0.5):
                if min_strand[i] + fragmentsize <= int(bin_start + binsize + 0.5):
                    c += 1
                num_reads += 1
                i += 1
            min_strand = min_strand[c:]

            i = 0
            c = 0
            while i < len(plus_strand) and plus_strand[i] <= int(bin_start + binsize + 0.5):
                if plus_strand[i] + fragmentsize <= int(bin_start + binsize + 0.5):
                    c += 1
                num_reads += 1
                i += 1
            plus_strand = plus_strand[c:]

            if rpkm:
                per_kb = num_reads * (1000.0 / binsize)
                row.append(per_kb / total_reads)
            else:
                row.append(num_reads)
            bin_start += binsize
        if feature.strand == "-":
            row = row[::-1]
        ret.append([feature.chrom, feature.start, feature.end] + row)
        count += 1
    track.close()
    del in_track
    if split:
        return ret
    else:
        return ["\t".join([str(x) for x in row]) for row in ret]


def load_heatmap_data(featurefile, datafile, bins=100, up=5000, down=5000, rmdup=True, rpkm=False, rmrepeats=True,
                      fragmentsize=None, dynam=False, guard=[]):
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
    result = fluff.fluffio.get_binned_stats(tmp.name, datafile, bins, rpkm=rpkm, rmdup=rmdup, rmrepeats=rmrepeats,
                                            fragmentsize=fragmentsize)
    # Retrieve original order
    r_data = numpy.array([[float(x) for x in row.split("\t")[3:]] for row in result])
    return os.path.basename(datafile), regions, r_data, guard  # [r_order]


def check_data(featurefile, up=5000, down=5000):
    guard = []
    for i, line in enumerate(open(featurefile)):
        if line.startswith("#") or line[:5] == "track":
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
        if start >= 0:
            guard.append(True)
        else:
            guard.append(False)
    return guard
