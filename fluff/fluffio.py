# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License
# 
import HTSeq
import pysam
import pybedtools
import sys
import os
from numpy import zeros,min,max
import tempfile
import fluff
import numpy

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

    if track_a.file_type != "bed" or  track_b.file_type != "bed":
        raise ValueError, "Need BED files\n"
    
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
                if len(features) == 1 and features[0][1:3] == ['-1','-1']:
                    yield last, []
                else:
                    yield last, features  
            features = []
            last = feature
        
        features.append(vals[field_len_a:]) 

    if len(features) == 1 and features[0][1:3] == ['-1','-1']:
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
        else:
            self.track = pybedtools.BedTool(fname)
            self.ftype = "bed"
            #self.tabix_track = self.track.tabix()     

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
                    #ailmnt.iv.length = almn.iv
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
                intervals.append(HTSeq.GenomicInterval(chrom, read.start, read.end, read.strand))
        
        return intervals            
    
    
    def count(self, rmdup=False, rmrepeats=False):
        if self.ftype == "bam":
            if (not rmdup and not rmrepeats):
                return self.track.mapped
            c = 0
            for read in self.track:
                if (not rmdup or not read.flag & 0x0400):
                    if (not rmrepeats or ('X0', 1) in read.tags or not 'X0' in [x[0] for x in read.tags]):
                        c += 1
            return c
        if self.ftype == "bed":
            if rmrepeats:
                sys.stderr.write("Warning: rmrepeats has no result on BED files!\n")
            if rmdup:
                sys.stderr.write("Warning: rmdup has no result on BED files! (yet...;))")

    def read_length(self):
        if self.ftype == "bam":
            for read in self.track.fetch():
                if read.alen:
                    return read.alen
        if self.ftype == "bed":
            for read in self.track:
                return read.end - read.start
    
    def _add_read_to_list(self, read, min_strand, plus_strand, rmrepeats=False):
        if (not rmrepeats) or (('X0',1) in read.tags or not 'X0' in [x[0] for x in read.tags]):
            if read.is_reverse:
                min_strand.append(read.pos)
            else:
                plus_strand.append(read.pos)

    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        """ Generator 
        """
   
        if self.ftype == "bed":
            if rmrepeats:
                sys.stderr.write("Warning: rmrepeats has no result on BED files!\n")

        if self.ftype == "bam":
            for feature in track:
                min_strand = []
                plus_strand = []
                if feature.start < 0:
                    feature.start = 0

                if feature.chrom in self.chroms:
                    self.track.fetch(feature.chrom, feature.start, feature.end, 
                                 callback=lambda x: self._add_read_to_list(x, 
                                                                           min_strand, 
                                                                           plus_strand, 
                                                                           rmrepeats
                                                                           )
                                )
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

    def close(self):
        if self.ftype == "bam":
            self.track.close()

def get_profile(interval, track, fragmentsize=200):
    chrom,start,end = interval
    profile = zeros(end - start, dtype="i")
    t = TrackWrapper(track)
    for iv in t[(chrom, start, end, ".")]:
        iv.length = fragmentsize
        profile[iv.start - start:iv.end - start] += 1
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
        cluster_data.setdefault(_convert_value(f.name), []).append("{0}:{1}-{2}".format(f.chrom, f.start, f.end))
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
        result =  [row.split("\t") for row in result]
        data[os.path.basename(datafile)] = dict([["{0}:{1}-{2}".format(vals[0], vals[1], vals[2]), [float(x) for x in vals[3:]]] for vals in result])
    return data

def load_profile(interval, tracks, fragmentsize=200):
    profiles = []
    for track_group in tracks:
        if type(track_group) == type([]):
            profile_group = []
            for track in track_group:
                profile = get_profile(interval, track, fragmentsize)
                profile_group.append(profile)
        else:
                track = track_group
                profile_group = get_profile(interval, track, fragmentsize)
        profiles.append(profile_group)
    return profiles

def get_free_track(overlap, start, end, max_end):
    for i,track in enumerate(overlap):
        if max(track[start:end]) == 0:
            track[start:end] += 1
            return overlap, i
    overlap.append(zeros(max_end, dtype="i"))
    overlap[-1][start:end] += 1
    return overlap, len(overlap) - 1

def load_annotation(interval, fname):
    genes = []
    chrom, start, end = interval
    for line in open(fname):
        if not line.startswith("#") and not line.startswith("track"):
            vals = line.strip().split("\t")
            if len(vals) != 12:
                sys.stderr.write("Need BED 12 format for annotation\n")
                sys.exit(1)
            for i in [1,2,6,7]:
                vals[i] = int(vals[i])
            if vals[0] == chrom:
                if (vals[1] > start and vals[1] < end) or (vals[2] > start and vals[2] < end):
                    sys.stderr.write("Adding {0}\n".format(vals[3]))
                    genes.append(vals)
    
    if len(genes) == 0:
        return {}

    min_start = min([gene[1] for gene in genes])
    max_end = max([gene[2] for gene in genes])
    overlap = []
    gene_tracks = {}
    for gene in genes:
        overlap,i = get_free_track(overlap, gene[1]  - min_start, gene[2] - min_start, max_end - min_start)    
        if gene_tracks.has_key(i):
            gene_tracks[i].append(gene)
        else:
            gene_tracks[i] = [gene]
    return gene_tracks

class SimpleFeature():
    def __init__(self, chrom, start, end, gene, value, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene = gene
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
            value = vals[3]
            gene = vals[4]
            return SimpleFeature(vals[0], start, end, gene, value, vals[5])
        else:
            self.f.close()
            raise StopIteration
        
def get_binned_stats(in_fname, data_fname, nbins, rpkm=False, rmdup=False, rmrepeats=False, fragmentsize=None):
    track = TrackWrapper(data_fname)
    
    readlength = track.read_length()
    if not fragmentsize:
        fragmentsize = readlength
    
    #sys.stderr.write("Using fragmentsize %s\n" % fragmentsize)

    total_reads = 1
    if rpkm:
        #sys.stderr.write("Counting reads...\n")
        total_reads = track.count(rmdup, rmrepeats) / 1000000.0
        #sys.stderr.write("Done.\n")

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
            while i < len(min_strand) and min_strand[i] <= int(bin_start + binsize + 0.5):
                num_reads += 1
                i += 1
            while len(min_strand) > 0 and min_strand[0] + fragmentsize <= int(bin_start + binsize + 0.5):
                min_strand.pop(0)
            i = 0
            while i < len(plus_strand) and plus_strand[i] <= int(bin_start + binsize + 0.5):
                num_reads += 1
                i += 1
            while len(plus_strand) > 0 and plus_strand[0] + fragmentsize <= int(bin_start + binsize + 0.5):
                plus_strand.pop(0)
        
            if rpkm:
                per_kb = num_reads * (1000.0 / binsize)
                row.append(per_kb / total_reads)
            else:
                row.append(num_reads)

            bin_start += binsize
        
        if feature.strand == "-":
            row = row[::-1]
        
        ret.append( [feature.chrom, feature.start, feature.end] + row)
        
        count += 1
        #if count % 1000 == 0:
        #    sys.stderr.write("%s processed\n" % count)
    
    track.close()
    del in_track

    return ["\t".join([str(x) for x in row]) for row in ret]

def load_heatmap_data(featurefile, datafile, bins=100, up=5000, down=5000, rmdup=True, rpkm=False, rmrepeats=True, fragmentsize=None):
    tmp = tempfile.NamedTemporaryFile(delete=False, prefix="fluff")
    regions = []
    order = {}
    count = 0
    for line in open(featurefile):
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
        start,end = middle, middle
        if strand == "+":
            start -= up
            end += down
        else:
            start -= down
            end += up
        if start >= 0:
            regions.append([vals[0], start, end, gene, strand])
            order["{0}:{1}-{2}:{3}".format(vals[0], start, end, gene)] = count
            count += 1
            tmp.write("{0}\t{1}\t{2}\t{3}\t0\t{4}\n".format(vals[0], start, end, gene, strand))
    tmp.flush()
    
    result = fluff.fluffio.get_binned_stats(tmp.name, datafile, bins, rpkm=rpkm, rmdup=rmdup, rmrepeats=rmrepeats, fragmentsize=fragmentsize)
    
    # Retrieve original order
    #r_regions = ["{}:{}-{}".format(*row.split("\t")[:3]) for row in result]
    #r_order = numpy.array([order[region] for region in r_regions]).argsort()[::-1]
    r_data = numpy.array([[float(x) for x in row.split("\t")[3:]] for row in result])

    return os.path.basename(datafile), regions, r_data#[r_order]

