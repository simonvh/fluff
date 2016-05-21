import os
import re
import sys
import tempfile
from collections import Counter
from warnings import warn

import HTSeq
import numpy as np
import pybedtools
import pysam
import pyBigWig

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

class Track(object):
    _registry = []
    _filetypes = []

    class __metaclass__(type):
        """ Keep track of "plugin" classes
        """
        def __init__(cls, name, bases, attrs):
            cls._registry.append((name,cls))

    def __init__(self, fname):
        """
        """

        raise NotImplementedError("please instantiate subclass")

    @classmethod
    def filetypes(cls):
        """
        Return all supported filetypes of the subclasses of Track
        
        Returns
        -------
        list
            List of supported filetypes
        """
        
        return list(set([ftype for _,t in cls._registry for ftype in t._filetypes]))

    def _get_interval(self, interval):
        """
        Translate interval to tuple of (chrom, start, end)
        
        Params
        ------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        Returns
        -------
        tuple
            (chrom, start, end)
        
        """
        try:
            chrom, start, end = interval
        except:
            chrom, start, end = re.split(r'[:-]', interval)
            start, end = int(start), int(end)

        return (chrom, start, end)

    @classmethod
    def load(self, fname, *args, **kwargs):
        _, ftype = os.path.splitext(fname)
        ftype = ftype.strip(".")
        for name, cls in self._registry:
            if ftype in cls._filetypes:
                return cls(fname, *args, **kwargs)
        raise ValueError, "can't guess type of file {}".format(fname)

class BinnedMixin(object):
    def binned_stats(self, in_fname, nbins, rpkm=False, split=False):
        readlength = self.read_length()
        fragmentsize = self.fragmentsize
        if not fragmentsize:
            fragmentsize = readlength
        total_reads = 1
        if rpkm:
            total_reads = self.count() / 1000000.0
        ret = []
        count = 1
        # Only use a BedTool if really necessary, as BedTools does not close open files
        # on object deletion
        if self.ftype == "bam":
            in_track = SimpleBed(in_fname)
        else:
            in_track = pybedtools.BedTool(in_fname)
        
        extend = fragmentsize - readlength
        for feature, min_strand, plus_strand in self.fetch_to_counts(in_track):
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
        
        del in_track
        if split:
            return ret
        else:
            return ["\t".join([str(x) for x in row]) for row in ret]

class BamTrack(BinnedMixin, Track):
    _filetypes = ["bam"] 
    
    def __init__(self, fname, **kwargs):
        """
        Track interface to a BAM file

        Parameters
        ----------
        fname: str
            filename of BAM file
        
        fragmentsize : int, optional
            Reads are extended to fragmentsize before summarizing the profile.
            If fragmentsize is None, the read length is used.
        
        rmdup : bool, optional
            Ignore duplicate reads if True, default False

        rmrepeats : bool, optional
            Ignore reads with mapping quality 0 (multi-mapped reads) if 
            True, default False

        """
        
        self.rmdup = kwargs.get("rmdup", False)
        self.rmrepeats = kwargs.get("rmrepeats", False)
        self.fragmentsize = kwargs.get("fragmentsize", None)

        if fname.endswith("bam"):
            self.track = pysam.AlignmentFile(fname, "rb")
            self.ftype = "bam"
            self.chroms = self.track.references
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def count(self):
        """
        Count total number of reads in file

        Returns
        -------
        int
            Number of reads
        """
        
        if (not self.rmdup and not self.rmrepeats):
            try:
                return self.track.mapped
            except:
                pass

        c = 0
        for read in self.track:
            # duplicates
            if (not self.rmdup or not read.flag & 0x0400):
                # multi-mappers / mapping quality 0
                if (not self.rmrepeats) or read.mapq > 0:
                    c += 1
        return c
      
    def read_length(self):
        """
        Return the read length

        This function returns the read length of the first read where it is defined

        Returns
        -------
        int
            Read length
        """
        it = self.track.head(100)
        lengths = [read.infer_query_length(always=False) for read in it]
        return Counter(lengths).most_common(1)[0][0]

    def fetch_to_counts(self, track):
        """ 
        Yields the number of reads for each feature in track

        Parameters
        ----------
        track : Object of <TODO: what exactly do we expect?>
            Track with features

        Yields
        ------
        tuple
            Return value consists of:
                feature
                number of reads on minus strand
                number of reads on plus strand
        """
        for feature in track:
            min_strand = []
            plus_strand = []
            if feature.start < 0:
                feature.start = 0
            if feature.chrom in self.chroms:
                for read in self.track.fetch(feature.chrom, feature.start, feature.end):
                    if (not self.rmrepeats) or read.mapq > 0:
                        if read.is_reverse:
                            min_strand.append(read.pos)
                        else:
                            plus_strand.append(read.pos)

                # Remove duplicates
                if self.rmdup:
                    min_strand = sorted(set(min_strand))
                    plus_strand = sorted(set(plus_strand))
                else:
                    min_strand = sorted(min_strand)
                    plus_strand = sorted(plus_strand)
            yield (feature, min_strand, plus_strand)
     
    def fetch_reads(self, interval, strand=None):
        warn("fetch_reads is deprecated, please use fetch", DeprecationWarning)    

    def fetch(self, interval, strand=None):
        """ 
        Retrieve all reads within a given window
        
        Parameters
        ----------
        window : list, tuple or str
            If window is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        strand : str, optional
            Either '+' or '-'. By default all reads are returned.

        Yields
        ------
        AlignedSegment
            Yields pysam AlignedSegment objects.
        """

        chrom, start, end = self._get_interval(interval)

        if chrom in self.chroms:
            for read in self.track.fetch(chrom, start, end):
                # duplicate reads
                if self.rmdup and (read.flag & 1024):
                    continue
                # multimappers / low mapping quality 
                if self.rmrepeats and read.mapping_quality < 10:
                    continue
                if strand:
                    if strand == "+" and read.is_reverse:
                        continue
                    elif strand == "-" and not read.is_reverse:
                        continue
                yield read

    def close(self):
        self.track.close()

    def get_profile(self, interval):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """
        chrom, start, end = self._get_interval(interval)
        profile = np.zeros(end - start, dtype="f")
        profile.fill(np.nan)

        strand = {True: "-", False: "+"}

        for read in self.fetch(interval):
            iv = HTSeq.GenomicInterval(
                    chrom, 
                    read.reference_start, 
                    read.reference_end, strand[read.is_reverse]
                    )
            if self.fragmentsize:
                iv.length = self.fragmentsize
            region = profile[iv.start - start:iv.end - start]
            region[np.isnan(region)] = 0
            region += 1

        return profile

class BedTrack(BinnedMixin, Track):
    _filetypes = ["bed"]
    
    def __init__(self, fname, **kwargs):
        """
        Parameters
        ----------
        
        fragmentsize : int, optional
            Reads are extended to fragmentsize before summarizing the profile.
            If fragmentsize is None, the read length is used.
        """ 
        
        self.fragmentsize = kwargs.get("fragmentsize", None)

        self.track = pybedtools.BedTool(fname)
        self.ftype = "bed"

    def count(self):
        """
        Count total number of features in file

        Returns
        -------
        int
            Number of features
        """
        
        return self.track.count()

    def read_length(self):
        if self.ftype == "bed":
            for read in self.track:
                return read.end - read.start

    def fetch_to_counts(self, track):
        """ 
        Yields the number of features for each feature in track

        Parameters
        ----------
        track : Object of <TODO: what exactly do we expect?>
            Track with features

        Yields
        ------
        tuple
            Return value consists of:
                feature
                number of reads on minus strand
                number of reads on plus strand
        """

        for feature, features in self._get_features_by_feature(track):
            min_strand = []
            plus_strand = []

            for f in features:
                if len(f) >= 6 and f[5] == "-":
                    min_strand.append(int(f[1]))
                else:
                    plus_strand.append(int(f[1]))
            yield (feature, min_strand, plus_strand)

    def fetch_reads(self, interval):
        warn("fetch_reads is deprecated, please use fetch", DeprecationWarning)    
  
    def _get_features_by_feature(self, track_a):
        """
        Return overlapping features 
        
        Parameters
        ----------
        
        track_a : BedTrack object 
        """

        track_b = self.track
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
    
    def _interval_bedtool(self, interval, strand=None):
        """
        Convert an interval to a BedTool

        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        strand : str, optional
            Either '+' or '-'. Default is no strand.
        
        Returns
        -------
        BedTool object
        """

        chrom, start, end = self._get_interval(interval)
        
        if strand == None:
            strand = "."
        
        if strand == ".":
            feature = pybedtools.BedTool(
                    "{0} {1} {2}".format(
                        chrom, start, end), 
                    from_string=True)
            s = False
        else:
            feature = pybedtools.BedTool(
                    "{0} {1} {2} 0 0 {3}".format(
                        chrom, start, end, strand), 
                    from_string=True)
            s = True
        return feature

    def fetch(self, interval, strand=None):
        """ 
        Retrieve all reads within a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        strand : str, optional
            Either '+' or '-'. By default all reads are returned.
        
        Yields
        ------
        GenomicInterval
            Yields HTSeq GenomicInterval objects.
        """
        
        feature = self._interval_bedtool(interval, strand=strand)
        chrom, start, end = self._get_interval(interval)
        for read in self.track.intersect(feature, u=True, stream=True, s=strand in ["+", "-"]):
            yield HTSeq.GenomicInterval(
                    chrom, 
                    read.start, 
                    read.end, 
                    str(read.strand))

    def close(self):
        pass

    def get_profile(self, interval):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        fragmentsize : int, optional
            Reads are extended to fragmentsize before summarizing the profile.
            If fragmentsize is None, the read length is used.
        
        rmdup : bool, optional
            Don't return duplicate reads if True, default False

        rmrepeats : bool, optional
            Don't return reads with mapping quality 0 (multi-mapped reads) if 
            True, default False

        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """
 
        chrom, start, end = self._get_interval(interval)
        profile = np.zeros(end - start, dtype="f")
        profile.fill(np.nan)

        for f in self.fetch(interval):
            iv = HTSeq.GenomicInterval(
                    chrom, 
                    f.start, 
                    f.end,
                    f.strand
                    )
            if self.fragmentsize:
                iv.length = self.fragmentsize
            region = profile[iv.start - start:iv.end - start]
            region[np.isnan(region)] = 0
            region += 1

        return profile

class WigTrack(Track):
    _filetypes = ["bg", "wig"]
    
    def __init__(self, fname, **kwargs):
        if fname.endswith("bg") or fname.endswith("wig"):
            self.htseq_track = HTSeq.WiggleReader(fname, verbose=True)
            self.ftype = "wig"
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """
        
        chrom, start, end = self._get_interval(interval)
        profile = np.zeros(end - start, dtype="f")
        profile.fill(np.nan)

        for iv, score in self.htseq_track:
            if iv.chrom == chrom:
                if iv.start <= end and iv.end >= start:
                    if iv.start < start:
                        iv.start = start
                    if iv.end > end:
                        iv.end = end

                    profile[iv.start - start:iv.end - start] = score

        return profile

class BigWigTrack(Track):
    _filetypes = ["bw"]
    
    def __init__(self, fname, **kwargs):
        if fname.endswith("bw"):
            self.track = pyBigWig.open(fname)
            self.ftype = "bw"
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """

        
        chrom, start, end = self._get_interval(interval)
        return np.array(self.track.values(chrom, start, end)) 


class TabixTrack(Track):
    def __init__(self, fname, **kwargs):
        if fname.endswith("gz") and os.path.exists(fname + ".tbi"):
            self.tabix_track = pysam.Tabixfile(fname)
            self.ftype = "tabix"
        else:
            raise NotImplementedError

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
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
