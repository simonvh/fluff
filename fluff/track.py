import os
import pyBigWig
import re
import tempfile
from collections import Counter
from warnings import warn

import HTSeq
import numpy as np
import pybedtools
import pysam
from scipy.stats import binned_statistic


class SimpleFeature(object):
    def __init__(self, chrom, start, end, value, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.value = value
        self.strand = strand

class SimpleBed(object):
    """
    BED file as a simple iterator
    """
    
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

    track_type = "profile"
    
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
        except Exception:
            chrom, start, end = re.split(r'[:-]', interval)
            start, end = int(start), int(end)

        return (chrom, start, end)

    @classmethod
    def load(self, fname, *args, **kwargs):
        """
        Load a track in one of the following formats:
            bam, bed, wig, bg, bw, tabix (bed.gz, bg.gz, wig.gz)
        The format is guessed by the file extension.

        Parameters
        ----------
        fname : str
            Filename

        Returns
        -------
        Track object of the specified type
        """
        
        _, ftype = os.path.splitext(fname)
        ftype = ftype.strip(".")
        for _, cls in self._registry:
            for filetype in cls._filetypes:
                if filetype.endswith(ftype):
                    return cls(fname, *args, **kwargs)
        raise ValueError("can't guess type of file {}".format(fname))

class BinnedMixin(object):
    def binned_stats(self, in_fname, nbins, split=False, **args):
        rpkm = args.get("rpkm", False)
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
        
        #extend = fragmentsize - readlength
        for feature, min_strand, plus_strand in self.fetch_to_counts(in_track):
            binsize = (feature.end - feature.start) / float(nbins)
            row = []
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
            return ["\t".join([str(x) for x in r]) for r in ret]

class BamTrack(BinnedMixin, Track):
    _filetypes = ["bam"] 
    track_type = "feature"

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

        if fname.split(".")[-1] in self._filetypes:
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
            except Exception:
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
     
    def fetch_reads(self, args, **kwargs):
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

    def get_profile(self, interval, **kwargs):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        scalepm : bool, optional
            Scale profile to per million reads
        
        scalefactor : float, optional
            Scale profile by this factor, default 1.0
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """
        scalefactor = kwargs.get("scalefactor", 1.0) 
        scalepm = kwargs.get("scalepm", False)
        if scalepm:
            scalefactor = scalefactor * 1e6 / float(self.count())
        
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
            else:
                for blockstart, blockend in read.get_blocks():
                    region = profile[blockstart - start:blockend - start]
                    region[np.isnan(region)] = 0
                    region += 1
        profile = profile * scalefactor
        return profile

class BedTrack(BinnedMixin, Track):
    _filetypes = ["bed"]
    track_type = "feature"
    
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

    def fetch_reads(self, *args, **kwargs):
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
            raise ValueError("Need BED files")
        for f in track_a:
            field_len_a = len(f.fields)
            break
        i = track_a.intersect(track_b, wao=True, stream=False)
        tmp = tempfile.NamedTemporaryFile(delete=False, prefix="fluff")
        _ = i.saveas(tmp.name)
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
        
        if strand is None:
            strand = "."
        
        if strand == ".":
            feature = pybedtools.BedTool(
                    "{0} {1} {2}".format(
                        chrom, start, end), 
                    from_string=True)
        else:
            feature = pybedtools.BedTool(
                    "{0} {1} {2} 0 0 {3}".format(
                        chrom, start, end, strand), 
                    from_string=True)
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

    def get_profile(self, interval, **kwargs):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        scalepm : bool, optional
            Scale profile to per million reads
        
        scalefactor : float, optional
            Scale profile by this factor, default 1.0
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array

        """
        scalefactor = kwargs.get("scalefactor", 1.0)
        scalepm = kwargs.get("scalepm", False)
        if scalepm:
            scalefactor = scalefactor * 1e6 / float(self.count())

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
    _filetypes = ["bg", "wig", "bdg", "bedGraph"]
    
    def __init__(self, fname, **kwargs):
        self.fname = fname
        
        if fname.split(".")[-1] in self._filetypes:
            self.track = pybedtools.BedTool(fname)
            self.ftype = "wig"
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def get_profile(self, interval, **kwargs):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        scalefactor : float, optional
            Scale profile by this factor, default 1.0
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array
        """
        scalefactor = kwargs.get("scalefactor", 1.0)
        
        chrom, start, end = self._get_interval(interval)
        int_bed = pybedtools.BedTool(
                "{} {} {}".format(chrom, start, end), 
                from_string=True)
        
        profile = np.zeros(end - start, dtype="f")
        profile.fill(np.nan)

        for f in self.track.intersect(int_bed, u=True):
            if f.chrom == chrom:
                if f.start <= end and f.end >= start:
                    if f.start < start:
                        f.start = start
                    if f.end > end:
                        f.end = end
                # in a wig file, 4th column is score
                profile[f.start - start:f.end - start] = float(f.name)
        
        profile = profile * scalefactor
        return profile

    def binned_stats(self, in_fname, nbins, split=False, **args):
        """
        Yields a binned statistic applied to the track values for
        every feature in in_fname.

        Parameters
        ----------
        in_fname : str
            BED file

        nbins : int
            number of bins

        statistic : str, optional
            Default is "mean", other options are "min", "max" and "std"
        """
        
        in_track = pybedtools.BedTool(in_fname)
       
        statistic = args.get("statistic", "mean")
        if statistic in ["min", "max", "std"]:
            statistic = eval(statistic)

        order = {}
        regions = []
        lens = []
        for i, f in enumerate(in_track):
            region = "{}:{}-{}".format(f.chrom, f.start, f.end)
            regions.append([f.chrom, f.start, f.end])
            order[region] = i
            lens.append(f.end - f.start)
        max_len = max(lens)

        profile = np.zeros((len(regions), max_len)) 
        for f in self.track.intersect(in_track, wo=True):
            start, end = [int(x) for x in f.fields[5:7]]
            region = "{}:{}-{}".format(*f.fields[4:7])
            pos = order[region]
           
            f_start, f_end = int(f[1]), int(f[2])

            if f_start < start:
                f_start = start
            if f_end > end:
                f_end = end
            
            profile[pos][f_start - start: f_end - start] = float(f[3])
        
        for l,region,row in zip(lens, regions, profile):
            h,_,_ = binned_statistic(np.arange(l), row, bins=nbins, statistic=statistic)
            yield region + list(h)
         

class BigWigTrack(Track):
    _filetypes = ["bw", "bigWig"]
    
    def __init__(self, fname, **kwargs):
        if fname.split(".")[-1] in self._filetypes:
            self.track = pyBigWig.open(fname)
            self.ftype = "bw"
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def get_profile(self, interval, **kwargs):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        scalefactor : float, optional
            Scale profile by this factor, default 1.0
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array
        """
        scalefactor = kwargs.get("scalefactor", 1.0)
        
        chrom, start, end = self._get_interval(interval)
        profile = np.array(self.track.values(chrom, start, end))
        profile = profile * scalefactor
        return profile

    def binned_stats(self, in_fname, nbins, split=False, **args):
        """
        Yields a binned statistic applied to the track values for
        every feature in in_fname.

        Parameters
        ----------
        in_fname : str
            BED file

        nbins : int
            number of bins

        statistic : str, optional
            Default is "mean", other options are "min", "max" and "std"
        """
 
        statistic = args.get("statistic", "mean")
        use_strand = args.get("use_strand", False)
        in_track = SimpleBed(in_fname)
        for f in in_track:
            try: 
                vals = self.track.stats(f.chrom, f.start, f.end, 
                        type=statistic, nBins=nbins)
                vals = np.array(vals, dtype="float")
                vals = np.nan_to_num(vals)
                if use_strand and f.strand == "-":
                    vals = vals[::-1]
                yield [f.chrom, f.start, f.end] + list(vals)
            except:
                yield [f.chrom, f.start, f.end] + [0.0] * nbins

class TabixTrack(Track):
    _filetypes = ["bg.gz", "wig.gz", "bed.gz", 
                    "bedGraph.gz", "bigWig.gz", "bdg.gz"]

    def __init__(self, fname, **kwargs):
        if fname.endswith("gz"):
            if not os.path.exists(fname + ".tbi"):
                raise ValueError("Can't find tabix index for {}".format(fname))
            for ftype in self._filetypes:
                if fname.endswith(ftype):
                    self.tabix_track = pysam.Tabixfile(fname)
                    self.ftype = "tabix"
                    return
            raise ValueError("Can't guess format of {}".format(fname))
        else:
            raise ValueError("Can only process bgzipped files.")

    def get_profile(self, interval, **kwargs):
        """
        Return summary profile in a given window
        
        Parameters
        ----------
        interval : list, tuple or str
            If interval is a list or tuple, it should contain chromosome (str), 
            start (int), end (int). If it is a string, it should be of the
            format chrom:start-end
        
        scalefactor : float, optional
            Scale profile by this factor, default 1.0
        
        Returns
        -------
        numpy array
            A summarized profile as a numpy array
        """
        scalefactor = kwargs.get("scalefactor", 1.0)
 
        chrom, start, end = self._get_interval(interval) 

        profile = np.zeros(end - start)
        profile.fill(np.nan)
        for f in self.tabix_track.fetch(chrom, start, end):
            f = f.split()
            fstart = int(f[1])
            fend = int(f[2])
            if fstart < start:
                fstart = start
            if fend > end:
                fend = end
            profile[fstart - start: fend - end] = float(f[3])
        
        profile = profile * scalefactor 
        return profile

    def binned_stats(self, in_fname, nbins, split=False, **args):
        """
        Yields a binned statistic applied to the track values for
        every feature in in_fname.

        Parameters
        ----------
        in_fname : str
            BED file

        nbins : int
            number of bins

        statistic : str, optional
            Default is "mean", other options are "min", "max" and "std"
        """
        
        statistic = args.get("statistic", "mean")
        in_track = SimpleBed(in_fname)
       
        if statistic in ["min", "max", "std"]:
            statistic = eval(statistic)

        for r in in_track: 
            profile = np.zeros(r.end - r.start)
            for f in self.tabix_track.fetch(r.chrom, r.start, r.end):
                f = f.split()
                start = int(f[1])
                end = int(f[2])
                if start < r.start:
                    start = r.start
                if end > r.end:
                    end = r.end
                profile[start - r.start: end - r.end] = float(f[3])
            h,_,_ = binned_statistic(
                    np.arange(r.end - r.start), 
                    profile, 
                    bins=nbins, 
                    statistic=statistic)
            yield [f[0], r.start, r.end] + list(h)
