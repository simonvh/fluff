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

class BamTrack(Track):
    _filetypes = ["bam"] 
    
    def __init__(self, fname):
        """
        Track interface to a BAM file

        Parameters
        ----------
        fname: str
            filename of BAM file

        """
        if fname.endswith("bam"):
            self.track = pysam.AlignmentFile(fname, "rb")
            self.ftype = "bam"
            self.chroms = self.track.references
        else:
            raise ValueError("filetype of {} is not supported".format(fname))

    def count(self, rmdup=False, rmrepeats=False):
        """
        Count total number of reads in file

        Parameters
        ----------
        rmdup : bool, optional
            Count duplicate reads only once if True, default False


        rmrepeats : bool, optional
            Count reads with mapping quality 0 (multi-mapped reads) only once if 
            True, default False

        Returns
        -------
        int
            Number of reads
        """
        
        if (not rmdup and not rmrepeats):
            try:
                return self.track.mapped
            except:
                pass

        c = 0
        for read in self.track:
            # duplicates
            if (not rmdup or not read.flag & 0x0400):
                # multi-mappers / mapping quality 0
                if (not rmrepeats) or read.mapq > 0:
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

    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        """ 
        Yields the number of reads for each feature in track

        Parameters
        ----------
        track : Object of <TODO: what exactly do we expect?>
            Track with features

        rmdup : bool, optional
            Count duplicate reads only once if True, default False

        rmrepeats : bool, optional
            Count reads with mapping quality 0 (multi-mapped reads) only once if 
            True, default False

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
     
    def fetch_reads(self, interval, strand=None, rmdup=False, rmrepeats=False):
        warn("fetch_reads is deprecated, please use fetch", DeprecationWarning)    

    def fetch(self, interval, strand=None, rmdup=False, rmrepeats=False):
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
        
        rmdup : bool, optional
            Don't return duplicate reads if True, default False

        rmrepeats : bool, optional
            Don't return reads with mapping quality 0 (multi-mapped reads) if 
            True, default False

        Yields
        _______
        AlignedSegment
            Yields pysam AlignedSegment objects.
        """

        chrom, start, end = self._get_interval(interval)

        if chrom in self.chroms:
            for read in self.track.fetch(chrom, start, end):
                # duplicate reads
                if rmdup and (read.flag & 1024):
                    continue
                # multimappers / low mapping quality 
                if rmrepeats and read.mapping_quality < 10:
                    continue
                if strand:
                    if strand == "+" and read.is_reverse:
                        continue
                    elif strand == "-" and not read.is_reverse:
                        continue
                yield read

    def close(self):
        self.track.close()

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
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
        _______
        numpy array
            A summarized profile as a numpy array

        """
        chrom, start, end = self._get_interval(interval)
        profile = np.zeros(end - start, dtype="f")
        profile.fill(np.nan)

        strand = {True: "-", False: "+"}

        for read in self.fetch(interval, rmdup=rmdup, rmrepeats=rmrepeats):
            iv = HTSeq.GenomicInterval(
                    chrom, 
                    read.reference_start, 
                    read.reference_end, strand[read.is_reverse]
                    )
            if fragmentsize:
                iv.length = fragmentsize
            region = profile[iv.start - start:iv.end - start]
            region[np.isnan(region)] = 0
            region += 1

        return profile

class BedTrack(Track):
    _filetypes = ["bed"]
    
    def __init__(self, fname):
       self.track = pybedtools.BedTool(fname)
       self.ftype = "bed"

    def count(self, rmdup=False, rmrepeats=False):
        if rmrepeats:
            sys.stderr.write("Warning: rmrepeats has no result on BED files!")
        if rmdup:
            sys.stderr.write("Warning: rmdup has no result on BED files! (yet...;))")

    def read_length(self):
        if self.ftype == "bed":
            for read in self.track:
                return read.end - read.start

    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        """ Generator 
        """
        if rmrepeats:
            sys.stderr.write("Warning: rmrepeats has no result on BED files!")
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
        raise NotImplementedError
        chrom, start, end, strand = window
        if strand == None:
            strand = "."
        intervals = []
        if strand == ".":
            feature = pybedtools.BedTool("{0} {1} {2}".format(chrom, start, end), from_string=True)
            s = False
        else:
            feature = pybedtools.BedTool("{0} {1} {2} 0 0 {3}".format(chrom, start, end, strand), from_string=True)
            s = True
        for read in self.track.intersect(feature, u=True, stream=True, s=s):
            intervals.append(HTSeq.GenomicInterval(chrom, read.start, read.end, str(read.strand)))
        return intervals

    def close(self):
        pass

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
        raise NotImplementedError

class WigTrack(Track):
    _filetypes = ["bg", "wig"]
    
    def __init__(self, fname):
        if fname.endswith("bg") or fname.endswith("wig"):
            self.htseq_track = HTSeq.WiggleReader(fname, verbose=True)
            self.ftype = "wig"
        else:
            raise NotImplementedError

    def __getitem__(self, window):
        """ 
        Retrieve all reads within a given window
        Arguments: window - a list or tuple containing chromosome, start, end and strand
        Returns a list of GenomicIntervals
        """
    
        raise NotImplementedError

    def count(self, rmdup=False, rmrepeats=False):
        raise NotImplementedError

    def read_length(self):
       raise NotImplementedError
    
    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        """ Generator 
        """
        raise NotImplementedError

    def fetch_reads(self, interval, rmdup=False, rmrepeats=False):
        """ Generator 
        """
        raise NotImplementedError

    def close(self):
        pass

    def get_profile(self, interval, fragmentsize=200, rmdup=False, rmrepeats=False):
        chrom, start, end = interval
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

class TabixTrack(Track):
    def __init__(self, fname):
        if fname.endswith("gz") and os.path.exists(fname + ".tbi"):
            self.tabix_track = pysam.Tabixfile(fname)
            self.ftype = "tabix"
        else:
            raise NotImplementedError

    def __getitem__(self, window):
        """ 
        Retrieve all reads within a given window
        Arguments: window - a list or tuple containing chromosome, start, end and strand
        Returns a list of GenomicIntervals
        """
        
        raise NotImplementedError

    def count(self, rmdup=False, rmrepeats=False):
        raise NotImplementedError

    def read_length(self):
        raise NotImplementedError
    
    def fetch_to_counts(self, track, rmdup=False, rmrepeats=False):
        raise NotImplementedError

    def fetch_reads(self, interval, rmdup=False, rmrepeats=False):
        raise NotImplementedError
    
    def close(self):
        pass

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

if __name__ == "__main__":
    bam = "/home/simon/git/fluff/tests/data/profile.bam"
    
    b = BamTrack(bam)
    print Track._registry
    print Track.filetypes()

    #b[["scaffold_1",44749422, 44750067, "+"]]
    print b.read_length()
    print len([i for i in b.fetch(["scaffold_1",44749422, 44750067])])
    print len([i for i in b.fetch("scaffold_1:44749422-44750067")])
    print len([i for i in b.fetch("scaffold_1:44749422-44750067", strand="+")])
    print len([i for i in b.fetch("scaffold_1:44749422-44750067", strand="-")])

    print b.get_profile(["scaffold_1", 44749422, 44750067], fragmentsize=None)
