Disclaimer
==========

Fluff is still under heavy development, mainly for personal use. No guarantees are given ;). However, if you find it useful and encounter problems, let me know and I'll try to help.

fluff
=====

Fluff is a Python package containing several scripts with the aim to produce pretty, publication-quality figures for next-generation sequencing experiments. I've tried to make sure the default settings produce figures that are ready to use.

It currently contains three scripts:
* fluff_heatmap.py
* fluff_profile.py
* fluff_bandplot.py

Plotting is handled by the excellent [matplotlib](http://matplotlib.sourceforge.net/) library, several image formats are supported (SVG, Postscript, PDF, PNG).

Fluff makes heavy use of [pysam](http://code.google.com/p/pysam/), [pybedtools](http://packages.python.org/pybedtools/), [HTseq](http://www-huber.embl.de/users/anders/HTSeq/) and indexed [BAM files](http://samtools.sourceforge.net/) for speed and easy of use. Due to the indexing, the fluff scripts can be quick, even when working with large files. Currently, [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) files are also supported, but performance will suffer.

Prerequisites
-------------
* pysam - http://code.google.com/p/pysam/
* pybedtools - http://packages.python.org/pybedtools/
* HTSeq - http://www-huber.embl.de/users/anders/HTSeq/
* matplotlib - http://matplotlib.sourceforge.net/
* numpy - http://numpy.scipy.org/
* scipy - http://www.scipy.org/

I'd recommend installing matplotlib, numpy and scipy using your preferred package manager. Ubuntu, Debian, Fedora, Gentoo, etc. all have packages providing these libraries.

Colors
------
One important feature of fluff is the use of color, as this is part of making your figures look good. I use the Set1 palette from the [ColorBrewer](http://colorbrewer2.org/) colors by default, however, all colors can be specified either by name, by R,G,B values or by Hex code.

The scripts
===========

fluff_heatmap.py
----------------
Produce a heatmap like [this example](add link). Features can be shown "as is", preserving the order in the input file, or can be clustered using hierarchical or k-means clustering. 

fluff_profile.py
----------------
Produces output like a Genome Browser screenshot. Currently only BAM profiles and gene annotation (in [BED12](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)) can be visualized.

fluff_bandplot.py
-----------------
While heatmaps can be very informative, sometimes you want to show the average profile. However, by just plotting the mean, you lose quite some information. This is my attempt to combine a boxplot with an average profile. The mean enrichment is shown using a black line. The 50th and 90th percentile are also visualized using a dark and light color respectively.
