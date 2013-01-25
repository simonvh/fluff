Disclaimer
==========

Fluff is still under heavy development, mainly for personal use. No guarantees are given ;). However, if you find it useful and encounter problems, let me know and I'll try to help. The gene annotation track that fluff_profile.py produces is still pretty ugly too.

fluff
=====

![heatmap](https://raw.github.com/simonvh/fluff/master/examples/heatmap.png) ![bandplot](https://raw.github.com/simonvh/fluff/master/examples/bandplot.png) ![profile](https://raw.github.com/simonvh/fluff/master/examples/profile.png) 

Fluff is a Python package containing several scripts with the aim to produce pretty, publication-quality figures for next-generation sequencing experiments. I've tried to make sure the default settings produce figures that are ready-to-use (resolution, font size, colors etc.)

It currently contains three scripts:
* fluff_heatmap.py
* fluff_bandplot.py
* fluff_profile.py

Plotting is handled by the excellent [matplotlib](http://matplotlib.sourceforge.net/) library, several image formats are supported (SVG, Postscript, PDF, PNG).

Fluff makes heavy use of [pysam](http://code.google.com/p/pysam/), [pybedtools](http://packages.python.org/pybedtools/), [HTseq](http://www-huber.embl.de/users/anders/HTSeq/) and indexed [BAM files](http://samtools.sourceforge.net/) for speed and ease of use. Due to the indexing, the fluff scripts can be quick, even when working with large files. Currently, [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) files are also supported, but performance will suffer.

Prerequisites
-------------
* pysam - http://code.google.com/p/pysam/
* pybedtools - http://packages.python.org/pybedtools/
* HTSeq - http://www-huber.embl.de/users/anders/HTSeq/
* matplotlib - http://matplotlib.sourceforge.net/
* numpy - http://numpy.scipy.org/
* scipy - http://www.scipy.org/
* colorbrewer - http://pypi.python.org/pypi/colorbrewer/

I'd recommend installing matplotlib, numpy and scipy using your preferred package manager. Ubuntu, Debian, Fedora, Gentoo, etc. all have packages providing these libraries.

Colors
------
One important feature of fluff is the use of color, as this is part of making your figures look good. There are three ways to specify colors: by name, palette or hex code. You can specify as few or as many colors as you want, seperated by commas. If fluff runs out of colors it will start at the beginning.

1. **Color names.**
Fluff knowns nine color names: red, blue, green, purple, orange, yellow, pink and grey. These are the colors from the Set1 palette (see below).

2. **Palettes.**
The palettes with all the [ColorBrewer](http://colorbrewer2.org/) colors can be specified (thanks to the Python [colorbrewer](http://pypi.python.org/pypi/colorbrewer/) package supplied by Michael Hoffman). If you specify the palette name only, the maximum number of colors will be used. Alternatively specify the number of colors with a colon. For example: *Set1:5,Set2:4* is a valid argument.

3. **Hex code.**
Finally, colors can be specified using a hex code. For instance, *ff0000,00FF00,00000FF* will give you red, blue and green. Indeed, case-insensitive.

All different color specification can be mixed and matched for extra fun.

The scripts
===========

fluff_heatmap.py
----------------
Produce a heatmap like [this example](add link). Features can be shown "as is", preserving the order in the input file, or can be clustered using hierarchical or k-means clustering. 

fluff_profile.py
----------------
Produces output like a Genome Browser screenshot. Currently only 1) profiles based on reads in BAM or BED format and 2) gene annotation (in [BED12](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)) can be visualized.

fluff_bandplot.py
-----------------
While heatmaps can be very informative, sometimes you want to show the average profile. However, by just plotting the mean, you lose quite some information. This is my attempt to combine a boxplot with an average profile. The mean enrichment is shown using a black line. The 50th and 90th percentile are also visualized using a dark and light color respectively.
