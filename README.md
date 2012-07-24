fluff
=====

Fluff is package that contains several scripts to produce pretty, publication-quality figures for next-generation sequencing experiments. I've tried to make sure the default settings produce figures that are ready to use.

It currently contains three scripts:
* fluff_heatmap.py
* fluff_profile.py
* fluff_bandplot.py

Plotting is handled by the excellent [matplotlib](http://matplotlib.sourceforge.net/) library, several image formats are supported (SVG, Postscript, PDF, PNG).

Fluff makes heavy use of [pysam](http://code.google.com/p/pysam/) and indexed [BAM files](http://samtools.sourceforge.net/) for speed and easy of use. Due to the indexing, the fluff scripts can be quick, even when working with large files.W hile [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) files can also be used, this will be much slower.

fluff_heatmap.py
----------------

fluff_profile.py
----------------

fluff_bandplot.py
