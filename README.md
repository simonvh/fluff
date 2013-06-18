Disclaimer
==========

Fluff is still under heavy development. No guarantees are given ;). However, if you find it useful and encounter problems, let me know and I'll try to help. The gene annotation track that fluff_profile.py produces is still pretty ugly too.

fluff
=====

![fluff examples](https://raw.github.com/simonvh/fluff/master/examples/test.png) 

Fluff is a Python package containing several scripts with the aim to produce pretty, publication-quality figures for next-generation sequencing experiments. The goal is to provide default settings which produce figures that are ready-to-use (resolution, font size, colors etc.)

It currently contains three scripts:
* fluff_heatmap.py
* fluff_bandplot.py
* fluff_profile.py

Plotting is handled by the excellent [matplotlib](http://matplotlib.sourceforge.net/) library, several image formats are supported (SVG, Postscript, PDF, PNG).

Fluff makes heavy use of [pysam](http://code.google.com/p/pysam/), [pybedtools](http://packages.python.org/pybedtools/), [HTseq](http://www-huber.embl.de/users/anders/HTSeq/) and indexed [BAM files](http://samtools.sourceforge.net/) for speed and ease of use. Due to the indexing, the fluff scripts can be quick, even when working with large files. Currently, [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) files are also supported, but performance will suffer.

Prerequisites
-------------
* tabix - http://sourceforge.net/projects/samtools/files/
* pysam - http://code.google.com/p/pysam/
* pybedtools - http://packages.python.org/pybedtools/
* HTSeq - http://www-huber.embl.de/users/anders/HTSeq/
* matplotlib - http://matplotlib.sourceforge.net/
* numpy - http://numpy.scipy.org/
* scipy - http://www.scipy.org/
* colorbrewer - http://pypi.python.org/pypi/colorbrewer/
* Pycluster - https://pypi.python.org/pypi/Pycluster
* pp - http://www.parallelpython.com/ (optional)

I'd recommend installing matplotlib, numpy and scipy using your preferred package manager. Ubuntu, Debian, Fedora, Gentoo, etc. all have packages providing these libraries.

Colors
------
One important feature of fluff is the use of color. There are three ways to specify colors: by name, palette or hex code. You can specify as few or as many colors as you want, seperated by commas. If fluff runs out of colors it will start at the beginning.

1. **Color names.**
Fluff knowns nine color names: red, blue, green, purple, orange, yellow, pink and grey. These are the colors from the Set1 palette (see below).

2. **Palettes.**
The palettes with all the [ColorBrewer](http://colorbrewer2.org/) colors can be specified (thanks to the Python [colorbrewer](http://pypi.python.org/pypi/colorbrewer/) package supplied by Michael Hoffman). If you specify the palette name only, the maximum number of colors will be used. Alternatively specify the number of colors with a colon. For example: *Set1:5,Set2:4* is a valid argument.

3. **Hex code.**
Finally, colors can be specified using a hex code. For instance, *ff0000,00FF00,00000FF* will give you red, blue and green. Indeed, case-insensitive.

All different color specification can be mixed and matched for extra fun.

General options
---------------

There are three options shared amongst all scripts:
* `-r` Enable this option to use RPKM values (read per kb per million reads) instead of read counts.
* `-R` When bam files are produced by bwa, reads mapping to non-duplicate regions of the genome are marked. By default fluff removes all these reads. When this option is used these reads will be kept.
* `-D` By default fluff removes all duplicate reads (when marked in the bam file). Enable the `-D` option to keep duplicates.

The scripts
===========

fluff_heatmap.py
----------------
Produce a heatmap like [this example](add link). Features can be shown "as is", preserving the order in the input file, or can be clustered using hierarchical or k-means clustering. This scripts creates two output files. One is an image that contains the heatmap (png by default), the other one contains the features in the same order as the heatmap including cluster number. This file can be used as input for `fluff_bandplot.py`. 
It is worth mentioning that the heatmap can be save in a vector-based format, however, the file quickly becomes prohibitively large and unwieldy. Be warned ;).

```
Usage: fluff_heatmap.py -f <bedfile> -d <file1>[,<file2>,...] -o <out> [options]

Options:
  --version     show program's version number and exit
  -h, --help    show this help message and exit
  -f FILE       BED file containing features
  -d FILE(S)    data files (reads in BAM or BED format)
  -o FILE       output file (type determined by extension)

  Optional:
    -C METHOD   kmeans, hierarchical or none
    -k INT      number of clusters
    -m          merge mirrored clusters (only with kmeans)
    -c NAME(S)  color(s) (name, colorbrewer profile or hex code)
    -B NAME(S)  background color(s) (name, colorbrewer profile or hex code)
    -e INT      extend (in bp)
    -b INT      bin size (default 100)
    -s SCALE    scale (absolute or percentage)
    -r          use RPKM instead of read counts
    -D          keep duplicate reads (removed by default)
    -R          keep repeats (removed by default, bwa only)
```

### Clustering ###
By default, fluff_heatmap.py will preserve the order of the features in the input BED file. This is equivalent to specifying `-C none`. Alternatively, one of two basic clustering methods can be specified using the `-C` parameter: `hierarchical` and `kmeans`. If `kmeans` is selected the number of clusters (`-k`) is mandatory. The profiles will be clustered using the Euclidian distance metric, based on normalized data (although the heatmap still shows the raw data). If the regions in the input BED file are not strand-specific (peaks, for instance), different clusters might describe the same strand-specific profile in two different orientations. In this case you can select the `-m` parameter to merge clusters which are similar, but mirrored around the center coordinate. Similarity is based on the chi-squared p-value of the mean profile per cluster. **All** tracks have to pass this similarity test for two clusters to be merged! Merging follows a iterative, greedy procedure where the two most similar clusters are merged after which the process is repeated.

### Region, bin size and scaling ###
The *center* of all features will be selected from the input BED file. A total region of 10kb will be centered at this coordinate (5kb up- and downstream). To change this, specify the `-e` parameter that controls the extension up- and downstream. For example,  `-e 10000` will select regions of 20kb.

The bin size is 100 bp, this can be controlled with the `-b` parameter.

The color scale of the heatmap profile can be adjusted with the `-s` parameter. It is difficult to set a generally applicable threshold, as this value might depend on what you'd like to highlight or analyze in your data. The `-s` argument can either be absolute (`-s 15` means a bin that contains 15 or more reads will have the darkest color) or a percentage (`-s 90%` sets the threshold at the 90th percentile of the values in all bins over all tracks and features).

### Background color ###
The default background color is white. This can be changed with the `-B` parameter, which works the same as the color specification.

fluff_bandplot.py
-----------------
While heatmaps can be very informative, sometimes you want to show the average profile. However, by just plotting the mean, you lose quite some information. This is my attempt to combine a boxplot with an average profile. The mean enrichment is shown using a black line. The 50th and 90th percentile are also visualized using a dark and light color respectively.

```
Usage: fluff_bandplot.py -c <bedfile> -d <file1>[,<file2>,...] -o <out> [options]

Options:
  --version     show program's version number and exit
  -h, --help    show this help message and exit
  -i FILE       BED file with cluster in 4th column
  -d FILE(S)    data files (reads in BAM or BED format)
  -o FILE       output file (type determined by extension)

  Optional:
    -c NAME(S)  color(s) (name, colorbrewer profile or hex code)
    -s GROUPS   scale groups
    -p INT,INT  Range of percentiles (default 50,90)
    -r          use RPKM instead of read counts
    -D          keep duplicate reads (removed by default)
    -R          keep repeats (removed by default, bwa only)
```

The input file should contain a cluster identifier in the fourth column. Any string or number can be used; all features with the same identifier will be grouped. The size of the regions is determined by the coordinates in the BED file, ie if you want to plot a region of 10kb all the entries in the BED file should span a region of 10kb. The `*_clusters.bed` file from `fluff_heatmap.py` can directly be used in this script.

### Scales and percentiles ###
To compare different datasets (for instance same antibody in a timeline or different conditions) the maximum of the Y-axis can be linked between tracks. As an example, `-s 1:3,4,5:7` will ensure that the first three rows share the same scale, the fourth track will be scaled to the maximum of only this single track and the last three rows will once again share a similar scale.

The `-p` parameter controls the percentile range. By default, the profile that contains 50% of the data will be shaded in a dark color and the profile that contains 90% of the data will be shaded in a lighter color. 

fluff_profile.py
----------------
Produces output like a Genome Browser screenshot. Currently only 1) profiles based on reads in BAM or BED format and 2) gene annotation (in [BED12](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)) can be visualized.
```
Usage: fluff_profile.py -i <loc1>[,<loc2>,...] -d <file1>[,<file2>,...] -o <out> [options]

Options:
  --version           show program's version number and exit
  -h, --help          show this help message and exit
  -i INTERVAL(S)      one or more genomic intervals (chrom:start-end)
  -d FILE(S)          data files (reads in BAM or BED format)
  -o FILE             output file name (type determined by extension)

  Optional:
    -a FILE           annotation in BED12 format
    -c NAME(S)        color(s) (name, colorbrewer profile or hex code)
    -t GROUPS         track groups
    -s GROUPS         scale groups
    -S SCALE          scale: 'auto' (default), 'off' or int for each track
    -b BACKGROUND     background color: white | color | stripes
    -f FRAGMENTSIZE   fragment length (default: 200)
```

###Scale and track groups###
Just as with `fluff_bandplot.py`, identical Y-axis scale can be used to compare tracks. In addition, different datasets can be overlayed on the same track wit the `-t` option. Use the `-S` option to manually set a the Y-axis upper limit, or to turn the label off (with the `off` argument).

###Colors###
The background color can be one of three different choices. The default is `white`. Specifying `color` will set the background of each track to a lighter version of the track color. The last option is `stripes` which will show an alternating pattern of white and light grey bands.

###Fragment length###
By default, all reads will be extended to 200bp before creating the profiles. Set an alternative fragment size with the `-f` option. Paired-end profiles are not yet supported, paired reads will be independently processed.





