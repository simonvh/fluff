Usage
=====

.. _quick-example:

fluff heatmap
-------------
::

    fluff heatmap -f <BED> -d <BAM> <BAM> -o <NAME>

Options
~~~~~~~

Required arguments:
~~~~~~~~~~~~~~~~~~~

-  ``-f`` FILE
BED file containing features

-  ``-d`` [FILE [FILE ...]]
data files (reads in BAM or BED format)

-  ``-o`` name
output file (type determined by extension)

Optional arguments:
~~~~~~~~~~~~~~~~~~~

-  ``-h``
show this help message and exit

-  ``-C`` METHOD
kmeans, hierarchical or none

-  ``-k`` INT
number of clusters

-  ``-M`` METHOD
Euclidean or Pearson (default: Euclidean)

-  ``-g``
Identify dynamics

-  ``-p`` PICK
pick specific data files to use for clustering

-  ``-e`` INT
extend (in bp. Default: 5000)

-  ``-b`` INT
bin size (default 100)

-  ``-F`` FRAGMENTSIZE
Fragment length (default: read length)

-  ``-r``
use RPKM instead of read counts

-  ``-D``
keep duplicate reads (removed by default)

-  ``-R``
keep repeats (removed by default, bwa only)

-  ``-m``
merge mirrored clusters (only with kmeans and without -g option)

-  ``-s`` SCALE
scale (absolute or percentage)

-  ``-c`` NAME(S)
color(s) (name, colorbrewer profile or hex code)

-  ``-B`` NAME(S)
background color(s) (name, colorbrewer profile or hex code)

-  ``-P`` INT
number of CPUs (default: 4)



fluff bandplot -f <BED> -d <BAM> <BAM> -o <NAME>
--------------

::

    fluff bandplot



Options
~~~~~~~

Required arguments:
~~~~~~~~~~~~~~~~~~~

-  ``-f`` FILE
BED file with cluster in 5th column
-  ``-d`` [FILE [FILE ...]]
data files (reads in BAM or BED format)
-  ``-counts`` FILE
Read Counts
-  ``-o`` name
output file (type determined by extension)

Optional arguments:
~~~~~~~~~~~~~~~~~~~

-  ``-h``
show this help message and exit
-  ``-S``
create summary graphs
-  ``-b`` INT
number of bins
-  ``-F`` FRAGMENTSIZE
fragment length (default: read length)
-  ``-r``
use RPKM instead of read counts
-  ``-D``
keep duplicate reads (removed by default)
-  ``-R``
keep repeats (removed by default, bwa only)
-  ``-s`` GROUPS
scale groups
-  ``-p`` INT,INT
range of percentiles (default 50,90)
-  ``-P`` INT
Percentile at which to extract score. Value should be in range [0,100] (default 90)
-  ``-c`` NAME(S)
color(s) (name, colorbrewer profile or hex code)





fluff profile
-------------

::

    fluff profile -i <GENOMIC LOCATION> -d <BAM> <BAM> -o <NAME>



Options
~~~~~~~

Required arguments:
~~~~~~~~~~~~~~~~~~~

-  ``i`` INTERVAL(S)
one or more genomic intervals (chrom:start-end)
-  ``d`` [FILE [FILE ...]]
data files (reads in BAM or BED format)
-  ``o`` name
output file (type determined by extension)

Optional arguments:
~~~~~~~~~~~~~~~~~~~

-  ``h``
show this help message and exit
-  ``a`` FILE
annotation in BED12 format
-  ``t`` GROUPS
track groups
-  ``s`` GROUPS
scale groups
-  ``S`` SCALE
scale: 'auto' (default), 'off' or int for each track
-  ``f`` FRAGMENTSIZE
fragment length (default: 200)
-  ``D``
keep duplicate reads (removed by default)
-  ``R``
keep repeats (removed by default, bwa only)
-  ``r``
reverse
-  ``c`` NAME(S)
color(s) (name, colorbrewer profile or hex code)
-  ``b`` BACKGROUND
background color: white | color | stripes
