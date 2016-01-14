Usage
=====

.. _quick-example:

fluff heatmap
---------------------------

::

    fluff_heatmap.py

Options
~~~~~~~~~~~~~~~~

-  ``-f``
    BED file containing features

-  ``-d``
    data files (reads in BAM or BED format)

-  ``-o`` FILE            output file (type determined by extension)

-  ``-p`` PICK          pick specific data files to use for clustering

-  ``-C`` METHOD        kmeans, hierarchical or none

-  ``-k`` INT           number of clusters

-  ``-m``               merge mirrored clusters (only with kmeans and without -g option)

-  ``-c`` NAME(S)       color(s) (name, colorbrewer profile or hex code)

-  ``-B`` NAME(S)       background color(s) (name, colorbrewer profile or hex code)

-  ``-e`` INT           extend (in bp. Default: 5000)

-  ``-b`` INT           bin size (default 100)

-  ``-s`` SCALE         scale (absolute or percentage)

-  ``-F`` FRAGMENTSIZE  Fragment length (default: read length)

-  ``-r``               use RPKM instead of read counts

-  ``-D``               keep duplicate reads (removed by default)

-  ``-R``               keep repeats (removed by default, bwa only)

-  ``-P`` INT           number of CPUs (default: 4)

-  ``-M`` METHOD        Euclidean or Pearson (default: Euclidean)

-  ``-g``               Identify dynamics by extending features 1kb up/down stream(just for clustering), cluster as 1 bin, diplay as original number of bins and with the default extend values


fluff bandplot
-------------------------

::

    fluff_bandplot.py

fluff profile
-------------------------

::

    fluff_profile.py

