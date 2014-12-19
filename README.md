tree2tax
========

Automatic taxonomy through consistent application of tree-based thresholding.

Phylogenetic tree insertion methods (such as [GraftM](https://github.com/geronimp/graftM) and 
[pplacer](http://matsen.fhcrc.org/pplacer/)) can place sequences into parts of the tree where taxonomists have not
yet come to agreed upon names. Here, taxonomy is automatically assigned based on tree distance thresholds, and then
nodes corresponding to these thresholds are given names taken from the usual taxonomy so they are meaningful to humans
and their formidable pattern finding skillz.

Tree2tax is inspired by, but not related to, the [tax2tree](https://github.com/biocore/tax2tree) software, which doesn't suggest new taxonomic clades.

Needs:
* scikit-bio
* an taxonomically annotated tree, tested on greengenes 99_otus.tree
* possibly other things

How:
```sh
git clone this repo
PYTHONPATH=$PYTHONPATH:/path/to/tree2tax /path/to/tree2tax/bin/tree2tax -h
```
