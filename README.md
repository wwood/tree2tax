tree2tax
========

[![Build Status](https://img.shields.io/travis/wwood/tree2tax.svg)](https://travis-ci.org/wwood/tree2tax)

Automatic taxonomy through consistent application of tree-based thresholding.

Phylogenetic tree insertion methods (such as [GraftM](https://github.com/geronimp/graftM) and 
[pplacer](http://matsen.fhcrc.org/pplacer/)) can place sequences into parts of the tree where taxonomists have not
yet come to agreed upon names. Here, 7 different granularities of operational taxonomic units (OTUs) are automatically assigned based on tree distance thresholds, and then
nodes corresponding to these thresholds are given names taken from the usual taxonomy so they are meaningful to humans
and their formidable pattern finding skillz.

Tree2tax is inspired by, but not related to, the [tax2tree](https://github.com/biocore/tax2tree) software, which doesn't suggest new taxonomic clades.

Pre-requisites:
* A taxonomically annotated tree, tested on greengenes 99_otus.tree
* some standard Unix tools

How:
```sh
pip install tree2tax
tree2tax -h
```

Auto-taxonomy Naming Convention
-----
The output OTU file names lineages according to the following naming convention. It is somewhat involved, and requires some explanation.

After a clade is defined, to name that clade tree2tax searches up the tree to find the closest named ancestral node (unless the node itself is already named). Some examples:
```
O__gHalococcus
```
An (approximately) order level grouping, where lowest ancestral named node is the genus Halococcus. The initial 'O' is capitalised to indicate
that this is not a true taxonomic level (a suggestion by Dr Donovan Parks (github @dparks1134)).

```
S__gHalobacteriaceae.1
```
An (approximately) species level grouping, where the lowest ancestral node is the genus Halobacteriaceae. The `.1` at the end indicates that >1 species level grouping has this same name, but the grouping above has the largest amount of tips (sequences) included: a `.2` would indicate the clade with the second largest number of tips, and so on.

```
C__cHalobacteria.oHalobacteriales
```
An approximately class level grouping where the ancestral named node has been annotated as both the class Halobacteria and the order Halobacteriales

```
P__kArchaea.pEuryarchaeota|cMethanomicrobia
```
An approximately phylum level grouping where the ancestral named node is Methanomicrobia. The `kArchaea.pEuryarchaeota` is included because the higher level parent kingdom grouping is `K__Root`, so the Archaea and Euryarchaeota labels would otherwise be missing.
