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

OTU Naming Convention
-----
The output taxonomy file names lineages according to the following naming convention. It is somewhat involved, and requires some explanation.

After a clade is defined, to name that clade tree2tax searches up the tree to find the closest named ancestral node (unless the node itself is already named). Some examples:
```
o__gHalococcus
```
An (approximately) order level grouping, where lowest ancestral named node is the genus Halococcus.

```
s__gHalobacteriaceae.1
```
An (approximately) species level grouping, where the lowest ancestral node is the genus Halobacteriaceae. The `.1` at the end indicates that >1 species level grouping has this same name, but the grouping above has the largest amount of tips included.

```
c__cHalobacteria.oHalobacteriales
```
An approximately class level grouping where the ancestral named node has been annotated as both the class Halobacteria and the order Halobacteriales

```
p__kArchaea.pEuryarchaeota|cMethanomicrobia
```
An approximately phylum level grouping where the ancestral named node is Methanomicrobia. The `kArchaea.pEuryarchaeota` is included because the higher level parent kingdom grouping is `k__Root`, so the Archaea and Euryarchaeota labels would otherwise be missing.
