#!/usr/bin/env python2.7

import logging
import os
import sys
import argparse
import subprocess
from tempdir import TempDir

parser = argparse.ArgumentParser(description='''--- autotax %s --- a pipeline creating a new taxonomy file with tree2tax suitable for use with taxtastic''' % tree2tax.__version__)
parser.add_argument('-t', '--tree', help='annotated newick format tree file to partition', required=True)
parser.add_argument('-o', '--output_directory', help='output directory for generated files', required=True)
parser.add_argument('--thresholds', nargs=7, help='tree distance thresholds to use for partitioning (one each for kingdom, phylum, class, order, family, genus, species)', type=float, default=[1.4,0.82,0.42,0.27,0.15,0.12,0.08])
parser.add_argument('--debug', help='output debug information', action="store_true")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
    

# Change directory to output directory
os.chdir(args.output_directory)

# Write thresholds file
with open('thresholds','w') as thresholds_fh:
    for thresh in args.thresholds:
        thresholds_fh.write(str(thresh))
        thresholds_fh.write("\n")
        
# cat thresholds |parallel ~/git/tree2tax/bin/tree2tax -t /srv/db/gg/2013_08/gg_13_8_otus/trees/97_otus.tree -d {} '|' sort -rn '|' cut -f2 '|' sed "'s/; /_/g'" '>' {}.taxonomies
cmd = "cat thresholds |parallel tree2tax -t %s -d {} '|' sort -rn '>' {}.two_column" % args.tree
logging.info("Running tree2tax")
logging.debug("Running cmd: %s" % cmd)
#subprocess.check_call(cmd, shell=True)

cmd = "cat thresholds |parallel cut -f2 {}.two_column '|' sed "'s/.__//g'" '|' sed \"'s/; /./g'\" '>' {}.taxonomies"
logging.debug("Running cmd: %s" % cmd)
subprocess.check_call(['bash','-c',cmd])

eg_taxonomy_file = "%s.two_column" % args.thresholds[0]
cmd = "cat %s |wc -l" % eg_taxonomy_file
logging.debug("Running cmd: %s" % cmd)
num_taxons = int(subprocess.check_output(cmd, shell=True))
logging.info("Found %i tips in the tree" % num_taxons)

one_letter_abbreviations = 'k p c o f g s'.split()
for i, abbrev in enumerate(one_letter_abbreviations):
    thresh = args.thresholds[i]
    # paste <(yes k__ |head -n %s) %s.taxonomies |sed 's/\t//' >k
    # paste <(yes p__ |head -n %s) %s.taxonomies |sed 's/\t//' >p
    # paste <(yes c__ |head -n %s) %s.taxonomies |sed 's/\t//' >c
    # paste <(yes o__ |head -n %s) %s.taxonomies |sed 's/\t//' >o
    # paste <(yes f__ |head -n %s) %s.taxonomies |sed 's/\t//' >f
    # paste <(yes g__ |head -n %s) %s.taxonomies |sed 's/\t//' >g
    # paste <(yes s__ |head -n %s) %s.taxonomies |sed 's/\t//' >s
    cmd = "paste <(yes %s__ |head -n %s) %s.taxonomies |sed 's/\t//' >%s" %(abbrev, num_taxons, thresh, abbrev)
    logging.debug("Running cmd: %s" % cmd)
    subprocess.check_call(['/bin/bash','-c',cmd])
    
# Create gg taxonomy file
cmd = "paste %s |sed 's/\t/; /g' |paste <(cut -f1 %s) - >taxonomy.gg.csv" % (' '.join(one_letter_abbreviations),
                                                                            eg_taxonomy_file)
logging.debug("Running cmd: %s" % cmd)
subprocess.check_call(['bash','-c',cmd])

# Create the taxtastic compatible files
cmd = "get_tax_n_seq2.py --taxonomy taxonomy.gg.csv --output_seqinfo seqinfo.taxtastic.csv --output_taxonomy taxonomy.taxtastic.csv"
logging.debug("Running cmd: %s" % cmd)
subprocess.check_call(['bash','-c',cmd])
    
    


# 
# num_taxons = cat thesholds[0].taxonomies |wc -l
# 
# paste <(yes k__ |head -n %s) %s.taxonomies |sed 's/\t//' >k
# paste <(yes p__ |head -n %s) %s.taxonomies |sed 's/\t//' >p
# paste <(yes c__ |head -n %s) %s.taxonomies |sed 's/\t//' >c
# paste <(yes o__ |head -n %s) %s.taxonomies |sed 's/\t//' >o
# paste <(yes f__ |head -n %s) %s.taxonomies |sed 's/\t//' >f
# paste <(yes g__ |head -n %s) %s.taxonomies |sed 's/\t//' >g
# paste <(yes s__ |head -n %s) %s.taxonomies |sed 's/\t//' >s
# 
# paste k p c o f g s |sed 's/\t/; /g' |paste <(cut -f1 thresholds[0].taxonomies) - >taxonomy.gg.csv
# 
# ~/git/Getaxnseq/bin/get_tax_n_seq2.py --taxonomy taxonomy.gg.csv --output_seqinfo seqinfo.taxtastic.csv --output_taxonomy taxonomy.taxtastic.csv
# logging.info 
# 
# logging.info("Reading tree file..")
# tree = TreeNode.read(args.tree)
# logging.info("Read in tree with %s tips" % tree.count(tips=True))
# 
# logging.info("Clustering..")
# clusters = Tree2Tax().named_clusters(tree, args.threshold)
# 
# logging.info("Found %i clusters" % len(clusters))
# # Output clade info for each tip
# clade_number = 1
# num_singleton_clades = 0
# for named_clade in clusters:
#     for tip in named_clade.tips:
#         print("\t".join([tip.name, named_clade.name()]))
#     clade_number += 1
#     if len(named_clade.tips) == 1: num_singleton_clades += 1
#     
# logging.info("Of these clusters, %s contained only a single sequence" % num_singleton_clades)
