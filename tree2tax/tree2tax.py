from skbio.tree import TreeNode
from sets import Set
import logging
import re
import IPython

class TaxonomyFunctions:
    @staticmethod
    def condense(taxonomy_string):
        '''e.g. if taxonomy is c__Halobacteria; o__Halobacteriales; f__MSP41
        then return 'cHalobacteria.oHalobacteriales.fMSP41' We cannot omit the 
        c, o or f because
        in some rare cases two taxonomic levels have the same name except for their
        level e.g. c__Gemmatimonadetes; o__Gemmatimonadetes in GreenGenes 
        2013_08'''
        splits = taxonomy_string.split('; ')
        
        # get rid of f__ prefixes etc.
        regex = re.compile(r'^.__')
        splits2 = []
        for s in splits:
            reg = regex.match(s)
            if reg:
                splits2.append(s[0]+s[reg.end():])
            else:
                logging.debug("Found unexpected form for taxonomy in %s", str(s))
                splits2.append(s)
        return '.'.join(splits2)
    
    @staticmethod
    def missing_taxonomy(tree, descendent_node, ancestral_node):
        '''given a tree, and two nodes where one is an ancestor of another,
        return a condensed taxonomy list representing the 
        taxonomic information that is contained in nodes between the descendent 
        and ancestral nodes, except for the first one encountered
        '''
        to_return = []
        
        current = descendent_node
        while current != ancestral_node:
            current = current.parent
            if current is None:
                raise Exception("descendent and ancestral nodes do not appear to be related as expected in #missing_taxonomy")
            
            if current != ancestral_node:
                tax = TaxonomyFunctions().taxonomy_from_node_name(current.name)
                if tax:
                    to_return.append(TaxonomyFunctions().condense(tax))

        # return in descending order
        return [r for r in reversed(to_return)]
            
    @staticmethod
    def taxonomy_from_node_name(node_name):
        '''return the taxonomy incorporated at a particular node, or None
        if it does not encode any taxonomy'''
        
        def isFloat(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
        
        if node_name is None:
            return None
        elif isFloat(node_name):
            # no name, just a bootstrap
            return None
        else:
            bootstrap_regex = re.compile(r'[\d\.]+:(.*)')
            reg = bootstrap_regex.match(node_name)
            if reg:
                # bootstrap in name
                return reg.groups(0)[0]
            else:
                # bootstrap not in name
                return node_name
    

class NamedCluster:
    def __init__(self, taxonomy, tips, lca_node):
        self.taxonomy = taxonomy
        self.tips = tips
        self.lca_node = lca_node
        self.cluster_number = None
        
    def name(self):
        if self.cluster_number:
            return "%s.%s" % (self.taxonomy, self.cluster_number)
        else:
            return self.taxonomy
        
    def __str__(self):
        return "%s %s: %s" % (self.taxonomy, self.cluster_number, self.tips)
    
    def condensed_name(self):
        '''e.g. if taxonomy is c__Halobacteria; o__Halobacteriales; f__MSP41
        then return 'cHalobacteria.oHalobacteriales.fMSP41' (cluster number
        is added too, if that is present. We cannot omit the c, o or f because
        in some rare cases two taxonomic levels have the same name except for their
        level e.g. c__Gemmatimonadetes; o__Gemmatimonadetes in GreenGenes 
        2013_08'''
        joined = TaxonomyFunctions.condense(self.taxonomy)
        if self.cluster_number:
            return "%s.%s" % (joined, self.cluster_number)
        else:
            return joined
        
        
class ThresholdAndClusters:
    def __init__(self, threshold, clusters):
        self.threshold = threshold
        self.clusters = clusters
        
    def tip_name_to_cluster(self, tip_name):
        key = tip_name
        try:
            # cached already?
            return self._tip_to_cluster[key]
        except AttributeError:
            # no dice. Have to create the hash
            self._tip_to_cluster = {}
            for cluster in self.clusters:
                for datip in cluster.tips:
                    if self._tip_to_cluster.has_key(datip.name):
                        logging.warn("Unexpectedly found multiple leaf nodes with the same name, undefined behaviour possibly imminent: %s" % tip_name)
                    else:
                        self._tip_to_cluster[datip.name] = cluster
        return self._tip_to_cluster[key]
        
    def tip_to_cluster(self, tip):
        '''return the NamedCluster to which the specified tip belongs'''
        return self.tip_name_to_cluster(tip.name)
        
                    
    def each_tip(self):
        '''iterate over all the tips in all clusters'''
        for cl in self.clusters:
            for tip in cl.tips:
                yield tip
            

class Tree2Tax:
    def named_clusters_for_several_thresholds(self, original_tree, thresholds):
        '''Given a list of thresholds, return a iterable of ThresholdAndClusters
        where the clustering has been done iteratively, providing a consistent
        taxonomic annotation scheme'''
        original_tree.assign_ids()
        tree = original_tree.copy() #this copy gets destructively pruned as part of the algorithm
        
        # sort from smallest to largest because smaller distance thresholds
        # need to be applied before larger thresholds
        sorted_thresholds = sorted(thresholds)
        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Found thresholds %s", str(sorted_thresholds))
        
        to_return = []

        # skbio seems to be quite slow running repeated find_by_id. Do some caching 
        # to speed things up
        id_to_node = {}
        for node in original_tree.postorder(include_self=True):
            id_to_node[node.id] = node
            
        clades_to_distances = {}
        for threshold in sorted_thresholds:
            if logging.getLogger().isEnabledFor(logging.INFO): logging.info("Clustering on threshold %s.." % str(threshold))
            
            clades_to_distances = self.destructively_cluster_tree(tree, threshold, clades_to_distances)
            
            clusters = []
            if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(tree.ascii_art())
            cluster_lcas = tree.tips()
            
            # So that several clades don't get named the same thing
            # it get a bit complex when several nodes are annotated as having the sam
            # taxonomy (it happens..., in gg at least). So group by taxonomy
            # rather than node ID
            taxonomy_to_named_nodes = {} 
            
            # when everything is in 1 cluster (doesn't happen in practice I suspect)
            # but there is a unit test..
            if len(tree) == 0:
                cl = NamedCluster('Root', original_tree.tips(), original_tree)
                cl.cluster_number = ''
                clusters = [cl]
    
            else:
                for lca in cluster_lcas:
                    orig = id_to_node[lca.id]
                    
                    # Find the closest ancestral named node (where ancestral does not include self)
                    current = orig
                    
                    #in a greengenes file, tips have names, but don't count these 
                    #because they are not taxonomy but rather prokMSA IDs
                    if current.is_tip(): current = current.parent
                    
                    # Keep proceeding up the tree until a node with taxonomy is found.
                    # If the tree has bootstraps, then current.name is a float (in string form)
                    # Ignore these floats because they aren't named taxonomy
                    while current.parent and (current.name is None or 
                         TaxonomyFunctions.taxonomy_from_node_name(current.name) is None):
                        current = current.parent
                        
                    if current.parent is None:
                        taxonomy = 'Root'
                    else:
                        taxonomy = TaxonomyFunctions.taxonomy_from_node_name(current.name)
                        
                    if orig.is_tip():
                        named_cluster = NamedCluster(taxonomy, [orig], orig)
                    else:
                        named_cluster = NamedCluster(taxonomy, list(orig.tips()), orig)
                        
                    clusters.append(named_cluster)
                    
                    try:
                        taxonomy_to_named_nodes[taxonomy].append(named_cluster)
                    except KeyError:
                        taxonomy_to_named_nodes[taxonomy] = [named_cluster]
                
                if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(taxonomy_to_named_nodes)
                # Assign cluster numbers by decreasing numbers of tips i.e. for each
                # taxonomy, the cluster with the most sequences is given the number 1,
                # second most abundant is number 2, etc.
                for taxonomy, named_clusters in taxonomy_to_named_nodes.items():
                    if len(named_clusters) > 1:
                        number = 1
                        # sort twice (stupid python). Consider number of tips
                        # to be more important than the minimum node ID. Need
                        # to sort repeatably for the purposes of testing
                        sorts1 = sorted(named_clusters, key = lambda c: min([tip.id for tip in c.tips]))
                        for clade in sorted(sorts1, reverse = True, key = lambda c: len(c.tips)):
                            clade.cluster_number = number
                            number += 1
                            
            threshold_and_clusters = ThresholdAndClusters(threshold, clusters) 
            to_return.append(threshold_and_clusters)
            
        return to_return
        
        
    def named_clusters(self, original_tree, threshold):
        '''cluster a tree given a threshold tree_distance'''
        array = self.named_clusters_for_several_thresholds(original_tree, [threshold])
        return array[0].clusters
    
    def destructively_cluster_tree(self, tree, threshold, clades_to_distances = {}):
        '''Given a tree, collapse tips using complete linkage such that 
        all sequences within the collapsed clade have at most the threshold 
        tree distance (and cluster as much as possible). Return hash of 
        tip.id => farthest distance found from that root that still obeyed
        the threshold'''
        
        tips_to_evaluate = Set(tree.tips())
        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Parsed tree with %s tips: %s" % (tree.count(tips=True), [t.name for t in tips_to_evaluate]))
        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(tree.ascii_art())
        
        while len(tips_to_evaluate) > 0:
            tip = tips_to_evaluate.pop() #raises KeyError when set is empty, jumping out of the while loop
            if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Evaluating tip %s/%s" % (tip.id, tip.name))
            
            # If this tip has only 1 sibling, then evaluate whether it and its sibling can be collapsed according to the threshold
            if tip.parent is not None: #if we aren't at the root
                bros = [n for n in tip.parent.children if n != tip and n.is_tip()]
                
                if len(bros) == 1:
                    # The sibling node is a tip
                    bro = bros[0]
                    distance_to_bro = tip.length + bro.length
                    if tip.id in clades_to_distances: distance_to_bro += clades_to_distances[tip.id]
                    if bro.id in clades_to_distances: distance_to_bro += clades_to_distances[bro.id]
                    if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("found bro %s/%s at distance %s" % (bro.id, bro.name, distance_to_bro))
                    if distance_to_bro <= threshold:
                        # cluster these two tips
                        # remove tip and bro from the tree
                        # add the parent node to the set so it can possibly can merge with further tips or parents that become tips
                        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Below threshold, removing these nodes")
                        if bro in tips_to_evaluate: tips_to_evaluate.remove(bro)
                        parental = tip.parent
                        tips_to_evaluate.add(parental)
                        if parental.remove(bro) is not True: raise Exception("Programming error")
                        if parental.remove(tip) is not True: raise Exception("Programming error")
                        
                        # Set the distance of this node to be the maximal distance to below
                        # a la complete linkage clustering
                        max_dist = 0.0
                        for n in [tip, bro]:
                            try:
                                # if n was originally an internal node,
                                # there should be a clades_to_distances entry 
                                # associated. Distance is then that distance
                                # plus the distance between n and parental
                                d = clades_to_distances[n.id] + n.length
                            except KeyError:
                                d = n.length
                            if d > max_dist: max_dist = d
                        clades_to_distances[parental.id] = max_dist
                        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Setting maximal distance at %s" % max_dist)

                        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(tree.ascii_art())
                    else:
                        # Don't merge, but remove bro node from the list to visit,
                        # if it isn't already removed
                        if bro in tips_to_evaluate: tips_to_evaluate.remove(bro)
        
                # else there is multiple siblings. Don't need to add this back into the set again because it will be rediscovered later.
                else:
                    if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("The sibling tree is not a single tip, ignoring this node")
                    
        return clades_to_distances