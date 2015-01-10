from skbio.tree import TreeNode
from sets import Set
import logging
import re
import IPython

class NamedCluster:
    def __init__(self, taxonomy, tips):
        self.taxonomy = taxonomy
        self.tips = tips
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
        then return 'Halobacteria.o__Halobacteriales.f__MSP41' (cluster number
        is added too, if that is present'''
        splits = self.taxonomy.split('; ')
        
        # get rid of f__ prefixes etc.
        regex = re.compile(r'^.__')
        splits2 = []
        for s in splits:
            reg = regex.match(s)
            if reg:
                splits2.append(s[reg.end():])
            else:
                logging.debug("Found unexpected form for taxonomy in %s", str(s))
                splits2.append(s)
                               
        joined = '.'.join(splits2)
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
                    if self._tip_to_cluster.has_key(key):
                        logging.warn("Unexpectedly found multiple leaf nodes with the same name, undefined behaviour possibly immenent: %s" % tip_name)
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
    def version(self):
        '0.0.1'
        
    def named_clusters_for_several_thresholds(self, original_tree, thresholds):
        '''Given a list of thresholds, return a iterable of ThresholdAndClusters
        where the clustering has been done iteratively, providing a consistent
        taxonomic annotation scheme'''
        original_tree2 = original_tree.copy()
        original_tree2.assign_ids()
        tree = original_tree2.copy() #this copy gets destructively pruned as part of the algorithm
        
        # sort from smallest to largest because smaller distance thresholds
        # need to be applied before larger thresholds
        sorted_thresholds = sorted(thresholds)
        if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug("Found thresholds %s", str(sorted_thresholds))
        
        to_return = []
        
        def isFloat(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
            
        bootstrap_regex = re.compile(r'^[\d\.]+:')

        # skbio seems to be quite slow running repeated find_by_id. Do some caching 
        # to speed things up
        id_to_node = {}
        for node in original_tree2.postorder(include_self=True):
            id_to_node[node.id] = node
            
        clades_to_distances = {}
        for threshold in sorted_thresholds:
            if logging.getLogger().isEnabledFor(logging.INFO): logging.info("Clustering on threshold %s.." % str(threshold))
            
            clades_to_distances = self.destructively_cluster_tree(tree, threshold, clades_to_distances)
            
            clusters = []
            if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(tree.ascii_art())
            cluster_lcas = tree.tips()
            
            # So that several clades don't get named the same thing
            node_id_to_named_nodes = {} 
            
            # when everything is in 1 cluster (doesn't happen in practice I suspect)
            # but there is a unit test..
            if len(tree) == 0:
                cl = NamedCluster('Root',original_tree2.tips())
                cl.cluster_number = ''
                clusters = [cl]
    
            else:
                for lca in cluster_lcas:
                    orig = id_to_node[lca.id]
                    
                    # Find the closest ancestral named node (where ancestral does not include self)
                    current = orig
                    #in a greengenes file, tips have names, but don't count these 
                    #because they are not taxonomy but rather prokMSA IDs
                    if orig.is_tip(): current = orig.parent
                    
                    # Keep proceeding up the tree until a node with taxonomy is found.
                    # If the tree has bootstraps, then current.name is a float (in string form)
                    # Ignore these floats because they aren't named taxonomy
                    while current.parent and (current.name is None or isFloat(current.name)):
                        current = current.parent
                        
                    if current.parent is None:
                        taxonomy = 'Root'
                    else:
                        taxonomy = current.name
                        reg = bootstrap_regex.match(taxonomy)
                        if reg:
                            taxonomy = taxonomy[reg.end():]
                        
                    if orig.is_tip():
                        named_cluster = NamedCluster(taxonomy, [orig])
                    else:
                        named_cluster = NamedCluster(taxonomy, list(orig.tips()))
                        
                    clusters.append(named_cluster)
                    
                    try:
                        node_id_to_named_nodes[current.id].append(named_cluster)
                    except KeyError:
                        node_id_to_named_nodes[current.id] = []
                        node_id_to_named_nodes[current.id].append(named_cluster)
                
                if logging.getLogger().isEnabledFor(logging.DEBUG): logging.debug(node_id_to_named_nodes)
                # Assign cluster numbers by decreasing numbers of tips i.e. for each
                # taxonomy, the cluster with the most sequences is given the number 1,
                # second most abundant is number 2, etc.
                for taxonomy, named_clusters in node_id_to_named_nodes.items():
                    if len(named_clusters) > 1:
                        number = 1
                        for clade in sorted(named_clusters, reverse = True, key = lambda c: len(c.tips)):
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