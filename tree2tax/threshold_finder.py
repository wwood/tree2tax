import re
import logging
from sets import Set



class CladeDistanceSet:
    def __init__(self):
        self.clade_distances = []
    
    def append(self, clade_distance):
        self.clade_distances.append(clade_distance)
        
    def average_distance(self):
        sum([cd.distance for cd in self.clade_distances]) / len(self.clade_distancess)

class CladeDistance:
    def __init__(self, parent_node, daughter_node1, daughter_node2, distance):
        self.parent_node = parent_node
        self.daughter_node1 = daughter_node1
        self.daughter_node2 = daughter_node2
        self.distance = distance
        
class Stack:
    def __init__(self):
        self.__storage = []

    def isEmpty(self):
        return len(self.__storage) == 0

    def push(self,p):
        self.__storage.append(p)

    def pop(self):
        return self.__storage.pop()
    
class ThresholdInconsistencyException(Exception): pass

class ThresholdFinder:
    def find_thresholds(self, annotated_tree, prefixes):
        '''Find tree distance thresholds that separate taxonomic ranks. Start
        by finding a list of examples of e.g. the distances between two genera
        in a family, and then return half-way between the medians of each 
        taxonomic level.
        
        Parameters
        ----------
        annotated_tree: TreeNode
            An tree with taxonomic annotation decorated on it, in the node
            names.
        prefixes: list of single characters strings
            A list of prefixes representing each taxonomic rank e.g.
            string.split('k p c o f g s')
        '''
        median_distances = []
        for i, level_prefix in enumerate(prefixes):
            if i==0: continue
            examples = self.find_examples(annotated_tree, prefixes[i-1], level_prefix)
            if len(examples) > 0:
                median_distances.append(self._median([e.distance for e in examples]))
                logging.info("Found distances for %i pairs for level prefix %s" % (len(examples), level_prefix))
            else:
                median_distances.append(None)
                logging.warn("No pairs were found for level prefix %s" % level_prefix)
        return self._find_thresholds_from_example_distances(median_distances)
            
    def _median(self, lst):
        '''Return median of a list of floats'''
        even = (0 if len(lst) % 2 else 1) + 1
        half = (len(lst) - 1) / 2
        return sum(sorted(lst)[half:half + even]) / float(even)
    
    def _find_thresholds_from_example_distances(self, example_distances):
        '''Given a list of floats or Nones, return a list which is halfway
        between each element of the list. However, deal gracefully with Nones
        which indicate unknown estimations, and require the distances to be
        monotonically decreasing (since tree distance thresholds must decrease
        as we progress down the taxonomic levels.
        
        Parameters
        ----------
        example_distances: list of float/None
            distances to find thresholds from
        '''
        # find first float
        # set Nones before to this top level
        # the last one might be big, set it to the next lowest
        # if None or a bigger one is detected in the middle, error out
        distances2 = []
        state = 'start'
        for i, dist in enumerate(example_distances):
            if state == 'start':
                if dist is not None:
                    for _ in range(0, i+1):
                        distances2.append(dist)
                        state = 'middle'
            elif state == 'middle':
                if dist is None or dist > distances2[i-1]:
                    state = 'possible_end'
                    distances2.append(distances2[i-1])
                else:
                    distances2.append(dist)
            elif state == 'possible_end':
                if dist is not None:
                    raise ThresholdInconsistencyException("Failed to sanitise taxonomy-wise distances %s" % str(example_distances))
                distances2.append(distances2[i-1])
            else:
                raise Exception("Programming error")
            
        if len(distances2) == 0:
            # All None
            raise ThresholdInconsistencyException("Failed to sanitise taxonomy-wise distances %s" % str(example_distances))
        
        thresholds = []
        last_dist = None
        for i, dist in enumerate(distances2):
            if i > 0:
                thresholds.append((dist+last_dist)/2.0)
            last_dist = dist
        
        return thresholds
        
        
    
    def find_examples(self, tree, upper_prefix, lower_prefix):
        '''return a CladeDistanceSet of pairs of clades from the lower_prefix
        rank that are sisters in the upper_prefix rank'''
        to_return = []
        
        # get a list of nodes that are at the upper threshold, by
        # descending the tree
        upper_prefix_regex = re.compile(r'(\d+:){0,1}.*(; ){0,1}%s__' % upper_prefix)
        upper_nodes = []
        for node in tree.non_tips():
            if node.name: 
                if upper_prefix_regex.match(node.name):
                    upper_nodes.append(node)
                else:
                    pass
                    #logging.debug("Node name does not match expected type: %s" % node.name)
        logging.debug("Found %s upper nodes" % len(upper_nodes))
        
        # For each of these upper class prefixes    
        # find all the lower level nodes
        lower_prefix_regex = re.compile(r'(\d+:){0,1}.*(; ){0,1}%s__' % lower_prefix)    
        for unode in upper_nodes:
            lower_nodes = []
            for node in unode.non_tips(True):
                if node.name and lower_prefix_regex.match(node.name):
                    lower_nodes.append(node)
            logging.debug("Lower nodes: %s" % lower_nodes)
                
            # if there is only 1, go to next upper because there are no pairs
            if len(lower_nodes) <= 1: continue
        
            # for each of the lower level nodes, find the maximal distance to the tips
            # special cases - if one of the lower nodes is the LCA,
            # then the max distance of that LCA cannot be derived from a 
            # tip in the other
            lower_max_distances = self.unique_best_distance(lower_nodes)
                                   
            # make/add to a list of all the pairwise distances,
            for i, first in enumerate(lower_nodes):
                for j, second in enumerate(lower_nodes):
                    if i <= j: continue #upper right triangle of distance matrices only
                    
                    # by taking each pair of lower nodes, 
                    # distance is distance of node1 to lca + distance of node2 to lca
                    # + max distance underneath of each node
                    
                    # there is currently (v 0.2.2) a bug in scikit-bio's implementation
                    # of distance - so workaround here
                    # https://github.com/biocore/scikit-bio/issues/807
                    lca = self.my_lca(tree, first, second)
                    distance = first.accumulate_to_ancestor(lca) + \
                        second.accumulate_to_ancestor(lca) + \
                        lower_max_distances[i] + \
                        lower_max_distances[j]
                        
                    to_return.append(CladeDistance(unode, first, second, distance))
        
        # return the list of distances
        return to_return
    
    def unique_best_distance(self, nodes):
        '''given a set of nodes, return an array of max distances to tips
        that are descendants of each node. However, the 'max tip' cannot be a
        descendant of any other node. Assumes no tips are counted as a clade'''
        to_return = []
        
        annotated_nodes = Set(nodes)
            
        for node in nodes:
            # Make a set for fast membership querying

            stack = Stack()
            max_distance = 0.0
            
            # initialise stack with first descendents
            for d in node.children:
                stack.push(d) 
                
            while True:
                try:
                    d = stack.pop()
                except IndexError:
                    break #stack empty, we are done
                
                if d not in annotated_nodes: # don't descend when we run into a different clade
                    if d.is_tip():
                        distance = d.accumulate_to_ancestor(node)
                        if distance > max_distance: 
                            max_distance = distance
                    else:
                        [stack.push(desc) for desc in d.children]

                
            to_return.append(max_distance)
        return to_return
    
    def my_lca(self, tree, first, second):
        '''v0.2.2 of skbio is buggy, here's a quick fix'''
        first_ancestors = Set([first] + first.ancestors())
        curr = second
        while curr:
            if curr in first_ancestors:
                return curr
            else:
                curr = curr.parent
                
        raise Exception("No common ancestors found, whack - are they from the same tree?") 
    
    
