from skbio.tree import TreeNode
from StringIO import StringIO
import sys
import os
import unittest

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..'))
from tree2tax.threshold_finder import ThresholdFinder,\
    ThresholdInconsistencyException

class Tests(unittest.TestCase):
    def assertSameCladeDistanceSet(self, expected, observed):
        self.assertEqual(len(expected), len(observed))
        for i, exp in enumerate(expected):
            self.assertEqual(exp[0], observed[i].parent_node.name)
            self.assertEqual(sorted([exp[1], exp[2]]), \
                          sorted([observed[i].daughter_node1.name, observed[i].daughter_node2.name]))
            self.assertEqual(exp[3], observed[i].distance)
    
    def testSimple(self):
        tree = TreeNode.read(StringIO("(((A:1, B:2)'g__genus1':3, (C:4, D:5)'g__genus2':6)'f__family':10)root;"))
        examples = ThresholdFinder().find_examples(tree, 'f', 'g')
        self.assertSameCladeDistanceSet([['f__family','g__genus1','g__genus2',16.0]],
                                         examples)
    
    def testTreeSubtree(self):
        '''one genus is a subtree of another'''
        tree = TreeNode.read(StringIO("((((A:1, B:2)'g__genus1':3, D:50)'g__genus2':6)'f__family':10)root;"))
        examples = ThresholdFinder().find_examples(tree, 'f', 'g')
        self.assertSameCladeDistanceSet([['f__family','g__genus1','g__genus2',55.0]],
                                         examples)
    def testTreeSubtree2(self):
        '''one genus is a subtree of another, and the longest branch is in both subtrees'''
        tree = TreeNode.read(StringIO("((((A:1, B:52)'g__genus1':3, D:50)'g__genus2':6)'f__family':10)root;"))
        examples = ThresholdFinder().find_examples(tree, 'f', 'g')
        self.assertSameCladeDistanceSet([['f__family','g__genus1','g__genus2',105.0]],
                                         examples)
        
    def testMultiplyNamedNode(self):
        tree = TreeNode.read(StringIO("(((A:1, B:2)'g__genus1':3, (C:4, D:5)'g__genus2; s__spec':6)'f__family':10)root;"))
        examples = ThresholdFinder().find_examples(tree, 'f', 'g')
        self.assertSameCladeDistanceSet([['f__family','g__genus1','g__genus2; s__spec',16.0]],
                                         examples)
        
    def testNoPairs(self):
        tree = TreeNode.read(StringIO("(((A:1, B:2):3, (C:4, D:5):6)'f__family; g__genoos':10)root;"))
        examples = ThresholdFinder().find_examples(tree, 'f', 'g')
        self.assertSameCladeDistanceSet([],
                                         examples)
        
    def test_find_thresholds_from_example_distances(self):
        finder = ThresholdFinder()
        self.assertEqual([2.5,1.5], finder._find_thresholds_from_example_distances([3,2,1]))
        self.assertEqual([3.0,2.5,1.5], finder._find_thresholds_from_example_distances([None,3,2,1]))
        self.assertEqual([2.5,1.5,1.0], finder._find_thresholds_from_example_distances([3,2,1,None]))
        self.assertEqual([2.5,1.5,1.0], finder._find_thresholds_from_example_distances([3,2,1,2]))
        with self.assertRaises(ThresholdInconsistencyException):
            finder._find_thresholds_from_example_distances([3,2,1,2,1])
        
if __name__ == "__main__":
    unittest.main()