from nose.tools import assert_equals, assert_true
from tree2tax.threshold_finder import ThresholdFinder
from skbio.tree import TreeNode
from StringIO import StringIO
from string import split as _

class TestThresholdFinder:
    def assertSameCladeDistanceSet(self, expected, observed):
        assert_equals(len(expected), len(observed))
        for i, exp in enumerate(expected):
            assert_equals(exp[0], observed[i].parent_node.name)
            assert_equals(sorted([exp[1], exp[2]]), \
                          sorted([observed[i].daughter_node1.name, observed[i].daughter_node2.name]))
            assert_equals(exp[3], observed[i].distance)
    
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
