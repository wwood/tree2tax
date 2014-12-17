from nose.tools import assert_equals, assert_true
from tree2tax.tree2tax import Tree2Tax
from skbio.tree import TreeNode
from StringIO import StringIO
from string import split as _
import IPython

class TestCoverageStats:
    def assertSameClusters(self, expected, observed):
        name_clusters = [sorted([n.name for n in named_cluster.tips]) for named_cluster in observed]
        expected_name_clusters = [sorted(node_names) for node_names in expected]
        #print "named clusters: %s" % name_clusters
        assert_equals(expected_name_clusters, name_clusters)
        
    def testSimple(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.25)
        self.assertSameClusters([['A','B'],['D']], clusters)
        
    def testNoClusering(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.05)
        self.assertSameClusters([['A'],['B'],['D']], clusters)
        assert_equals(_('C1 C2 Root1'), [c.name() for c in clusters])
         
    def testClusterEverything(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.5)
        self.assertSameClusters([['A','B','D']], clusters)
        assert_equals('Root1',clusters[0].name())
         
    def testClusterOnInternalNode(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12)C:10, D:9)E:20, F:20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([_('A B D'), ['F']], clusters)
         
    def testClusterOnTwoInternalNodes(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12)C:10, (H:8, D:9)I:3)E:20, F:20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([_('A B D H'), ['F']], clusters)
         
    def testClusterIntoThree(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12)C:10, (H:8, D:9)I:3)E:20, F:20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 25)
        self.assertSameClusters([_('A B'), _('D H'), ['F']], clusters)
        
    def testClusterNamingOnTwoInternalNodes(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12):10, (H:8, D:9):3):20, F:20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([_('A B D H'), ['F']], clusters)
        assert_equals(_('G1 G2'), [c.name() for c in clusters])
        
    def testClusterNamingOnTwoInternalNodesReverseOrder(self):
        tree = TreeNode.read(StringIO('((F:20, ((A:11, B:12):10, (H:8, D:9):3):20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([['F'], _('A B D H')], clusters)
        assert_equals(_('G2 G1'), [c.name() for c in clusters])
        