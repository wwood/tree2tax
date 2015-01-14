from nose.tools import assert_equals, assert_true
from tree2tax.tree2tax import Tree2Tax, NamedCluster, TaxonomyFunctions
from skbio.tree import TreeNode
from StringIO import StringIO
from string import split as _

class TestTaxonomyFunctions:
    def test_taxonomy_from_node_name(self):
        assert_equals('C', TaxonomyFunctions().taxonomy_from_node_name('C'))
        assert_equals('C', TaxonomyFunctions().taxonomy_from_node_name('0.997:C'))
        assert_equals(None, TaxonomyFunctions().taxonomy_from_node_name('0.1'))
        
    def test_missing_taxonomy(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12)C:10, D:9)E:20, F:20)G:30)root;'))
        assert_equals(['C'], TaxonomyFunctions().missing_taxonomy(tree, tree.find('A'), tree.find('E')))
        assert_equals([], TaxonomyFunctions().missing_taxonomy(tree, tree.find('A'), tree.find('A')))
        assert_equals(['E','C'], TaxonomyFunctions().missing_taxonomy(tree, tree.find('A'), tree.find('G')))

class TestTree2TaxNamedClusters:
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
        assert_equals(_('C.1 C.2 Root'), [c.name() for c in clusters])
          
    def testClusterEverything(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.5)
        self.assertSameClusters([['A','B','D']], clusters)
        assert_equals('Root',clusters[0].name())
          
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
        assert_equals(_('G.1 G.2'), [c.name() for c in clusters])
         
    def testClusterNamingOnTwoInternalNodesReverseOrder(self):
        tree = TreeNode.read(StringIO('((F:20, ((A:11, B:12):10, (H:8, D:9):3):20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([['F'], _('A B D H')], clusters)
        assert_equals(_('G.2 G.1'), [c.name() for c in clusters])
         
    def testNamingWithBootstraps(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)0.091:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.05)
        self.assertSameClusters([['A'],['B'],['D']], clusters)
        assert_equals(_('Root.1 Root.2 Root.3'), [c.name() for c in clusters])
         
    def testClusterNamingWithBootstraps(self):
        tree = TreeNode.read(StringIO("((F:20, ((A:11, B:12):10, (H:8, D:9):3):20)'0.7:G':30)root;"))
        clusters = Tree2Tax().named_clusters(tree, 40)
        self.assertSameClusters([['F'], _('A B D H')], clusters)
        assert_equals(_('G.2 G.1'), [c.name() for c in clusters])
         
    def testClusterNamingConventionsWithSomeUnnamed(self):
        tree = TreeNode.read(StringIO('((((A:11, B:12):10, D:9):20, F:20)G:30)root;'))
        clusters = Tree2Tax().named_clusters(tree, 0.05) #i.e. everything is a separate cluster
        self.assertSameClusters([['A'],['B'],['D'],['F']], clusters)
        assert_equals(['G.1', 'G.2', 'G.3', 'G.4'], [c.name() for c in clusters])
 
         
class TestTree2TaxNamedClusterSets:
    def assertSameClusterSets(self, expected, observed):
        assert_equals(len(expected), len(observed))
        for i, e in enumerate(expected):
            assert_equals(e[0], observed[i].threshold) #same thresholds
            self.assertSameClusters(e[1], observed[i].clusters)
     
    def assertSameClusters(self, expected, observed):
        name_clusters = [sorted([n.name for n in named_cluster.tips]) for named_cluster in observed]
        expected_name_clusters = [sorted(node_names) for node_names in expected]
        #print "named clusters: %s" % name_clusters
        assert_equals(expected_name_clusters, name_clusters)
         
    def testSimple(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [0.25])
        self.assertSameClusterSets([[0.25,[['A','B'],['D']]]], clusters)
 
    def testSimpleTwice(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [0.25, 0.25])
        self.assertSameClusterSets([[0.25,[['A','B'],['D']]], [0.25,[['A','B'],['D']]]], clusters)
         
    def testNormalType(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [0.25, 0.05])
        self.assertSameClusterSets([[0.05,[['A'],['B'],['D']]], [0.25,[['A','B'],['D']]]], clusters)
        assert_equals(_('C.1 C.2 Root'), [c.name() for c in clusters[0].clusters])
         
    def testOppositeSorting(self):
        tree = TreeNode.read(StringIO('((A:0.11, B:0.12)C:0.1, D:0.2)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [0.05, 0.25])
        self.assertSameClusterSets([[0.05,[['A'],['B'],['D']]], [0.25,[['A','B'],['D']]]], clusters)
        assert_equals(_('C.1 C.2 Root'), [c.name() for c in clusters[0].clusters])
         
    def testNaming(self):
        tree = TreeNode.read(StringIO('((F:20, ((A:11, B:12):10, (H:8, D:9):3):20)G:30)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [40, 25])
        self.assertSameClusterSets([[25,[['F'], _('A B'), _('D H')]], [40,[['F'], _('A B D H')]]], clusters)
        assert_equals(_('G.3 G.1 G.2'), [c.name() for c in clusters[0].clusters])
        assert_equals(_('G.2 G.1'), [c.name() for c in clusters[1].clusters])
 
    def testTipToCluster(self):
        tree = TreeNode.read(StringIO('((F:20, ((A:11, B:12):10, (H:8, D:9):3):20)G:30)root;'))
        clusters = Tree2Tax().named_clusters_for_several_thresholds(tree, [40, 25])
        self.assertSameClusterSets([[25,[['F'], _('A B'), _('D H')]], [40,[['F'], _('A B D H')]]], clusters)
        assert_equals(_('G.3 G.1 G.2'), [c.name() for c in clusters[0].clusters])
        assert_equals(_('G.2 G.1'), [c.name() for c in clusters[1].clusters])
         
        tip = tree.find('F')
        assert_equals('G.3', clusters[0].tip_to_cluster(tip).name())
        assert_equals('G.2', clusters[1].tip_to_cluster(tip).name())
 
        tip = tree.find('D')
        assert_equals('G.2', clusters[0].tip_to_cluster(tip).name())
        assert_equals('G.1', clusters[1].tip_to_cluster(tip).name())
 
class TestNamedCluster:
    def testCondensedName(self):
        nc = NamedCluster('c__Halo; o__fu', ['notips'], None)
        assert_equals('cHalo.ofu', nc.condensed_name())
         
        nc.cluster_number = 6
        assert_equals('cHalo.ofu.6', nc.condensed_name())
         
        nc.taxonomy = 'c__Halo'
        assert_equals('cHalo.6', nc.condensed_name())
     