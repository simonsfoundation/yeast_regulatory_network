
from jp_gene_viz import dGraph, dNetwork, dLayout, LExpression, dExpression
from jp_gene_viz import multiple_network, getData, HMap
from jp_gene_viz import grid_forest
import os

dNetwork.load_javascript_support()

TF = "TF.GN"
TARG = "Targ.GN"
TFOLN = "TF.OLN"
TARGOLN = "Targ.OLN"
PRECISION = "precision"
GS = 'in.GS'
YD = 'in.YD'
YI =  'in.YI'
K =  'in.K'
SGD = 'in.SGD'
COND = 'cond.clust'
WHICH = 'which.clust'
CLUSTER = 'gene.clust'
SIGN = "sign"
IN_NETWORKS = (GS, YD, YI, SGD, K)

class YeastData(object):

    def __init__(self, filename="data/compare-old-vs-new.tsv", expression_fn="data/expression-reduced-200.tsv", limit=2000):
        self.filename = filename
        self.expression_fn = expression_fn
        self.file = open(filename)
        self.headers = ["counter"] + self.get_line_data(strict=True)
        #print "headers", self.headers
        line_dicts = []
        oln_to_gn = {}
        while len(line_dicts) < limit:
            line_dict = self.get_line_dict()
            if line_dict is None:
                break
            #print "line_dict", line_dict
            line_dicts.append(line_dict)
            try:
                oln_to_gn[line_dict[TFOLN]] = line_dict[TF]
                oln_to_gn[line_dict[TARGOLN]] = line_dict[TARG]
            except:
                print len(line_dicts), "bad line_dict", line_dict
                raise
        self.oln_to_gn = oln_to_gn
        self.line_dicts = line_dicts
        self.gene_clusters = set(line_dict[CLUSTER] for line_dict in line_dicts)
        self.cond_clusters = set(line_dict[COND] for line_dict in line_dicts)
        G = self.precision_graph = self.graph()
        layout_path = filename + "." + str(limit)+ ".layout.json"
        layout_exists = os.path.exists(layout_path)
        if layout_exists:
            print ("Loading saved layout", layout_path)
            self.layout = dLayout.load(layout_path)
        else:
            print("computing network layout.")
            #self.layout = dLayout.group_layout(self.precision_graph)[0]
            (layout, rectangles) = grid_forest.forest_layout(G)
            self.layout = layout
            print ("Saving layout", layout_path)
            dLayout.dump(self.layout, layout_path)

    def pairs_in_cluster(self, cluster_value, cluster_field=CLUSTER):
        for line_dict in self.line_dicts:
            if line_dict[cluster_field] != cluster_value:
                continue
            tf = line_dict[TF]
            targ = line_dict[TARG]
            value = float(line_dict[SIGN])
            yield (tf, targ, value)

    def pairs_in_which(self, index):
        for line_dict in self.line_dicts:
            which = line_dict[WHICH]
            if len(which) <= index:
                print "bad which indicator string ", repr(which)
                continue
            if which[index] != "1":
                continue
            tf = line_dict[TF]
            targ = line_dict[TARG]
            value = float(line_dict[SIGN])
            yield (tf, targ, value)

    def which_cluster_multi(self):
        L = [self.network(self.pairs_in_which(i), "which:" + repr(i)) for i in range(4)]
        L = [self.expression_network()] + L
        return self.two_column_multi(L)

    def precisions(self):
        for line_dict in self.line_dicts:
            tf = line_dict[TF]
            targ = line_dict[TARG]
            value = float(line_dict[PRECISION])
            yield (tf, targ, value)

    def cond_clusters_multi(self):
        L = [self.network(self.pairs_in_cluster(v, COND), v) for v in self.cond_clusters]
        L = [self.expression_network()] + L
        return self.two_column_multi(L)

    def gene_clusters_multi(self):
        L = [self.network(self.pairs_in_cluster(v), v) for v in self.gene_clusters]
        L = [self.expression_network()] + L
        return self.two_column_multi(L)

    def in_networks_multi(self):
        L = [self.expression_network(), self.get_network("in.any")] + [self.get_network(n) for n in IN_NETWORKS]
        return self.two_column_multi(L)

    def two_column_multi(self, networks):
        L = []
        networks = networks[:]
        while networks:
            L.append(networks[:2])
            networks = networks[2:]
        return multiple_network.MultipleNetworks(L)

    def in_network(self, field=GS):
        for line_dict in self.line_dicts:
            test = int(line_dict[field])
            if test != 0: 
                tf = line_dict[TF] * test
                targ = line_dict[TARG]
                value = float(line_dict[PRECISION]) * test
                yield (tf, targ, value)

    def in_any_network(self):
        for line_dict in self.line_dicts:
            test = 0
            for name in IN_NETWORKS:
                test += int(line_dict[name])
            if test != 0: 
                tf = line_dict[TF] * test
                targ = line_dict[TARG]
                value = float(line_dict[PRECISION])
                yield (tf, targ, value)

    def graph(self, iterator=None):
        if iterator is None:
            iterator = self.precisions()
        G = dGraph.WGraph()
        for (tf, targ, value) in iterator:
            G.add_edge(tf, targ, value)
        return G

    def get_network(self, name=PRECISION):
        if name == PRECISION:
            return self.network()
        elif name in IN_NETWORKS:
            return self.network(self.in_network(name), name)
        elif name == "in.any":
            return self.network(self.in_any_network(), name)
        else:
            raise ValueError("Unknown network name " + repr(name))

    def network(self, iterator=None, name=PRECISION, threshold=0, N=None):
        if iterator is None:
            iterator = self.precisions()
        G = self.graph(iterator)
        if N is None:
            N = dNetwork.NetworkDisplay()
        N.layout_dropdown.value = dNetwork.SPOKE
        N.load_data(G, self.layout, draw=False)
        N.threshhold_slider.value = threshold
        N.title_html.value = name
        N.restore_click()
        N.container_dropdown.value = dNetwork.CANVAS
        return N

    def expression_network(self, side_length=500):
        filename = self.expression_fn
        L = LExpression.LinkedExpressionNetwork()
        N = L.network
        N = self.network(N=N)
        E = L.expression
        #dExpression.display_heat_map(filename, E)
        H = HMap.HeatMap()
        (rows, cols, data) = getData.read_tsv(filename)
        rows = [unquote(r) for r in rows]
        rows = [self.oln_to_gn.get(r, r) for r in rows]
        cols = [unquote(c) for c in cols]
        H.set_data(rows, cols, data)
        E.load_data(H, side_length=side_length)
        return L

    def get_line_dict(self):
        line_data = self.get_line_data()
        if line_data is None:
            return None
        return {h: d for (h, d) in zip(self.headers, line_data)}

    def get_line_data(self, strict=False):
        line_data = []
        line = self.file.readline()
        if not line:
            if strict:
                raise ValueError("read line past eof")
            else:
                return None
        line = line.strip()
        for datum in line.split("\t"):
            line_data.append(unquote(datum))
        return line_data

def unquote(s):
    if s and s[0] == s[-1] and s[0] in "\"'":
        return s[1:-1]
    return s
