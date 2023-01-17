#!/usr/bin/env python3
#################################
# This is the script for tree analysis
# Establish Raman tree and visualize it
# This script includes three class
# - Raman_create_tree : generate a newick file for Raman df (n x d)
# - Raman_read_tree : read a newick file
# - Raman_analyze_tree : dictionary comparison and Mantel Test are included here. More will be integrated after tests in R.
#################################

from load import *

class create_tree:

    def __init__(self,data):
        """
        The class is to create a newick file for a Raman dataframe: cell x feature
        :param data: DataFrame
        """
        self.data = data
        self.node = None
        self.labels = [i.replace(' ','_') if 'O157:H7' not in i else \
                       i.replace(' ','_').replace('O157:H7','O157_H7')\
                       for i in self.data.index.values.tolist()]


    def getNewick(self,node, newick, parentdist, leaf_names):
        #     print(len(leaf_names))
        if node.is_leaf():
            #         print(node.id)
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = self.getNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick

    def get_plot_dendrogram(self,model='AC', fs=(10, 10), **kwargs):
        if model == 'AC':
            model = sklearn.cluster.AgglomerativeClustering(n_clusters=None,
                                                     distance_threshold=0,
                                                     affinity='euclidean',
                                                     linkage='ward')
            model = model.fit(self.data)
        # Create linkage matrix and then plot the dendrogram

        # create the counts of samples under each node
        counts = np.zeros(model.children_.shape[0])
        n_samples = len(model.labels_)
        for i, merge in enumerate(model.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  # leaf node
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count

        linkage_matrix = np.column_stack(
            [model.children_, model.distances_, counts]
        ).astype(float)
        fig, ax = plt.subplots(figsize=fs)
        # Plot the corresponding dendrogram
        scipy.cluster.hierarchy.dendrogram(linkage_matrix, ax=ax, labels=self.labels,
                                           truncate_mode="level",
                                           p=0,
                                            orientation='left',
                                            leaf_rotation=0,
                                            **kwargs)
        tree = scipy.cluster.hierarchy.to_tree(linkage_matrix)
        self.node = tree#save the node to attribute
        nw = self.getNewick(self.node, "", tree.dist, self.labels).replace('Escherichia coli serotype O157:H7', \
                                                        'Escherichia_coli_O157_H7')
        nw = nw.replace(' ','_')
        return nw, fig

class read_tree:

    def __init__(self,in_path,type=None,title=None,formatete3=1,format='newick',quoted_node_names=False):
        """
        :param ts: tree style
                in_path: path for newick file or newick string
        """
        # circular
        # ts.mode = "c"
        # ts.arc_start = -180 # 0 degrees = 3 o'clock
        # ts.arc_span = 180

        if type=='16S':
            t = "16S Phylogenetic Tree"
            self.title = t
        elif type=='Raman':
            t = "SCRS Phenotyping Tree"
            self.title = t
        if title != None:
            self.title = title

        tree_g = ete3.Tree(in_path,formatete3,quoted_node_names)  # the E coli name is not correct without underline
        ts = ete3.TreeStyle()
        ts.mode = 'c'#circular
        ts.arc_span = 360
        ts.show_leaf_name = True
        ts.title.add_face(ete3.TextFace(self.title, fsize=20), column=0)  # add legend and title

        self.ts = ts
        self.tree_g = tree_g

        #get distance matrix
        tree_g2 = Bio.Phylo.read(in_path,format)#for the distance matrix calculation
        self.tree_g2 = tree_g2
        dm = self.get_distance_matrix()
        normed_dm = Raman_preprocess.normalize(dm)
        self.dm = dm
        self.normed_dm = normed_dm# the normed one here is applied for Mantel test

    def get_distance_matrix(self):
        t = self.tree_g2
        d = {}
        for x, y in itertools.combinations(t.get_terminals(), 2):
            v = t.distance(x, y) # Add and return the sum of the branchlengths between two nodes.
            d[x.name] = d.get(x.name, {})
            d[x.name][y.name] = v
            d[y.name] = d.get(y.name, {})
            d[y.name][x.name] = v
        for x in t.get_terminals():
            d[x.name][x.name] = 0
        dm = pd.DataFrame(d).sort_index().sort_index(axis=1)
        return dm

    def set_color_range(self, genus_to_color):
        def get_node_color(tree_g, genus_to_color):
            """
            tree_g: ETE tree object
            genus_to_color: dictionary obtained by genus_to_color
            """
            node_name = []
            node_color = []
            for node in self.tree_g.traverse("postorder"):
                temp_key = node.name.replace(' ','_').split('_')[0]
                if len(temp_key) >= 1:  # in case the empty string
                    node_name.append(node.name)
                    try:
                        temp_c = genus_to_color[temp_key]
                    except KeyError:
                        try:
                            temp_c = genus_to_color["'"+temp_key.replace("'","")+"'"]
                        except KeyError:
                            print(node.name)
                    node_color.append(temp_c)
            node_genus_color = dict(zip(node_name, node_color))
            return node_genus_color

        node_genus_color = get_node_color(self.tree_g, genus_to_color)
        for genus, color in node_genus_color.items():
            temp_g = self.tree_g & genus
            temp_g.img_style["bgcolor"] = color
        self.tree_g.set_style(ete3.NodeStyle())
        self.tree_g.img_style["bgcolor"] = "white"
        self.tree_g.img_style["size"] = 5
        self.tree_g.img_style["fgcolor"] = 'black'

        return self.tree_g, self.ts

class analyze_tree:
    """
    This method takes in a dictionary: tree_dict[tree_name] = Raman_read_tree object
    """
    def __init__(self,tree_dict):
        self.tree_dict = tree_dict
        combine = itertools.combinations(self.tree_dict.keys(), 2)
        self.combine  = combine # the combination of any two of the dictionary keys

    def compare_tree(self):
        tree_dict = self.tree_dict
        combine = self.combine
        compare_tree = {}
        for i in combine:
            tree1 = tree_dict[i[0]]
            tree2 = tree_dict[i[1]]
            temp_retuls = ete3.Tree.compare(tree1.tree_g, tree2.tree_g, unrooted=True)

            a = temp_retuls.pop('common_edges')
            temp_retuls.pop('source_edges')
            temp_retuls.pop('ref_edges')
            compare_tree[i] = temp_retuls
        return compare_tree

    def Mantel_analysis(self,title='Mantel Analysis',normed=True,method='pearson',**kwargs):
        tree_dict = self.tree_dict
        combine = self.combine
        Mantel_analysis = {}

        for i in combine:
            tree1 = tree_dict[i[0]]
            tree2 = tree_dict[i[1]]
            if normed == True:
                # print(tree1.normed_dm)
                temp_retuls = skbio.stats.distance.mantel(tree1.normed_dm,tree2.normed_dm,method=method,\
                                                          **kwargs) # corr_coeff, p_value, n
            elif normed == False:
                temp_retuls = skbio.stats.distance.mantel(tree1.dm, tree2.dm,method=method, \
                                                          **kwargs)  # corr_coeff, p_value, n
            Mantel_analysis[i] = list(temp_retuls)
        Mantel_df = pd.DataFrame.from_dict(Mantel_analysis)
        Mantel_df.columns = [i[0]+'_'+i[1] for i in Mantel_analysis.keys()]
        Mantel_df.index = ['correlation coefficient','p_value','n']
        Mantel_df = Mantel_df.T
        fig,ax = plt.subplots(figsize=(5,3))
        ax = Mantel_df.plot.bar(y='correlation coefficient',ax=ax)
        xticks = ax.get_xticks()
        #add significance star
        for i in range(len(xticks)):
            p_value = Mantel_df.iloc[i,1]
            y = Mantel_df.iloc[i,0]
            if p_value<0.05:
                ax.text(xticks[i],1.008*y,'*',)
        ax.set_ylabel(method, fontsize = 'x-large')
        ax.set_title(title)
        ax.legend(frameon=False)
        return Mantel_analysis, fig

