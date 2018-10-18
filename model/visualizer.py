import math
import networkx as nx
import matplotlib.pyplot as plt


class Visualizer:

    @staticmethod
    def _get_color_by_label(label):
        color = ''
        if label == 'Header':
            color = 'cyan'
        elif label == 'Data':
            color = 'green'
        elif label == 'Attributes':
            color = 'magenta'
        elif label == 'Derived':
            color = 'red'
        elif label == 'Metadata':
            color = 'darkorange'  # darkorange
        return color

    @staticmethod
    def _cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    @staticmethod
    def draw_fragmented_graph(graph, fragments, draw_direction=False):
        subgraph_nodes = set()
        edges_to_keep = set()
        for frag in fragments:
            subgraph_nodes.update(frag)
            frag_edges = graph.edges(frag)
            for edge in frag_edges:
                if edge[0] in frag and edge[1] in frag:
                    edges_to_keep.add(edge)

        partitioned = graph.subgraph(subgraph_nodes)
        to_omit = set(partitioned.edges()) - edges_to_keep
        partitioned.remove_edges_from(to_omit)

        Visualizer.draw_worksheet_graph(partitioned, draw_direction)

    @staticmethod
    def draw_worksheet_graph(graph, draw_direction=False):

        node_x_pos = dict()
        node_y_pos = dict()
        node_pos = dict()
        node_sizes = list()
        node_colors = list()
        for node in graph.nodes_iter(data=True):

            address = node[0]
            label = node[1]['label']
            xcenter = node[1]['xcenter']
            ycenter = node[1]['ycenter']
            area = node[1]['area']

            # we color the nodes according to their label
            color = Visualizer._get_color_by_label(label)
            node_colors.append(color)

            # the size of the node is adjusted based on its area
            lognorm_size = math.log(area, 2)
            node_sizes.append((300 + lognorm_size * 150))

            # the location of the node is determined
            # by the coordinates of its centroid
            node_pos[address] = []
            if xcenter in node_x_pos.keys():
                node_x_pos[xcenter].append(address)
            else:
                node_x_pos[xcenter] = [address]

            if ycenter in node_y_pos.keys():
                node_y_pos[ycenter].append(address)
            else:
                node_y_pos[ycenter] = [address]

        xcount = 0
        for x_pos in sorted(node_x_pos.keys()):
            for node_adrs in node_x_pos[x_pos]:
                node_pos[node_adrs].append(xcount)
            xcount += 1

        ycount = 0
        for y_pos in sorted(node_y_pos.keys()):
            for node_adrs in node_y_pos[y_pos]:
                node_pos[node_adrs].append(ycount)
            ycount -= 1

        # prepare edge labels
        edge_labels = dict()
        # prepare edge line format
        same_table_edges = list()
        multi_table_edges = list()
        diff_table_edges = list()
        tables_attributes = nx.get_node_attributes(graph, 'tables')
        for e in graph.edges_iter(data=True):

            # currently we use the weight as edge label
            w = float(e[2]['weight'])
            d = e[2]['direction']
            if draw_direction:
                edge_labels[(e[0], e[1])] = "%s\n%s" % (str(round(w, 3)), d[0])
            else:
                edge_labels[(e[0], e[1])] = "%s" % str(round(w, 3))

            # get table id for nodes in this edge
            tables_1 = set(ta[0] for ta in tables_attributes[e[0]])
            tables_2 = set(ta[0] for ta in tables_attributes[e[1]])

            # check if nodes have at least a table_id in common
            intersection = set.intersection(tables_1, tables_2)
            if len(intersection) > 0:
                difference = set.union(set.difference(tables_1, tables_2), set.difference(tables_2, tables_1))
                if len(difference) > 0:
                    multi_table_edges.append((e[0], e[1]))
                else:
                    same_table_edges.append((e[0], e[1]))
            else:
                diff_table_edges.append((e[0], e[1]))

        # set plot figure size
        xcount = xcount
        x = (xcount * 4) if (xcount * 4) <= 52 else 52
        ycount = ycount
        y = (math.fabs(ycount) * 4) if (math.fabs(ycount) * 4) <= 80 else 80
        plt.figure(figsize=Visualizer._cm2inch(x, y))

        # draw the graph
        nx.draw_networkx_nodes(graph, pos=node_pos, node_size=node_sizes,
                               node_color=node_colors, node_shape='s', alpha=0.25)
        nx.draw_networkx_labels(graph, pos=node_pos, labels=None, font_size=8, alpha=0.4)
        nx.draw_networkx_edges(graph, pos=node_pos, edgelist=same_table_edges, style='solid', alpha=0.3, arrows=False)
        nx.draw_networkx_edges(graph, pos=node_pos, edgelist=diff_table_edges, style='dotted', alpha=0.3, arrows=False)
        nx.draw_networkx_edges(graph, pos=node_pos, edgelist=multi_table_edges, style='dashdot', alpha=0.4,
                               arrows=False)
        nx.draw_networkx_edge_labels(graph, pos=node_pos, edge_labels=edge_labels)

        # display the graph
        plt.axis('off')
        plt.show()
