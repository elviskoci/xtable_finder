from random import randint, shuffle
import networkx as nx


class Generator:

    @staticmethod
    def get_random_fragmentation(G, variable_edges=None, fixed_edges=None, as_subgraphs=False):

        if G.number_of_nodes() < 2:  # no fragmentation is possible
            return list()

        if variable_edges:
            edges_list = list(variable_edges)
        else:
            edges_list = list()
            for edge in G.edges():
                if (edge[1], edge[0]) not in edges_list:
                    edges_list.append((edge[0], edge[1]))

        shuffle(edges_list)
        index = randint(0, len(edges_list) - 1)  # avoid picking all edges
        selected = edges_list[:index]
        if fixed_edges:
            selected.extend(fixed_edges)

        temp = nx.Graph()
        temp.add_nodes_from(G)
        temp.add_edges_from(selected)
        if not as_subgraphs:
            return list(nx.connected_components(temp))
        else:
            return [G.subgraph(comp) for comp in list(nx.connected_components(temp))]

    @staticmethod
    def get_random_fragmentations(G, num, as_subgraphs=False):
        # determine what is the total number of possible fragmentations for the given graph
        # this is an estimation, i.e., an upper bound
        # the presence of cycles reduces the possible number of fragmentations
        unique_undirected = set()
        if G.number_of_edges() == 0 and G.number_of_nodes() > 1:
            possible_fragmentations = 1
        elif G.number_of_edges() == 0:
            possible_fragmentations = 0
        else:
            for edge in G.edges():
                if (edge[1], edge[0]) not in unique_undirected:
                    unique_undirected.add((edge[0], edge[1]))
            possible_fragmentations = 2 ** len(unique_undirected) - 1

        # update the target number of fragmentations to generate
        if num > possible_fragmentations:
            threshold = possible_fragmentations
        else:
            threshold = num

        # use this to terminate the loop, in case the threshold is not achieved
        iteration = 1
        max_iteration = threshold + int(threshold * 2.0)
        if max_iteration < 10:
            max_iteration = 10

        fragmentations = set()
        while len(fragmentations) < threshold:
            frag = Generator.get_random_fragmentation(G, variable_edges=unique_undirected, as_subgraphs=False)
            frozen_frag = frozenset(frozenset(part) for part in frag)

            if len(frozen_frag) > 1 and frozen_frag not in fragmentations:
                fragmentations.add(frozen_frag)

            if iteration == max_iteration:
                break
            iteration += 1

        output = list()
        if as_subgraphs:
            for frag in fragmentations:
                output.append(list(G.subgraph(node_set) for node_set in frag))
        else:
            output.extend(list(frag) for frag in fragmentations)
        return output
