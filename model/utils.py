import networkx as nx


class Utils:
    @staticmethod
    def remove_inview_edges(G, is_strict=False):
        updated = G.copy()
        edge_ttypes = nx.get_edge_attributes(G, 'target_type')
        for item in edge_ttypes.items():
            if item[1] == 'inview':
                if updated.has_edge(item[0][1], item[0][0], key=item[0][2]):
                    if is_strict or edge_ttypes[(item[0][1], item[0][0], item[0][2])] == 'inview':
                        updated.remove_edge(item[0][0], item[0][1])
                        updated.remove_edge(item[0][1], item[0][0])
        return updated

    @staticmethod
    def get_node_region_coordinates(nodes_attributes):

        half_width = nodes_attributes['width'] / 2
        xcenter = nodes_attributes['xcenter']
        xmin = xcenter - half_width
        xmax = half_width + xcenter

        half_height = nodes_attributes['height'] / 2
        ycenter = nodes_attributes['ycenter']
        ymin = ycenter - half_height
        ymax = half_height + ycenter

        return dict([('xmin', int(xmin)), ('ymin', int(ymin)), ('xmax', int(xmax)), ('ymax', int(ymax))])

    @staticmethod
    def group_node_area_tuples_by_tableId(graph, one_tbl_rule=False):

        tables_per_node = nx.get_node_attributes(graph, 'tables')
        nodes_per_table = dict()
        for node_tables in tables_per_node.items():
            node = node_tables[0]
            tables = node_tables[1]
            if one_tbl_rule:
                for table in sorted(tables, key=lambda item: item[1], reverse=True):
                    table_id = table[0]
                    area = table[1]
                    if table_id in nodes_per_table.keys():
                        nodes_per_table[table_id][node] = area
                    else:
                        nodes_per_table[table_id] = {node: area}
                    break
            else:
                for table in tables:
                    table_id = table[0]
                    area = table[1]
                    if table_id in nodes_per_table.keys():
                        nodes_per_table[table_id][node] = area
                    else:
                        nodes_per_table[table_id] = {node: area}

        return nodes_per_table

    @staticmethod
    def calculate_weighted_jaccard(dictA, dictB):

        keysA = set(dictA.keys())
        keysB = set(dictB.keys())

        intersection = set.intersection(keysA, keysB)
        difference = set.union(set.difference(keysA, keysB), set.difference(keysB, keysA))

        min_sum_intersection = 0
        max_sum_intersection = 0
        for key in intersection:
            min_sum_intersection += min(dictA[key], dictB[key])
            max_sum_intersection += max(dictA[key], dictB[key])

        sum_difference = 0
        for key in difference:
            if key in keysA:
                sum_difference += dictA[key]
            else:
                sum_difference += dictB[key]

        weighted_jaccard = min_sum_intersection / (max_sum_intersection + sum_difference)
        return weighted_jaccard

    @staticmethod
    def evaluate_identified_table(sheet_graph, subgraph):

        node_area_tuples_per_table = Utils.group_node_area_tuples_by_tableId(sheet_graph)
        jaccard_per_table_id = dict()
        subgraph_node_area_tuples = nx.get_node_attributes(subgraph, 'area')
        for tid in node_area_tuples_per_table.keys():
            table_node_area_tuples = node_area_tuples_per_table[tid]
            weighted_jaccard_index = Utils.calculate_weighted_jaccard(table_node_area_tuples, subgraph_node_area_tuples)
            if weighted_jaccard_index > 0.0:
                jaccard_per_table_id[tid] = weighted_jaccard_index

        return jaccard_per_table_id
