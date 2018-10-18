import networkx as nx


class Builder:

    @staticmethod
    def create_graph_for_sheet(worksheet, weight_function=None,  kwargs={}):

        # initial graph is directed and allows multiple edges between pairs of nodes.
        DG = nx.MultiDiGraph(corpus=worksheet.corpus_name, file=worksheet.file_name, sheet=worksheet.sheet_name)

        if worksheet.row_heights and worksheet.column_widths:
            DG.graph['row heights'] = worksheet.row_heights
            DG.graph['column widths'] = worksheet.column_widths
        else:
            DG.graph['row heights'] = None
            DG.graph['column widths'] = None

        # create graph nodes
        for reg in worksheet.regions:
            DG.add_node(reg.address, label=reg.label, xcenter=reg.xcenter, ycenter=reg.ycenter, area=reg.area,
                        width=reg.width, height=reg.height, tables=reg.tables)

        # create graph edges
        for rel in worksheet.relations:
            DG.add_edge(rel.source, rel.target, direction=rel.direction, distance=rel.distance,
                        overlap_ratio=rel.overlap_ratio, target_type=rel.target_type)

        if weight_function is None:
            weight_function = Builder.relevance_as_weight  # currently, relevance is used as default edge weight

        return weight_function(DG, **kwargs)

    @staticmethod
    def distance_as_weight(G):
        updated = G.copy()
        for edge in updated.edges_iter(data=True):
            distance = edge[2]['distance']
            edge[2]['weight'] = distance
        return updated

    @staticmethod
    def inverted_distance_as_weight(G, dist_exp=2):
        updated = G.copy()
        for edge in updated.edges_iter(data=True):
            distance = edge[2]['distance']
            weight = 1 / ((distance + 1) ** dist_exp)
            edge[2]['weight'] = weight
        return updated

    @staticmethod
    def relevance_as_weight(G, dist_exp=2):
        updated = G.copy()
        for edge in updated.edges_iter(data=True):
            distance = edge[2]['distance']
            overlap_ratio = edge[2]['overlap_ratio']
            weight = overlap_ratio * (1 / ((distance + 1) ** dist_exp))
            edge[2]['weight'] = weight

        return updated
