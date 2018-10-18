from model import *
from collections import OrderedDict
import matplotlib.pyplot as plt
import networkx as nx


def spit_horizontally(sheet_graph):
    copy_graph = sheet_graph.copy()
    to_remove = list()
    for edge in sheet_graph.edges_iter(data=True):
        direction = edge[2]['direction']
        distance = edge[2]['distance']
        if (direction == 'LEFT' or direction == 'RIGHT') and distance > 0:
            to_remove.append((edge[0], edge[1]))
    copy_graph.remove_edges_from(to_remove)
    return copy_graph


def split_vertically(sheet_graph, reject_single_col=False):
    # distinguish header nodes from the rest
    headers, all_nodes = dict(), dict()
    for node, data in sheet_graph.nodes_iter(data=True):
        label = data['label']
        if label == 'Header':
            # get the region coordinates
            headers[node] = Utils.get_node_region_coordinates(data)
            all_nodes[node] = Utils.get_node_region_coordinates(data)
        else:
            all_nodes[node] = Utils.get_node_region_coordinates(data)

    # sort header nodes in descending order of their maximum row, followed by ascending order of minimum row
    ordered_headers = OrderedDict(sorted(headers.items(), key=lambda item: (-item[1]['ymax'], item[1]['ymin'])))
    # sort graph nodes in descending order of their minimum row
    ordered_all = OrderedDict(sorted(all_nodes.items(), key=lambda item: (-item[1]['ymin'])))

    # identify tables arranged vertically in each connected component of the graph
    proposed_tables = list()
    undetermined = list()
    for comp in nx.connected_components(sheet_graph.to_undirected()):
        comp_headers = OrderedDict(item for item in ordered_headers.items() if item[0] in comp)
        comp_all = OrderedDict(item for item in ordered_all.items() if item[0] in comp)
        frontier = set()  # Holds node set of last proposed table

        # analyze each header node and determine if it can form a table
        # with the other nodes from the connected component
        for hnode in list(comp_headers.keys()):
            if hnode in frontier:
                continue
            else:
                # if header node is not paired already
                hcoord = comp_headers[hnode]
                candidate_nodes = set()
                for n in comp_all.keys():
                    # collect nodes starting at the same row or lower than the header
                    # this includes the header node itself
                    if hcoord['ymin'] <= comp_all[n]['ymin']:
                        candidate_nodes.add(n)
                    else:
                        break

                # a header node can form a table only with nodes that are connected to it
                connected_nodes = set()
                temp_graph = sheet_graph.subgraph(candidate_nodes)
                for node_set in nx.connected_components(temp_graph.to_undirected()):
                    if hnode in node_set:
                        connected_nodes.update(node_set)
                        break

                # identify header nodes at the left or right of the selected header
                # bring them all together to form a group header
                header_group = set()
                connected_headers = (set(comp_headers.keys()) & connected_nodes)
                for ch in connected_headers:
                    if (comp_all[ch]['ymin'] < hcoord['ymax']) and (comp_all[ch]['ymin'] >= hcoord['ymin']):
                        header_group.add(ch)

                # measure the alignment of the header group with the other nodes
                hgroup_columns, other_columns = set(), set()
                for cn in connected_nodes:
                    coord = comp_all[cn]
                    if cn in header_group:
                        hgroup_columns.update([j for j in range(coord['xmin'], coord['xmax'])])
                    else:
                        other_columns.update([j for j in range(coord['xmin'], coord['xmax'])])

                # the alignment ratio must be more than the threshold
                common_columns = hgroup_columns & other_columns
                alignment_ratio = len(common_columns) / len(other_columns) if len(other_columns) > 0 else 0.0
                if alignment_ratio >= 0.5 and len(hgroup_columns) > 1:
                    # save the previous candidate table's node set
                    # it is considered complete by this point
                    if len(frontier) > 0:
                        if reject_single_col:
                            frontier_headers = (set(ordered_headers.keys()) & frontier)
                            fh_cols = set()
                            for fh in frontier_headers:
                                fh_cols.update([j for j in
                                                range(ordered_headers[fh]['xmin'], ordered_headers[fh]['xmax'])])
                            if len(fh_cols) > 1:
                                proposed_tables.append(frontier)
                            else:
                                undetermined.extend(frontier)
                        else:
                            proposed_tables.append(frontier)

                    frontier = set(connected_nodes)
                    for cn in connected_nodes:
                        del comp_all[cn]
                else:
                    # if only headers, try to append them to the node set
                    # of the last identified table (from this component)
                    if len(frontier) > 0 and len(other_columns) == 0:
                            frnt_cols = set()
                            for fn in frontier:
                                frnt_cols.update([j for j in range(ordered_all[fn]['xmin'], ordered_all[fn]['xmax'])])
                            # append only those header that overlap horizontally with the frontier
                            # which means they have columns in common
                            if len(hgroup_columns & frnt_cols) > 0:
                                frontier.update(connected_headers)
                                for cnode in connected_headers:
                                    del comp_all[cnode]

                    # hgroup_rows, other_rows = set(), set()
                    # for cn in connected_nodes:
                    #     coord = comp_all[cn]
                    #     if cn in header_group:
                    #         hgroup_rows.update([j for j in range(coord['ymin'], coord['ymax'])])
                    #     else:
                    #         other_rows.update([j for j in range(coord['ymin'], coord['ymax'])])
                    #
                    # max_hgrow = max(hgroup_rows)
                    # max_orow = max(other_rows) if len(other_rows) > 1 else -1.0
                    #
                    # if len(frontier) > 0 and max_hgrow >= max_orow and len(hgroup_columns) > len(other_columns):
                    #     frnt_cols = set()
                    #     for fn in frontier:
                    #         frnt_cols.update([j for j in range(ordered_all[fn]['xmin'], ordered_all[fn]['xmax'])])
                    #     # append only those header that overlap horizontally with the frontier
                    #     # which means they have columns in common
                    #
                    #     #for ch in connected_headers:
                    #     #    ch_cols = set([j for j in range(comp_all[ch]['xmin'], comp_all[ch]['xmax'])])
                    #     if len((hgroup_columns | other_columns) & frnt_cols) > 0:
                    #         frontier.update(connected_nodes)
                    #         for cnode in connected_nodes:
                    #             del comp_all[cnode]

        # in case there are unpaired nodes (i.e., not part of a proposed table)
        if len(comp_all) > 0:
            undetermined.extend(list(comp_all.keys()))
        # append the last frontier
        if len(frontier) > 0:
            if reject_single_col:
                frontier_headers = (set(comp_headers.keys()) & frontier)
                fh_cols = set()
                for fh in frontier_headers:
                    fh_cols.update([j for j in range(comp_headers[fh]['xmin'], comp_headers[fh]['xmax'])])
                if len(fh_cols) > 1:
                    proposed_tables.append(frontier)
                else:
                    undetermined.extend(frontier)
            else:
                proposed_tables.append(frontier)

    return proposed_tables, undetermined


def pair_left_right(sheet_graph, candidate_tables, unpaired, dist_th=-1):
    paired = set(sheet_graph.nodes()) - set(unpaired)
    for u in list(unpaired):
        neighbors = set(sheet_graph[u].keys())
        if len(neighbors & paired) > 0:
            neighboring_tables = dict()
            for tix in range(0, len(candidate_tables)):
                intersection = candidate_tables[tix] & neighbors
                if len(intersection) > 0:
                    for v in intersection:
                        edge = sheet_graph[u][v][0]
                        if (edge['direction'] == 'LEFT') or (edge['direction'] == 'RIGHT'):
                            if edge['distance'] not in neighboring_tables.keys():
                                neighboring_tables[edge['distance']] = set()
                            neighboring_tables[edge['distance']].add(tix)

            if len(neighboring_tables) > 0:
                min_dist = min(list(neighboring_tables.keys()))
                if len(neighboring_tables[min_dist]) == 1 and (dist_th == -1 or min_dist <= dist_th):
                    tix = neighboring_tables[min_dist].pop()
                    candidate_tables[tix].add(u)
                    unpaired.remove(u)

    return candidate_tables, unpaired


def handle_overlapping_tables(sheet_graph, candidate_tables):

    nodes_coords = dict()
    nodes_labels = dict()
    for node, data in sheet_graph.nodes_iter(data=True):
        nodes_labels[node] = data['label']
        nodes_coords[node] = Utils.get_node_region_coordinates(data)

    header_rows, header_cols, other_rows, other_cols = list(), list(), list(), list()
    for table in candidate_tables:
        hcols, hrows, ocols, orows = set(), set(), set(), set()
        for node in table:
            coord = nodes_coords[node]
            # TODO: Consider only valid header
            if nodes_labels[node] == 'Header':
                hcols.update([x for x in range(coord['xmin'], coord['xmax'])])
                hrows.update([y for y in range(coord['ymin'], coord['ymax'])])
            else:
                ocols.update([x for x in range(coord['xmin'], coord['xmax'])])
                orows.update([y for y in range(coord['ymin'], coord['ymax'])])
        header_cols.append(hcols), header_rows.append(hrows)
        other_cols.append(ocols), other_rows.append(orows)

    overlap_indices = set()
    to_merge = list()
    undecided = list()
    updated_candidates = list()
    for i in range(0, len(candidate_tables)-1):
        t1_rows, t1_cols = (header_rows[i] | other_rows[i]), (header_cols[i] | other_cols[i])
        for j in range(i+1, len(candidate_tables)):
            t2_rows, t2_cols = (header_rows[j] | other_rows[j]), (header_cols[j] | other_cols[j])
            if len(t1_rows & t2_rows) > 0 and len(t1_cols & t2_cols) > 0:
                overlap_indices.update({i, j})
                comb_header_rows = (header_rows[i] | header_rows[j])
                comb_other_rows = (other_rows[i] | other_rows[j])
                comb_other_above = set(y for y in comb_other_rows if y <= max(comb_header_rows))
                if len(comb_other_above - comb_header_rows) == 0:
                    to_merge.append({i, j})
                else:
                    top_tix, bottom_tix = i, j
                    if min(header_rows[j]) < min(header_rows[i]):
                        top_tix, bottom_tix = j, i
                    cut_row = min(header_rows[bottom_tix])

                    top_cand = set(candidate_tables[top_tix])
                    for node in candidate_tables[top_tix]:
                        coord = nodes_coords[node]
                        if coord['ymin'] >= cut_row:
                            undecided.append(node)
                            top_cand.remove(node)
                        elif coord['ymin'] < cut_row < coord['ymax']:
                            # TODO: Split the node and update the graph
                            print("Split top for %s" % sheet_graph.graph)
                            diff = cut_row - coord['ymin']
                            if diff < ((coord['ymax'] - coord['ymin']) * 0.5):
                                undecided.append(node)
                                top_cand.remove(node)

                    bottom_cand = set(candidate_tables[bottom_tix])
                    for node in candidate_tables[bottom_tix]:
                        coord = nodes_coords[node]
                        if coord['ymax'] <= cut_row:
                            undecided.append(node)
                            bottom_cand.remove(node)
                        elif coord['ymin'] < cut_row < coord['ymax']:
                            # TODO: Split the node and update the graph
                            print("Split bottom for %s" % sheet_graph.graph)
                            diff = cut_row - coord['ymin']
                            if diff < ((coord['ymax'] - coord['ymin']) * 0.5):
                                undecided.append(node)
                                bottom_cand.remove(node)

                    # add updated candidates
                    updated_candidates.append(top_cand)
                    updated_candidates.append(bottom_cand)

    # collect tables that do not overlap with any other table
    not_overlap = list(candidate_tables[tix] for tix in range(0, len(candidate_tables)) if tix not in overlap_indices)
    updated_candidates.extend(not_overlap)  # add non overlapping

    # identify merged chains
    if len(to_merge) > 1:
        updated = to_merge[1:]
        for i in range(0, len(to_merge)):
            temp = list()
            merged = to_merge[i]
            for j in range(0, len(updated)):
                pair_j = updated[j]
                if len(merged & pair_j) > 0:
                    merged = merged | pair_j
                else:
                    temp.append(pair_j)
            temp.append(merged)
            updated = list(temp)
        to_merge = updated

    # add merged candidates
    for entry in to_merge:
        merged_candidate = set()
        for tix in entry:
            merged_candidate.update(candidate_tables[tix])
            updated_candidates.append(merged_candidate)
    return updated_candidates, undecided


def display_evaluation_results(scores_per_sheet, worksheets, print_threshold=0.9):
    # prepare evaluation results for display
    all_scores = list()
    all_multi_scores = list()
    all_single_scores = list()

    start_threshold = 80
    increase_by = 5
    end_threshold = 100
    score_stats = dict()

    print("\n")
    for sheet_index in scores_per_sheet.keys():
        scores = scores_per_sheet[sheet_index]
        is_multi_table = False
        if len(scores.keys()) > 1:
            if -1 in scores.keys():
                if len(scores.keys()) > 2:
                    is_multi_table = True
            else:
                is_multi_table = True

        for table_id, table_score in scores.items():
            if table_id != -1:
                all_scores.append(table_score)

                if table_score < print_threshold:
                    print("index=%d, info=\"%s\", tableId=%d, score=%f" % (
                        sheet_index, worksheets[sheet_index].getId(), table_id, table_score))

                if is_multi_table:
                    all_multi_scores.append(table_score)
                else:
                    all_single_scores.append(table_score)

                for t in range(start_threshold, end_threshold, increase_by):

                    threshold = t / 100
                    if threshold not in score_stats.keys():
                        score_stats[threshold] = list([0, 0, 0])

                    if table_score < threshold:
                        score_stats[threshold][0] += 1
                        if is_multi_table:
                            score_stats[threshold][1] += 1
                        else:
                            score_stats[threshold][2] += 1

    sum_scores = sum(all_scores)
    nr_of_tables = len(all_scores)

    multi_sum_scores = sum(all_multi_scores)
    multi_nr_of_tables = len(all_multi_scores)

    single_sum_scores = sum(all_single_scores)
    single_nr_of_tables = len(all_single_scores)

    # display scores box plot
    print("\n")
    data = list()
    data.append(all_scores)
    data.append(all_single_scores)
    data.append(all_multi_scores)

    results = plt.boxplot(data, labels=["All", "Single", "Multi"], meanline=True, showmeans=True)
    plt.show()

    print("\nSum scores = %f " % sum_scores)
    print("Number of tables = %f " % nr_of_tables)
    print("Median Score = %f" % results['medians'][0].get_ydata(orig=True)[0])
    print("Average Score = %f" % (sum_scores / nr_of_tables))
    print("Lower Whisker = %f" % min(results['whiskers'][0].get_ydata(orig=True)))
    print("Min Score = %f" % min(all_scores))

    if single_nr_of_tables > 0:
        print("\nSingle Sum scores = %f " % single_sum_scores)
        print("Single number of tables = %f " % single_nr_of_tables)
        print("Single Median Score = %f" % results['medians'][1].get_ydata(orig=True)[0])
        print("Single Average Score = %f" % (single_sum_scores / single_nr_of_tables))
        print("Single Lower Whisker = %f" % min(results['whiskers'][2].get_ydata(orig=True)))
        print("Single Min Score = %f" % min(all_single_scores))

    if multi_nr_of_tables > 0:
        print("\nMulti Sum scores = %f " % multi_sum_scores)
        print("Multi number of tables = %f " % multi_nr_of_tables)
        print("Multi Median Score = %f" % results['medians'][2].get_ydata(orig=True)[0])
        print("Multi Average Score = %f" % (multi_sum_scores / multi_nr_of_tables))
        print("Multi Lower Whisker = %f" % min(results['whiskers'][4].get_ydata(orig=True)))
        print("Multi Min Score = %f" % min(all_multi_scores))

    for key in score_stats.keys():
        print("\nStats for threshold: %.2f" % key)
        all_above_th = (100 - ((score_stats[key][0] / nr_of_tables) * 100))
        print("Tables above threshold: %.2f%%" % all_above_th)
        multi_above_th = -1.0
        single_above_th = -1.0
        if multi_nr_of_tables > 0:
            multi_above_th = (100 - ((score_stats[key][1] / multi_nr_of_tables) * 100))
            print("Multi tables above threshold: %.2f%%" % multi_above_th)
        if single_nr_of_tables > 0:
            single_above_th = (100 - ((score_stats[key][2] / single_nr_of_tables) * 100))
            print("Single tables above threshold: %.2f%%" % single_above_th)
        print("%.2f, %.2f, %.2f" % (all_above_th, single_above_th, multi_above_th))


def summarize_evaluation_results(scores_per_sheet):
    # prepare evaluation results for display
    all_scores = list()
    all_multi_scores = list()
    all_single_scores = list()

    start_threshold = 80
    increase_by = 5
    end_threshold = 100
    score_stats = dict()

    for sheet_index in scores_per_sheet.keys():
        scores = scores_per_sheet[sheet_index]
        is_multi_table = False
        if len(scores.keys()) > 1:
            if -1 in scores.keys():
                if len(scores.keys()) > 2:
                    is_multi_table = True
            else:
                is_multi_table = True

        for table_id, table_score in scores.items():
            if table_id != -1:
                all_scores.append(table_score)
                if is_multi_table:
                    all_multi_scores.append(table_score)
                else:
                    all_single_scores.append(table_score)

                for t in range(start_threshold, end_threshold, increase_by):
                    threshold = t / 100
                    if threshold not in score_stats.keys():
                        score_stats[threshold] = list([0, 0, 0])
                    if table_score < threshold:
                        score_stats[threshold][0] += 1
                        if is_multi_table:
                            score_stats[threshold][1] += 1
                        else:
                            score_stats[threshold][2] += 1

    nr_of_tables = len(all_scores)
    multi_nr_of_tables = len(all_multi_scores)
    single_nr_of_tables = len(all_single_scores)

    results_per_threshold = dict()
    for key in score_stats.keys():
        all_above_th = (100 - ((score_stats[key][0] / nr_of_tables) * 100))
        multi_above_th = -1.0
        single_above_th = -1.0
        if multi_nr_of_tables > 0:
            multi_above_th = (100 - ((score_stats[key][1] / multi_nr_of_tables) * 100))
        if single_nr_of_tables > 0:
            single_above_th = (100 - ((score_stats[key][2] / single_nr_of_tables) * 100))

        results_per_threshold[key] = list([all_above_th, single_above_th, multi_above_th])
    return results_per_threshold


def identify(graph, reject_single_col=False, handle_overlaps=False, dist_th=-1):
    # handle horizontally aligned tables
    horizontal_graph = spit_horizontally(graph)

    # handle vertically aligned tables
    predicted, undecided = split_vertically(horizontal_graph, reject_single_col)
    # correct discrepancies resulting out of overlapping tables
    if handle_overlaps:
        predicted, cut_out = handle_overlapping_tables(graph, predicted)
        undecided.extend(cut_out)
    # pair with the nearest left or right table
    predicted, undecided = pair_left_right(graph, predicted, undecided, dist_th)

    # graph_of_undecided = graph.subgraph(undecided)
    # undecided_fragments = list(nx.connected_components(graph_of_undecided.to_undirected()))

    undecided_fragments = [{node} for node in undecided]
    return predicted, undecided_fragments


def identify_and_evaluate(graphs, reject_single_col=False, handle_overlaps=False, dist_th=-1):
    scores_per_graph = dict()
    for gix in range(0, len(graphs)):

        # handle horizontally aligned tables
        horizontal_graph = spit_horizontally(graphs[gix])

        # handle vertically aligned tables
        predicted, undecided = split_vertically(horizontal_graph, reject_single_col)
        # correct discrepancies resulting out of overlapping tables
        if handle_overlaps:
            predicted, cut_out = handle_overlapping_tables(graphs[gix], predicted)
            undecided.extend(cut_out)
        # pair with the nearest left or right table
        predicted, undecided = pair_left_right(graphs[gix], predicted, undecided, dist_th)

        # evaluate the predicted tables
        max_jaccard_per_table = dict()
        predicted_fragments = list(predicted)
        graph_of_undecided = graphs[gix].subgraph(undecided)
        undecided_fragments = list(nx.connected_components(graph_of_undecided.to_undirected()))
        # predicted_fragments.extend(undecided_fragments)

        tables_node_areas = Utils.group_node_area_tuples_by_tableId(graphs[gix])
        table_ids = set(tables_node_areas.keys())

        if len(predicted_fragments) > 0:
            for frag in predicted_fragments:
                jaccard_per_table = Utils.evaluate_identified_table(graphs[gix], graphs[gix].subgraph(frag))
                for tid in jaccard_per_table.keys():
                    if tid != -1:
                        if tid in max_jaccard_per_table.keys():
                            if max_jaccard_per_table[tid] < jaccard_per_table[tid]:
                                max_jaccard_per_table[tid] = jaccard_per_table[tid]
                        else:
                            max_jaccard_per_table[tid] = jaccard_per_table[tid]
        else:
            for tbl_id in table_ids:
                max_jaccard_per_table[tbl_id] = 0.0

        if (-1 in table_ids) and undecided_fragments:
            for frag in undecided_fragments:
                jaccard_per_table = Utils.evaluate_identified_table(graphs[gix], graphs[gix].subgraph(frag))
                if -1 in jaccard_per_table.keys():
                    if -1 in max_jaccard_per_table.keys():
                        if max_jaccard_per_table[-1] < jaccard_per_table[-1]:
                            max_jaccard_per_table[-1] = jaccard_per_table[-1]
                    else:
                        max_jaccard_per_table[-1] = jaccard_per_table[-1]

        # for tables that can not be mapped to any of the predicted fragments
        remaining = table_ids - set(max_jaccard_per_table.keys())
        for tbl_id in remaining:
            max_jaccard_per_table[tbl_id] = 0.0

        scores_per_graph[gix] = max_jaccard_per_table

    return scores_per_graph


def inspect(index, graph, reject_single_col=False, handle_overlaps=False, dist_th=-1, visualize_steps=False):
    # the original graph
    print("%d. %s" % (index, graph.graph))
    if visualize_steps:
        print("\nDraw original graph")
        Visualizer.draw_worksheet_graph(graph)

    # handle horizontally aligned tables
    horz_graph = spit_horizontally(graph)
    if visualize_steps:
        print("Draw horizontal graph")
        Visualizer.draw_worksheet_graph(horz_graph)

    # handle vertically aligned tables
    pt, un = split_vertically(horz_graph, reject_single_col)
    fragments = list(pt)
    fragments.extend([{entry} for entry in un])
    if visualize_steps:
        print("Draw vertical graph")
        Visualizer.draw_fragmented_graph(horz_graph, fragments)

    if handle_overlaps:
        # correct discrepancies resulting out of overlapping tables
        pt, cut = handle_overlapping_tables(graph, pt)
        un.extend(cut)
        fragments = list(pt)
        fragments.extend([{entry} for entry in un])
        if visualize_steps:
            print("Draw graph after overlaps")
            Visualizer.draw_fragmented_graph(graph, fragments)

    # attempt to pair remaining with the nearest left/right (identified) table
    pt, un = pair_left_right(graph, pt, un, dist_th)
    fragments = list(pt)
    fragments.extend([{entry} for entry in un])
    if visualize_steps:
        print("Draw graph after distance pairing")
        Visualizer.draw_fragmented_graph(graph, fragments)

    # evaluate the predicted tables
    tables_node_areas = Utils.group_node_area_tuples_by_tableId(graph)
    table_ids = set(tables_node_areas.keys())

    # group the fragments by type
    predicted_fragments = list(pt)
    graph_of_undecided = graph.subgraph(un)
    undecided_fragments = list(nx.connected_components(graph_of_undecided.to_undirected()))

    max_jaccard_per_table = dict()
    if len(predicted_fragments) > 0:
        for frag in predicted_fragments:
            jaccard_per_table = Utils.evaluate_identified_table(graph, graph.subgraph(frag))
            for tid in jaccard_per_table.keys():
                if tid != -1:
                    if tid in max_jaccard_per_table.keys():
                        if max_jaccard_per_table[tid] < jaccard_per_table[tid]:
                            max_jaccard_per_table[tid] = jaccard_per_table[tid]
                    else:
                        max_jaccard_per_table[tid] = jaccard_per_table[tid]
    else:
        for tbl_id in table_ids:
            max_jaccard_per_table[tbl_id] = 0.0

    if (-1 in table_ids) and undecided_fragments:
        for frag in undecided_fragments:
            jaccard_per_table = Utils.evaluate_identified_table(graph, graph.subgraph(frag))
            if -1 in jaccard_per_table.keys():
                if -1 in max_jaccard_per_table.keys():
                    if max_jaccard_per_table[-1] < jaccard_per_table[-1]:
                        max_jaccard_per_table[-1] = jaccard_per_table[-1]
                else:
                    max_jaccard_per_table[-1] = jaccard_per_table[-1]

    # for tables that can not be mapped to any of the predicted fragments
    remaining = table_ids - set(max_jaccard_per_table.keys())
    for tbl_id in remaining:
        max_jaccard_per_table[tbl_id] = 0.0
    print(max_jaccard_per_table)


if __name__ == '__main__':

    gold_standard_heuristics = "./datasets/gold_standard_heuristics.json"
    selected_path = gold_standard_heuristics
    sheets = Reader.read_json_dataset(selected_path)

    sheet_graphs = list()
    for sh_ix in range(len(sheets)):
        # create graphs
        multi_digraph = Builder.create_graph_for_sheet(sheets[sh_ix])
        # remove in-view edges
        multi_digraph = Utils.remove_inview_edges(multi_digraph, is_strict=False)
        # store graph for later use
        sheet_graphs.append(multi_digraph)
    print("Number of graphs in the dataset: %s" % len(sheet_graphs))

    scores_per_sheet_graph = identify_and_evaluate(sheet_graphs, reject_single_col=True,
                                                   handle_overlaps=True, dist_th=2)
    display_evaluation_results(scores_per_sheet_graph, sheets, print_threshold=0.9)
