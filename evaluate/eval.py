import matplotlib.pyplot as plt
import networkx as nx

from model import Utils


def get_max_jaccard_per_target_fragment(graph, predicted_fragments):
    max_jaccard_per_target = dict()

    nodes_labels = nx.get_node_attributes(graph, 'label')
    data_nodes, head_nodes = set(), set()
    for node, label in nodes_labels.items():
        if label == 'Header':
            head_nodes.add(node)
        else:
            data_nodes.add(node)

    # if it contains a header and a data, than it is a candidate table
    candidate_tables, other_fragments = list(), list()
    for frag in predicted_fragments:
        if (frag & head_nodes) and (frag & data_nodes):
            candidate_tables.append(frag)
        else:
            other_fragments.append(frag)

    tables_node_areas = Utils.group_node_area_tuples_by_tableId(graph)
    table_ids = set(tables_node_areas.keys())
    if len(candidate_tables) > 0:
        for frag in candidate_tables:
            jaccard_per_fragment = Utils.evaluate_identified_table(graph, graph.subgraph(frag))
            for tbl_id in jaccard_per_fragment.keys():
                if tbl_id != -1:  # -1 is a special id reserved for regions outside the annotated tables
                    if tbl_id in max_jaccard_per_target.keys():
                        if max_jaccard_per_target[tbl_id] < jaccard_per_fragment[tbl_id]:
                            max_jaccard_per_target[tbl_id] = jaccard_per_fragment[tbl_id]
                    else:
                        max_jaccard_per_target[tbl_id] = jaccard_per_fragment[tbl_id]
    else:
        for tbl_id in table_ids:
            max_jaccard_per_target[tbl_id] = 0.0

    if (-1 in table_ids) and other_fragments:
        for frag in other_fragments:
            jaccard_per_fragment = Utils.evaluate_identified_table(graph, graph.subgraph(frag))
            if -1 in jaccard_per_fragment.keys():
                if -1 in max_jaccard_per_target.keys():
                    if max_jaccard_per_target[-1] < jaccard_per_fragment[-1]:
                        max_jaccard_per_target[-1] = jaccard_per_fragment[-1]
                else:
                    max_jaccard_per_target[-1] = jaccard_per_fragment[-1]

    # when there are tables that do not overlap with any of the predicted candidate_tables
    remaining = table_ids - set(max_jaccard_per_target.keys())
    for tbl_id in remaining:
        max_jaccard_per_target[tbl_id] = 0.0
    return max_jaccard_per_target


def display_evaluation_results(scores_per_sheet, sheet_graphs, print_threshold=0.9):
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

                if table_score < print_threshold:  # (start_threshold / 100)
                    print("index=%d, info=\"%s\", tableId=%d, score=%f" % (
                        sheet_index, str((sheet_graphs[sheet_index].graph['file'],
                                         sheet_graphs[sheet_index].graph['sheet'])), table_id, table_score))

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
