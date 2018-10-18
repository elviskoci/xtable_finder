from model import *
from collections import namedtuple
import random as rand
import math

# A Bucket contains a list of sheet graphs
Bucket = namedtuple('Bucket', 'sheet_graphs count_tables count_multi')


def distribute_graphs_to_buckets(sheet_graphs, nr_of_buckets, seed=-1):

    if seed > -1:
        rand.seed(seed)

    # for each bucket collect sheet-graphs and other info
    buckets_graphs = [list() for _ in range(nr_of_buckets)]
    buckets_table_counts = [0] * nr_of_buckets
    buckets_multi_counts = [0] * nr_of_buckets

    for shg in sheet_graphs:
        bix = rand.randint(0, nr_of_buckets-1)

        # append the sheet graph to the selected bucket
        buckets_graphs[bix].append(shg)

        # get the number of tables for this sheet graph
        tables_nodes_areas = Utils.group_node_area_tuples_by_tableId(shg)
        table_ids = list(tables_nodes_areas.keys())
        if -1 in table_ids:  # -1 is a reserved special value for independent regions
            table_ids.remove(-1)  # -1 stands for "no table"
        buckets_table_counts[bix] += len(table_ids)
        if len(table_ids) > 1:  # count sheet graphs containing multiple table
            buckets_multi_counts[bix] += 1
    return [Bucket(buckets_graphs[j], buckets_table_counts[j], buckets_multi_counts[j]) for j in range(nr_of_buckets)]


# This function distributes the sheet graphs into buckets.
# it first groups the sheets graphs by file name, and then
# it randomly distributes the groups, making sure that
# the number of items per bucket is uniform
def distribute_graphs_to_buckets_by_file(sheet_graphs, nr_of_buckets, balance_multis=False, seed=-1):

    if seed > -1:
        rand.seed(seed)

    # sheets from the same file tend to have similar characteristics
    graphs_by_file = dict()  # group graphs by file name
    files_table_counts = dict()  # record the number of tables per file
    files_multis_counts = dict()  # record the number of sheets with multiple tables per file
    for shg in sheet_graphs:
        key_ = (shg.graph['file'], shg.graph['corpus'])
        if key_ not in graphs_by_file.keys():
            graphs_by_file[key_] = list()
            files_table_counts[key_] = 0
        # append to group
        graphs_by_file[key_].append(shg)
        # get the number of tables for this sheet graph
        tables_nodes_areas = Utils.group_node_area_tuples_by_tableId(shg)
        table_ids = list(tables_nodes_areas.keys())
        if -1 in table_ids:  # -1 is a reserved special value for independent regions
            table_ids.remove(-1)  # -1 stands for "no table"
        files_table_counts[key_] += len(table_ids)
        if len(table_ids) > 1:  # count sheet graphs containing multiple table
            if key_ not in files_multis_counts.keys():
                files_multis_counts[key_] = 0
            files_multis_counts[key_] += 1

    # for each bucket collect sheet-graphs and other info
    buckets_graphs = [list() for _ in range(nr_of_buckets)]
    buckets_table_counts = [0] * nr_of_buckets
    buckets_multi_counts = [0] * nr_of_buckets

    # we handle the multi-table sheets first
    # to ensure that they are distributed evenly
    if balance_multis:
        # calculate the lower bound threshold
        # all buckets must have at least this many multi-table sheets
        total_multi_sheets = sum(files_multis_counts.values())
        mth = math.floor(total_multi_sheets/nr_of_buckets)

        for key_, n_multis in sorted(files_multis_counts.items(), key=lambda kv: kv[1], reverse=True):
            indices = list(range(0, nr_of_buckets))
            bix = rand.choice(indices)  # randomly choose one of the buckets
            # if bucket is full choose another one
            while buckets_multi_counts[bix] >= mth:
                indices.remove(bix)
                if not indices:  # all buckets achieved the lower bound
                    indices = list(range(0, nr_of_buckets))
                    # update the threshold, to reflect the upper bound
                    mth = math.ceil(total_multi_sheets/nr_of_buckets)
                bix = rand.choice(indices)

            buckets_table_counts[bix] += files_table_counts[key_]
            buckets_multi_counts[bix] += n_multis
            buckets_graphs[bix].extend(graphs_by_file[key_])
            del graphs_by_file[key_]

    # determine the maximum number of sheets graphs per bucket
    # first initialize this threshold using the lower bound (i.e., floor)
    gth = math.floor(len(sheet_graphs) / nr_of_buckets)
    for key_, graphs_list in sorted(graphs_by_file.items(),  key=lambda kv: len(kv[1]), reverse=True):
        indices = list(range(0, nr_of_buckets))
        bix = rand.choice(indices)  # randomly choose one of the buckets

        # if bucket is full choose another one
        while len(buckets_graphs[bix]) >= gth:
            indices.remove(bix)
            if not indices:  # all buckets achieved the lower bound
                indices = list(range(0, nr_of_buckets))
                # update the threshold, to reflect the upper bound
                gth = math.ceil(len(sheet_graphs) / nr_of_buckets)
            bix = rand.choice(indices)

        buckets_multi_counts[bix] += files_multis_counts[key_] if key_ in files_multis_counts.keys() else 0
        buckets_table_counts[bix] += files_table_counts[key_]
        buckets_graphs[bix].extend(graphs_list)
    return [Bucket(buckets_graphs[j], buckets_table_counts[j], buckets_multi_counts[j]) for j in range(nr_of_buckets)]


if __name__ == '__main__':
    test_dataset_path = "C:/Users/Elvis/Documents/Python Scripts/sptblrec/evolve/datasets/gold_standard_graphs.json"
    error_dataset_path = "C:/Users/Elvis/Documents/Python Scripts/sptblrec/evolve/datasets/graph_dataset_error_1p.json"

    path = test_dataset_path
    sheets = Reader.read_json_dataset(path)

    worksheet_graphs = list()
    for i in range(len(sheets)):
        # create graph using one of the weight functions
        multi_digraph = Builder.create_graph_for_sheet(sheets[i], weight_function=Builder.distance_as_weight)
        # remove in-view edges
        multi_digraph = Utils.remove_inview_edges(multi_digraph, is_strict=False)
        # store graph for later use
        worksheet_graphs.append(multi_digraph)

    graph_buckets_rand = distribute_graphs_to_buckets_by_file(worksheet_graphs, 10, balance_multis=True, seed=3)

    for gb in graph_buckets_rand:
        print("\n\nn_graphs=%d, n_tables=%d, n_multi=%d" % (len(gb.sheet_graphs), gb.count_tables, gb.count_multi))
        for grf in gb.sheet_graphs:
            print(str(grf.graph))
