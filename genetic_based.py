import numpy as np
# from sklearn.preprocessing import RobustScaler, MinMaxScaler

from evaluate import eval
from evaluate.bucket import distribute_graphs_to_buckets_by_file
from model import Reader, Builder, Utils
from optimz.tuner import get_sample_for_optimization, get_optimal_weights
from optimz import Metrics
from evolve import reduced_search


def calculate_error_rate(tbl_values, alt_values, opt_weights):
    # multiply the values with the weights
    multi_wt = np.multiply(opt_weights, tbl_values)
    multi_wa = np.multiply(opt_weights, alt_values)
    # calculate weighted sum to determine the fitness of each instance
    sum_wt = np.sum(multi_wt, axis=1)
    sum_wa = np.sum(multi_wa, axis=1)
    differences = sum_wt - sum_wa

    bad_points = [val for val in differences if val >= 0]
    return len(bad_points) / len(differences)


def fit_transformer(transformer, graphs, metrics_function, sample_per_sheet=200, sample_per_table=10):
    t_measurements, p_measurements, t_counts, g_indices = \
        get_sample_for_optimization(graphs, metrics_function, max_n_frgm=sample_per_sheet,
                                    n_frgm_per_table=sample_per_table)

    all_measurements = list()
    for k in range(len(t_measurements)):
        all_measurements.append(list(t_measurements[k].values()))
        all_measurements.append(list(p_measurements[k].values()))

    return transformer.fit(np.array(all_measurements))


def main(test_dataset_path, train_dataset_path, metrics_calc_fun, weight_bounds, n_buckets=10, seed=1,
         n_opt_runs=10, n_searches=10, max_n_frgm=200, n_frgm_per_table=10, tf=None):

    # read json datasets
    print("\nUsing training dataset = %s" % train_dataset_path)
    print("Using test dataset = %s" % test_dataset_path)
    test_sheets = Reader.read_json_dataset(test_dataset_path)
    error_sheets = Reader.read_json_dataset(train_dataset_path)
    test_graphs = list()
    error_graphs = list()
    for i in range(len(error_sheets)):
        # create graphs
        er_multi_digraph = Builder.create_graph_for_sheet(error_sheets[i])
        ts_multi_digraph = Builder.create_graph_for_sheet(test_sheets[i])
        # remove in-view edges
        er_multi_digraph = Utils.remove_inview_edges(er_multi_digraph, is_strict=False)
        ts_multi_digraph = Utils.remove_inview_edges(ts_multi_digraph, is_strict=False)
        # store graph for later use
        error_graphs.append(er_multi_digraph)
        test_graphs.append(ts_multi_digraph)

    # distribute the graphs into buckets for the cross validation
    print("Distribute into %d buckets using seed=%d" % (n_buckets, seed))
    error_buckets = distribute_graphs_to_buckets_by_file(error_graphs, n_buckets, balance_multis=True, seed=seed)
    test_buckets = distribute_graphs_to_buckets_by_file(test_graphs, n_buckets, balance_multis=True, seed=seed)

    print("Execute %d times the genetic search" % n_searches)
    print("Perform %d weight optimization rounds" % n_opt_runs)
    # fit the transformer for scaling the metrics' values
    print("Using max_n_frag=%d, n_frag_per_table=%d" % (max_n_frgm, n_frgm_per_table))
    if not tf:
        tf_fit = None
        print("Working with unscaled data")
    else:
        tf_fit = fit_transformer(tf, error_graphs, metrics_calc_fun, max_n_frgm * n_opt_runs,
                                 (n_frgm_per_table * n_opt_runs) if n_frgm_per_table > 0 else -1)
        print("Scale the data using '%s'" % tf.__class__.__name__)

    # collect results
    print("Using the function '%s' to calculate fragmentation metrics" % str(metrics_calc_fun.__name__))

    scores_per_graph = list()
    target_scores, predicted_scores = list(), list()
    # for each fold train and test
    for j in range(len(error_buckets)):
        print("\n\nFold Nr.%d" % j)
        train_buckets = list(error_buckets[:j])
        train_buckets.extend(error_buckets[j + 1:])
        # the sheet graphs from the other buckets
        train_graphs = list()
        for bk in train_buckets:
            train_graphs.extend(bk.sheet_graphs)

        # identify the optimal weights based on examples from the training set.
        # for better accuracy, we perform the optimization several times
        # and take the weighted average of the results
        error_rates = list()
        weights = list()
        for _ in range(n_opt_runs):
            # create an optimization sample
            tables_metrics, part_metrics, tables_counts, graph_indices = \
                get_sample_for_optimization(train_graphs, metrics_calc_fun,
                                            max_n_frgm=max_n_frgm, n_frgm_per_table=n_frgm_per_table)

            tables_values, part_values = list(), list()
            for index in range(len(tables_metrics)):
                tables_values.append(list(tables_metrics[index].values()))
                part_values.append(list(part_metrics[index].values()))

            if tf_fit:
                t_values_scaled = tf_fit.transform(np.array(tables_values))
                p_values_scaled = tf_fit.transform(np.array(part_values))
            else:
                t_values_scaled = tables_values
                p_values_scaled = part_values

            optw = get_optimal_weights(t_values_scaled, p_values_scaled, weight_bounds,
                                       tbl_counts=None, method='slsqp_jac')
            weights.append(optw)
            error = calculate_error_rate(t_values_scaled, p_values_scaled, optw)
            error_rates.append(error)

        # calculate the weighted average of all optimization iterations
        # we use the error rate to weight the instances (i.e., the optimal weights)
        inverse_error_rates = [1 - er for er in error_rates]  # the smaller the error the better
        weights_t = np.transpose(weights)  # transpose to bring to the right shape for multiplication
        multi_ew = np.multiply(inverse_error_rates, weights_t)
        sum_ew = np.sum(multi_ew, axis=1)
        inverse_error_sum = sum(inverse_error_rates)
        avg_opt_weights = list(np.divide(sum_ew, inverse_error_sum))
        print("Optimal Weights: %s" % str(avg_opt_weights))
        print("Inverse Errors: %s" % str(inverse_error_rates))

        # the averaged weights for each metric
        metrics_weights = dict()
        # metrics_weights["count_empty_cols"] = avg_opt_weights[0]
        # metrics_weights["count_empty_rows"] = avg_opt_weights[1]
        metrics_weights["avg_adj_emt_width"] = avg_opt_weights[0]
        metrics_weights["avg_adj_emt_height"] = avg_opt_weights[1]
        metrics_weights["data_above_ratio"] = avg_opt_weights[2]
        metrics_weights["neg_head_align_ratio"] = avg_opt_weights[3]
        metrics_weights["neg_data_align_ratio"] = avg_opt_weights[4]
        metrics_weights["count_other_hvalid"] = avg_opt_weights[5]
        metrics_weights["all_in_one_column"] = avg_opt_weights[6]
        # metrics_weights["count_only_data"] = avg_opt_weights[7]
        # metrics_weights["count_only_head"] = avg_opt_weights[8]
        metrics_weights["overlap_ratio"] = avg_opt_weights[7]

        # identify the tables for each sheet graph in the test bucket
        for bk_gix in range(len(test_buckets[j].sheet_graphs)):
            bk_graph = test_buckets[j].sheet_graphs[bk_gix]

            # print("\nGraph={file='%s', sheet='%s'}" % (bk_graph.graph['file'], bk_graph.graph['sheet']))
            jaccards_per_table = dict()
            for _ in range(n_searches):
                # search for best fragmentation. search can be either exhaustive or genetic, depending on the graph size
                predicted, p_score, true_tables, t_score = \
                    reduced_search.identify_tables(bk_graph, metrics_weights, metrics_calc_fun,
                                                   genetic_algorithm='varOr', scaler=tf_fit, use_seed=True,
                                                   multi_process=False, verbose=False, show_stats=False)
                predicted_scores.append(p_score)
                target_scores.append(t_score)

                # evaluate the predictions, by finding per true table the best match from the predicted fragments
                max_jaccards = eval.get_max_jaccard_per_target_fragment(bk_graph, predicted)
                for tid in max_jaccards.keys():
                    if tid not in jaccards_per_table:
                        jaccards_per_table[tid] = list()
                    jaccards_per_table[tid].append(max_jaccards[tid])

            # calculate the average jaccard per table.
            avg_jaccard_per_table = dict()
            for tid in jaccards_per_table.keys():
                avg_jaccard_per_table[tid] = sum(jaccards_per_table[tid]) / n_searches
            scores_per_graph.append(avg_jaccard_per_table)
            # print("Avg Jaccard per Table = %s" % str(avg_jaccard_per_table))

    # display the evaluation results
    ordered_sheet_graphs = list()
    for bk in test_buckets:
        ordered_sheet_graphs.extend(bk.sheet_graphs)

    dict_scores_per_graph = dict()
    for ix in range(len(ordered_sheet_graphs)):
        dict_scores_per_graph[ix] = scores_per_graph[ix]
    print("\n\n\nOverall Results")
    eval.display_evaluation_results(dict_scores_per_graph, ordered_sheet_graphs, print_threshold=0.9)


if __name__ == '__main__':
    gold_standard = "./datasets/gold_standard.json"
    induced_error = "./datasets/induced_error.json"
    metrics_fun = Metrics.calculate_fragmentation_metrics_opt3
    bnds = ((0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000), (0, 1000))

    # robust = RobustScaler(quantile_range=(25, 75))
    # min_max = MinMaxScaler()
    main(gold_standard, induced_error, metrics_fun, bnds, n_buckets=10, seed=1, n_opt_runs=2,
         n_searches=2, max_n_frgm=200, n_frgm_per_table=10, tf=None)
