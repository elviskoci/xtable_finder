import numpy as np
from scipy.optimize import minimize

from optimz import Generator
from model import Utils


def quality_ratio(w, tc, pc, tn):
    multis_wtc = np.multiply(w, tc)
    sums_wtc = np.sum(multis_wtc, axis=1) + 1
    multis_wpc = np.multiply(w, pc)
    sums_wpc = np.sum(multis_wpc, axis=1) + 1

    divisions = np.divide(sums_wtc, sums_wpc)
    pow_div = np.power(divisions, 2)
    multi_div = np.multiply(tn, pow_div)
    return np.sum(multi_div)


def quality_ratio_constrained(w, tc, pc, tn):
    multis_wtc = np.multiply(w, tc)
    sums_wtc = np.sum(multis_wtc, axis=1) + 1
    multis_wpc = np.multiply(w, pc)
    sums_wpc = np.sum(multis_wpc, axis=1) + 1

    divisions = np.divide(sums_wtc, sums_wpc)
    multi_div = np.multiply(tn, divisions)
    return np.sum(multi_div)


def quality_ratio_constrained_derivative(w, tc, pc, tn):
    # using quotient
    # dfdx * g(x) - dgdx * f(x) / pow(g(x),2)
    # f(x) is the numerator and g(x) is the denominator
    multis_wpc = np.multiply(w, pc)
    gx = np.sum(multis_wpc, axis=1) + 1
    gxpow2 = np.power(gx, 2)

    multis_wtc = np.multiply(w, tc)
    sums_wtc = np.sum(multis_wtc, axis=1) + 1
    fx = np.multiply(sums_wtc, tn)

    derivatives = list()
    for i in range(0, len(w)):
        dfdxi = np.multiply(tc[:, i], tn)
        dgdxi = pc[:, i]
        fx_dgdxi = np.multiply(dgdxi, fx)
        gx_dfdxi = np.multiply(dfdxi, gx)

        numerator = np.add(gx_dfdxi, np.multiply(fx_dgdxi, -1))
        derivatives.append(np.sum(np.divide(numerator, gxpow2)))

    return np.array(derivatives)


def get_sample_for_optimization(graphs, metrics_function, max_n_frgm=10, n_frgm_per_table=-1, verbose=False):
    annotated_tables, other_fragmentations = list(), list()
    table_counts, count_single, count_multi = list(), list(), list()
    graph_indices = list()
    for gix in range(len(graphs)):
        sheet_graph = graphs[gix]

        # collect and prepare the attributes for each node in this graph
        nodes_attributes = dict()
        for node, data in sheet_graph.nodes_iter(data=True):
            selected_data = dict()
            coordinates = Utils.get_node_region_coordinates(data)
            selected_data['area'] = data['area']
            selected_data['coordinates'] = coordinates
            selected_data['columns'] = list(range(coordinates['xmin'], coordinates['xmax']))
            selected_data['rows'] = list(range(coordinates['ymin'], coordinates['ymax']))
            selected_data['label'] = data['label']
            nodes_attributes[node] = selected_data

        # identify the tables in this graph
        tables_node_areas = Utils.group_node_area_tuples_by_tableId(sheet_graph, one_tbl_rule=True)
        # the ids of the tables in the graph
        table_ids = list(tables_node_areas.keys())
        n_tables = len(table_ids) if -1 not in table_ids else len(table_ids) - 1

        # the fragments that represent true tables (a.k.a, the target solution)
        tbl_fragments = list()
        for tbl_id in table_ids:
            table_nodes = set(tables_node_areas[tbl_id].keys())
            if tbl_id == -1:  # a reserved value for independent regions (i.e., located outside the annotated tables)
                # treat this independent regions individually (separately)
                for node_region in table_nodes:
                    tbl_fragments.append({node_region})
            else:
                tbl_fragments.append(table_nodes)

        tmx = metrics_function(nodes_attributes, tbl_fragments, sheet_graph.graph['column widths'],
                               sheet_graph.graph['row heights'])

        # determine the number of random subgraphs to generate
        if n_frgm_per_table > 0:
            n_frgm = len(table_ids) * n_frgm_per_table
            if n_frgm > max_n_frgm:
                n_frgm = max_n_frgm
        else:
            n_frgm = max_n_frgm

        # generate random subgraphs.
        fragmentations = Generator.get_random_fragmentations(sheet_graph, n_frgm, as_subgraphs=False)
        for frag in fragmentations:
            fmx = metrics_function(nodes_attributes, frag, sheet_graph.graph['column widths'],
                                   sheet_graph.graph['row heights'])
            annotated_tables.append(tmx)
            other_fragmentations.append(fmx)
            table_counts.append(n_tables)
            graph_indices.append(gix)

        if n_tables > 1:
            count_multi.append(len(fragmentations))
        else:
            count_single.append(len(fragmentations))

    if verbose:
        print("\nSingle contributed %d" % sum(count_single))
        print("Count single %d" % len(count_single))
        print("Multi contributed %d" % sum(count_multi))
        print("Count multi %d" % len(count_multi))
    return annotated_tables, other_fragmentations, table_counts, graph_indices


def get_optimal_weights(tbl_values, alt_values, bnds, tbl_counts=None, scaler=None, method="slsqp_jac"):

    if scaler:
        t_coefs = scaler.transform(np.array(tbl_values))
        p_coefs = scaler.transform(np.array(alt_values))
    else:
        t_coefs = np.array(tbl_values)
        p_coefs = np.array(alt_values)

    if tbl_counts:
        t_counts = np.array(tbl_counts)
    else:
        t_counts = np.ones(len(tbl_values))
    initial_weights = np.array([1.0] * len(bnds))

    if method == 'slsqp_jac':
        res_slsqp_jac = minimize(quality_ratio_constrained, initial_weights, args=(t_coefs, p_coefs, t_counts),
                                 method='SLSQP', jac=quality_ratio_constrained_derivative, bounds=bnds,
                                 options={'disp': False, 'maxiter': 10000})
        return res_slsqp_jac.x
    elif method == 'slsqp':
        res_slsqp = minimize(quality_ratio_constrained, initial_weights, args=(t_coefs, p_coefs, t_counts),
                             method='SLSQP', bounds=bnds, options={'disp': False, 'maxiter': 10000})
        return res_slsqp.x
    elif method == 'bfgs_jac':
        res_bfgs_jac = minimize(quality_ratio_constrained, initial_weights, args=(t_coefs, p_coefs, t_counts),
                                method='L-BFGS-B', bounds=bnds, jac=quality_ratio_constrained_derivative,
                                options={'disp': False, 'maxiter': 10000})
        return res_bfgs_jac.x
    elif method == 'bfgs':
        res_bfgs = minimize(quality_ratio_constrained, initial_weights, args=(t_coefs, p_coefs, t_counts),
                            method='L-BFGS-B', bounds=bnds, options={'disp': False, 'maxiter': 10000})
        return res_bfgs.x
