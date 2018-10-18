import networkx as nx
import numpy as np
import random as rand
import matplotlib.pyplot as plt
import multiprocessing
import array
from collections import OrderedDict
from deap import base, creator, tools, algorithms

from model import *
from optimz import Metrics
import heuristics_based

# define types used for genetic algorithm
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", array.array, typecode="B", fitness=creator.FitnessMin)


def evaluate(individual, basic_fragments, fragments_indices, fragment_connections, **kwargs):
    selected_connections = list()
    for j in range(len(individual)):
        if individual[j]:
            selected_connections.append(fragment_connections[j])

    # selected_connections = fragment_connections[individual]
    partitions = get_partitions(basic_fragments, fragments_indices, selected_connections)

    metrics_function = kwargs.get('metrics_function')
    nodes_attrb = kwargs.get('node_features')
    column_widths = kwargs.get('column_widths')
    row_heights = kwargs.get('row_heights')
    weights = kwargs.get('weights')
    scaler = kwargs.get('scaler')

    measurements = metrics_function(nodes_attrb, partitions, column_widths, row_heights)
    fitness = Metrics.calculate_quality(measurements, weights, scaler)

    return fitness,


def get_partitions(basic_fragments, fragments_indices, selected_connections):
    g = nx.Graph()
    g.add_nodes_from(fragments_indices)
    g.add_edges_from(selected_connections)

    partitions = list()
    for comp in nx.connected_components(g):
        part = list()
        for ix in comp:
            part.extend(basic_fragments[ix])
        partitions.append(set(part))
    return partitions


def _spit_horizontally(sheet_graph):
    copy_graph = sheet_graph.copy()
    to_remove = list()
    for edge in sheet_graph.edges_iter(data=True):
        direction = edge[2]['direction']
        distance = edge[2]['distance']
        if (direction == 'LEFT' or direction == 'RIGHT') and distance > 0:
            to_remove.append((edge[0], edge[1]))
    copy_graph.remove_edges_from(to_remove)
    return copy_graph


def _spit_vertically(sheet_graph):
    copy_graph = sheet_graph.copy()
    to_remove = list()
    for edge in sheet_graph.edges_iter(data=True):
        direction = edge[2]['direction']
        distance = edge[2]['distance']
        if (direction == 'TOP' or direction == 'BOTTOM') and distance > 0:
            to_remove.append((edge[0], edge[1]))
    copy_graph.remove_edges_from(to_remove)
    return copy_graph


def spit_by_empty_rows(nodes_attrib, horizontal_fragments):
    fragments_by_empty_rows = list()
    for fragment in horizontal_fragments:
        union_rows = set()
        frag_node_attrib = dict()
        for node in fragment:
            union_rows.update(nodes_attrib[node]['rows'])
            frag_node_attrib[node] = nodes_attrib[node]['coordinates']

        min_row = min(union_rows)
        horizontal_cuts = set(r - 1 for r in (union_rows - {min_row})) - union_rows
        horizontal_cuts.add(min_row)

        ordered_frag_nodes = [item[0] for item in sorted(frag_node_attrib.items(),
                                                         key=lambda kv: (-kv[1]['ymin'], -kv[1]['xmin']))]
        continue_at = 0
        for hcut in sorted(horizontal_cuts, reverse=True):
            part = set()
            for ix in range(continue_at, len(ordered_frag_nodes)):
                node = ordered_frag_nodes[ix]
                if hcut <= frag_node_attrib[node]['ymin']:
                    part.add(node)
                else:
                    continue_at = ix
                    break
            fragments_by_empty_rows.append(part)
    return fragments_by_empty_rows


def get_basic_fragments_with_heuristics(sheet_graph):
    predicted, undecided = heuristics_based.identify(sheet_graph, reject_single_col=True, handle_overlaps=True, dist_th=2)
    heuristic_parts = predicted + undecided

    basic_fragments = list()
    hp_index = 0
    parent_indices = list()
    for hp in heuristic_parts:
        hsplit_graph = _spit_horizontally(sheet_graph.subgraph(hp))
        vsplit_graph = _spit_vertically(hsplit_graph)
        vsplit_components = list(nx.connected_components(vsplit_graph.to_undirected()))
        basic_fragments.extend(vsplit_components)
        parent_indices.extend([hp_index] * len(vsplit_components))
        hp_index += 1

    return basic_fragments, parent_indices


def get_basic_fragments(sheet_graph, nodes_attrib, vertical_split=False):
    undirected_copy = sheet_graph.to_undirected()
    hsplit_graph = _spit_horizontally(undirected_copy)

    # distinguish header nodes from the rest
    headers, all_nodes = dict(), dict()
    for node, attrib in nodes_attrib.items():
        label = attrib['label']
        if label == 'Header':
            # get the region coordinates
            headers[node] = attrib['coordinates']
            all_nodes[node] = attrib['coordinates']
        else:
            all_nodes[node] = attrib['coordinates']

    # sort header nodes in descending order of their maximum row, followed by ascending order of minimum row
    ordered_headers = OrderedDict(sorted(headers.items(), key=lambda item: (-item[1]['ymax'], item[1]['ymin'])))
    # sort graph nodes in descending order of their minimum row
    ordered_all = OrderedDict(sorted(all_nodes.items(), key=lambda item: (-item[1]['ymin'])))

    # identify tables arranged vertically in each connected component of the graph
    basic_fragments = list()
    head_comp_index = 0
    parent_indices = list()
    for comp in nx.connected_components(hsplit_graph):
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
                graph_of_candidates = hsplit_graph.subgraph(candidate_nodes)
                for node_set in nx.connected_components(graph_of_candidates):
                    if hnode in node_set:
                        connected_nodes.update(node_set)
                        break

                connected_headers = (set(comp_headers.keys()) & connected_nodes)
                if len(connected_nodes) > len(connected_headers):
                    if frontier:
                        if vertical_split:
                            subgraph_frontier = hsplit_graph.subgraph(frontier)
                            vsplit_graph = _spit_vertically(subgraph_frontier)
                            vsplit_components = list(nx.connected_components(vsplit_graph))
                            basic_fragments.extend(vsplit_components)
                            parent_indices.extend([head_comp_index]*len(vsplit_components))
                        else:
                            basic_fragments.append(connected_nodes)
                            parent_indices.append(head_comp_index)
                        head_comp_index += 1
                    frontier = set(connected_nodes)

                elif len(connected_nodes) == len(connected_headers):
                    frontier_columns = set()
                    for fnd in frontier:
                        frontier_columns.update(list(range(ordered_all[fnd]['xmin'], ordered_all[fnd]['xmax'])))

                    header_columns = set()
                    for hnd in connected_headers:
                        header_columns.update(list(range(ordered_all[hnd]['xmin'], ordered_all[hnd]['xmax'])))

                    if len(header_columns & frontier_columns) > 0:
                        frontier.update(connected_headers)
                    else:
                        if vertical_split:
                            subgraph_frontier = hsplit_graph.subgraph(frontier)
                            vsplit_graph = _spit_vertically(subgraph_frontier)
                            vsplit_components = list(nx.connected_components(vsplit_graph))
                            basic_fragments.extend(vsplit_components)
                            parent_indices.extend([head_comp_index] * len(vsplit_components))
                        else:
                            basic_fragments.append(connected_nodes)
                            parent_indices.append(head_comp_index)
                        head_comp_index += 1
                        frontier = set(connected_nodes)

                for cn in connected_nodes:
                    del comp_all[cn]

        if frontier:
            if vertical_split:
                subgraph_frontier = hsplit_graph.subgraph(frontier)
                vsplit_graph = _spit_vertically(subgraph_frontier)
                vsplit_components = list(nx.connected_components(vsplit_graph))
                basic_fragments.extend(vsplit_components)
                parent_indices.extend([head_comp_index] * len(vsplit_components))
            else:
                basic_fragments.append(frontier)
                parent_indices.append(head_comp_index)
            head_comp_index += 1

        # in case there are still nodes
        if len(comp_all) > 0:
            remaining_nodes = set(comp_all.keys())
            if vertical_split:
                subgraph_remaining = hsplit_graph.subgraph(remaining_nodes)
                vsplit_graph = _spit_vertically(subgraph_remaining)
                vsplit_components = list(nx.connected_components(vsplit_graph))
                basic_fragments.extend(vsplit_components)
                parent_indices.extend([head_comp_index] * len(vsplit_components))
            else:
                basic_fragments.append(remaining_nodes)
                parent_indices.append(head_comp_index)

    return basic_fragments, parent_indices


def get_fragment_connections(complete_graph, fragments):
    connections = list()
    for a in range(len(fragments)-1):
        frag_a = fragments[a]
        frag_edges = complete_graph.edges(frag_a)
        outside_nodes = set()
        for edge in frag_edges:
            if edge[0] not in frag_a:
                outside_nodes.add(edge[0])
            elif edge[1] not in frag_a:
                outside_nodes.add(edge[1])

        if outside_nodes:
            for b in range(a+1, len(fragments)):
                frag_b = fragments[b]
                intersection = outside_nodes & frag_b
                if intersection:
                    outside_nodes.difference_update(intersection)
                    connections.append((a, b))
                if not outside_nodes:
                    break
    return connections


def get_seed_individual(basic_connections, parent_indices):
    parents_basic_frags = dict()
    # basic fragment index (denoted as bix)
    # parent fragment index (denoted as pix)
    for bix in range(len(parent_indices)):
        pix = parent_indices[bix]
        if pix not in parents_basic_frags.keys():
            parents_basic_frags[pix] = list()
        parents_basic_frags[pix].append(bix)

    parent_connections = set()
    for pix in parents_basic_frags.keys():
        basic_indices = parents_basic_frags[pix]
        for m in range(len(basic_indices)-1):
            for n in range(m+1, len(basic_indices)):
                parent_connections.add((basic_indices[m], basic_indices[n]))

    seed_individ = list()
    for bc in basic_connections:
        if bc in parent_connections:
            seed_individ.append(1)
        else:
            seed_individ.append(0)
    return seed_individ


def get_population_from_seed(toolbox, seed_individ, pop_size, indpb=0.2):
    pop = list()
    for _ in range(pop_size):
        ind = toolbox.individual()
        for j in range(len(seed_individ)):
            if rand.random() < indpb:
                ind[j] = not seed_individ[j]
            else:
                ind[j] = seed_individ[j]
        pop.append(ind)
    return pop


def search_genetic(basic_fragments, basic_connections, gen_alg='varOr', seed_individ=None,
                   multi_process=False, show_stats=False, verbose=False, **kwargs):
    frags_indices = list(range(len(basic_fragments)))

    # register initialization methods
    toolbox = base.Toolbox()
    toolbox.register("attribute", rand.randint, 0, 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attribute, n=len(basic_connections))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # register genetic operators
    toolbox.register("mate", tools.cxUniform, indpb=0.5)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.1)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # register evaluation function, to measure the fitness
    toolbox.register("evaluate", evaluate, basic_fragments=basic_fragments, fragments_indices=frags_indices,
                     fragment_connections=basic_connections, **kwargs)
    pool = None
    if multi_process:
        pool = multiprocessing.Pool()
        toolbox.register("map", pool.map)
        if verbose:
            print("\nUsing a multiprocessing pool")

    pop_size = min(500, 100 + len(basic_connections) * 2)
    if verbose:
        print("Population Size = %d" % pop_size)

    if seed_individ:
        rand_pop_size = pop_size - int(pop_size/2)
        pop = toolbox.population(n=rand_pop_size)
        seed_pop = get_population_from_seed(toolbox, seed_individ, (pop_size - rand_pop_size - 1), indpb=0.1)
        pop.extend(seed_pop)

        the_seed = toolbox.individual()
        for vix in range(len(seed_individ)):
            the_seed[vix] = seed_individ[vix]
        pop.append(the_seed)
        if verbose:
            seed_fitness = evaluate(the_seed, basic_fragments, frags_indices, basic_connections, **kwargs)
            print("The seed fitness is: %s" % seed_fitness[0])
    else:
        pop = toolbox.population(n=pop_size)

    # create  hall of fame list
    hof = tools.HallOfFame(1)

    if show_stats:
        # specify which statistics to collect
        fit_stats = tools.Statistics(lambda ind: ind.fitness.values)
        fit_stats.register("mean", np.mean, axis=0)
        fit_stats.register("std", np.std, axis=0)
        fit_stats.register("min", np.min, axis=0)
        fit_stats.register("max", np.max, axis=0)
    else:
        fit_stats = None

    # specify the parameter values
    cxpb, mutpb, ngen = 0.5, 0.5, 100
    if gen_alg == 'varOr':
        # execute genetic varOr search
        mu_es, lampda_es = pop_size, pop_size
        if verbose:
            print("execute varOr search")
            print("VarOr: cxpb=%s, mutpb=%s, ngen=%s, mu_es=%s, lampda_es=%s" % (cxpb, mutpb, ngen, mu_es, lampda_es))
        result, log = algorithms.eaMuPlusLambda(pop, toolbox, mu=mu_es, lambda_=lampda_es, cxpb=cxpb, mutpb=mutpb,
                                                ngen=ngen, halloffame=hof, verbose=False, stats=fit_stats)
        best_ind = hof[0]
    elif gen_alg == 'varAnd':
        if verbose:
            print("execute varAnd search")
            print("VarAnd: cxpb=%s, mutpb=%s, ngen=%s" % (cxpb, mutpb, ngen))
        # execute simple genetic varAnd search
        result, log = algorithms.eaSimple(pop, toolbox, cxpb=cxpb, mutpb=mutpb,
                                          ngen=ngen, halloffame=hof, verbose=False, stats=fit_stats)
        best_ind = hof[0]
    else:
        raise ValueError('Value specified for genetic algorithm is not recognized! gen_alg=%s' % gen_alg)

    if multi_process and pool:
        pool.close()
        pool.join()
        if verbose:
            print("\nClosed the multiprocessing pool")

    if show_stats:
        plt.figure(figsize=(11, 4))
        plots = plt.plot(log.select('min'), 'c-', log.select('std'), 'y-', log.select('mean'), 'r-')
        plt.legend(plots, ('Minimum fitness', 'Std fitness', 'Mean fitness'), frameon=True)
        plt.ylabel('Fitness')
        plt.xlabel('Iterations')
        plt.show()

    selected_connections = list()
    for j in range(len(best_ind)):
        if best_ind[j]:
            selected_connections.append(basic_connections[j])
    best_ind_parts = get_partitions(basic_fragments, frags_indices, np.array(selected_connections))
    return best_ind_parts, best_ind.fitness.values[0]


def search_exhaustive(basic_fragments, fragment_connections, **kwargs):

    fragments_indices = list(range(len(basic_fragments)))

    # exhaustively generate all possible individuals
    pop = list()
    format_spec = "{0:0%db}" % len(fragment_connections)
    for num in range(0, 2 ** len(fragment_connections)):
        pop.append([int(digit) > 0.5 for digit in format_spec.format(num)])

    # identify the best individual
    best_fitness = evaluate(pop[0], basic_fragments, fragments_indices, fragment_connections, **kwargs)
    best_ind = pop[0]
    for ind in pop[1:]:
        fitness = evaluate(ind, basic_fragments, fragments_indices, fragment_connections, **kwargs)
        if fitness < best_fitness:
            best_fitness = fitness
            best_ind = ind

    selected_connections = list()
    for j in range(len(best_ind)):
        if best_ind[j]:
            selected_connections.append(fragment_connections[j])
    return get_partitions(basic_fragments, fragments_indices, selected_connections), best_fitness[0]


def identify_tables(sheet_graph, weights, metrics_function, scaler=None, genetic_algorithm='varOr', use_seed=False,
                    multi_process=False, verbose=False, show_stats=False):
    # collect node (region) coordinates and labels
    nodes_features = dict()
    for node, features in sheet_graph.nodes_iter(data=True):
        selected_features = dict()
        coordinates = Utils.get_node_region_coordinates(features)
        # coordinates follow the Cartesian system, where the origin (0,0)
        # is placed at the top-left corner of the worksheet.
        # the y-axis (row-wise) is inverted (i.e., values increase as we move downwards)
        selected_features['coordinates'] = coordinates
        # besides the coordinates, we maintain the column and row numbers (for the node regions).
        # the latter follows the R0C0 reference system, which starts enumerating rows and columns from 0.
        # note, the R1C1 reference system (standard Excel/VBA), starts from 1.
        selected_features['columns'] = list(range(coordinates['xmin'], coordinates['xmax']))
        selected_features['rows'] = list(range(coordinates['ymin'], coordinates['ymax']))
        selected_features['area'] = features['area']
        selected_features['label'] = features['label']
        nodes_features[node] = selected_features

    # the basic fragment are a result of horizontal split, followed by header split, and finally vertical split
    # basic_fragments, parent_indices = get_basic_fragments_with_heuristics(sheet_graph)
    basic_fragments, parent_indices = get_basic_fragments(sheet_graph, nodes_features, vertical_split=True)
    fragment_connections = get_fragment_connections(sheet_graph, basic_fragments)

    # determine the fitness of the target fragmentation
    target_fragments = list()
    for table_id, node_areas in Utils.group_node_area_tuples_by_tableId(sheet_graph, one_tbl_rule=True).items():
        node_keys = set(list(node_areas.keys()))
        if table_id == -1:  # these fragments do not form a table, they are independent (i.e., outside the tables)
            # treat this independent regions individually
            for node in node_keys:
                target_fragments.append({node})
        else:
            target_fragments.append(node_keys)

    target_measurements = metrics_function(nodes_features, target_fragments,
                                           sheet_graph.graph['column widths'], sheet_graph.graph['row heights'])
    target_fitness = Metrics.calculate_quality(target_measurements, weights, scaler=scaler)

    # prepare arguments needed for evaluation
    eval_args = dict()
    eval_args['metrics_function'] = metrics_function
    eval_args['node_features'] = nodes_features
    eval_args['column_widths'] = sheet_graph.graph['column widths']
    eval_args['row_heights'] = sheet_graph.graph['row heights']
    eval_args['weights'] = weights
    eval_args['scaler'] = scaler

    # search for tables using one of the available methods
    if len(fragment_connections) == 0:
        if verbose:
            print("\nWarning: graph has only one basic fragment!!!")
        return target_fragments, target_fitness, target_fragments, target_fitness
    elif len(fragment_connections) <= 7:
        if verbose:
            print("\nPerforming exhaustive search, fragment_connections = %d" % len(fragment_connections))
        fittest_part, fitness = search_exhaustive(basic_fragments, fragment_connections, **eval_args)
        return fittest_part, fitness, target_fragments, target_fitness

    if verbose:
        print("\nPerforming genetic search, fragment_connections = %d" % len(fragment_connections))
    if use_seed:
        seed_individ = get_seed_individual(fragment_connections, parent_indices)
        if verbose:
            print("Search using initial seed!")
    else:
        if verbose:
            print("Search without seed!")
        seed_individ = None

    fittest_part, fitness = search_genetic(basic_fragments, fragment_connections, gen_alg=genetic_algorithm,
                                           seed_individ=seed_individ, multi_process=multi_process, verbose=verbose,
                                           show_stats=show_stats, **eval_args)

    return fittest_part, fitness, target_fragments, target_fitness
