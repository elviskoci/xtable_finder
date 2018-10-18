
from collections import OrderedDict, namedtuple
import numpy as np

Coordinates = namedtuple('Coordinates', 'xmin ymin xmax ymax')


class Metrics:

    @staticmethod
    def _group_regions_by_exact_alignment(regions_coordinates, orientation='H'):
        groups = dict()
        for reg, crds in regions_coordinates.items():

            if orientation == 'H':
                start = crds['xmin']
                end = crds['xmax']
            elif orientation == 'V':
                start = crds['ymin']
                end = crds['ymax']
            else:
                raise Exception('Unrecognized orientation specification!')

            interval = (start, end)
            if interval not in groups.keys():
                groups[interval] = list()
            groups[interval].append(reg)

        return groups

    @staticmethod
    def _group_regions_by_overlaps(regions_coordinates, orientation='H', strict_overlap=False):
        region_groups = Metrics._group_regions_by_exact_alignment(regions_coordinates, orientation)
        if len(region_groups.keys()) > 1:
            ordered_groups = OrderedDict(
                sorted(region_groups.items(), key=lambda item: (item[0][0], -item[0][1])))
            intervals = list(ordered_groups.keys())
            tmp = list(intervals[1:])
            merged_groups = dict()
            i = 0
            while i < len(intervals):
                intva = intervals[i]
                start_a = intva[0]
                end_a = intva[1]
                interval_nodes = list()
                is_overlap = False
                for intvb in tmp:
                    start_b = intvb[0]
                    end_b = intvb[1]
                    if (not (end_a <= start_b) and strict_overlap) or (not (end_a < start_b) and not strict_overlap):
                        start_a = min(start_a, start_b)
                        end_a = max(end_a, end_b)
                        interval_nodes.extend(region_groups[intvb])
                        is_overlap = True
                        del region_groups[intvb]
                        i = i + 1
                    else:
                        break
                if is_overlap:
                    interval_nodes.extend(region_groups[intva])
                    del region_groups[intva]
                    merged_groups[(start_a, end_a)] = interval_nodes
                i = i + 1
                tmp = list(intervals[(i + 1):])
            region_groups.update(merged_groups)
        return OrderedDict(sorted(region_groups.items(), key=lambda item: (item[0][0], -item[0][1])))

    @staticmethod
    def _get_top_most_headers(headers_coordinates):
        header_groups = Metrics._group_regions_by_overlaps(headers_coordinates, orientation='V', strict_overlap=True)
        intervals = list(header_groups.keys())
        if len(intervals) > 1:
            top, next_ = intervals[0], intervals[1]
            frontier = header_groups[top]
            top_headers_list = list(header_groups[top])
            i = 2
            while (next_ is not None) and (next_[0] - top[1] == 0):
                are_aligned = False
                for ht in frontier:
                    coord_ht = headers_coordinates[ht]
                    start_a, end_a = coord_ht['xmin'], coord_ht['xmax']
                    for hn in header_groups[next_]:
                        coord_hn = headers_coordinates[hn]
                        start_b, end_b = coord_hn['xmin'], coord_hn['xmax']
                        if not(end_a < start_b) and not(end_b < start_a):
                            are_aligned = True
                            break
                    if are_aligned:
                        top = (top[0], next_[1])
                        top_headers_list.extend(header_groups[next_])
                        frontier = header_groups[next_]
                        break
                next_ = intervals[i] if i < len(intervals) else None
                i = i + 1
            return top, top_headers_list
        else:
            return intervals[0], header_groups[intervals[0]]

    @staticmethod
    def _get_grouped_headers(headers_coordinates):
        vheader_groups = Metrics._group_regions_by_overlaps(headers_coordinates, orientation='V', strict_overlap=True)
        intervals = list(vheader_groups.keys())
        if len(intervals) > 1:
            updated_hgroups = dict()
            top_intv = intervals[0]
            top_headers = list(vheader_groups[top_intv])
            for i in range(1, len(intervals)+1):
                if (i < len(intervals)) and (intervals[i][0] - top_intv[1]) == 0:
                    are_aligned = False
                    for ht in vheader_groups[intervals[i-1]]:
                        coord_ht = headers_coordinates[ht]
                        start_a, end_a = coord_ht['xmin'], coord_ht['xmax']
                        for hb in vheader_groups[intervals[i]]:
                            coord_hn = headers_coordinates[hb]
                            start_b, end_b = coord_hn['xmin'], coord_hn['xmax']
                            if not (end_a <= start_b) and not (end_b <= start_a):
                                are_aligned = True
                                break
                        if are_aligned:
                            break

                    if not are_aligned:
                        updated_hgroups[top_intv] = list(top_headers)
                        top_intv = tuple(intervals[i])
                        top_headers = list(vheader_groups[intervals[i]])
                    else:
                        top_intv = (top_intv[0], intervals[i][1])
                        top_headers.extend(vheader_groups[intervals[i]])
                else:
                    updated_hgroups[top_intv] = list(top_headers)
                    if i < len(intervals):
                        top_intv = tuple(intervals[i])
                        top_headers = list(vheader_groups[intervals[i]])

            return OrderedDict(sorted(updated_hgroups.items(), key=lambda item: (item[0][0], -item[0][1])))
        else:
            return vheader_groups

    @staticmethod
    def calculate_quality(metrics, weights=None, scaler=None):
        """
         :param metrics a dictionary associating metrics (their names) to calculated values
         :param weights a dictionary associating metrics (their names) to weights. if weight not specified for a metric,
         assume 1.0
         :param scaler one of the sklearn pre-processing transformer scalers
        """
        metrics_values, weights_values = list(), list()
        for key in metrics.keys():
            if weights and (key in weights.keys()):
                metrics_values.append(metrics[key])
                weights_values.append(weights[key])
            else:
                if weights and key not in weights.keys():
                    print("Warning: no weight specified for metric '%s'!" % key)
                metrics_values.append(metrics[key])
                weights_values.append(1.0)

        if scaler:
            scaled_values = scaler.transform(np.array(metrics_values).reshape((1, -1)))[0]
            quality = sum([scaled_values[i] * weights_values[i] for i in range(len(scaled_values))])
        else:
            quality = sum([metrics_values[i] * weights_values[i] for i in range(len(metrics_values))])

        return quality

    @staticmethod
    def calculate_fragmentation_metrics_opt2(regions, fragments, column_widths, row_heights):

        metrics_values = OrderedDict()
        metrics_values["count_empty_cols"] = 0.0
        metrics_values["count_empty_rows"] = 0.0
        metrics_values["data_above_ratio"] = 0.0
        metrics_values["neg_head_align_ratio"] = 0.0
        metrics_values["neg_data_align_ratio"] = 0.0
        metrics_values["count_other_hvalid"] = 0.0
        metrics_values["all_in_one_column"] = 0.0
        metrics_values["overlap_ratio"] = 0.0

        # collect the coordinates to determine if fragments overlap with each other
        fragment_coordinates = list()
        total_used_area = 0.0
        # identify unique empty columns and rows
        # count once empty rows/cols in common among the fragments
        empty_cols = set()
        empty_rows = set()

        used_cols = set()
        used_rows = set()

        # fragments consisting solely of Data or Header regions
        only_data_frag = 0.0
        only_head_frag = 0.0
        for frag in fragments:

            dcols, drows = set(), set()
            hcols, hrows = set(), set()
            total_dcells, total_hcells = 0, 0
            headers_coordinates, data_coordinates = dict(), dict()
            # distinguish among Data and Header, and gather info
            for reg in frag:
                if regions[reg]['label'] == 'Header':
                    headers_coordinates[reg] = regions[reg]['coordinates']
                    hcols.update(regions[reg]['columns'])
                    hrows.update(regions[reg]['rows'])
                    total_hcells += regions[reg]['area']
                else:
                    data_coordinates[reg] = regions[reg]['coordinates']
                    dcols.update(regions[reg]['columns'])
                    drows.update(regions[reg]['rows'])
                    total_dcells += regions[reg]['area']

            frag_cols = hcols | dcols
            frag_rows = hrows | drows
            # determine the fragment coordinates
            frag_coord = (min(frag_cols), min(frag_rows), max(frag_cols) + 1, max(frag_rows) + 1)
            fragment_coordinates.append(frag_coord)
            total_used_area += (frag_coord[3] - frag_coord[1]) * (frag_coord[2] - frag_coord[0])
            # gather empty rows and columns from the fragment
            empty_cols.update((set(range(frag_coord[0] + 1, frag_coord[2])) - frag_cols))
            empty_rows.update((set(range(frag_coord[1] + 1, frag_coord[3])) - frag_rows))

            used_cols.update(frag_cols)
            used_rows.update(frag_rows)

            dcells_above_ratio = 0.0
            alignment_jaccard = 0.0
            header_align_ratio = 1.0
            data_align_ratio = 0.0
            count_other_hvalid = 0.0
            if total_hcells > 0:
                if total_dcells > 0:
                    # create header groups based on their vertical alignment (i.e., row-wise alignment)
                    ordered_hgroups = Metrics._get_grouped_headers(headers_coordinates)
                    hg_intervals = list(ordered_hgroups.keys())
                    ordered_data = sorted(data_coordinates.items(), key=lambda kv: (kv[1]['ymin'], kv[1]['xmin']))

                    # identify data cells that are higher than the top header group
                    top_intv = hg_intervals[0]
                    dcells_above = 0
                    for j in range(len(ordered_data)):
                        dcoord = ordered_data[j][1]  # the coordinates of the data region
                        if dcoord['ymin'] < top_intv[0]:
                            # calculate the number of cells, considering that this data region
                            # might be only partially higher than the header group
                            dcells_above += (min(top_intv[0], dcoord['ymax']) - dcoord['ymin']) * \
                                            (dcoord['xmax'] - dcoord['xmin'])
                        else:
                            break
                    dcells_above_ratio = dcells_above / total_dcells

                    # measure the alignment of the top header-group with the data regions
                    top_hcols = set()
                    for h in ordered_hgroups[top_intv]:
                        top_hcols.update(regions[h]['columns'])
                    count_common = len(top_hcols & dcols)  # consider all data columns
                    alignment_jaccard = count_common / len(frag_cols)
                    data_align_ratio = count_common / len(dcols)
                    header_align_ratio = count_common / len(top_hcols)

                    # determine if there are other header groups,
                    # which could potentially play the role of a table header.
                    # in other words, here we check if the fragment carries more than one table.
                    if len(hg_intervals) > 1:
                        for j in range(1, len(hg_intervals)):
                            hg_cols = set()
                            for h in ordered_hgroups[hg_intervals[j]]:
                                hg_cols.update(regions[h]['columns'])
                            if len(hg_cols) >= 2:  # TODO: Change to 2
                                count_other_hvalid += 1
                else:
                    data_align_ratio = 1.0
                    header_align_ratio = 1.0
                    only_head_frag += 1
            else:
                only_data_frag += 1

            # aggregate metrics values
            metrics_values['data_above_ratio'] += dcells_above_ratio
            metrics_values['neg_head_align_ratio'] += 1.0 - header_align_ratio
            metrics_values['neg_data_align_ratio'] += 1.0 - data_align_ratio
            metrics_values['all_in_one_column'] += 1.0 if (len(hcols) == 1) and (len(dcols) == 1) and \
                                                          (alignment_jaccard == 1) else 0
            metrics_values['count_other_hvalid'] += count_other_hvalid

        # identify overlaps and measure their size
        total_used_area = len(used_cols) * len(used_rows)
        sum_overlap_ratio = 0.0
        fragment_coordinates.sort(key=lambda coord: (coord[1], coord[0]))
        for i in range(len(fragment_coordinates) - 1):
            f1_coord = fragment_coordinates[i]
            for j in range(i + 1, len(fragment_coordinates)):
                f2_coord = fragment_coordinates[j]
                if f2_coord[1] < f1_coord[3]:
                    if not (f1_coord[2] <= f2_coord[0] or f1_coord[0] >= f2_coord[2]):
                        # identify the coordinates of the overlap
                        min_ymax = min(f1_coord[3], f2_coord[3])
                        max_ymin = f2_coord[1]
                        min_xmax = min(f1_coord[2], f2_coord[2])
                        max_xmin = max(f1_coord[0], f2_coord[0])

                        overlap_size = (min_ymax - max_ymin) * (min_xmax - max_xmin)
                        overlap_ratio = overlap_size / total_used_area
                        sum_overlap_ratio += overlap_ratio
                else:
                    break

        metrics_values["overlap_ratio"] = sum_overlap_ratio
        metrics_values['count_empty_cols'] = len(empty_cols)
        metrics_values['count_empty_rows'] = len(empty_rows)
        return metrics_values

    @staticmethod
    def calculate_fragmentation_metrics_opt3(regions, fragments, column_widths, row_heights):

        metrics_values = OrderedDict()
        metrics_values["avg_adj_emt_width"] = 0.0
        metrics_values["avg_adj_emt_height"] = 0.0
        metrics_values["data_above_ratio"] = 0.0
        metrics_values["neg_head_align_ratio"] = 0.0
        metrics_values["neg_data_align_ratio"] = 0.0
        metrics_values["count_other_hvalid"] = 0.0
        metrics_values["all_in_one_column"] = 0.0
        # metrics_values["count_only_data"] = 0.0
        # metrics_values["count_only_head"] = 0.0
        metrics_values["overlap_ratio"] = 0.0
        # metrics_values["n_fragments"] = 0.0

        # collect the coordinates to determine if fragments overlap with each other
        fragment_coordinates = list()

        # identify unique empty columns and rows
        # count once empty rows/cols in common among the fragments
        empty_cols, empty_rows = set(), set()
        empty_adjacent_columns = set()
        empty_adjacent_rows = set()
        used_cols, used_rows = set(), set()

        # fragments consisting solely of Data or Header regions
        only_data_frag = 0.0
        only_head_frag = 0.0
        for frag in fragments:

            dcols, drows = set(), set()
            hcols, hrows = set(), set()
            total_dcells, total_hcells = 0, 0
            headers_coordinates, data_coordinates = dict(), dict()
            # distinguish among Data and Header, and gather info
            for reg in frag:
                if regions[reg]['label'] == 'Header':
                    headers_coordinates[reg] = regions[reg]['coordinates']
                    hcols.update(regions[reg]['columns'])
                    hrows.update(regions[reg]['rows'])
                    total_hcells += regions[reg]['area']
                else:
                    data_coordinates[reg] = regions[reg]['coordinates']
                    dcols.update(regions[reg]['columns'])
                    drows.update(regions[reg]['rows'])
                    total_dcells += regions[reg]['area']

            frag_cols = hcols | dcols
            frag_rows = hrows | drows
            # determine the fragment coordinates
            frag_coord = (min(frag_cols), min(frag_rows), max(frag_cols) + 1, max(frag_rows) + 1)
            fragment_coordinates.append(frag_coord)
            # gather empty rows and columns from the fragment
            frag_empty_cols = (set(range(frag_coord[0] + 1, frag_coord[2])) - frag_cols)
            frag_empty_rows = (set(range(frag_coord[1] + 1, frag_coord[3])) - frag_rows)
            empty_cols.update(frag_empty_cols)
            empty_rows.update(frag_empty_rows)
            used_cols.update(frag_cols)
            used_rows.update(frag_rows)

            # group adjacent empty rows for this fragment
            frag_adj_empty_rows = set()
            if frag_empty_rows:
                ordered_empty_rows = sorted(frag_empty_rows)
                prev_row = ordered_empty_rows[0]
                adj_rows = {prev_row}
                for row in ordered_empty_rows[1:]:
                    if row - prev_row == 1:
                        adj_rows.add(row)
                    else:
                        frag_adj_empty_rows.add(frozenset(adj_rows))
                        adj_rows = {row}
                    prev_row = row
                frag_adj_empty_rows.add(frozenset(adj_rows))
            empty_adjacent_rows.update(frag_adj_empty_rows)

            # group adjacent empty columns for this fragment
            frag_adj_empty_cols = set()
            if frag_empty_cols:
                ordered_empty_cols = sorted(frag_empty_cols)
                prev_col = ordered_empty_cols[0]
                adj_cols = {prev_col}
                for col in ordered_empty_cols[1:]:
                    if col - prev_col == 1:
                        adj_cols.add(col)
                    else:
                        frag_adj_empty_cols.add(frozenset(adj_cols))
                        adj_cols = {col}
                    prev_col = col
                frag_adj_empty_cols.add(frozenset(adj_cols))
            empty_adjacent_columns.update(frag_adj_empty_cols)

            # calculate alignment related metrics
            dcells_above_ratio = 0.0
            alignment_jaccard = 0.0
            header_align_ratio = 1.0
            data_align_ratio = 0.0
            count_other_hvalid = 0.0
            if total_hcells > 0:
                if total_dcells > 0:

                    # create header groups based on their vertical alignment (i.e., row-wise alignment)
                    ordered_hgroups = Metrics._get_grouped_headers(headers_coordinates)
                    hg_intervals = list(ordered_hgroups.keys())
                    ordered_data = sorted(data_coordinates.items(), key=lambda kv: (kv[1]['ymin'], kv[1]['xmin']))

                    # identify data cells that are higher than the top header group
                    top_intv = hg_intervals[0]
                    dcells_above = 0
                    for j in range(len(ordered_data)):
                        dcoord = ordered_data[j][1]  # the coordinates of the data region
                        if dcoord['ymin'] < top_intv[0]:
                            # calculate the number of cells, considering that this data region
                            # might be only partially higher than the header group
                            dcells_above += (min(top_intv[0], dcoord['ymax']) - dcoord['ymin']) * \
                                            (dcoord['xmax'] - dcoord['xmin'])
                        else:
                            break
                    dcells_above_ratio = dcells_above / total_dcells

                    # measure the alignment of the top header-group with the data regions
                    top_hcols = set()
                    for h in ordered_hgroups[top_intv]:
                        top_hcols.update(regions[h]['columns'])
                    count_common = len(top_hcols & dcols)  # consider all data columns
                    alignment_jaccard = count_common / len(frag_cols)
                    data_align_ratio = count_common / len(dcols)
                    header_align_ratio = count_common / len(top_hcols)

                    # determine if there are other header groups,
                    # which could potentially play the role of a table header.
                    # in other words, here we check if the fragment carries more than one table.
                    if len(hg_intervals) > 1:
                        for j in range(1, len(hg_intervals)):
                            hg_cols = set()
                            for h in ordered_hgroups[hg_intervals[j]]:
                                hg_cols.update(regions[h]['columns'])
                            if len(hg_cols) >= 2:  # TODO: try to set to 2
                                count_other_hvalid += 1
                else:
                    data_align_ratio = 1.0
                    header_align_ratio = 0.0
                    only_head_frag += 1
            else:
                # header_align_ratio = 1.0
                # data_align_ratio = 1.0
                only_data_frag += 1

            # aggregate metrics values
            metrics_values['data_above_ratio'] += dcells_above_ratio
            metrics_values['neg_head_align_ratio'] += 1.0 - header_align_ratio
            metrics_values['neg_data_align_ratio'] += 1.0 - data_align_ratio
            metrics_values['all_in_one_column'] += 1.0 if (len(hcols) == 1) and (len(dcols) == 1) and \
                                                          (alignment_jaccard == 1) else 0
            metrics_values['count_other_hvalid'] += count_other_hvalid

        # identify overlaps and measure their size
        total_used_area = len(used_cols) * len(used_rows)
        sum_ov_ratio = 0.0
        sum_overlaps = 0.0
        fragment_coordinates.sort(key=lambda coord: (coord[1], coord[0]))
        for i in range(len(fragment_coordinates) - 1):
            f1_coord = fragment_coordinates[i]
            for j in range(i + 1, len(fragment_coordinates)):
                f2_coord = fragment_coordinates[j]
                if f2_coord[1] < f1_coord[3]:
                    if not (f1_coord[2] <= f2_coord[0] or f1_coord[0] >= f2_coord[2]):
                        # identify the coordinates of the overlap
                        min_ymax = min(f1_coord[3], f2_coord[3])
                        max_ymin = f2_coord[1]
                        min_xmax = min(f1_coord[2], f2_coord[2])
                        max_xmin = max(f1_coord[0], f2_coord[0])

                        overlap_size = (min_ymax - max_ymin) * (min_xmax - max_xmin)
                        overlap_ratio = overlap_size / total_used_area
                        sum_overlaps += overlap_size
                        sum_ov_ratio += overlap_ratio
                else:
                    break

        # calculate the average width/height for consecutive empty columns/rows
        avg_adj_emt_width = 0.0
        # max_adj_emt_width = 0.0
        if empty_adjacent_columns:
            adj_empty_widths = list()
            for col_group in empty_adjacent_columns:
                col_group_sum = 0.0
                for c in col_group:
                    col_group_sum += column_widths[c]
                adj_empty_widths.append(col_group_sum)
            avg_adj_emt_width = sum(adj_empty_widths)/len(adj_empty_widths)
            # max_adj_emt_width = max(adj_empty_widths)

        avg_adj_emt_height = 0.0
        # max_adj_emt_height = 0.0
        if empty_adjacent_rows:
            adj_empty_heights = list()
            for row_group in empty_adjacent_rows:
                row_group_sum = 0.0
                for r in row_group:
                    row_group_sum += row_heights[r]
                adj_empty_heights.append(row_group_sum)
            avg_adj_emt_height = sum(adj_empty_heights)/len(adj_empty_heights)
            # max_adj_emt_height = max(adj_empty_heights)

        metrics_values["overlap_ratio"] = sum_ov_ratio
        metrics_values['avg_adj_emt_width'] = avg_adj_emt_width
        metrics_values['avg_adj_emt_height'] = avg_adj_emt_height
        # metrics_values["n_fragments"] = len(fragments)/len(regions.keys())
        # metrics_values["count_only_data"] = only_data_frag
        # metrics_values["count_only_head"] = only_head_frag
        return metrics_values
