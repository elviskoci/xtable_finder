import json


class Worksheet:
    def __init__(self, corpus_name, file_name, sheet_name):
        self._corpus_name = corpus_name
        self._file_name = file_name
        self._sheet_name = sheet_name
        self._regions = list()
        self._relations = list()
        self._column_widths = dict()
        self._row_heights = dict()

    def getStrRepr(self):
        return "Corpus : \"%s\", File: \"%s\", Sheet: \"%s\"" % (self._corpus_name, self._file_name, self._sheet_name)

    def __str__(self):
        return self.getStrRepr()

    @property
    def corpus_name(self):
        return self._corpus_name

    @property
    def file_name(self):
        return self._file_name

    @property
    def sheet_name(self):
        return self._sheet_name

    @property
    def regions(self):
        return self._regions

    @regions.setter
    def regions(self, value):
        self._regions = value

    @property
    def relations(self):
        return self._relations

    @relations.setter
    def relations(self, value):
        self._relations = value

    @property
    def column_widths(self):
        return self._column_widths

    @column_widths.setter
    def column_widths(self, columns_to_widths):
        self._column_widths = columns_to_widths

    @property
    def row_heights(self):
        return self._row_heights

    @row_heights.setter
    def row_heights(self, rows_to_heights):
        self._row_heights = rows_to_heights

    def getId(self):
        return list([self._corpus_name, self._file_name, self._sheet_name])


class Region:
    def __init__(self, label, address, xcenter, ycenter, area, width, height, tables):
        self._label = label
        self._address = address
        self._xcenter = float(xcenter)
        self._ycenter = float(ycenter)
        self._area = float(area)
        self._width = float(width)
        self._height = float(height)
        self._tables = list(tables)

    def getStrRepr(self):
        return "Label : %s, Address: %s " % (self.label, self.address)

    def __str__(self):
        return self.getStrRepr()

    @property
    def label(self):
        return self._label

    @property
    def address(self):
        return self._address

    @property
    def xcenter(self):
        return self._xcenter

    @property
    def ycenter(self):
        return self._ycenter

    @property
    def area(self):
        return self._area

    @property
    def height(self):
        return self._height

    @property
    def width(self):
        return self._width

    @property
    def tables(self):
        return self._tables


class Relation:
    def __init__(self, source, target, relevance, attraction, distance, overlap_ratio, direction, target_type):
        self._source = source
        self._target = target
        self._relevance = float(relevance)
        self._attraction = float(attraction)
        self._distance = distance
        self._overlap_ratio = overlap_ratio
        self._direction = direction
        self._target_type = target_type

    def getStrRepr(self):
        return "Source: %s, Target: %s, Direction: %s" % (self._source, self._target, self._direction)

    def __str__(self):
        return self.getStrRepr()

    @property
    def source(self):
        return self._source

    @property
    def target(self):
        return self._target

    @property
    def relevance(self):
        return self._relevance

    @property
    def attraction(self):
        return self._attraction

    @property
    def distance(self):
        return self._distance

    @property
    def overlap_ratio(self):
        return self._overlap_ratio

    @property
    def direction(self):
        return self._direction

    @property
    def target_type(self):
        return self._target_type


class Reader:
    @staticmethod
    def read_json_dataset(file_path):
        file = open(file_path, 'r')
        contents = file.read()
        file.close()

        obj = json.loads(contents)

        worksheets = list()
        for sheet_obj in obj['dataset']:

            corpus_name = sheet_obj['corpus']
            file_name = sheet_obj['file']
            sheet_name = sheet_obj['sheet']

            ws = Worksheet(corpus_name, file_name, sheet_name)

            regions = list()
            relations = list()
            for region_obj in sheet_obj['regions']:

                label = region_obj['label']
                address = region_obj['address']
                xcenter = region_obj['xcenter']
                ycenter = region_obj['ycenter']
                area = region_obj['area']
                width = region_obj['width']
                height = region_obj['height']

                tables_obj = region_obj['tables']
                tables = list()
                for entry in tables_obj:
                    table_id = entry['table']
                    shared_area = entry['area']
                    tables.append((table_id, shared_area))

                reg = Region(label, address, xcenter, ycenter, area, width, height, tables)
                regions.append(reg)

                relations_obj = region_obj['relations']
                for entry in relations_obj:
                    source_address = address
                    target_address = entry['address']
                    relevance = entry['relevance']
                    attraction = entry['attraction']
                    distance = entry['distance']
                    overlap_ratio = entry['overlap_ratio']
                    direction = entry['direction']
                    target_type = entry['type']

                    rel = Relation(source_address, target_address, relevance, attraction, distance, overlap_ratio,
                                   direction, target_type)
                    relations.append(rel)

            column_widths = dict()
            if 'column_widths' in sheet_obj:
                for colNum in sheet_obj['column_widths'].keys():
                    column_widths[int(colNum)] = sheet_obj['column_widths'][colNum]

            row_heights = dict()
            if 'row_heights' in sheet_obj:
                for rowNum in sheet_obj['row_heights'].keys():
                    row_heights[int(rowNum)] = sheet_obj['row_heights'][rowNum]

            ws.regions = regions
            ws.relations = relations
            ws.column_widths = column_widths
            ws.row_heights = row_heights
            worksheets.append(ws)

        return worksheets


if __name__ == '__main__':
    # sheets = Reader.read_json_dataset(
    #     "C:/Users/Elvis/Documents/PyNotebooks/cell_classification/3_labels/gold_standard_graphs.json")

    test = "C:/Users/Elvis/Documents/Workspace/XPostProcessor/test.json"
    sheets = Reader.read_json_dataset(test)
    # print(len(sheets))

    for sheet in sheets:
        print(sheet.getId())

