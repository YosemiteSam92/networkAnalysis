#!/usr/bin/env python3
import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))
if __name__ == "__main__":

    from streaming.dump_reader import (
        SingleRecordReader,
        AggregateRecordReader,
        MultipleRecordReader,
    )

    data_path1 = "dumps/dump.1_atomId_type_scaledCoords_nptCuring_gelation"
    data_path2 = "dumps/dump.1_btype_batom1_batom2_nptCuring_gelation"

    data_file1 = open(data_path1, "r")
    data_file2 = open(data_path2, "r")

    data = None

    reader1 = SingleRecordReader(data)
    reader2 = SingleRecordReader(data)

    aggregate_reader = AggregateRecordReader(data, [reader1, reader2])

    multi_reader = MultipleRecordReader(data, aggregate_reader)

    aggr = multi_reader.readRecord([data_file1, data_file2])

    aggr.print_info()