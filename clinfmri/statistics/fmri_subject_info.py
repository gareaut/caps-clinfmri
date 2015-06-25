#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import csv
from collections import defaultdict


def get_onsets(csvfile, condition_name, onset_name, duration_name,
               delimiter="\t", start=0):
    """ Load .csv behavioral data
    """

    # read the csv
    onsets = defaultdict(list)
    with open(csvfile, 'rb') as ocsv:
        # move the pointer to the header
        for i in range(start):
            next(ocsv)
        reader = csv.DictReader(ocsv, delimiter=delimiter, quotechar='|')
        for row in reader:
            for (k, v) in row.items():
                onsets[k].append(v)

    # python 2.7 only (not working on Centos 6.4)
    # conditions = { x : {}
    #              for x in set([x for x in onsets[condition_name] if x!=""]) }
    conditions = {}
    for x in onsets[condition_name]:
        if x != "":
            conditions[x] = {}

    for condition in conditions.keys():
        indices = [i for i, x in enumerate(onsets[condition_name])
                     if x == condition]
        conditions[condition]["onsets"] = [onsets[onset_name][i]
                                           for i in indices]
        conditions[condition]["durations"] = [onsets[duration_name][i]
                                              for i in indices]

    return conditions

if __name__ == "__main__":
    get_onsets("/volatile/nsap/1level/mid.csv",
               condition_name="Trial Category",
               onset_name="Trial Start Time (Onset)",
               duration_name="Target Phase Duration",
               start=1)
