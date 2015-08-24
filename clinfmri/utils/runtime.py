#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import sys


def progress_bar(ratio, title="", bar_length=40):
    """ Display a custom progress bar.

    Parameters
    ----------
    ratio: float (mandatory)
        float describing the current processing status: 0 < ratio < 1.
    title: str (optional, default '')
        a title to identify the progress bar.
    bar_length: int (optional, default 40)
        the length of the bar that will be displayed.
    """
    progress = int(ratio * 100.)
    block = int(round(bar_length * ratio))
    text = "\r{2} in Progress: [{0}] {1}%".format(
        "=" * block + " " * (bar_length - block), progress, title)
    sys.stdout.write(text)
    sys.stdout.flush()
