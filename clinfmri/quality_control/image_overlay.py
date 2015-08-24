#! /usr/bin/env python
##########################################################################
# Nsap - Neurospin - Berkeley - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os

# Nilearn import
from nilearn import datasets
from nilearn import plotting, image


def edges_overlay(input_file, template_file, prefix="e", output_directory=None):
    """ Plot image outline on top of another image (useful for checking
    registration)

    <process>
        <return name="edges_file" type="File" desc="A snap with two overlayed
            images."/>
        <input name="input_file" type="File" desc="An image to display."/>
        <input name="template_file" type="File" desc="The target image to
            extract the edges from."/>
        <input name="prefix" type="String" desc="the prefix of the result
            file."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder." optional="True"/>
    </process>
    """
    # Check the input images exist on the file system
    for in_file in [input_file, template_file]:
        if not os.path.isfile(in_file):
            raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the outdir is valid
    if output_directory is not None:
        if not os.path.isdir(output_directory):
            raise ValueError("'{0}' is not a valid directory.".format(
                output_directory))
    else:
        output_directory = os.path.dirname(input_file)

    # Create the plot
    display = plotting.plot_anat(input_file, title="Overlay")
    display.add_edges(template_file)
    edges_file = os.path.join(
        output_directory,
        prefix + os.path.basename(input_file).split(".")[0] + ".pdf")
    display.savefig(edges_file)
    display.close()

    return edges_file
