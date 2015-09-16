#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


# CAPS import
from fmri_subject_info import get_onsets

# system import
import os
import json


def spm_model_specification(behavioral_data, fmri_sessions, onset_name,
                            condition_name, duration_name, time_repetition,
                            realignment_parameters, delimiter, start,
                            concatenate_runs, high_pass_filter_cutoff,
                            output_directory):
    """ Specify the SPM model used in the GLM and estimate the design matrix.

    .. note::

        * `fmri_sessions` and `behavioral_data` must have the same number
          of elements.
        * `onsets` and `durations` values must have the same units as the
          TR used in the processings (ie. seconds).

    <unit>
        <input name="behavioral_data" type="List" content="File" desc="list of
            .csv session behavioral data." />
        <input name="fmri_sessions" type="List" content="File" desc="list of
            path to fMRI sessions." />
        <input name="onset_name" type="String" desc="the name of the column
            in the `behavioral_data` file containing the onsets."/>
        <input name="condition_name" type="String" desc="the name of the
            column in the `behavioral_data` file containing the conditions."/>
        <input name="duration_name" type="String" desc="the name of the column
            in the `behavioral_data` file containing the condition durations.
            "/>
        <input name="time_repetition" type="Float" desc="the repetition time
            in seconds (in seconds)."/>
        <input name="realignment_parameters" type="File" desc="path to the SPM
            realign output parameters."/>
        <input name="delimiter" type="String" desc="separator used to split
            the `behavioral_data` file."/>
        <input name="start" type="Int" desc="line from which we start reading
            the `behavioral_data` file."/>
        <input name="concatenate_runs" type="Bool" desc="concatenate all runs
            to look like a single session."/>
        <input name="high_pass_filter_cutoff" type="Float" desc="high-pass
            filter cutoff in secs."/>
        <input name="output_directory" type="Directory" desc="Where to store
            the output file"/>
        <output name="session_info" type="Any" desc="session info to leverage
            the first level design."/>
        <output name="model_specifications" type="File" desc="file containing
            all model specifications" />
    </unit>
    """
    # Local imports
    from nipype.interfaces.base import Bunch
    from nipype.algorithms.modelgen import SpecifySPMModel

    # Assert that we have one behavioral data per session
    if len(behavioral_data) != len(fmri_sessions):
        raise ValueError("One behavioral data per session is required, "
                         "got {0} behaviral data and {1} session.".format(
                             len(behavioral_data), len(fmri_sessions)))

    # Get each session acquisition conditions
    info = []
    for csvfile in behavioral_data:

        # Parse the behavioural file
        all_onsets = get_onsets(csvfile,
                                condition_name, onset_name, duration_name,
                                delimiter, start)

        # Create a nipype Bunch (dictionary-like) structure
        conditions = []
        onsets = []
        durations = []
        for condition_name, item in all_onsets.items():
            conditions.append(condition_name)
            onsets.append([float(x) for x in item["onsets"]])
            durations.append([float(x) for x in item["durations"]])
        info.append(
            Bunch(conditions=conditions, onsets=onsets, durations=durations))

    # Make a model specification compatible with spm designer
    spec_interface = SpecifySPMModel(
        concatenate_runs=concatenate_runs,
        input_units="secs",
        output_units="secs",
        time_repetition=time_repetition,
        high_pass_filter_cutoff=high_pass_filter_cutoff,
        functional_runs=fmri_sessions,
        subject_info=info,
        realignment_parameters=realignment_parameters)
    spec_interface.run()

    # The previous interface use numpy in double precision. In order to be
    # python-json compliant need to cast expicitely all float items
    def cast_to_float(obj):
        """ Recursive method that cast numpy.double items.

        Parameters
        ----------
        obj: object
            a generic python object.

        Returns
        -------
        out: object
            the float-casted input object.
        """
        # Deal with dictionary
        if isinstance(obj, dict):
            out = {}
            for key, val in obj.items():
                out[key] = cast_to_float(val)

        # Deal with tuple and list
        elif isinstance(obj, (list, tuple)):
            out = []
            for val in obj:
                out.append(cast_to_float(val))
            if isinstance(obj, tuple):
                out = tuple(out)

        # Otherwise cast if it is a numpy.double
        else:
            out = obj
            if isinstance(obj, float):
                out = float(obj)

        return out

    session_info = cast_to_float(spec_interface.aggregate_outputs().get()[
        "session_info"])

    model_specifications = os.path.join(output_directory,
                                        "model_specifications.json")

    # save the design parameters
    with open(model_specifications, "w") as _file:
        json.dump(session_info, _file, indent=4)

    return session_info, model_specifications
