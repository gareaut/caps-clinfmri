#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# CAPSUL import
try:
    from capsul.utils.pilot import pilotfunction
except:
    def pilotfunction(func):
        return func


@pilotfunction
def pilot_smoothing():
    """ 
    Smoothing
    =========
    """
    # Pilot imports
    import os
    from caps.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from capsul.process import get_process_instance

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/pclinfmri/spmsmoothing"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "MatlabConfig", "SPMConfig",
                 "NipypeConfig"],
        use_smart_caching=True,
        matlab_exec="/neurospin/local/bin/matlab",
        use_matlab=True,
        spm_directory="/i2bm/local/spm8",
        use_spm=True,
        output_directory=working_dir)

    """
    Load the toy dataset
    --------------------

    To do so, we use the get_sample_data function to download the toy
    dataset on the local file system (here localizer data):
    """
    toy_dataset = get_sample_data("localizer")

    """
    The toy_dataset is an Enum structure with some specific elements of
    interest:

        * **??**: ??.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <pclinfmri.preproc.pipeline.SliceTiming>` that
    define the different step of the processings:
    """
    pipeline = get_process_instance("pclinfmri.utils.spm_smoothing.xml")
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.image_file = toy_dataset.fmri

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=True, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


@pilotfunction
def pilot_bet():
    """ 
    BET
    ===
    """
    # Pilot imports
    import os
    from caps.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from pclinfmri.utils.pipeline import FslBet

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/pclinfmri/fslbet"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=True,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,        
        output_directory=working_dir)

    """
    Load the toy dataset
    --------------------

    To do so, we use the get_sample_data function to download the toy
    dataset on the local file system (here localizer data):
    """
    toy_dataset = get_sample_data("localizer")

    """
    The toy_dataset is an Enum structure with some specific elements of
    interest:

        * **??**: ??.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <pclinfmri.preproc.pipeline.SliceTiming>` that
    define the different step of the processings:
    """
    pipeline = FslBet()
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.input_image_file = toy_dataset.anat
    pipeline.generate_binary_mask = True
    pipeline.bet_threshold = 0.5

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=True, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


@pilotfunction
def pilot_newsegment():
    """ 
    New Segment
    ===========
    """
    # Pilot imports
    import os
    from caps.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from pclinfmri.utils.pipeline import SpmNewSegment

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/pclinfmri/spmnewsegment"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "MatlabConfig", "SPMConfig",
                 "NipypeConfig"],
        use_smart_caching=True,
        matlab_exec="/neurospin/local/bin/matlab",
        use_matlab=True,
        spm_directory="/i2bm/local/spm8",
        use_spm=True,
        output_directory=working_dir)

    """
    Load the toy dataset
    --------------------

    To do so, we use the get_sample_data function to download the toy
    dataset on the local file system (here localizer data):
    """
    toy_dataset = get_sample_data("localizer")

    """
    The toy_dataset is an Enum structure with some specific elements of
    interest:

        * **??**: ??.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <pclinfmri.preproc.pipeline.SliceTiming>` that
    define the different step of the processings:
    """
    pipeline = SpmNewSegment()
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.coregistered_struct_file = toy_dataset.mean

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=True, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


if __name__ == "__main__":
    pilot_smoothing()
    #pilot_bet()
    #pilot_newsegment()
