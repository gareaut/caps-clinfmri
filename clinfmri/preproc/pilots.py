#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2015
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
def pilot_slice_timing():
    """
    Slice Timing Correction
    =======================
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/slicetiming"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance("clinfmri.preproc.slice_timing.xml")

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "MatlabConfig", "SPMConfig",
                 "FSLConfig",
                 "NipypeConfig"],
        use_smart_caching=True,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
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
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.fmri_file = toy_dataset.fmri
    pipeline.select_slicer = "fsl"
    pipeline.force_repetition_time = toy_dataset.TR
    pipeline.force_slice_orders = [index + 1 for index in range(40)]

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
def pilot_realignement():
    """
    Realignement
    ============
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/realignement"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance("clinfmri.preproc.spm_realignement.xml")

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
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.fmri_file = toy_dataset.fmri
    pipeline.register_to_mean = True

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
def pilot_coregistration():
    """
    Coregistration
    ==============
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/coregistration"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance("clinfmri.preproc.spm_coregistration.xml")

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
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.reference_image = toy_dataset.mean
    pipeline.moving_image = toy_dataset.anat
    pipeline.fwhm = [7, 7]
    pipeline.jobtype = "estwrite"

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
def pilot_normalization():
    """
    Normalization
    =============
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/normalization"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance("clinfmri.preproc.spm_normalization.xml")

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
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.reference_image = toy_dataset.mean
    pipeline.coregistered_struct_file = toy_dataset.mean
    pipeline.fmri_file = toy_dataset.fmri

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
def pilot_template_registration():
    """
    Template Registration
    =====================
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/template_registration"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance(
        "clinfmri.preproc.pipeline.spm_template_registration.xml")

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
    template_dataset = get_sample_data("mni_1mm")

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
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.template_file = template_dataset.brain
    pipeline.coregistered_struct_file = toy_dataset.mean
    pipeline.fmri_file = toy_dataset.fmri

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
def pilot_preproc():
    """
    FMRI preprocessings
    ===================
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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/fmri_preproc"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Now we get the pipeline from its definition (xml file)
    """
    pipeline = get_process_instance(
        "clinfmri.preproc.pipeline.fmri_preproc.xml")

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "MatlabConfig", "SPMConfig",
                 "FSLConfig",
                 "NipypeConfig"],
        use_smart_caching=True,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
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
    template_dataset = get_sample_data("mni_1mm")

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
    pipeline = get_process_instance("pclinfmri.preproc.fmri_preproc.xml")
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.fmri_file = toy_dataset.fmri
    pipeline.structural_file = toy_dataset.anat
    pipeline.realign_register_to_mean = True
    pipeline.select_slicer = "none"
    pipeline.select_registration = "template"
    pipeline.template_file = template_dataset.brain
    pipeline.force_repetition_time = toy_dataset.TR
    pipeline.force_slice_orders = [index + 1 for index in range(40)]

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
    #pilot_slice_timing()
    #pilot_realignement()
    #pilot_coregistration()
    #pilot_normalization()
    pilot_template_registration()
    #pilot_preproc()
