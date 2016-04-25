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
def pilot_bet(enable_display=False):
    """
    Brain extractio Tool
    ====================
    """
    import os
    from mmutils.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from capsul.api import get_process_instance

    working_dir = "/volatile/nsap/catalogue/pclinfmri/fmri_bet"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Then define the study configuration:
    """
    study_config = StudyConfig(
        modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=False,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
        matlab_exec="/neurospin/local/bin/matlab",
        use_matlab=False,
        spm_directory="/i2bm/local/spm8",
        use_spm=False,
        output_directory=working_dir)

    """
    Load the toy dataset
    --------------------

    To do so, we use the get_sample_data function to download the toy
    dataset on the local file system (here localizer data):
    """
    template_dataset = get_sample_data("mni_1mm")

    """
    Processing definition
    ---------------------
    """
    pipeline = get_process_instance("clinfmri.utils.converted_fsl_bet")
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.input_image_file = template_dataset.brain
    pipeline.generate_binary_mask = True

    """
    It is possible to display the pipeline.
    """
    if enable_display:
        import sys
        from PySide import QtGui
        from capsul.qt_gui.widgets import PipelineDevelopperView

        app = QtGui.QApplication(sys.argv)
        view = PipelineDevelopperView(pipeline)
        view.show()
        app.exec_()

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=False, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


@pilotfunction
def pilot_preproc_spm_fmri(enable_display=False):
    """
    FMRI preprocessings
    ===================

    Preprocessing with the SPM slice timing and a normalization to a given
    template.

    Start to import required modules:
    """
    import os
    from mmutils.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from capsul.api import get_process_instance

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/catalogue/pclinfmri/fmri_preproc_spm_fmri"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Then define the study configuration:
    """
    study_config = StudyConfig(
        modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=False,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
        matlab_exec="/neurospin/local/bin/matlab",
        use_matlab=True,
        spm_directory="/i2bm/local/spm8",
        use_spm=True,
        output_directory=working_dir,
        number_of_cpus=1,
        generate_logging=True,
        use_scheduler=True,)

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

        * fmri: the functional volume.
        * anat: the structural volume.
        * TR: the repetition time.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <clinfmri.preproc.FmriPreproc>` that
    define the different step of the processings:
    """
    pipeline = get_process_instance("clinfmri.preproc.converted_fmri_preproc")
    print pipeline.get_input_spec()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.fmri_file = toy_dataset.fmri
    pipeline.structural_file = toy_dataset.anat
    pipeline.realign_register_to_mean = True
    pipeline.select_slicer = "spm"
    pipeline.select_normalization = "fmri"
    pipeline.template_file = template_dataset.brain
    pipeline.force_repetition_time = toy_dataset.TR
    pipeline.force_slice_orders = [index + 1 for index in range(40)]

    """
    It is possible to display the pipeline.
    """
    if enable_display:
        import sys
        from PySide import QtGui
        from capsul.qt_gui.widgets import PipelineDevelopperView

        app = QtGui.QApplication(sys.argv)
        view = PipelineDevelopperView(pipeline)
        view.show()
        app.exec_()

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=False, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


@pilotfunction
def pilot_preproc_fsl_anat(enable_display=False):
    """
    FMRI preprocessings
    ===================

    Preprocessing with the FSL slice timing and the unified registration/class
    segmentation to the MNI template.

    Start to import required modules:
    """
    import os
    from mmutils.toy_datasets import get_sample_data
    from capsul.study_config import StudyConfig
    from capsul.process import get_process_instance

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/catalogue/pclinfmri/fmri_preproc_fsl_anat"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    Then define the study configuration:
    """
    if 0:
        study_config = StudyConfig(
            modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
            use_smart_caching=False,
            fsl_config="/etc/fsl/4.1/fsl.sh",
            use_fsl=True,
            matlab_exec="/neurospin/local/bin/matlab",
            use_matlab=True,
            spm_directory="/i2bm/local/spm8",
            use_spm=True,
            output_directory=working_dir,
            number_of_cpus=1,
            generate_logging=True,
            use_scheduler=True)

    study_config = StudyConfig(
        modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=False,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
        use_matlab=False,
        use_spm=False,
        spm_exec="/i2bm/local/bin/spm12",
        spm_standalone=True,
        use_nipype=True,
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

        * fmri: the functional volume.
        * anat: the structural volume.
        * TR: the repetition time.

    Processing definition
    ---------------------

    First create the
    :ref:`slice timing pipeline <clinfmri.preproc.FmriPreproc>` that
    define the different step of the processings:
    """
    pipeline = get_process_instance("clinfmri.preproc.fmri_preproc")
    print pipeline.get_input_spec()

    """
    It is possible to display the pipeline.
    """
    if enable_display:
        import sys
        from PySide import QtGui
        from capsul.qt_gui.widgets import PipelineDevelopperView

        app = QtGui.QApplication(sys.argv)
        view = PipelineDevelopperView(pipeline)
        view.show()
        app.exec_()

    """
    Now we need now to parametrize this pipeline:
    """
    pipeline.fmri_file = toy_dataset.fmri
    pipeline.structural_file = toy_dataset.anat
    pipeline.realign_register_to_mean = True
    pipeline.select_slicer = "fsl"
    pipeline.select_normalization = "anat"
    pipeline.force_repetition_time = toy_dataset.TR
    pipeline.force_slice_orders = [index + 1 for index in range(40)]

    """
    The pipeline is now ready to be run:
    """
    study_config.run(pipeline, executer_qc_nodes=False, verbose=1)

    """
    Results
    -------

    Finally, we print the pipeline outputs:
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in pipeline.get_outputs().items():
        print("{0}: {1}".format(trait_name, trait_value))


if __name__ == "__main__":
    #pilot_bet(enable_display=True)
    pilot_preproc_spm_fmri(enable_display=True)
    #pilot_preproc_fsl_anat(enable_display=True)
