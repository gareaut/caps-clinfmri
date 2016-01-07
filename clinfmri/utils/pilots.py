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
def pilot_bet(enable_display=False):
    """ 
    BET
    ===

    Brain extraction with FSL. 

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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/fsl_bet"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration:
    """
    study_config = StudyConfig(
        modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=False,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,
        output_directory=working_dir,
        number_of_cpus=1,
        generate_logging=True,
        use_scheduler=True)

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
    :ref:`slice timing pipeline <clinfmri.preproc.FslBet>` that
    define the different step of the processings:
    """
    pipeline = get_process_instance("clinfmri.utils.fsl_bet.xml")
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
    pipeline.input_image_file = toy_dataset.anat
    pipeline.generate_binary_mask = True
    pipeline.bet_threshold = 0.5

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
def pilot_new_segment(enable_display=False):
    """ 
    New Segment
    ===========

    Unifed SPM segmentation: segments, bias corrects and spatially normalises. 

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
    working_dir = "/volatile/nsap/catalogue/pclinfmri/spm_newsegment"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration:
    """
    study_config = StudyConfig(
        modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
        use_smart_caching=False,
        matlab_exec="/neurospin/local/bin/matlab",
        use_matlab=True,
        spm_directory="/i2bm/local/spm8",
        use_spm=True,
        output_directory=working_dir,
        number_of_cpus=1,
        generate_logging=True,
        use_scheduler=True)

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
    :ref:`slice timing pipeline <clinfmri.utils.SpmNewSegment>`
    that define the different step of the processings:
    """
    pipeline = get_process_instance("clinfmri.utils.spm_new_segment.xml")
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
    pipeline.channel_files = [toy_dataset.mean]
    pipeline.reference_volume = template_dataset.brain

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
    pilot_bet(enable_display=True)
    pilot_new_segment(enable_display=True)
