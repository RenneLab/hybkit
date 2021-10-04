
hybkit Toolkit
==================================

    The hybkit toolkit contains command-line scripts for analysis and manipulation of 
    hyb and fold files. Scripts are implemented in Python3, and the hybkit module must be on the 
    user's PYTHONPATH for script execution.

    Effort has been taken to document each of the command-line options and flags using the 
    Python3 argparse module. Relevant settings pertaining to specific hybkit classes are accessible
    via command-line flags, as demonstrated in the "shell_analysis" implementations in the 
    :ref:`Example Pipelines`. This 

    This version of hybkit includes the following executables:


        =================================== ===========================================================
        Utility                             Description
        =================================== ===========================================================
        :ref:`hyb_check`                    Parse a hyb file and check for errors
        :ref:`hyb_eval`                     Evaluate hyb records to identify segment types and miRNAs
        :ref:`hyb_filter`                   Filter a hyb file to a specific subset of sequences
        :ref:`hyb_analyze`                  Perform a type, miRNA, summary, or target analysis 
                                            on a hyb file
        :ref:`hyb_exclude_fold`             Filter a fold file using an exclusion table created by
                                            hyb_filter
        :ref:`hyb_fold_analyze`             Perform a fold analysis on a hyb and a RNA secondary
                                            structure (fold) file in ".vienna" or ".ct" format
                                            on a hyb file
        =================================== ==========================================================

    Detailed descriptions and usage information are avilable at each respective script page.

.. toctree::
   :maxdepth: 4
   :caption: Hybkit Toolkit Executables:

   toolkit/hyb_check
   toolkit/hyb_filter
   toolkit/hyb_eval
   toolkit/hyb_analyze
   toolkit/hyb_fold_analyze
   toolkit/hyb_exclude_fold


