
hybkit Toolkit
==================================

    The hybkit toolkit contains command-line scripts for analysis and manipulation of
    hyb and fold files. Scripts are implemented in Python3, and the hybkit module must be on the
    user's PYTHONPATH for script execution.

    The command-line options and flags are generated with the
    Python3 argparse module. Relevant settings pertaining to specific hybkit classes are accessible
    via command-line flags, as demonstrated in the "shell_analysis" implementations in the
    :ref:`Example Analyses`.

    This version of hybkit includes the following executables:


        =================================== ===========================================================
        Utility                             Description
        =================================== ===========================================================
        :ref:`hyb_check`                    Parse a hyb (/fold) file and check for errors
        :ref:`hyb_eval`                     Evaluate hyb (/fold) records to identify segment types and miRNAs
        :ref:`hyb_filter`                   Filter a hyb (/fold) file to a specific subset of sequences
        :ref:`hyb_analyze`                  Perform a type, miRNA, summary, or target analysis
                                            on a hyb (/fold) file
        =================================== ===========================================================

    Detailed descriptions and usage information are available at each respective script page.

.. toctree::
   :maxdepth: 4
   :caption: Hybkit Toolkit Executables:

   toolkit/hyb_check
   toolkit/hyb_filter
   toolkit/hyb_eval
   toolkit/hyb_analyze


