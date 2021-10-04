******
hybkit
******
.. image:: https://img.shields.io/github/v/release/RenneLab/hybkit?include_prereleases&logo=github
   :target: https://github.com/RenneLab/hybkit/releases
   :alt: GitHub release (latest by date including pre-releases)
.. image:: https://img.shields.io/pypi/v/hybkit?logo=pypi&logoColor=white
   :target: https://pypi.org/project/hybkit/
   :alt: PyPI Package Version
.. image:: https://img.shields.io/readthedocs/hybkit?logo=read-the-docs
   :target: https://hybkit.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/pypi/pyversions/hybkit?logo=python&logoColor=white
   :target: https://pypi.org/project/hybkit/
   :alt: PyPI - Python Version
.. image:: https://img.shields.io/badge/License-GPLv3+-blue?logo=GNU
   :target: https://www.gnu.org/licenses/gpl-3.0.en.html
   :alt: GNU GPLv3+ License

| Welcome to *hybkit*, a toolkit for analysis of ".hyb" (hyb-file) format chimeric 
  (hybrid) RNA sequence data defined with the Hyb software package by |Travis2014|_.
  This genomic data-type is generated from RNA proximitiy-ligation and ribonomics 
  techniques such as Crosslinking, Ligation, and 
  Sequencing of Hybrids (CLASH; |Helwak2013|_) and Quick CLASH (qCLASH; |Gay2018|_). 
| This software is available via Github, at http://www.github.com/RenneLab/hybkit .
| Full project documentation is available at |docs_link|_.

Project components:
    #. hybkit toolkit of command-line utilities for manipulating,
       analyzing, and plotting hyb-format data.
    #. The hybkit python API, an extendable documented codebase
       for creation of custom analyses of hyb-format data.
    #. Example analysis pipelines for analysis of publically available qCLASH hybrid 
       sequence data implemented either with toolkit scripts or the hybkit Python API.

Hybkit Toolkit:
    The hybkit toolkit includes several command-line utilities 
    for manipulation of ".hyb" format data:

        =================================== ===========================================================
        Utility                             Description
        =================================== ===========================================================
        hyb_check                           Parse a hyb file and check for errors
        hyb_eval                            Evaluate hyb records to identify segment types and miRNAs
        hyb_filter                          Filter a hyb file to a specific subset of sequences
        hyb_analyze                         Perform a type, miRNA, summary, or target analysis 
                                            on a hyb file
        hyb_exclude_fold                    Filter a fold file using an exclusion table created by
                                            hyb_filter
        hyb_fold_analyze                    Perform a fold analysis on a hyb and a RNA secondary
                                            structure (fold) file in ".vienna" or ".ct" format
                                            on a hyb file
        =================================== ===========================================================
        
    These scripts are used on the command line with hyb files. For example, to filter a 
    hyb file to contain only hybrids with a sequence identifier containing the string "kshv"

    Example:

        ::

            $ hyb_filter -i my_hyb_file.hyb --filter any_seg_contains kshv

    Further detail on the usage of each script is provided in 
    the |hybkit Toolkit| section of |docs_link|_.


Hybkit API:
    Hybkit provides a Python3 module with a documented API for interacting with 
    records in hyb files. 
    This capability was inspired by the `BioPython Project <https://biopython.org/>`_. 
    The primary utility is provided by classes for hyb records (HybRecord), classes
    for fold records (FoldRecord, DynamicFoldRecord), and file-iterator classes 
    (HybFile, ViennaFile, CTFile).
    Record attributes can be analyzed, set, and evaluated using included class methods.

    For example, a workflow to print the identifiers of only sequences within a ".hyb" file
    that contain a miRNA can be performed as such::

        #!/usr/bin/env python3
        import hybkit
        in_file = '/path/to/my_hyb_file.hyb'

        # Open a hyb file as a HybFile Object:
        with hybkit.HybFile.open(in_file, 'r') as hyb_file:

            # Return each line in a hyb file as a HybRecord object
            for hyb_record in hyb_file:

                # Analyze each record to assign segment types
                hyb_record.eval_types()

                # If the record contains a long noncoding RNA type, print the record identifier.
                if hyb_record.has_prop('any_seg_type_contains', 'lncRNA')
                    print(hyb_record.id)

    Further documentation on the hybkit API can be found in the 
    |hybkit API| section of |docs_link|_.

Pipelines:
    Hybkit provides several example pipelines for analysis of "hyb" data using the 
    utilities provided in the toolkit. These include:
    
        ============================= ===========================================================
        Pipeline                      Description
        ============================= ===========================================================
        Summary Analysis              Quantify the sequence and miRNA types in a hyb file
        Target Analysis               Analyze targets of a set of miRNAs from a single 
                                      experiment
        Grouped Target Analysis       Analyze and plot targets of a set of miRNAs from 
                                      pooled experimental replicates
        Fold Analysis                 Analyze and plot predicted miRNA folding patterns in
                                      miRNA-containing hybrids
        ============================= ===========================================================

    These pipelines provide analysis results in both tabular and graph form.
    As an illustration, the example summary analysis includes the return of 
    the contained hybrid sequence types as both a csv table and as a pie chart:

        `CSV Output <https://raw.githubusercontent.com/RenneLab/hybkit/master/sample_01_summary_analysis/example_output/combined_analysis_type_hybrids.csv>`_

        |sample_01_image|

    Further detail on each provided pipeline can be found in 
    the |Example Pipelines| section of |docs_link|_.

Installation:
    
    Dependencies:
        * Python3.6+
        * `matplotlib <https://matplotlib.org/>`_ 
        * `BioPython <https://biopython.org/>`_ 
    
    The recommended installation method is via hybkit's 
    `PyPI Package Index <https://pypi.org/project/hybkit/>`_ using 
    `python3 pip <https://pip.pypa.io/en/stable/>`_, which will 
    automatically handle version control and dependency installation::
        
        $ pip install hybkit
    
    Acquisition of the package can also be performed by cloning the project's Github repository::
    
        $ git clone git://github.com/RenneLab/hybkit
    
    Or by downloading the zipped package::
    
        $ curl -OL https://github.com/dstrib/hybkit/archive/master.zip
        $ unzip master.zip
    
    Followed by installation using python's setuptools::
    
        $ python setup.py install
    
    Further documentation on hybkit usage can be found in |docs_link|_.

.. |hybkit Toolkit| replace:: *hybkit Toolkit*
.. |Example Pipelines| replace:: *Example Pipelines*
.. |hybkit API| replace:: *hybkit API*
.. |docs_link| replace:: hybkit's ReadTheDocs
.. _docs_link: https://hybkit.readthedocs.io#
.. |Helwak2013| replace:: *Helwak et al. (Cell 2013)*
.. _Helwak2013: https://doi.org/10.1016/j.cell.2013.03.043
.. |Travis2014| replace:: *Travis et al. (Methods 2014)*
.. _Travis2014: https://doi.org/10.1016/j.ymeth.2013.10.015 
.. |Gay2018| replace:: *Gay et al. (J. Virol. 2018)*
.. _Gay2018: https://doi.org/10.1128/JVI.02138-17
.. |sample_01_image| image:: sample_01_summary_analysis/example_output/combined_analysis_type_hybrids.png

.. include:: docs_readme_format.rst
