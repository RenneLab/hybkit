******
hybkit
******
.. image:: https://img.shields.io/github/v/release/RenneLab/hybkit?include_prereleases&logo=github
   :target: https://github.com/RenneLab/hybkit/releases
   :alt: GitHub release (latest by date including pre-releases)
.. image:: https://img.shields.io/pypi/v/hybkit?logo=pypi&logoColor=white
   :target: https://pypi.org/project/hybkit/
   :alt: PyPI Package Version
.. image:: https://img.shields.io/conda/vn/bioconda/hybkit?logo=anaconda
   :target: http://bioconda.github.io/recipes/hybkit/README.html
   :alt: Bioconda Release
.. image:: https://img.shields.io/conda/dn/bioconda/hybkit?logo=Anaconda
   :target: http://bioconda.github.io/recipes/hybkit/README.html
   :alt: Hybkit Downloads on Bioconda
.. image:: https://img.shields.io/conda/vn/bioconda/hybkit?color=lightgrey&label=Image%20%28quay.io%29&logo=docker
   :target: https://quay.io/repository/biocontainers/hybkit?tab=tags
   :alt: Docker Image Version
.. image:: https://img.shields.io/circleci/build/github/RenneLab/hybkit?label=CircleCI&logo=circleci
   :target: https://app.circleci.com/pipelines/github/RenneLab/hybkit
   :alt: Circle-CI Build Status
.. image:: https://img.shields.io/readthedocs/hybkit?logo=read-the-docs
   :target: https://hybkit.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/coveralls/github/RenneLab/hybkit?logo=coveralls
   :target: https://coveralls.io/github/RenneLab/hybkit
   :alt: Coveralls Coverage
.. image:: https://img.shields.io/pypi/pyversions/hybkit?logo=python&logoColor=white
   :target: https://pypi.org/project/hybkit/
   :alt: PyPI - Python Version
.. image:: https://img.shields.io/badge/License-GPLv3+-blue?logo=GNU
   :target: https://www.gnu.org/licenses/gpl-3.0.en.html
   :alt: GNU GPLv3+ License

| Welcome to *hybkit*, a toolkit for analysis of hyb-format chimeric
  (hybrid) RNA sequence data defined with the Hyb software package by |Travis2014|_.
  This genomic data-type is generated from RNA proximity-ligation and ribonomics
  techniques such as Crosslinking, Ligation, and
  Sequencing of Hybrids (CLASH; |Helwak2013|_) and quick CLASH (qCLASH; |Gay2018|_).
| This software is available via Github, at http://www.github.com/RenneLab/hybkit .
| Full project documentation is available at |docs_link|_.

Project components:
    #. hybkit toolkit of command-line utilities for manipulating,
       analyzing, and plotting hyb-format data.
    #. The hybkit python API, an extendable documented codebase
       for creation of custom analyses of hyb-format data.
    #. Integrated analysis of predicted secondary structure (fold) information for
       the API and command-line utilities.
    #. Example analyses for publicly available qCLASH hybrid
       sequence data implemented in each of the command-line scripts and hybkit Python API.

Hybkit Toolkit:
    The hybkit toolkit includes several command-line utilities
    for manipulation of hyb-format data:

        =================================== ===========================================================
        Utility                             Description
        =================================== ===========================================================
        hyb_check                           Parse hyb (and fold) files and check for errors
        hyb_eval                            Evaluate hyb (and fold) records to identify / assign
                                            segment types and miRNAs using custom criteria
        hyb_filter                          Filter hyb (and fold) records to a specific
                                            custom subset
        hyb_analyze                         Perform an energy, type, miRNA, target, or fold analysis
                                            on hyb (and fold) files and plot results
        =================================== ===========================================================

    These scripts are used on the command line with hyb (and associated "vienna" or "CT") files.
    For example, to filter a
    hyb and corresponding vienna file to contain only hybrids with
    a sequence identifier containing the string "kshv":

    Example:

        ::

            $ hyb_filter -i my_hyb_file.hyb -f my_hyb_file.vienna --filter any_seg_contains kshv

    Further detail on the usage of each script is provided in
    the |hybkit Toolkit| section of |docs_link|_.


Hybkit API:
    Hybkit provides a Python3 module with a documented API for interacting with
    records in hyb files and associated vienna or CT files.
    This capability was inspired by the `BioPython Project <https://biopython.org/>`_.
    The primary utility is provided by a class for hyb records (HybRecord), a class
    for fold records (FoldRecord), and file-iterator classes
    (HybFile, ViennaFile, CTFile, HybFoldIter).
    Record attributes can be analyzed, set, and evaluated using included class methods.

    For example, a workflow to print the identifiers of only sequences within a ".hyb" file
    that contain a miRNA can be performed as such:

    .. code-block:: Python

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

Example Analyses:
    Hybkit provides several example analyses for hyb data using the
    utilities provided in the toolkit. These include:

        ============================= ===========================================================
        Analysis                      Description
        ============================= ===========================================================
        Type/miRNA Analysis           Quantify sequence types and miRNA types in a hyb file
        Target Analysis               Analyze targets of a set of miRNAs from a single
                                      experimental replicate
        Grouped Target Analysis       Analyze and plot targets of a set of miRNAs from
                                      pooled experimental replicates
        Fold Analysis                 Analyze and plot predicted miRNA folding patterns in
                                      miRNA-containing hybrids
        ============================= ===========================================================

    These analyses provide analysis results in both tabular and graph form.
    As an illustration, the example summary analysis includes the return of
    the contained hybrid sequence types as both a csv table and as a pie chart:

        `CSV Output <https://raw.githubusercontent.com/RenneLab/hybkit/master/example_01_type_mirna_analysis/example_output/combined_analysis_type_hybrid_types.csv>`_

        |example_01_image|

    Further detail on each provided analysis can be found in
    the |Example Analyses| section of |docs_link|_.

Installation:
    Dependencies:
        * Python3.8+
        * `matplotlib <https://matplotlib.org/>`_ >= 3.7.1 (|Hunter2007|_)
        * `BioPython <https://biopython.org/>`_ >= 1.79 (|Cock2009|_)
        * `typing_extensions <https://pypi.org/project/typing-extensions/>` >= 4.8.0

    Via PyPI / Python PIP:
        |PipVersion|

        The recommended installation method is via hybkit's
        `PyPI Package Index <https://pypi.org/project/hybkit/>`_ using
        `python3 pip <https://pip.pypa.io/en/stable/>`_, which will
        automatically handle version control and dependency installation:

        .. code-block:: bash

            $ python3 -m pip install hybkit

    Via Conda:
        |CondaVersion| |InstallBioconda|

        For users of conda, the hybkit package and dependencies are hosted on the
        the `Bioconda <https://bioconda.github.io/>`_ channel, and can be installed
        using conda:

        .. code-block:: bash

            $ conda install -c bioconda hybkit

    Via Docker/Singularity:
        |DockerVersion|

        The hybkit package is also available as a `Docker <https://www.docker.com/>`_
        image and `Singularity <https://sylabs.io/singularity/>`_ container, hosted
        via the `BioContainers <https://biocontainers.pro/>`_ project on
        `quay.io <https://quay.io/repository/biocontainers/hybkit?tab=tags>`_.
        To pull the image via docker:

        .. code-block:: bash

            $ docker pull quay.io/biocontainers/hybkit:0.3.3--pyhdfd78af_0

        To pull the image via singularity:

        .. code-block:: bash

            $ singularity pull docker://quay.io/biocontainers/hybkit:0.3.3--pyhdfd78af_0

    Manually Download and Install:
        |GithubVersion|

        Use git to clone the project's Github repository:

        .. code-block:: bash

            $ git clone git://github.com/RenneLab/hybkit

        *OR* download the zipped package:

        .. code-block:: bash

            $ curl -OL https://github.com/RenneLab/hybkit/archive/master.zip
            $ unzip master.zip

        Then install using python setuptools:

        .. code-block:: bash

            $ python setup.py install

    Further documentation on hybkit usage can be found in |docs_link|_.

Setup Testing:
    Hybkit provides a suite of unit tests to maintain stability of the API and script
    functionalities. To run the API test suite, install pytest and run the tests from the
    root directory of the hybkit package:

    .. code-block:: bash

        $ pip install pytest
        $ pytest

    Command-line scripts can be tested by running the auto_test.sh script in
    the auto_tests directory:

    .. code-block:: bash

        $ ./auto_tests/auto_test.sh


Copyright:
    | hybkit is a free, sharable, open-source project.
    | All source code and executable scripts contained within this package are considered
        part of the "hybkit" project and are distributed without any warranty or implied warranty
        under the GNU General Public License v3.0 or any later version, described in the "LICENSE"
        file.

.. |Helwak2013| replace:: *Helwak et al. (Cell 2013)*
.. _Helwak2013: https://doi.org/10.1016/j.cell.2013.03.043
.. |Travis2014| replace:: *Travis et al. (Methods 2014)*
.. _Travis2014: https://doi.org/10.1016/j.ymeth.2013.10.015
.. |Gay2018| replace:: *Gay et al. (J. Virol. 2018)*
.. _Gay2018: https://doi.org/10.1128/JVI.02138-17
.. |Hunter2007| replace:: *Hunter JD. (Computing in Science & Engineering 2007)*
.. _Hunter2007: https://doi.org/10.1109/MCSE.2007.55
.. |Cock2009| replace:: *Cock et al. (Bioinformatics 2009)*
.. _Cock2009: https://doi.org/10.1093/bioinformatics/btp163
.. |PipVersion| image:: https://img.shields.io/pypi/v/hybkit?logo=pypi&logoColor=white
   :target: https://pypi.org/project/hybkit/
   :alt: PyPI Package Version
.. |InstallBioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&logo=anaconda
   :target: http://bioconda.github.io/recipes/hybkit/README.html
   :alt: Install with Bioconda
.. |CondaVersion| image:: https://img.shields.io/conda/vn/bioconda/hybkit?logo=anaconda
   :target: http://bioconda.github.io/recipes/hybkit/README.html
   :alt: Bioconda Release
.. |DockerVersion| image:: https://img.shields.io/conda/vn/bioconda/hybkit?color=lightgrey&label=Image%20%28quay.io%29&logo=docker
   :target: https://quay.io/repository/biocontainers/hybkit?tab=tags
   :alt: Docker Image Version
.. |GithubVersion| image:: https://img.shields.io/github/v/release/RenneLab/hybkit?include_prereleases&logo=github
   :target: https://github.com/RenneLab/hybkit/releases
   :alt: GitHub release (latest by date including pre-releases)

.. include:: ./index_suffix.rst
