
About
=====

Renne Lab
---------
    | Principal Investigator: Rolf Renne
    | Henry E. Innes Professor of Cancer Research
    | University of Florida
    | UF Health Cancer Center
    | UF Department of Molecular Genetics and Microbiology
    | UF Genetics Institute
    | http://www.rennelab.com

Lead Developer
--------------
    * | Daniel Stribling <ds@ufl.edu>
      | https://www.danthescienceman.com
      | https://orcid.org/0000-0002-0649-9506
      | University of Florida, Renne Lab

Changelog
---------

    * 0.3.5 (2024-06) Changes include:

      * Misc Bugfixes

    * 0.3.4 (2023-11) Changes include:

      * Misc Bugfixes and Refinements
      * Switch code linting to Ruff
      * Add hybkit.errors module and HybkitError classes
      * Moved printing of warnings to python logging module
      * Add option for direct passage of file-like object for construction of
        HybFile and ViennaFile
      * Add HybRecord.to_csv_header(), HybRecord.to_fields(),
        and HybRecord.to_fields_header() methods
      * Refine description of HybFile.open() constructor method
      * Add typing_extensions dependency
      * Add Python3.8-compatible type hints to API
      * Documentation Updates

    * 0.3.3 (2023-09) Changes include:

      * Misc Bugfixes and Refinements
      * Update integer bar-plot functions

    * 0.3.2 (2023-08) Changes include:

      * Misc Bugfixes and Refinements
      * Add duplicate hybrid filtration (by HybRecord.id) options to hyb_filter
      * Add duplicate hybrid filtration to example analyses

    * 0.3.1 (2023-08) Changes include:

      * Misc Bugfixes and Refinements
      * Add --version flag to scripts
      * Change move scripts output file description to argparse epilog
      * Refine plot functions
      * Change default plot colors to the Bang Wong scheme [Wong2011]_ for
        colorblind accessibility
      * Documentation corrections
      * Spellcheck

    * 0.3.0 (2023-04) Major Codebase And API Overhaul. Changes include:

      * Simplifying HybRecord API
      * Simplifying FoldRecord API
      * Unifying settings information for argparse and modules
      * Removing Support for ViennaD format
      * Moving identifier-parsing code to module type_finder
      * Moving target region analysis code to module region_finder
      * Moving code for settings into a "settings" module
      * Renamed HybRecord type_analysis and mirna_analysis to
        eval_types and eval_mirna, respectively
        to differentiate from analysis module functions
      * Reimplemented analyses methods within a single Analysis class
      * Added error checking / catching to HybFoldIter
      * Removed Target-Region Analysis and associated files
        due to lack of archival database information,
        pending future development
      * Added "dynamic" seq_type to FoldRecord for non-identical fold/hybrid sequence handling
      * Added shell implementation to all example analyses
      * Remove support for Python3.6, Python3.7
      * Migrate to CircleCI for CI/CD
      * Added Pytest unit testing integrated with CircleCI
      * Other Misc. Improvements / Bugfixes

    * 0.2.0  (2020-03) Added Command-line Toolkit. Code Refinements.

    * 0.1.9  (2020-03) Fix for Module Path Finding for Python > 3.6

    * 0.1.8  (2020-03) Streamlining, PyPI / PIP Initial Release

    * 0.1.0  (2020-01) Initial Implementation


References
==========

    .. [ViennaFormat]
         | `ViennaRNA Vienna File Format Description <https://www.tbi.univie.ac.at/RNA/tutorial/#sec2_7>`_
         | `UNAFold Vienna File Format Description <http://www.unafold.org/doc/formats.php#VIENNA>`_

    .. [CTFormat]
          | `UNAFold CT Format Description <http://www.unafold.org/doc/formats.php#CT>`_
          | `RNAStructure CT Format Description
            <https://rna.urmc.rochester.edu/Text/File_Formats.html#CT>`_
    .. [Zuker2003] Zuker M. Mfold web server for nucleic acid folding and hybridization
          prediction. Nucleic Acids Res. 2003 Jul 1;31(13):3406-15.
          doi: `10.1093/nar/gkg595 <https://doi.org/10.1093/nar/gkg595>`_.
          PMID: 12824337; PMCID: PMC169194.
    .. [Hunter2007] J. Hunter, "Matplotlib: A 2D Graphics Environment" in Computing in
           Science & Engineering, vol. 9, no. 03, pp. 90-95, 2007.
           doi: `10.1109/MCSE.2007.55 <https://doi.org/10.1109/MCSE.2007.55>`_
    .. [Cock2009] Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I,
           Hamelryck T, Kauff F, Wilczynski B, de Hoon MJ. Biopython: freely available
           Python tools for computational molecular biology and bioinformatics. Bioinformatics.
           2009 Jun 1;25(11):1422-3. doi:
           `10.1093/bioinformatics/btp163 <https://doi.org/10.1093/bioinformatics/btp163>`_.
           Epub 2009 Mar 20.
           PMID: 19304878; PMCID: PMC2682512.
    .. [Lorenz2011] Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C. et al.
           ViennaRNA Package 2.0. Algorithms Mol Biol 6, 26 (2011).
           doi: `10.1186/1748-7188-6-26 <https://doi.org/10.1186/1748-7188-6-26>`_
    .. [Wong2011] Wong, B. Points of view: Color blindness. Nat Methods 8, 441 (2011).
           doi: `10.1038/nmeth.1618 <https://doi.org/10.1038/nmeth.1618>`_
    .. [Helwak2013] Helwak A, Kudla G, Dudnakova T, Tollervey D. Mapping the human miRNA
           interactome by CLASH reveals frequent noncanonical binding. Cell. 2013
           Apr 25;153(3):654-65. doi:
           `10.1016/j.cell.2013.03.043 <https://doi.org/10.1016/j.cell.2013.03.043>`_.
           PMID: 23622248; PMCID: PMC3650559.
    .. [Travis2014] Travis AJ, et al. Hyb: a bioinformatics pipeline for the analysis of
           CLASH (crosslinking, ligation and sequencing of hybrids) data.
           Methods. 2014 Feb;65(3):263-73.
           doi: `10.1016/j.ymeth.2013.10.015 <https://doi.org/10.1016/j.ymeth.2013.10.015>`_.
    .. [Gay2018] Gay LA, Sethuraman S, Thomas M, Turner PC, Renne R. Modified Cross-Linking,
           Ligation, and Sequencing of Hybrids (qCLASH) Identifies Kaposi's
           Sarcoma-Associated Herpesvirus MicroRNA Targets in Endothelial Cells.
           J Virol. 2018 Mar 28;92(8):e02138-17.
           doi: `10.1128/JVI.02138-17 <https://doi.org/10.1128/JVI.02138-17>`_.
           PMID: 29386283; PMCID: PMC5874430.


    * [ViennaFormat]_
    * [CTFormat]_
    * [Zuker2003]_
    * [Hunter2007]_
    * [Cock2009]_
    * [Lorenz2011]_
    * [Wong2011]_
    * [Helwak2013]_
    * [Travis2014]_
    * [Gay2018]_





