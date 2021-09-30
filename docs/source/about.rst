
References
==========

    #. "Travis, Anthony J., et al. "Hyb: a bioinformatics pipeline for the analysis of CLASH
       (crosslinking, ligation and sequencing of hybrids) data."
       Methods 65.3 (2014): 263-273."
    #. Gay, Lauren A., et al. "Modified cross-linking, ligation, and sequencing of
       hybrids (qCLASH) identifies Kaposi's Sarcoma-associated herpesvirus microRNA
       targets in endothelial cells." Journal of virology 92.8 (2018): e02138-17.
    #. The Vienna File Format: http://unafold.rna.albany.edu/doc/formats.php#VIENNA


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

    * | 0.3.0b (dev)     Further Refactoring:
    * | 0.3.0a (2021-09) Major Codebase And API Overhaul. Changes include:
      | Simplifying HybRecord API
      | Simplifying FoldRecord API
      | Unifying settings information for argparse and module
      | Removing Support for ViennaD format
      | Addition of "NON" value to "target_reg" flag in specification
      | Moving identifier-parsing code to module type_finder
      | Moving target region analysis code to module region_finder
      | Moving code for settings into a "settings" module.
      | Renamed HybRecord type_analysis, mirna_analysis, and target_analysis to 
        eval_types, eval_mirna, eval_target, respectively
        to differentiate from analysis module functions
      | Reimplemented analyses methods as classes.
        
    * 0.2.0  (2020-03) Added Command-line Toolkit. Code Refinements.
    * 0.1.9  (2020-03) Fix for Module Path Finding for Python > 3.6
    * 0.1.8  (2020-03) Streamlining, PyPI / PIP Initial Release
    * 0.1.0  (2020-01) Initial Implementation




