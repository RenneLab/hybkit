# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

Example Target Analysis
=======================

This directory contains a sample analysis of Hyb-format data, published in 
the quick Crosslinking and Sequencing of Hybrids (qCLASH) experiment described in:
Gay, Lauren A., et al. "Modified cross-linking, ligation, and sequencing of hybrids 
(qCLASH) identifies Kaposi's Sarcoma-associated 
herpesvirus microRNA targets in endothelial cells." 
Journal of virology 92.8 (2018): e02138-17.

This analysis specifically investigates and characterizes miRNA arising from 
Kaposi's Sarcoma Herpesvirus, which are given the type name "KSHV_miRNA".
Both individual and summary output files are produced.
 
Hybrid sequence information created by the Hyb program  information is 
available at NCBI Gene Expression Ombnibus (GEO) GSE101978, at:

    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978

The data files can be downloaded and uncompressed by using the command:
  
    "sh ./download_data.sh"

The unpacked hyb data-file require ~130 Mb of space.
The completed output of the analysis requires ~20 Mb of space.

Target Analysis Example Output
------------------------------

.. image:: /../sample_02_target_analysis/example_output/GSM2720020_WT_BR1_KSHV_only_kshv-miR-K12-5.png

.. image:: /../sample_02_target_analysis/example_output/GSM2720020_WT_BR1_KSHV_only_kshv-miR-K12-5_types.png
