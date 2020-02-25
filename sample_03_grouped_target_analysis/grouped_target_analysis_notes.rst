# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

Example Grouped Target Analysis
===============================

This directory contains a sample analysis of Hyb-format data, published in 
the quick Crosslinking and Sequencing of Hybrids (qCLASH) experiment described in:
Gay, Lauren A., et al. "Modified cross-linking, ligation, and sequencing of hybrids 
(qCLASH) identifies Kaposi's Sarcoma-associated 
herpesvirus microRNA targets in endothelial cells." 
Journal of virology 92.8 (2018): e02138-17.

This analysis specifically investigates and characterizes miRNA arising from 
six experimental replicates from two conditions with cells infected with 
Kaposi's Sarcoma Herpesvirus, which are given the type name "KSHV_miRNA". 
The hybrid reads from KSHV miRNA are grouped and analyzed toghether.
Both individual and summary output files are produced.
 
Hybrid sequence information created by the Hyb program  information is 
available at NCBI Gene Expression Ombnibus (GEO) GSE101978, at:

    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978

The data files can be downloaded and uncompressed by using the command:

    "sh ./download_data.sh"

The unpacked hyb data-file require ~1.3 Gb of space.
The completed output of the analysis requires ~40 Mb of space.

Grouped Target Analysis Example Output
--------------------------------------

.. image:: /../sample_03_grouped_target_analysis/example_output/KSHV_Hyb_Combined_kshv-miR-K12-1star.png

.. image:: /../sample_03_grouped_target_analysis/example_output/KSHV_Hyb_Combined_kshv-miR-K12-1star_types.png
