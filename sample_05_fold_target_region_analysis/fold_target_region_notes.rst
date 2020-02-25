# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

Example Target Region Analysis
==============================

This directory contains a sample analysis of Hyb-format and Viennad-format data, published in 
the quick Crosslinking and Sequencing of Hybrids (qCLASH) experiment described in:
Gay, Lauren A., et al. "Modified cross-linking, ligation, and sequencing of hybrids 
(qCLASH) identifies Kaposi's Sarcoma-associated 
herpesvirus microRNA targets in endothelial cells." 
Journal of virology 92.8 (2018): e02138-17.

This analysis investigates the folding pattern of miRNA 
from an experimental replicate infected with 
Kaposi's Sarcoma Herpesvirus, which are given the type name "KSHV_miRNA". 
Data from the predicted folding pattern for each hybrid record produced 
by the "Hyb" program are analyzed. Coding targets are identified using a database
prepared from the Ensembl biomart, and included in the databases section. 
The region of the coding target is then determined. 
Records containing both an miRNA and a coding target are sorted by whether 
the miRNA is cellular or from KSHV, and then by coding target region. Fold pattern
outputs are then created for each category, and for combined coding and cellular miRNA.
 
Hybrid sequence information created by the Hyb program and the fold output are
provided with the hybkit package in the databases directory. They were created 
created by the downstream analysis of files 
available at NCBI Gene Expression Ombnibus (GEO) GSE101978, at:
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978

The data files can be downloaded and uncompressed by using the command:
  "sh ./prepare_data.sh"

The unpacked data-files require ~0.2 Gb of space.
The completed output of the analysis requires ~50 Mb of space.

