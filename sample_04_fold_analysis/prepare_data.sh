#!/usr/bin/env sh
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

# usage: "sh ./prepare_data.sh"

echo -e "\nPreparing data for sample fold analysis."

echo -e "\nCopying data to analysis directory."
cp -v ../databases/WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb.gz \
      ../databases/WT_BR1_comp_hOH7_KSHV_hybrids_ua.viennad.gz \
      ./
  

# Unzip (gunzip) individual hyb-format data-files.
echo -e "\nUnpacking gz-compressed hyb data files...\n"
gunzip WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb.gz \
       WT_BR1_comp_hOH7_KSHV_hybrids_ua.viennad.gz

echo -e "Done.\n"
