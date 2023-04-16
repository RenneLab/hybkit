#!/usr/bin/env sh
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "sh ./download_data.sh"

if [ -f WT_BR1_comp_hOH7_KSHV_hybrids_ua.vienna ]; then
    printf "Data already downloaded.\n"
else
    printf "\nPreparing data for example analysis."

      printf "\nCopying data to analysis directory."
      cp -v ../ref_data/WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb.gz \
            ../ref_data/WT_BR1_comp_hOH7_KSHV_hybrids_ua.vienna.gz \
            ./


      # Unzip (gunzip) individual hyb-format data-files.
      printf "\nUnpacking gz-compressed hyb data files...\n"
      gunzip WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb.gz \
            WT_BR1_comp_hOH7_KSHV_hybrids_ua.vienna.gz
fi

printf "Done.\n"
