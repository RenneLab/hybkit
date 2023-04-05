#!/usr/bin/env sh
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

# usage: "sh ./download_data.sh"

if [ -f GSM2720020_WT_BR1.hyb ]; then
    printf "Data already downloaded.\n"
else
    printf "\nDownloading Data for example data analysis."

    # Download tar-compressed hyb-format genomic sequence datafile GSM2720020_WT_BR1.hyb.txt.gz from NCBI Gene Expression Ombnibus (GEO) GSE101978, at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978 , with the specific file available at ascension GSM2720020 at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2720020
    printf "Downloading gz-compressed data-file..."
    curl -o GSM2720020_WT_BR1.hyb.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2720020&format=file&file=GSM2720020%5FWT%5FBR1%2Ehyb%2Etxt%2Egz"

    # Unzip (gunzip) hyb-format data-file.
    printf "\nUnpacking gz-compressed hyb data file...\n"
    gunzip GSM2720020_WT_BR1.hyb.txt.gz

    printf "Removing '.txt' suffix from hyb files.\n"
    mv GSM2720020_WT_BR1.hyb.txt GSM2720020_WT_BR1.hyb
fi
printf "Done.\n"
