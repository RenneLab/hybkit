#!/usr/bin/env sh
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "sh ./download_data.sh"

echo -e "\nDownloading Data for sample grouped target analysis."

# Download tar-compressed hyb-format genomic sequence datafile GSE101978_RAW.tar from NCBI Gene Expression Ombnibus (GEO) GSE101978, at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978
echo -e "Downloading tar-compressed data-file..."
#wget -O GSE101978_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101978&format=file"
curl -o GSE101978_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101978&format=file"

# Unpack tar-compressed file:
echo -e "Unpacking tar-compressed data-file...\n"
tar -xf GSE101978_RAW.tar
rm GSE101978_RAW.tar

# Remove unneded (non-KSHV-relevant) hyb files.
rm GSM2720017_UI_BR1.hyb.txt.gz \
   GSM2720018_UI_BR2.hyb.txt.gz \
   GSM2720019_UI_BR3.hyb.txt.gz

# Unzip (gunzip) individual hyb-format data-files.
echo -e "Unpacking gz-compressed hyb data files...\n"
gunzip GSM2720020_WT_BR1.hyb.txt.gz \
       GSM2720021_WT_BR2.hyb.txt.gz \
       GSM2720022_WT_BR3.hyb.txt.gz \
       GSM2720023_D11_BR1.hyb.txt.gz \
       GSM2720024_D11_BR2.hyb.txt.gz \
       GSM2720025_D11_BR3.hyb.txt.gz

echo -e "Removing '.txt' suffix from hyb files.\n"
mv GSM2720020_WT_BR1.hyb.txt GSM2720020_WT_BR1.hyb
mv GSM2720021_WT_BR2.hyb.txt GSM2720021_WT_BR1.hyb
mv GSM2720022_WT_BR3.hyb.txt GSM2720022_WT_BR1.hyb
mv GSM2720023_D11_BR1.hyb.txt GSM2720023_D11_BR1.hyb
mv GSM2720024_D11_BR2.hyb.txt GSM2720024_D11_BR1.hyb
mv GSM2720025_D11_BR3.hyb.txt GSM2720025_D11_BR1.hyb

echo -e "Done.\n"
