#!/usr/bin/env sh
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

echo -e "\nDownloading Data for sample data analysis."

# Download tar-compressed hyb-format genomic sequence datafile GSE101978_RAW.tar from NCBI Gene Expression Ombnibus (GEO) GSE101978, at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101978
echo -e "Downloading tar-compressed data-file..."
wget -O GSE101978_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101978&format=file"

# Unpack tar-compressed file:
echo -e "Unpacking tar-compressed data-file...\n"
tar -xf GSE101978_RAW.tar
rm GSE101978_RAW.tar 

# Unzip (gunzip) individual hyb-format data-files.
echo -e "Unpacking gz-compressed hyb data files...\n"
gunzip GSM2720017_UI_BR1.hyb.txt.gz \
       GSM2720018_UI_BR2.hyb.txt.gz \
       GSM2720019_UI_BR3.hyb.txt.gz \
       GSM2720020_WT_BR1.hyb.txt.gz \
       GSM2720021_WT_BR2.hyb.txt.gz \
       GSM2720022_WT_BR3.hyb.txt.gz \
       GSM2720023_D11_BR1.hyb.txt.gz \
       GSM2720024_D11_BR2.hyb.txt.gz \
       GSM2720025_D11_BR3.hyb.txt.gz

echo -e "Removing '.txt' suffix from hyb files.\n"
mv GSM2720017_UI_BR1.hyb.txt GSM2720017_UI_BR1.hyb
mv GSM2720018_UI_BR2.hyb.txt GSM2720018_UI_BR1.hyb
mv GSM2720019_UI_BR3.hyb.txt GSM2720019_UI_BR1.hyb
mv GSM2720020_WT_BR1.hyb.txt GSM2720020_WT_BR1.hyb
mv GSM2720021_WT_BR2.hyb.txt GSM2720021_WT_BR1.hyb
mv GSM2720022_WT_BR3.hyb.txt GSM2720022_WT_BR1.hyb
mv GSM2720023_D11_BR1.hyb.txt GSM2720023_D11_BR1.hyb
mv GSM2720024_D11_BR2.hyb.txt GSM2720024_D11_BR1.hyb
mv GSM2720025_D11_BR3.hyb.txt GSM2720025_D11_BR1.hyb

echo -e "Done.\n"
