#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "bash ./analysis_shell.sh"

NOTES="""
Grouped target analysis pipeline performed using shell scripts.

Provided as an example of usage of hybkit shell executable scripts.
File names are hardcoded, and functions are accessed directly.
This will produce identical output to analysis_python.py version,
though that implementation is more efficeint.

See: 'summary_analysis_notes.rst' for more information.
"""
echo -e """${NOTES}"""

# Stop if error 
set -e -u -o pipefail

OUT_DIR="output_shell"

if ! [ -d ${OUT_DIR} ]; then
  mkdir ${OUT_DIR}
fi

IN_FILES[0]="GSM2720020_WT_BR1.hyb"
IN_FILES[1]="GSM2720021_WT_BR2.hyb"
IN_FILES[2]="GSM2720022_WT_BR3.hyb"
IN_FILES[3]="GSM2720023_D11_BR1.hyb"
IN_FILES[4]="GSM2720024_D11_BR2.hyb"
IN_FILES[5]="GSM2720025_D11_BR3.hyb"
ALL_QC_FILES=""
COMBINED_FILE="${OUT_DIR}/kshv_combined.hyb"
STRING_MATCH_LEGEND_FILE="string_match_legend.csv"

echo "Analyzing Files:"
for fn in ${IN_FILES[*]}; do
  echo "    ${fn}"
done

set -v 
# Check input files for errors.
hyb_check -i ${IN_FILES[*]} --verbose

# Evaluate segment types and miRNA information and add to record flags
hyb_eval -i ${IN_FILES[*]} --verbose \
         --out_dir ${OUT_DIR} \
         --eval_types type mirna \
         --type_method string_match \
         --type_params_file string_match_legend.csv \
         --mirna_types miRNA KSHV-miRNA \
         --hybformat_id True \
         --set_dataset  

# Generate file names to match python analysis
EVAL_FILES=( $(ls ${OUT_DIR}/*evaluated*.hyb) )
QC_FILES=""
for EVAL_FILE in ${EVAL_FILES[*]}; do
  QC_FILES+="${EVAL_FILE/evaluated/KSHV_only} "
done

# Filter records to only those where any reference contains the string "kshv"
# Also Filter record by removing undesirable record types
hyb_filter -i ${EVAL_FILES[*]} --verbose \
           -o ${QC_FILES} \
           --filter has_mirna \
           --filter_2 mirna_contains kshv \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \

# Cleanup intermediate eval files
rm -v ${OUT_DIR}/*evaluated*.hyb 

# Create combined hyb file
echo -e "\nCreating Combined Analysis .hyb file:"
echo    "    ${COMBINED_FILE}"

cat ${QC_FILES} > ${COMBINED_FILE}

# Cleanup intermediate qc files
rm -v ${QC_FILES}

hyb_analyze -i ${COMBINED_FILE} --verbose \
            --out_dir ${OUT_DIR} \
            -u '' \
            --analysis_type target \
            --analysis_name "KSHV Combined" \
            --allow_mirna_dimers \
            --write_individual
set +v


echo -e "\nDone.\n"
