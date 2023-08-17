#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "bash ./analysis_shell.sh"

NOTES="""
Analysis for type/mirna analysis performed using shell scripts.

Provided as an example of usage of hybkit shell executable scripts.
This will produce identical output to the analysis_python.py version,
though that implementation is more efficient.

See: 'README.rst' for this analysis for more information.
"""

printf """${NOTES}"""

# Stop if error
set -e -u -o pipefail

OUT_DIR="output_shell"

if ! [ -d ${OUT_DIR} ]; then
  mkdir ${OUT_DIR}
fi

#IN_FILES[0]="GSM2720017_UI_BR1.hyb"
#IN_FILES[1]="GSM2720018_UI_BR2.hyb"
#IN_FILES[2]="GSM2720019_UI_BR3.hyb"
IN_FILES[0]="GSM2720020_WT_BR1.hyb"
IN_FILES[1]="GSM2720021_WT_BR2.hyb"
IN_FILES[2]="GSM2720022_WT_BR3.hyb"
#IN_FILES[6]="GSM2720023_D11_BR1.hyb"
#IN_FILES[7]="GSM2720024_D11_BR2.hyb"
#IN_FILES[8]="GSM2720025_D11_BR3.hyb"
COMBINED_FILE="${OUT_DIR}/combined_analysis.hyb"
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
         --type_params_file ${STRING_MATCH_LEGEND_FILE} \
         --mirna_types miRNA KSHV-miRNA \
         --hybformat_id True \
         --set_dataset

# Generate file names to match python analysis
EVAL_FILES=( $(ls ${OUT_DIR}/*evaluated*.hyb) )
QC_FILES=""
for EVAL_FILE in ${EVAL_FILES[*]}; do
  QC_FILES+="${EVAL_FILE/evaluated/qc} "
done

# Filter records by removing undesirable record types
hyb_filter -i ${EVAL_FILES[*]} --verbose \
           -o ${QC_FILES} \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \
           --skip_dup_id_after \

# Cleanup intermediate files
rm -v ${OUT_DIR}/*evaluated*.hyb

# Run Summary analysis on individual files
for fn in ${QC_FILES[*]}; do
  hyb_analyze -i ${fn} --verbose \
              --out_dir ${OUT_DIR} \
              -u '' \
              --analysis_types "type" "mirna" \
              --analysis_name "$(basename ${fn/_qc.hyb/})"
done
set +v


# Create combined hyb file
echo -e "\nCreating Combined Analysis .hyb file:"
echo    "    ${COMBINED_FILE}"

set -v
cat ${QC_FILES} > ${COMBINED_FILE}

# Run Summary analysis on all records
hyb_analyze -i ${COMBINED_FILE} --verbose \
              -o ${COMBINED_FILE/.hyb/} \
              -u '' \
              --analysis_types "type" "mirna" \
              --analysis_name "Combined Analysis"

rm -v ${COMBINED_FILE}
set +v
echo -e "\nDone.\n"

