#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "bash ./analysis_shell.sh"

NOTES="""
Example fold analysis performed using shell scripts.

Provided as an example of usage of hybkit shell executable scripts.
This will produce identical output to analysis_python.py version,
though that implementation is more efficient.

See: 'README.rst' for this analysis for more information.
"""
echo -e """${NOTES}"""

# Stop if error
set -e -u -o pipefail

OUT_DIR="output_shell"

if ! [ -d ${OUT_DIR} ]; then
  mkdir ${OUT_DIR}
fi

IN_HYB_FILES[0]="WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb"
IN_FOLD_FILES[0]="WT_BR1_comp_hOH7_KSHV_hybrids_ua.vienna"
STRING_MATCH_LEGEND_FILE="string_match_legend.csv"

echo "Analyzing Files:"
for fn_i in "0"; do
  echo "    ${IN_HYB_FILES[${fn_i}]}"
  echo "    ${IN_FOLD_FILES[${fn_i}]}"
done

set -v
# Check input hyb files for errors.
hyb_check -i ${IN_HYB_FILES[*]} -f ${IN_FOLD_FILES[*]} \
    --verbose \
    --seq_type dynamic \
    --allowed_mismatches 3 \
    --iter_error_mode 'skip'

# Evaluate segment types and miRNA information and add to record flags
hyb_eval -i ${IN_HYB_FILES[*]} \
         -f ${IN_FOLD_FILES[*]} \
         --out_dir ${OUT_DIR} \
         --eval_types type mirna \
         --type_method string_match \
         --type_params_file "${STRING_MATCH_LEGEND_FILE}" \
         --mirna_types miRNA KSHV-miRNA \
         --hybformat_id True \
         --seq_type dynamic \
         --set_dataset \
         --verbose \


# Generate file names to match python analysis
EVAL_HYB_FILES=( $(ls ${OUT_DIR}/*evaluated*.hyb) )
EVAL_FOLD_FILES=( $(ls ${OUT_DIR}/*evaluated*.vienna) )
QC_FILES=""
for EVAL_FILE in ${EVAL_HYB_FILES[*]}; do
  QC_FILES+="${EVAL_FILE/evaluated/qc} "
done

# Filter records to only those where any reference contains the string "kshv"
# Also Filter record by removing undesirable record types
for fn_i in "0"; do
  hyb_filter -i ${EVAL_HYB_FILES[${fn_i}]} \
             -f ${EVAL_FOLD_FILES[${fn_i}]} \
             --verbose \
             -o ${QC_FILES[${fn_i}]} \
             --out_dir ${OUT_DIR} \
             --seq_type dynamic \
             --filter mirna_not_dimer \
             --exclude any_seg_type_is rRNA \
             --exclude_2 any_seg_type_is mitoch-rRNA \


  # OUT_FOLD="${OUT_DIR}/${FOLD_BASENAME/.vienna/_qc.vienna}"
  # echo ${OUT_FOLD}
done

ls ${OUT_DIR}

# Cleanup intermediate eval files
rm -v ${OUT_DIR}/*evaluated*.hyb
rm -v ${OUT_DIR}/*evaluated*.vienna

ALL_FOLD_QC_FILES=( $(ls ${OUT_DIR}/*qc*.vienna) )

for fn_i in "0"; do
  hyb_analyze -i ${QC_FILES[${fn_i}]} --verbose \
              -f ${ALL_FOLD_QC_FILES[${fn_i}]} \
              --out_dir ${OUT_DIR} \
              --seq_type dynamic \
              -u '' \
              --analysis_types fold energy \
              --analysis_name "WT_BR1" \
              --seq_type dynamic \
              --allowed_mismatches 3
done
# Cleanup intermediate qc files
rm -v ${QC_FILES[*]} ${ALL_FOLD_QC_FILES[*]}

set +v


echo -e "\nDone.\n"

