#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "bash ./analysis_shell.sh"

NOTES="""
Example target analysis performed using shell scripts.

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

IN_FILES[0]="GSM2720020_WT_BR1.hyb"
STRING_MATCH_LEGEND_FILE="string_match_legend.csv"

printf "Analyzing Files:"
for fn in ${IN_FILES[*]}; do
  printf "    %s" "${fn}"
done

set -v

# Check input files for errors.
hyb_check -i "${IN_FILES[*]}" --verbose

# Evaluate segment types and miRNA information and add to record flags
hyb_eval -i "${IN_FILES[*]}" --verbose \
         --out_dir "${OUT_DIR}" \
         --eval_types type mirna \
         --type_method string_match \
         --type_params_file "${STRING_MATCH_LEGEND_FILE}" \
         --mirna_types miRNA KSHV-miRNA \
         --hybformat_id True \
         --set_dataset 

# Generate file names to match python analysis
EVAL_FILES=( $(ls ${OUT_DIR}/*evaluated*.hyb) )
OUT_FILES=""
for EVAL_FILE in ${EVAL_FILES[*]}; do
  OUT_FILES+="${EVAL_FILE/evaluated/kshv-miR-K12-5_only} "
done

# Filter records to only those where any reference contains the string "kshv"
# Also Filter record by removing undesirable record types
hyb_filter -i "${EVAL_FILES[*]}" --verbose \
           -o ${OUT_FILES} \
           --filter has_mirna \
           --filter_2 mirna_contains kshv-miR-K12-5 \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \
           --skip_dup_id_after

# Cleanup intermediate files
rm -v ${OUT_DIR}/*evaluated*.hyb

# Run target analysis on individual files
for fn in ${OUT_FILES[*]} ; do
  hyb_analyze -i ${fn} --verbose \
              --out_dir ${OUT_DIR} \
              --out_suffix '' \
              --analysis_types target \
              --analysis_name "$(basename ${fn/_kshv-miR-K12-5_only.hyb/})" \

done
set +v

printf "\nDone.\n"
