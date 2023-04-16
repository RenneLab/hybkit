#!/usr/bin/env bash
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

# usage: "./auto_test.sh"

NOTES="""
Automated testing for hybkit scripts.
"""
echo -e """${NOTES}"""

# Stop if error
set -e -u -o pipefail

OUT_DIR="output_autotest"

if ! [ -d "${OUT_DIR}" ]; then
  mkdir "${OUT_DIR}"
fi

IN_HYB="test_hybrid.hyb"
FULL_IN_HYB="test_data_files/${IN_HYB}"
IN_VIENNA="test_hybrid.vienna"
FULL_IN_VIENNA="test_data_files/${IN_VIENNA}"

# Run tests
hyb_check -i "${FULL_IN_HYB}" --verbose
hyb_check -i "${FULL_IN_HYB}" -f "${FULL_IN_VIENNA}" --verbose

hyb_eval -i "${FULL_IN_HYB}" --verbose \
         --out_dir "${OUT_DIR} "\
         --eval_types type mirna \
         --hybformat_id True \
         --set_dataset

hyb_eval -i "${FULL_IN_HYB}" \
         -f "${FULL_IN_VIENNA}" \
         --verbose \
         --out_dir "${OUT_DIR} "\
         --eval_types type mirna \
         --hybformat_id True \
         --seq_type dynamic \
         --set_dataset

hyb_filter -i "${OUT_DIR}/${IN_HYB/.hyb/_evaluated.hyb}" --verbose \
           --out_dir "${OUT_DIR}" \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \

hyb_filter -i "${OUT_DIR}/${IN_HYB/.hyb/_evaluated.hyb}" \
           -f "${OUT_DIR}/${IN_HYB/.hyb/_evaluated.vienna}" \
           --verbose \
           --seq_type dynamic \
           --out_dir "${OUT_DIR}" \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \

for mode in "energy" "type" "mirna" "target" "fold" "energy type mirna target fold"; do

hyb_analyze -i "${OUT_DIR}/${IN_HYB/.hyb/_evaluated_filtered.hyb}" --verbose \
            -f "${OUT_DIR}/${IN_HYB/.hyb/_evaluated_filtered.vienna}" \
            --out_dir "${OUT_DIR}" \
            --analysis_types ${mode} \
            --analysis_name "TEST_FOLD" \
            --seq_type dynamic \
            --allowed_mismatches 0

done

set +v
echo -e "\nDone with Autotests\n"

