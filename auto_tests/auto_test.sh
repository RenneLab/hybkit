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

if ! [ -d ${OUT_DIR} ]; then
  mkdir ${OUT_DIR}
fi

IN_HYB="test_hybrid.hyb"
FULL_IN_HYB="test_data_files/${IN_HYB}"
IN_VIENNA="test_hybrid.vienna"
FULL_IN_VIENNA="test_data_files/${IN_VIENNA}"

# Run tests
hyb_check -i ${FULL_IN_HYB} --verbose

hyb_eval -i ${FULL_IN_HYB} --verbose \
         --out_dir ${OUT_DIR} \
         --eval_types type mirna \
         --hybformat_id True \
         --set_dataset

ls ${OUT_DIR}

hyb_filter -i ${OUT_DIR}/${IN_HYB/.hyb/_evaluated.hyb} --verbose \
           --out_dir ${OUT_DIR} \
           --exclusion_table \
           --exclude any_seg_type_is rRNA \
           --exclude_2 any_seg_type_is mitoch-rRNA \

ls ${OUT_DIR}


hyb_exclude_fold -f ${FULL_IN_VIENNA} --verbose \
                 -e ${OUT_DIR}/${IN_HYB/.hyb/_evaluated_exclude.csv} \
                 --out_dir ${OUT_DIR} \
                 --error_mode warn_return \

ls ${OUT_DIR}

for mode in $(echo "type mirna summary target"); do
  hyb_analyze -i ${OUT_DIR}/${IN_HYB/.hyb/_evaluated_filtered.hyb} --verbose \
              --out_dir ${OUT_DIR} \
              --analysis_type ${mode} \
              --analysis_name "TEST"
done
set +v

ls ${OUT_DIR}


hyb_fold_analyze -i ${OUT_DIR}/${IN_HYB/.hyb/_evaluated_filtered.hyb} --verbose \
                 -f ${OUT_DIR}/${IN_VIENNA/.vienna/_filtered.vienna} \
                 --out_dir ${OUT_DIR} \
                 --analysis_type fold \
                 --analysis_name "TEST_FOLD" \
                 --foldrecord_type dynamic \
                 --allowed_mismatches 0

ls ${OUT_DIR}

echo -e "\nDone with Autotests\n"

