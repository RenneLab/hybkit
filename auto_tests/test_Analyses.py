#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""

import os
import sys
import copy
from contextlib import nullcontext as does_not_raise
import argparse
import pytest
import hybkit


# ----- Import Testing Helper Data -----
from auto_tests.test_helper_data import *
# Includes the following variables:
# TEST_HYBID_STR, TEST_SEQ_STR, TEST_FOLD_STR, TEST_ENERGY_STR
# ART_HYB_PROPS_1, ART_HYB_PROPS_ALL, ART_BAD_HYB_STRS
# ID_ALLOWED_TYPES, SEQ_ALLOWED_TYPES, FOLD_ALLOWED_TYPES, ENERGY_ALLOWED_TYPES
# test_out_dir, hyb_autotest_file_name, hyb_file_name


# ----- Import Testing Helper Functions -----
from auto_tests.test_helper_functions import *
# Includes the following functions:
# get_expected_result_string(is_allowed=False)
# get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError))

# ----- Begin old tests -----
def old_test_analyses():
    analysis_classes = {
        'type': hybkit.analysis.TypeAnalysis,
        'mirna': hybkit.analysis.MirnaAnalysis,
        'summary': hybkit.analysis.SummaryAnalysis,
        'target': hybkit.analysis.TargetAnalysis,
    }
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    with pytest.raises(RuntimeError):
        hybkit.analysis.BaseAnalysis()
    last_analysis = None
    for analysis_type in analysis_classes:
        Analysis = analysis_classes[analysis_type]
        analysis_basename = out_basename + '_' + analysis_type
        analysis = Analysis(name='Test')
        analysis2 = Analysis()
        analysis2.plot(analysis_basename)
        analysis2.add(hyb_record)
        analysis.update(analysis2)
        if last_analysis is not None:
            with pytest.raises(RuntimeError):
                analysis.update(last_analysis)
        print(analysis.results())
        print(analysis.results(None, True))
        analysis.write(analysis_basename)
        analysis.plot(analysis_basename)
        if hasattr(analysis, 'write_individual'):
            analysis.write_individual(analysis_basename)
            analysis.plot_individual(analysis_basename)
        last_analysis = analysis

    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    type_analysis = hybkit.analysis.TypeAnalysis()
    mirna_analysis = hybkit.analysis.MirnaAnalysis()
    target_analysis = hybkit.analysis.TargetAnalysis()
    type_analysis.add(hyb_record)
    hyb_record.set_flag('miRNA_seg', '3p')
    type_analysis.add(hyb_record)
    mirna_analysis.add(hyb_record)
    hyb_record.set_flag('miRNA_seg', 'B')
    type_analysis.add(hyb_record)
    mirna_analysis.add(hyb_record)
    init_val = target_analysis.settings['allow_mirna_dimers']
    target_analysis.settings['allow_mirna_dimers'] = True
    target_analysis.add(hyb_record)
    target_analysis.settings['allow_mirna_dimers'] = init_val
    hyb_record.set_flag('miRNA_seg', 'N')
    mirna_analysis.add(hyb_record)


def old_test_fold_analysis():
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)
    analysis = hybkit.analysis.FoldAnalysis(name='Test')
    analysis2 = hybkit.analysis.FoldAnalysis()
    hyb_record.set_fold_record(fold_record)
    analysis_basename = out_basename + '_fold'
    analysis2.plot(analysis_basename)
    analysis2.add(hyb_record)
    analysis.update(analysis2)
    print(analysis.results())
    print(analysis.results(None, True))
    analysis.write(analysis_basename)
    analysis.plot(analysis_basename)

    hyb_record.fold_record.fold = ')))' + hyb_record.fold_record.fold[3:]
    hyb_record.set_flag('miRNA_seg', 'B')
    init_val = analysis.settings['allow_mirna_dimers']
    analysis.settings['allow_mirna_dimers'] = True
    analysis.add(hyb_record)
    analysis.settings['allow_mirna_dimers'] = init_val



# if __name__ == '__main__':
#     test_type_finder()
