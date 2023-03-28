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
def old_test_viennafile():
    fold_record = hybkit.FoldRecord.from_vienna_string(VIENNA_STR_1)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_direct('testval')
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        with pytest.raises(RuntimeError):
            vienna_file.write_record('NotRecord')
        vienna_file.write_record(fold_record)
        vienna_file.write_records([fold_record, fold_record])
    assert hybkit.util.vienna_exists(vienna_file_name)
    assert hybkit.util.fold_exists(vienna_file_name)
    with hybkit.FoldFile.open(vienna_file_name, 'r') as fold_file:
        with pytest.raises(RuntimeError):
            fold_file.read_record()
        with pytest.raises(RuntimeError):
            fold_file._to_record_string(fold_record, False)
    with hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        fold_record_read = vienna_file.read_records()[0]

    assert fold_record == fold_record_read
    for test_prop in ['id', 'seq', 'fold', 'energy']:
        print('Comparing:', test_prop,
              getattr(fold_record, test_prop), getattr(fold_record_read, test_prop))
        assert getattr(fold_record, test_prop) == getattr(fold_record_read, test_prop)


def old_test_hybfolditer():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        hyb_file.write_records([hyb_record] * 30)
    fold_record = hybkit.FoldRecord.from_vienna_string(VIENNA_STR_1)
    dyn_fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)
    dyn_fold_record_mis = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1_MIS)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_records([fold_record] * 30)

    # Compare FoldRecord and HybRecord
    assert not fold_record.matches_hyb_record(hyb_record)
    fold_record.count_hyb_record_mismatches(hyb_record)
    with pytest.raises(RuntimeError):
        fold_record.ensure_matches_hyb_record(hyb_record)

    # Compare DynamicFoldRecord and HybRecord
    assert dyn_fold_record.matches_hyb_record(hyb_record)
    dyn_fold_record.count_hyb_record_mismatches(hyb_record)
    dyn_fold_record.ensure_matches_hyb_record(hyb_record)

    # Compare DynamicFoldRecord With Mismatch and HybRecord
    assert not dyn_fold_record_mis.matches_hyb_record(hyb_record)
    dyn_fold_record_mis.count_hyb_record_mismatches(hyb_record)
    with pytest.raises(RuntimeError):
        dyn_fold_record_mis.ensure_matches_hyb_record(hyb_record)

    hybkit.settings.FoldFile_settings['foldrecord_type'] = 'strict'
    hybkit.settings.HybFoldIter_settings['error_mode'] = 'raise'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        with pytest.raises(RuntimeError):
            for ret_item in hf_iter:
                print(ret_item)

    hybkit.settings.HybFoldIter_settings['error_mode'] = 'warn_skip'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        with pytest.raises(RuntimeError):
            for ret_item in hf_iter:
                print(ret_item)

    hybkit.settings.HybFoldIter_settings['error_mode'] = 'warn_return'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        for ret_item in hf_iter:
            print(ret_item)

    hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
    for combine in [True, False]:
        with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
             hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
            hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file, combine=combine)
            for ret_item in hf_iter:
                print(ret_item)
            hf_iter.report()
            hf_iter.print_report()

    hyb_record, fold_record = ret_item


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


def old_test_util():
    original_abspath = hybkit.settings._USE_ABSPATH
    hybkit.settings._USE_ABSPATH = True
    assert hybkit.util._bool_from_string(True)
    assert hybkit.util._bool_from_string('yes')
    assert not hybkit.util._bool_from_string('no')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util._bool_from_string('invalid')
    assert hybkit.util.dir_exists('~')
    assert hybkit.util.dir_exists('${PWD}')
    assert hybkit.util.file_exists(__file__)
    assert hybkit.util.out_path_exists(hyb_file_name)

    hybkit.util.make_out_file_name(hyb_file_name, name_suffix='out',
                                   in_suffix='.hyb', out_suffix='.new',
                                   out_dir='', seg_sep='_'
                                   )
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.dir_exists('nonexistent_dir')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.file_exists('nonexistent_file')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.hyb_exists(match_legend_file_name)
    use_namespace = argparse.Namespace()
    setattr(use_namespace, 'hybformat_id', True)
    hybkit.util.set_settings(use_namespace, verbose=True)
    hybkit.settings._USE_ABSPATH = original_abspath

    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.out_hybs_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', hyb_file_name, hyb_file_name, '-o', hyb_file_name]
    )
    with pytest.raises(SystemExit):
        hybkit.util.validate_args(args, script_parser)
    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.in_folds_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', hyb_file_name, hyb_file_name, '-f', vienna_file_name]
    )
    with pytest.raises(SystemExit):
        hybkit.util.validate_args(args, script_parser)


def old_test_type_finder():
    # Generic Tests
    def do_nothing(*args, **kwargs):
        pass
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    hybkit.type_finder.TypeFinder.set_custom_method(do_nothing)

    # Defualt Hybformat
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hybkit.type_finder.TypeFinder.set_method('hybformat')
    hyb_record.eval_types()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    # Non-Defualt String-Match
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params('badfile')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad1_match_legend_file_name)
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad2_match_legend_file_name)
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(
        bad3_match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.settings['check_complete_seg_types'] = True
    with pytest.raises(RuntimeError):
        hyb_record.eval_types()
    # match_params = {'startswith':[('MIMAT', 'MIMAT')], 'endswith':[('microRNA', 'miRNA')]}
    # hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    # with pytest.raises(RuntimeError):
    #     hyb_record.eval_types()
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    # Non-Defualt ID-Map
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params('wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params(type_file_pairs='wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params([bad1_id_map_legend_file_name])
    with pytest.raises(RuntimeError):
        id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
            type_file_pairs=[('seq1type', id_map_legend_file_name),
                             ('seq2type', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
        type_file_pairs=[('seqtype', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params([id_map_legend_file_name])
    hybkit.type_finder.TypeFinder.set_method('id_map', id_map_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.seg1_props['ref_name'] = 'not_real_name'
    with pytest.raises(RuntimeError):
        hyb_record.eval_types(allow_unknown=False)

# if __name__ == '__main__':
#     test_type_finder()
