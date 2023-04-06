#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""


import copy
from contextlib import nullcontext as does_not_raise
import pytest
# mport hybkit

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


# ----- Begin Analysis Tests -----
# ----- Test Analysis Class Standard energy, type, and mirna Analyses -----
test_parameters = [
    ('individual_add', True),
    ('add_hyb_records', False),
]


@pytest.mark.parametrize("test_name,individual_add", [*test_parameters])
def test_analysis_hyb(test_name, individual_add, tmp_path):
    hyb_autotest_file_path = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    test_hyb_strs = [props['hyb_str'] for props in ART_HYB_PROPS_ALL]
    test_hyb_strs_dup = []
    for test_hyb_str in test_hyb_strs:
        test_hyb_strs_dup += ([test_hyb_str] * 15)
    record_count = len(test_hyb_strs_dup)
    combined_hyb_strs = ''.join(test_hyb_strs_dup)
    with open(hyb_autotest_file_path, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(combined_hyb_strs)

    hyb_analysis = hybkit.analysis.Analysis(
        analysis_types=['energy', 'type', 'mirna'],
        name='test_analysis',
    )

    with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file:
        if individual_add:
            for hyb_record in hyb_file:
                hyb_record.eval_types()
                hyb_record.eval_mirna()
                hyb_analysis.add_hyb_record(hyb_record)
        else:
            hyb_analysis.add_hyb_records(hyb_file, eval_types=True, eval_mirna=True)

    all_results = hyb_analysis.get_all_results()
    for key in ['energy', 'type', 'mirna']:
        assert key in all_results

    # Check energy results
    assert all_results['energy']['energy_analysis_count'] == record_count
    assert all_results['energy']['has_energy_val'] == record_count
    assert all_results['energy']['no_energy_val'] == 0
    assert all_results['energy']['energy_min'] == -10.0
    assert all_results['energy']['energy_max'] == -10.0
    assert all_results['energy']['energy_mean'] == -10.0
    assert all_results['energy']['energy_std'] == 0.0
    assert all_results['energy']['binned_energy_vals'][-10] == record_count

    # Check type results
    assert all_results['type']['types_analysis_count'] == record_count
    assert all_results['type']['hybrid_types'] == {
        ('microRNA', 'microRNA'): record_count / 4,
        ('microRNA', 'mRNA'): record_count / 4,
        ('mRNA', 'microRNA'): record_count / 4,
        ('mRNA', 'mRNA'): record_count / 4,
    }
    assert all_results['type']['reordered_hybrid_types'] == {
        ('microRNA', 'microRNA'): record_count / 4,
        ('mRNA', 'microRNA'): record_count / 2,
        ('mRNA', 'mRNA'): record_count / 4,
    }
    assert all_results['type']['mirna_hybrid_types'] == {
        ('microRNA', 'microRNA'): record_count / 4,
        ('microRNA', 'mRNA'): record_count / 2,
        ('mRNA', 'mRNA'): record_count / 4,
    }
    assert all_results['type']['seg1_types'] == {
        'microRNA': record_count / 2,
        'mRNA': record_count / 2,
    }
    assert all_results['type']['seg2_types'] == {
        'microRNA': record_count / 2,
        'mRNA': record_count / 2,
    }
    assert all_results['type']['all_seg_types'] == {
        'microRNA': record_count,
        'mRNA': record_count,
    }

    # Check mirna results
    assert all_results['mirna']['mirna_analysis_count'] == record_count
    assert all_results['mirna']['has_mirna'] == record_count * 3 / 4
    assert all_results['mirna']['non_mirna'] == record_count / 4
    assert all_results['mirna']['mirnas_5p'] == record_count / 4
    assert all_results['mirna']['mirnas_3p'] == record_count / 4
    assert all_results['mirna']['mirna_dimers'] == record_count / 4

    # Check analysis results fetching:
    mirna_results = hyb_analysis.get_analysis_results('mirna')
    assert mirna_results == all_results['mirna']

    # Check specific results fetching:
    mirnas_5p_results = hyb_analysis.get_specific_result('mirnas_5p')
    assert mirnas_5p_results == all_results['mirna']['mirnas_5p']

    # Get Analysis Delim String:
    out_analysis_file_name = os.path.join(tmp_path, 'analysis_delim_str.csv')
    mirna_analysis_delim_str = hyb_analysis.get_analysis_delim_str(analysis='mirna')
    analysis_delim_str = hyb_analysis.get_analysis_delim_str()
    hyb_analysis.write_analysis_delim_str(analysis='mirna', out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_delim_str(analysis=['mirna'], out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_delim_str(out_file_name=out_analysis_file_name)

    out_special_file_base = os.path.join(tmp_path, 'special_analysis')
    hyb_analysis.write_analysis_delim_str(analysis=['mirna'], out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base,
                                                analysis=['mirna'])
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base,
                                                analysis='mirna')
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base)

    # Check erroring on bad results requests:
    with pytest.raises(RuntimeError):
        hyb_analysis.get_analysis_results('bad_analysis')
    with pytest.raises(RuntimeError):
        hyb_analysis.get_specific_result('bad_result')


# ----- Test Analysis Class Standard energy, type, and mirna Analyses -----
def test_analysis_problems(tmp_path):
    hyb_autotest_file_path = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    test_hyb_strs = [props['hyb_str'] for props in ART_HYB_PROPS_ALL]
    test_hyb_strs_dup = []
    for test_hyb_str in test_hyb_strs:
        test_hyb_strs_dup += ([test_hyb_str] * 15)
    record_count = len(test_hyb_strs_dup)
    combined_hyb_strs = ''.join(test_hyb_strs_dup)
    with open(hyb_autotest_file_path, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(combined_hyb_strs)

    hyb_analysis = hybkit.analysis.Analysis(
        analysis_types=['energy', 'type', 'mirna'],
        name='test_analysis',
    )

    # Test raise error if no energy values
    with pytest.raises(RuntimeError):
        with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file:
            for hyb_record in hyb_file:
                hyb_record.energy = None
                hyb_analysis.add_hyb_record(hyb_record)

    # Test raise error if no type values
    with pytest.raises(RuntimeError):
        with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file:
            for hyb_record in hyb_file:
                hyb_analysis.add_hyb_record(hyb_record)

    # Test raise error if no mirna values
    with pytest.raises(RuntimeError):
        with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file:
            for hyb_record in hyb_file:
                hyb_record.eval_types()
                hyb_analysis.add_hyb_record(hyb_record)

    # Test raise error if no fold_record values
    hyb_analysis = hybkit.analysis.Analysis(
        analysis_types=['mirna', 'fold'],
        name='test_analysis',
    )
    with pytest.raises(RuntimeError):
        with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file:
            for hyb_record in hyb_file:
                hyb_record.eval_types()
                hyb_record.eval_mirna()
                hyb_analysis.add_hyb_record(hyb_record)


# ----- Test Analysis Class Standard fold Analysis -----
def test_analysis_fold(tmp_path):
    hyb_autotest_file_path = os.path.join(tmp_path, 'hyb_autotest_file.hyb')
    vienna_autotest_file_path = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    all_hyb_vienna_props = [ART_HYB_VIENNA_PROPS_1, ART_HYB_VIENNA_PROPS_2]
    test_hyb_strs = [props['hyb_str'] for props in all_hyb_vienna_props]
    test_vienna_strs = [props['vienna_str'] for props in all_hyb_vienna_props]
    test_hyb_strs_dup = []
    test_veinna_strs_dup = []
    for test_hyb_str in test_hyb_strs:
        test_hyb_strs_dup += ([test_hyb_str] * 15)
    for test_vienna_str in test_vienna_strs:
        test_veinna_strs_dup += ([test_vienna_str] * 15)
    record_count = len(test_hyb_strs_dup)
    combined_hyb_strs = ''.join(test_hyb_strs_dup)
    combined_vienna_strs = ''.join(test_veinna_strs_dup)
    with open(hyb_autotest_file_path, 'w') as hyb_autotest_file:
        hyb_autotest_file.write(combined_hyb_strs)
    with open(vienna_autotest_file_path, 'w') as vienna_autotest_file:
        vienna_autotest_file.write(combined_vienna_strs)

    hyb_analysis = hybkit.analysis.Analysis(
        analysis_types=['fold'],
        name='test_analysis',
    )
    with open(hyb_autotest_file_path, 'r') as hyb_autotest_file, \
         open(vienna_autotest_file_path, 'r') as vienna_autotest_file:
        assert len(hyb_autotest_file.readlines()) == record_count
        assert len(vienna_autotest_file.readlines()) == record_count * 3

    with hybkit.HybFile(hyb_autotest_file_path, 'r') as hyb_file, \
         hybkit.ViennaFile(
             vienna_autotest_file_path, 'r', seq_type='dynamic') as vienna_file:
        hyb_fold_iter = hybkit.HybFoldIter(
            hyb_file,
            vienna_file,
            combine=True,
            iter_error_mode='warn_skip',
        )
        for hyb_record in hyb_fold_iter:
            hyb_record.eval_types()
            hyb_record.eval_mirna()
            hyb_analysis.add_hyb_record(hyb_record)

        # assert not hyb_fold_iter.report()

    all_results = hyb_analysis.get_all_results()
    for key in ['fold']:
        assert key in all_results

    fold_results = all_results['fold']
    assert fold_results['fold_analysis_count'] == record_count
    assert fold_results['folds_recorded'] == record_count
    assert fold_results['fold_match_counts'] == {
        14: (record_count / 2),
        15: (record_count / 2),
    }
    assert fold_results['mirna_nt_fold_counts'][4] == record_count
    assert fold_results['mirna_nt_fold_counts'][24] == (record_count / 2)
    assert fold_results['mirna_nt_fold_props'][4] == 1.0
    assert fold_results['mirna_nt_fold_props'][24] == 0.5

    # Get Analysis Delim String:
    out_analysis_file_name = os.path.join(tmp_path, 'analysis_delim_str.csv')
    fold_analysis_delim_str = hyb_analysis.get_analysis_delim_str(analysis='fold')
    analysis_delim_str = hyb_analysis.get_analysis_delim_str()
    hyb_analysis.write_analysis_delim_str(analysis='fold', out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_delim_str(analysis=['fold'], out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_delim_str(out_file_name=out_analysis_file_name)

    out_special_file_base = os.path.join(tmp_path, 'special_analysis')
    hyb_analysis.write_analysis_delim_str(analysis=['fold'], out_file_name=out_analysis_file_name)
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base,
                                                analysis=['fold'])
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base,
                                                analysis='fold')
    hyb_analysis.write_analysis_results_special(out_file_basename=out_special_file_base)

# ----- Begin old tests -----
def old_test_analyses():
    out_basename = None
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
    out_basename = None
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
