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
import argparse
import pytest
import hybkit

hyb_str_1 = ('695_804	ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC	.	'
             'MIMAT0000078_MirBase_miR-23a_microRNA	1	21	1	21	0.0027	'
             'ENSG00000188229_ENST00000340384_TUBB2C_mRNA	23	49	1181	1207	1.2e-06	dataset=test;'
             )
vie_str_1 = ('>695_804_MIMAT0000078_MirBase_miR-23a_microRNA_1_21-'
             '695_804_ENSG00000188229_ENST00000340384_TUBB2C_mRNA_1181_1207\n'
             'ATCACATTGCCAGGGATTTCCATCCCCAACAATGTGAAAACGGCTGTC\n'
             '.((((((((...(((((....)))))...))))))))...........	(-15)\n'
             )
vie_str_1_mis = ('>695_804_MIMAT0000078_MirBase_miR-23a_microRNA_1_21-'
                 '695_804_ENSG00000188229_ENST00000340384_TUBB2C_mRNA_1181_1207\n'
                 'ATCACATAGCCAGGGATTTCCATCCCCAACAATGTGAAAACGGCTGTC\n'
                 '.((((((((...(((((....)))))...))))))))...........	(-15)\n'
                 )


# ----- Set Testing Variables -----
auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
test_out_dir = os.path.join(auto_tests_dir, 'output_autotest')
hyb_file_name = os.path.join(test_out_dir, 'test_hybrid_py.hyb')
vienna_file_name = os.path.join(test_out_dir, 'test_hybrid_py.vienna')
out_basename = hyb_file_name.replace('.hyb', '')
match_legend_file_name = os.path.join(auto_tests_dir, 'test_string_match.csv')
bad1_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_1.csv')
bad2_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_2.csv')
bad3_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_3.csv')
id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_id_map.csv')
bad1_id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_id_map_1.csv')

# hybkit.settings.HybFile_settings['hybformat_id'] = True
# hybkit.settings.HybFile_settings['hybformat_ref'] = True
# hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
# hybkit.settings.Analysis_settings['allow_mirna_dimers'] = True

# Function for generating HybRecord objects for testing.
def def_hyb_records():
    """Generate two HybRecord objects for tests"""
    hyb_record_1 = hybkit.HybRecord.from_line(hyb_str_1,
                                              hybformat_id=True,
                                              hybformat_ref=True)
    fold_record = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1)
    hyb_record_2 = hybkit.HybRecord(id=hyb_record_1.id, seq=hyb_record_1.seq,
                                    seg1_props=copy.deepcopy(hyb_record_1.seg1_props),
                                    seg2_props=copy.deepcopy(hyb_record_1.seg2_props),
                                    fold_record=fold_record,
                                    )
    return hyb_record_1, hyb_record_2, fold_record

# ----- Set Testing Variables -----
def test_hybrecord():
    """Test Functions of HybRecord Class"""

    # Setup test variables
    # Test Record construction
    hyb_record_1, hyb_record_2, fold_record = def_hyb_records()  # Reset test vars
    test_id = hyb_record_1.id
    test_seq = hyb_record_1.seq
    
    # Test execution of record construction with seg1_props_naming
    hyb_record_2 = hybkit.HybRecord(id=test_id, seq=test_seq,
                                    seg1_props={'ref_name': 'noname1'},
                                    seg2_props={'ref_name': 'noname2'},
                                    read_count=4,
                                    )
    # Error: disallowed (bad) flag in __init__ constructor.
    with pytest.raises(RuntimeError):
        bad_hyb_record = hybkit.HybRecord(id=test_id, seq=test_seq,
                                          read_count=4, flags={'badflag': True},
                                          fold_record=fold_record)
    # Error: ???
    with pytest.raises(RuntimeError):
        bad_hyb_record = hybkit.HybRecord(id=test_id, seq=test_seq,
                                          read_count=4, flags={'read_count': 4})
    # Error: ???
    with pytest.raises(ValueError):
        bad_hyb_record = hybkit.HybRecord(id=test_id, seq=test_seq,
                                          seg1_props={'read_start': 'not_int'},
                                          read_count=4, flags={'read_count': 4})
    # Error: bad seg1_type flag 
    with pytest.raises(RuntimeError):
        bad_hyb_record = hybkit.HybRecord.from_line(hyb_str_1 + 'seg1_type=badtype', 
                                                    hybformat_ref=True)
    # Error: "Ensure" attribute "energy" not set.
    with pytest.raises(RuntimeError):
        hyb_record_2._ensure_set('energy')

    # Test execution of evaluation functions and detail functions. 
    hyb_record_1, hyb_record_2, fold_record = def_hyb_records()  # Reset test vars
    hyb_record_2.eval_types(allow_unknown=True)
    hyb_record_2._flagset = None
    hyb_record_2._make_flags_dict({})
    hyb_record_2._flagset = None
    hyb_record_2._get_ordered_flag_keys()
    hyb_record_2._flagset = None 
    
    # Test execution of eval_mirna() method with multiple segment permutations
    hyb_record_1, hyb_record_2, fold_record = def_hyb_records()  # Reset test vars
    hyb_record_2.set_flag('seg1_type', 'microRNA')
    hyb_record_2.set_flag('seg2_type', 'microRNA')
    hyb_record_2.eval_mirna()
    del(hyb_record_2.flags['miRNA_seg'])
    hyb_record_2.set_flag('seg1_type', 'mRNA')
    hyb_record_2.set_flag('seg2_type', 'microRNA')
    hyb_record_2.eval_mirna()
    del(hyb_record_2.flags['miRNA_seg'])
    hyb_record_2.set_flag('seg1_type', 'mRNA')
    hyb_record_2.set_flag('seg2_type', 'mRNA')
    hyb_record_2.eval_mirna()

    hyb_record_1, hyb_record_2, fold_record = def_hyb_records()  # Reset test vars
    hyb_record.set_flag('miRNA_seg', 'B')
    # Error: Execute miRNA-detail with dimer without specifically allowing.
    with pytest.raises(RuntimeError):
        hyb_record_2.mirna_detail('all') 
    
    # Class Detail-Retreival Testing
    hyb_record.mirna_detail('all', allow_mirna_dimers=True)
    hyb_record.set_flag('miRNA_seg', '3p')
    hyb_record.mirna_detail('all', allow_mirna_dimers=True)
    hyb_record.to_fasta_record('mirna')
    hyb_record.to_fasta_record('target')
    hyb_record.set_flag('miRNA_seg', 'N')
    with pytest.raises(RuntimeError):
        hyb_record.to_fasta_record('mirna')
    hyb_record.set_flag('miRNA_seg', 'INVALID')
    with pytest.raises(RuntimeError):
        hyb_record.mirna_detail('all')

    hyb_record_2 = hybkit.HybRecord(id=test_id, seq=test_seq,
                                    seg1_props={'ref_name': 'noname1'},
                                    seg2_props={},
                                    )
    with pytest.raises(RuntimeError):
        hyb_record_2.has_prop('seg2_is', 'blah')
    hyb_record_2.to_line()

    hyb_record = hybkit.HybRecord.from_line(hyb_str_1,
                                            hybformat_id=True,
                                            hybformat_ref=True)
    print(str(hyb_record))
    assert hyb_record == hyb_record
    assert not (hyb_record != hyb_record)
    assert bool(hyb_record)
    hash(hyb_record)
    print(len(hyb_record))
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    hyb_record.set_flag('count_total', 1)
    print(hyb_record._format_seg_props(hyb_record.seg1_props))

    test_functions = [
        (hyb_record.get_seg1_type, 'microRNA', [], {}),
        (hyb_record.get_seg1_type, 'microRNA', [True], {}),
        (hyb_record.get_seg2_type, 'mRNA', [], {}),
        (hyb_record.get_seg2_type, 'mRNA', [True], {}),
        (hyb_record.get_seg_types, ('microRNA', 'mRNA'), [], {}),
        (hyb_record.get_seg_types, ('microRNA', 'mRNA'), [True], {}),
        (hyb_record.get_read_count, 804, [], {}),
        (hyb_record.get_read_count, 804, [True], {}),
        (hyb_record.get_record_count, 1, [], {}),
        (hyb_record.get_record_count, 1, [True], {}),
        (hyb_record.get_count, 804, ['read'], {}),
        (hyb_record.get_count, 804, ['read', True], {}),
        (hyb_record.get_count, 1, ['record'], {}),
        (hyb_record.get_count, 1, ['record', True], {}),
    ]
    for fn, expected, args, kwargs in test_functions:
        output = fn(*args, **kwargs)
        print(fn, expected, output, args, kwargs)
        assert output == expected

    test_conversions = [
        (hyb_record.to_line, [], {}),
        (hyb_record.to_line, [True], {}),
        (hyb_record.to_csv, [], {}),
        (hyb_record.to_csv, [True], {}),
        (hyb_record.to_fasta_record, ['hybrid', False], {}),
        (hyb_record.to_fasta_record, ['hybrid', True], {}),
        (hyb_record.to_fasta_record, ['seg1', False], {}),
        (hyb_record.to_fasta_record, ['seg1', True], {}),
        (hyb_record.to_fasta_record, ['seg2', False], {}),
        (hyb_record.to_fasta_record, ['seg2', True], {}),
        (hyb_record.to_fasta_record, ['mirna', False], {}),
        (hyb_record.to_fasta_record, ['mirna', True], {}),
        (hyb_record.to_fasta_record, ['target', False], {}),
        (hyb_record.to_fasta_record, ['target', True], {}),
        (hyb_record.to_fasta_str, ['target', True], {}),
    ]

    for fn, args, kwargs in test_conversions:
        print(fn, args, kwargs)
        output = fn(*args, **kwargs)

    with pytest.raises(RuntimeError):
        hyb_record.to_fasta_record('notallowed')
    hyb_record.flags['miRNA_seg'] = None
    with pytest.raises(RuntimeError):
        hyb_record.to_fasta_record('mirna')
    hyb_record.flags['miRNA_seg'] = 'B'
    with pytest.raises(RuntimeError):
        hyb_record.to_fasta_record('mirna')
    hyb_record.flags['miRNA_seg'] = '5p'

    fail_functions = [
        (hyb_record.get_count, None, ['blah'], {}),
    ]
    for fn, expected, args, kwargs in fail_functions:
        try:
            print(fn, expected, args, kwargs)
            output = fn(*args, **kwargs)
        except RuntimeError as e:
            print(e)

    with pytest.raises(RuntimeError):
        hyb_record.is_set('badprop')
    with pytest.raises(RuntimeError):
        hyb_record.has_prop('badprop')
    with pytest.raises(RuntimeError):
        hyb_record.has_prop('any_seg_type_contains', None)

    assert not hyb_record.has_prop('has_indels')

    for prop in hybkit.HybRecord.STR_PROPS + hybkit.HybRecord.MIRNA_STR_PROPS:
        print('Prop:', prop, hyb_record.has_prop(prop, 'AAAAAAA'))
        assert not (hyb_record.has_prop(prop, 'AAAAAAA'))

    for prop in hybkit.HybRecord.MIRNA_PROPS:
        print('Prop:', prop, hyb_record.has_prop(prop))
        if prop in {'has_mirna', 'mirna_not_dimer', '5p_mirna'}:
            assert hyb_record.has_prop(prop)
        else:
            assert not hyb_record.has_prop(prop)

    hyb_record.set_flag('target_reg', 'U')
    for prop in hybkit.HybRecord.TARGET_PROPS:
        print('Prop:', prop, hyb_record.has_prop(prop))
        if prop == 'target_unknown':
            assert hyb_record.has_prop(prop)
        else:
            assert not hyb_record.has_prop(prop)

    with pytest.raises(RuntimeError):
        hyb_record.set_fold_record(None)
    with pytest.raises(RuntimeError):
        hyb_record.set_fold_record('not_fold_record')
    with pytest.raises(RuntimeError):
        hyb_record.mirna_detail('disallowed_detail')

    # Test Private Methods
    hyb_record.seg1_props['read_start'] = None
    with pytest.raises(RuntimeError):
        hyb_record._ensure_props_read_start_end()
    with pytest.raises(RuntimeError):
        hyb_record._get_seg_seq(hyb_record.seg1_props)
    with pytest.raises(RuntimeError):
        hyb_record._get_flag('fake_flag')
    with pytest.raises(RuntimeError):
        hyb_record._make_flags_dict('not_dict')
    with pytest.raises(RuntimeError):
        hyb_record._make_flags_dict({'bad_flag': True})
    with pytest.raises(RuntimeError):
        hyb_record._parse_hybformat_id('bad_id_name_continues_on')
    with pytest.raises(RuntimeError):
        hyb_record._parse_hybformat_ref('bad_ref_name_continues_on')
    with pytest.raises(RuntimeError):
        hyb_record._read_flags('bad_flag=B;')
    hyb_record._get_flag_keys(reorder_flags=False)
    hyb_record._read_flags('miRNA_seg=B;')
    hyb_record._read_flags('bad_flag=B;', allow_undefined_flags=True)

    hyb_record = hybkit.HybRecord.from_line(hyb_str_1,
                                            hybformat_id=True,
                                            hybformat_ref=True)
    with pytest.raises(RuntimeError):
        hyb_record._ensure_set('type_eval')


def test_hybfile():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)

    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        with pytest.raises(RuntimeError):
            hyb_file.write_record('NotRecord')
        hyb_file.write_record(hyb_record)
        hyb_file.write_records([hyb_record, hyb_record])
    assert hybkit.util.hyb_exists(hyb_file_name)
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file:
        hyb_record_read = hyb_file.read_records()[0]
    assert hyb_record == hyb_record_read
    for test_prop in ['id', 'seq', 'energy', 'seg1_props', 'seg2_props', 'flags']:
        print('Comparing:', test_prop,
              getattr(hyb_record, test_prop), getattr(hyb_record_read, test_prop))
        assert getattr(hyb_record, test_prop) == getattr(hyb_record_read, test_prop)


def test_foldrecord():
    for TestClass in [hybkit.FoldRecord, hybkit.DynamicFoldRecord]:
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_string(vie_str_1, 'bad_mode')
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['', '', ''], 'bad_mode')
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['', ''])
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', 'singleval'])
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', '(99'])
        fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', '(99'],
                                                          error_mode='warn_return')
        fold_record = TestClass.from_vienna_string(vie_str_1)
        print(str(fold_record))
        assert fold_record == fold_record
        assert not (fold_record != fold_record)
        assert bool(fold_record)
        hash(fold_record)
        print(len(fold_record))
        if not isinstance(fold_record, hybkit.DynamicFoldRecord):
            print(fold_record._get_seg_fold({'read_start': 1, 'read_end': 10}))
        print(fold_record.to_vienna_lines())
        print(fold_record.to_vienna_lines(True))
        print(fold_record.to_vienna_string())
        print(fold_record.to_vienna_string(True))


def test_viennafile():
    fold_record = hybkit.FoldRecord.from_vienna_string(vie_str_1)
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


def test_hybfolditer():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        hyb_file.write_records([hyb_record] * 30)
    fold_record = hybkit.FoldRecord.from_vienna_string(vie_str_1)
    dyn_fold_record = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1)
    dyn_fold_record_mis = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1_mis)
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


def test_analyses():
    analysis_classes = {
        'type': hybkit.analysis.TypeAnalysis,
        'mirna': hybkit.analysis.MirnaAnalysis,
        'summary': hybkit.analysis.SummaryAnalysis,
        'target': hybkit.analysis.TargetAnalysis,
    }
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
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

    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
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


def test_fold_analysis():
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    fold_record = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1)
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


def test_util():
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


def test_type_finder():
    # Generic Tests
    def do_nothing(*args, **kwargs):
        pass
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    hybkit.type_finder.TypeFinder.set_custom_method(do_nothing)

    # Defualt Hybformat
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
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
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.settings['check_complete_seg_types'] = True
    with pytest.raises(RuntimeError):
        hyb_record.eval_types()
    # match_params = {'startswith':[('MIMAT', 'MIMAT')], 'endswith':[('microRNA', 'miRNA')]}
    # hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    # with pytest.raises(RuntimeError):
    #     hyb_record.eval_types()
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
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
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.seg1_props['ref_name'] = 'not_real_name'
    with pytest.raises(RuntimeError):
        hyb_record.eval_types(allow_unknown=False)

# if __name__ == '__main__':
#     test_type_finder()
