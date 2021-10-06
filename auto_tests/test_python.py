#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""

import os
import sys
import argparse
import pytest
import hybkit

hyb_str_1 = ('695_804	ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC	.	'
             'MIMAT0000078_MirBase_miR-23a_microRNA	1	21	1	21	0.0027	'
             'ENSG00000188229_ENST00000340384_TUBB2C_mRNA	23	49	1181	1207	1.2e-06	'
            )
vie_str_1 = ('>695_804_MIMAT0000078_MirBase_miR-23a_microRNA_1_21-'
             '695_804_ENSG00000188229_ENST00000340384_TUBB2C_mRNA_1181_1207\n'
             'ATCACATTGCCAGGGATTTCCATCCCCAACAATGTGAAAACGGCTGTC\n'
             '.((((((((...(((((....)))))...))))))))...........	(-15)\n'
            )


auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
test_out_dir = os.path.join(auto_tests_dir, 'output_autotest')
hyb_file_name = os.path.join(test_out_dir, 'test_hybrid_py.hyb')
vienna_file_name = os.path.join(test_out_dir, 'test_hybrid_py.vienna')
out_basename = hyb_file_name.replace('.hyb', '')
match_legend_file_name = os.path.join(auto_tests_dir, 'string_match_legend.csv')
bad1_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_1.csv')
bad2_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_2.csv')
bad3_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_3.csv')
id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_id_map.csv')
bad1_id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_id_map_1.csv')


#remove_types = ['rRNA', 'mitoch-rRNA']

#hybkit.settings.HybFile_settings['hybformat_id'] = True
#hybkit.settings.HybFile_settings['hybformat_ref'] = True
#hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
#hybkit.settings.Analysis_settings['allow_mirna_dimers'] = True

def test_hybrecord():
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1, 
                                            hybformat_id=True, 
                                            hybformat_ref=True)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    hyb_record.set_flag('count_total', 1)
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

    fail_functions = [
        (hyb_record.get_count, None, ['blah'], {}),
    ]
    for fn, expected, args, kwargs in fail_functions:
        try:
            print(fn, expected, args, kwargs)
            output = fn(*args, **kwargs)
        except Exception as e:
            print(e)

    with pytest.raises(Exception):
        hyb_record.has_prop('badprop')

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


def test_hybfile():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)

    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        hyb_file.write_record(hyb_record)
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
        fold_record = TestClass.from_vienna_string(vie_str_1)
        print(str(fold_record))
        assert fold_record == fold_record
        assert not (fold_record != fold_record)
        assert bool(fold_record)
        hash(fold_record)
        print(len(fold_record))
        #print(fold_record._get_seg_fold({'read_start':1, 'read_end':10}))
        print(fold_record.to_vienna_lines())
        print(fold_record.to_vienna_lines(True))
        print(fold_record.to_vienna_string())
        print(fold_record.to_vienna_string(True))

       
def test_viennafile():           
    fold_record = hybkit.FoldRecord.from_vienna_string(vie_str_1)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_record(fold_record)               
    assert hybkit.util.vienna_exists(vienna_file_name)
    assert hybkit.util.fold_exists(vienna_file_name)
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
        hyb_file.write_record(hyb_record)
    fold_record = hybkit.FoldRecord.from_vienna_string(vie_str_1)
    dyn_fold_record = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_record(fold_record)               

    #Compare FoldRecord and HybRecord
    assert not fold_record.matches_hyb_record(hyb_record) 
    fold_record.count_hyb_record_mismatches(hyb_record) 
    with pytest.raises(Exception):
        fold_record.ensure_matches_hyb_record(hyb_record) 
  
    #Compare DynamicFoldRecord and HybRecord
    assert dyn_fold_record.matches_hyb_record(hyb_record) 
    dyn_fold_record.count_hyb_record_mismatches(hyb_record) 
    dyn_fold_record.ensure_matches_hyb_record(hyb_record) 

    hybkit.settings.HybFoldIter_settings['error_mode'] = 'raise'
    hybkit.settings.FoldFile_settings['foldrecord_type'] = 'strict'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
            hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
            with pytest.raises(Exception):
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
    with pytest.raises(Exception):
        BaseAnalysis()
    for analysis_type in analysis_classes:
        Analysis = analysis_classes[analysis_type]
        analysis_basename = out_basename + '_' + analysis_type
        analysis = Analysis(name='Test')
        analysis2 = Analysis()
        #with pytest.raises(Exception):
        #    analysis2.plot(analysis_basename)
        analysis2.plot(analysis_basename)
        analysis2.add(hyb_record)
        analysis.update(analysis2)
        with pytest.raises(Exception):
            analysis.update(last_analysis)
        print(analysis.results())
        print(analysis.results(None, True))
        analysis.write(analysis_basename) 
        analysis.plot(analysis_basename)
        if hasattr(analysis, 'write_individual'):
            analysis.write_individual(analysis_basename)
            analysis.plot_individual(analysis_basename)
        last_analysis = analysis

    analysis = hybkit.analysis.TypeAnalysis()
       
 
def test_fold_analysis():
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    fold_record = hybkit.DynamicFoldRecord.from_vienna_string(vie_str_1)
    analysis = hybkit.analysis.FoldAnalysis(name='Test')
    analysis2 = hybkit.analysis.FoldAnalysis()
    hyb_record.set_fold_record(fold_record)
    analysis_basename = out_basename + '_fold'
    with pytest.raises(Exception):
        analysis2.plot(analysis_basename)
    analysis2.add(hyb_record)
    analysis.update(analysis2)
    print(analysis.results())
    print(analysis.results(None, True))
    analysis.write(analysis_basename)
    analysis.plot(analysis_basename)


def test_util():
    original_abspath = hybkit.settings._USE_ABSPATH
    hybkit.settings._USE_ABSPATH = True
    assert hybkit.util._bool_from_string(True)
    assert hybkit.util._bool_from_string('yes')
    assert not hybkit.util._bool_from_string('no')
    with pytest.raises(Exception):
        hybkit.util._bool_from_string('invalid')
    assert hybkit.util.dir_exists('~')
    assert hybkit.util.dir_exists('${PWD}')
    assert hybkit.util.file_exists(__file__)
    assert hybkit.util.out_path_exists(hyb_file_name)

    hybkit.util.make_out_file_name(hyb_file_name, name_suffix='out', 
        in_suffix='.hyb', out_suffix='.new', out_dir='', seg_sep='_')
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

def test_type_finder():
    # Defualt Hybformat
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
    # Non-Defualt String-Match
    with pytest.raises(Exception):
        hybkit.HybRecord.TypeFinder.make_string_match_params('badfile')
        hybkit.HybRecord.TypeFinder.make_string_match_params(bad1_match_legend_file_name)
        hybkit.HybRecord.TypeFinder.make_string_match_params(bad2_match_legend_file_name)
    match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(
        bad3_legend_file_name)
    hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.settings['check_complete_seg_types'] = True
    with pytest.raises(Exception):
        hyb_record.eval_types()
    match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(match_legend_file_name)
    hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
    # Non-Defualt ID-Map
    with pytest.raises(Exception):
        hybkit.HybRecord.TypeFinder.make_id_map_params()
        hybkit.HybRecord.TypeFinder.make_id_map_params(id_map_legend_file_name)
        hybkit.HybRecord.TypeFinder.make_id_map_params(type_file_pairs=id_map_legend_file_name)
        hybkit.HybRecord.TypeFinder.make_id_map_params(bad1_id_map_legend_file_name)
    id_map_params = hybkit.HybRecord.TypeFinder.make_id_map_params(
        type_file_pairs=[('seqtype', id_map_legend_file_name)])
    id_map_params = hybkit.HybRecord.TypeFinder.make_id_map_params([id_map_legend_file_name])
    hybkit.HybRecord.TypeFinder.set_method('id_map', id_map_params)
    hyb_record = hybkit.HybRecord.from_line(hyb_str_1)
    hyb_record.eval_types()
 

