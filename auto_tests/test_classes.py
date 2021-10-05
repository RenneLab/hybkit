#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""

import os
import sys
import hybkit

hybkit.settings.HybFile_settings['hybformat_id'] = True
hybkit.settings.HybFile_settings['hybformat_ref'] = True
hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
hybkit.settings.Analysis_settings['allow_mirna_dimers'] = True

auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
hyb_file = os.path.join(auto_tests_dir, 'test_hybrid.hyb')
vienna_file = os.path.join(auto_tests_dir, 'test_hybrid.vienna')
out_dir = os.path.join(auto_tests_dir, 'output_autotest')
out_hyb_file = os.path.join(out_dir, hyb_file.replace('.hyb', '_py.hyb'))
out_analysis_basename = out_hyb_file.replace('.hyb', '')
match_legend_file = os.path.join(auto_tests_dir, 'string_match_legend.csv')
with hybkit.HybFile.open(hyb_file, 'r') as hyb_file_obj:
    for HYB_RECORD in hyb_file_obj:
        pass
HYB_RECORD.eval_types()
HYB_RECORD.eval_mirna()
with hybkit.ViennaFile.open(vienna_file, 'r') as vienna_file_obj:
    for FOLD_RECORD in vienna_file_obj:
        pass

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
#count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)


remove_types = ['rRNA', 'mitoch-rRNA']
match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(match_legend_file)
hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)

def test_hybrecord():
    with hybkit.HybFile.open(hyb_file, 'r') as input_hyb:
        for hyb_record in input_hyb:
            assert hyb_record == HYB_RECORD
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
            ]

            for fn, args, kwargs in test_conversions:
                print(fn, args, kwargs)
                output = fn(*args, **kwargs)

            fail_functions = [
                (hyb_record.get_count, None, ['blah'], {}),
            ]
            for fn, expected, args, kwargs in test_functions:
                try:
                    print(fn, expected, args, kwargs)
                    output = fn(*args, **kwargs)
                except Exception as e:
                    print(e)
            
            #test_props = []
            #for prefix in ['


def test_foldrecord():           
    for foldrecord_type in ['strict', 'dynamic']:
        hybkit.settings.FoldFile_settings['foldrecord_type'] = foldrecord_type
        for combine in [True, False]:
            with hybkit.HybFile.open(hyb_file, 'r') as input_hyb, \
                 hybkit.ViennaFile.open(vienna_file, 'r') as input_vienna, \
                 hybkit.HybFile.open(out_hyb_file, 'w') as out_hyb:
                hf_iter = hybkit.HybFoldIter(input_hyb, input_vienna, combine=combine)
                for ret_item in hf_iter:
                    pass
                hf_iter.report()
                hf_iter.print_report()
                

def test_analyses():
    analysis_classes = {
        'type': hybkit.analysis.TypeAnalysis, 
        'mirna': hybkit.analysis.MirnaAnalysis, 
        'summary': hybkit.analysis.SummaryAnalysis, 
        'target': hybkit.analysis.TargetAnalysis, 
    }
    for analysis_type in analysis_classes:
        Analysis = analysis_classes[analysis_type]
        analysis = Analysis(name='Test')
        analysis2 = Analysis()
        analysis2.add(HYB_RECORD)
        analysis.update(analysis2)
        print(analysis.results())
        analysis.write(out_dir + '/py_coverage_' + analysis_type) 
        analysis.plot(out_dir + '/py_coverage_' + analysis_type) 
        if hasattr(analysis, 'write_individual'):
            analysis.write_individual(out_dir + '/py_coverage_' + analysis_type) 
            analysis.plot_individual(out_dir + '/py_coverage_' + analysis_type) 
     
def test_fold_analysis():
    analysis = hybkit.analysis.FoldAnalysis(name='Test')
    analysis2 = hybkit.analysis.FoldAnalysis()
    HYB_RECORD.set_fold_record(FOLD_RECORD)
    #analysis2.add(HYB_RECORD)
    #analysis.update(analysis2)
    #analysis.write(out_dir + '/py_coverage_fold') 
    #analysis.plot(out_dir + '/py_coverage_fold') 

 
