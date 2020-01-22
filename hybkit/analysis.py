#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Methods for analyzing HybRecord and FoldRecord Objects.
'''

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __depreciated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import hybkit


# Public Methods : HybRecord Analysis Preparation
def type_analysis_dict():
    'Create a dictionary with keys for running type analyses.'
    ret_dict = {
                'hybrid_type_counts':{},
                'seg1_types':{},
                'seg2_types':{},
                'all_seg_types':{},
               } 
    return ret_dict

# Public Methods : HybRecord Analysis Preparation
def mirna_analysis_dict():
    'Create a dictionary with keys for running miRNA analyses.'
    ret_dict = type_analysis_dict()
    ret_dict.update({
                     '5p_mirna_count':0,
                     '3p_mirna_count':0,
                     'mirna_dimer_count':0,
                     'all_mirna_count':0,
                     'no_mirna_count':0,
                     'mirna_fold_count':0,
                     'mirna_fold_details':{i:0 for i in range(1,26)},
                     'mirna_kmers':{},
                     })
    return ret_dict

# Public Methods : HybRecord Analysis
def running_type_analysis(record, analysis_dict):
    '''
    Add information regarding various properties from the HybRecord object provided in 
    "record" to the dictionary provided in the analysis_dict argument.
    '''
    if not record.has_property('has_seg_types'): 
        message = 'seg_type flag is required for record analysis.'
        print(message)
        raise Exception(message)

    hybrid_type = '---'.join(record.seg_types_sorted())
    seg1_type = record.seg1_type()
    seg2_type = record.seg2_type()
   
    _add_count(analysis_dict['hybrid_type_counts'], hybrid_type)
    _add_count(analysis_dict['seg1_types'], seg1_type)
    _add_count(analysis_dict['seg2_types'], seg2_type)
    for seg_type in seg1_type, seg2_type:
        _add_count(analysis_dict['all_seg_types'], seg_type)


# Public Methods : HybRecord Analysis
def running_mirna_analysis(record, analysis_dict):
    '''
    Add information regarding various properties to the dictionary provided in
    the analysis_dict argument.
    '''
    record._ensure_mirna_analysis()

    # Perform type-analysis on record 
    running_type_analysis(record, analysis_dict)

    # Add mirna-analysis specific details
    if record.has_property('has_mirna_dimer'):
        _add_count(analysis_dict, 'mirna_dimer_count')
        _add_count(analysis_dict, 'all_mirna_count')
    elif record.has_property('3p_mirna'):
        _add_count(analysis_dict, '3p_mirna_count')
        _add_count(analysis_dict, 'all_mirna_count')
    elif record.has_property('5p_mirna'):
        _add_count(analysis_dict, '5p_mirna_count')
        _add_count(analysis_dict, 'all_mirna_count')
    else:
        _add_count(analysis_dict, 'no_mirna_count')

    if record.has_property('has_mirna_fold'):
        _add_count(analysis_dict, 'mirna_fold_count')
        mirna_fold = record.mirna_details['mirna_fold']
        for i in range(1, (len(mirna_fold) + 1)):
            _add_count(analysis_dict['mirna_fold_details'], i)
 

# Public Methods : HybRecord Type Analysis Parsing
def format_type_analysis_hybrid_type_counts(analysis_dict, sep=','):
    'Return the results of hybrid_type_counts in a list of sep-delimited lines.'
    # Sort by count in descending order
    ret_lines = ['hybrid_type' + sep + 'count'] 
    sorted_pairs = sorted(analysis_dict['hybrid_type_counts'].items(), 
                          key=lambda item: item[1], reverse=True)
    ret_lines += ['%s%s%i' % (key, sep, count) for (key, count) in sorted_pairs]
    return ret_lines

# Public Methods : HybRecord Type Analysis Parsing
def format_type_analysis_all_seg_types(analysis_dict, sep=','):
    'Return the results of all_seg_types in a list of sep-delimited lines.'
    print(analysis_dict)
    # Sort by count in descending order
    sorted_pairs = sorted(analysis_dict['all_seg_types'].items(), 
                          key=lambda item: item[1], reverse=True)
    ret_lines = ['seg_type' + sep + 'count'] 
    ret_lines += ['%s%s%i' % (key, sep, count) for (key, count) in sorted_pairs]
    return ret_lines

# Public Methods : HybRecord Type Analysis Parsing
def format_type_analysis(analysis_dict, sep=','):
    'Return the results of a type_analysis in a list of sep-delimited lines.'
    ret_lines = []
    ret_lines += format_type_analysis_hybrid_type_counts(analysis_dict, sep)
    ret_lines.append('')
    ret_lines += format_type_analysis_all_seg_types(analysis_dict, sep)
    return ret_lines

# Public Methods : HybRecord Type Analysis Writing
def write_type_analysis_file(file_name, analysis_dict, sep=','):
    'Write the results of the type-analysis to the file provided in file_name.'
    with open(file_name, 'w') as out_file:
        out_file.write('\n'.join(format_type_analysis(analysis_dict, sep)))

# Public Methods : HybRecord Type Analysis Writing
def write_type_analysis_multi_files(file_name_base, analysis_dict, sep=',', file_suffix='.csv'):
    '''
    Write the results of the type-analysis to a series of files with names based
    on file_name_base.
    '''
    analyses = [
                ('hybrid_types', format_type_analysis_hybrid_type_counts),
                ('seg_types', format_type_analysis_all_seg_types),
               ]
    for analysis_name, analysis_method in analyses:
        analysis_file_name = file_name_base + '_' + analysis_name + file_suffix
        with open(analysis_file_name, 'w') as out_file:
            out_file.write('\n'.join(analysis_method(analysis_dict, sep)))

# Public Methods : HybRecord miRNA Analysis Parsing
def format_mirna_analysis_counts(analysis_dict, sep=','):
    'Return the results of mirna analysis in a list of sep-delimited lines.'
    ret_lines = ['miRNA_type' + sep + 'count']
    for key in ['5p_mirna_count', '3p_mirna_count', 'mirna_dimer_count', 'all_mirna_count',
                'no_mirna_count', ]:
        ret_lines.append('%s%s%i' % (key, sep, analysis_dict[key]))
    return ret_lines

# Public Methods : HybRecord miRNA Analysis Parsing
def format_mirna_analysis(analysis_dict, sep=','):
    'Return the results of a mirna_analysis in a list of sep-delimited lines.'
    ret_lines = []
    ret_lines += format_type_analysis_hybrid_type_counts(analysis_dict, sep)
    ret_lines.append('')
    ret_lines += format_type_analysis_all_seg_types(analysis_dict, sep)
    ret_lines.append('')
    ret_lines += format_mirna_analysis_counts(analysis_dict, sep)
    return ret_lines

# Public Methods : HybRecord miRNA Analysis Writing
def write_mirna_analysis_file(file_name, analysis_dict, sep=','):
    'Write the results of the mirna-analysis to the file provided in file_name.'
    with open(file_name, 'w') as out_file:
        out_file.write('\n'.join(format_mirna_analysis(analysis_dict, sep)))

# Public Methods : HybRecord Type Analysis Writing
def write_mirna_analysis_multi_files(file_name_base, analysis_dict, sep=',', file_suffix='.csv'):
    '''
    Write the results of the mirna_analysis to a series of files with names based
    on file_name_base.
    '''
    analyses = [
                ('hybrid_types', format_type_analysis_hybrid_type_counts),
                ('seg_types', format_type_analysis_all_seg_types),
                ('mirna', format_mirna_analysis_counts),
               ]
    for analysis_name, analysis_method in analyses:
        analysis_file_name = file_name_base + '_' + analysis_name + file_suffix
        with open(analysis_file_name, 'w') as out_file:
            out_file.write('\n'.join(analysis_method(analysis_dict, sep)))


# Private Methods : Utility
def _add_count(count_dict, key):
    if key not in count_dict:
        count_dict[key] = 0
    count_dict[key] += 1 
