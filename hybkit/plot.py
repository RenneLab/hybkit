#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Methods for plotting analyses of HybRecord and FoldRecord Objects.
'''

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import hybkit



# Public Methods : HybRecord Type Analysis Plotting
def plot_type_analysis_hybrid_type_counts(analysis_dict, sep=','):
    'Plot the results of hybrid_type_counts in a list of sep-delimited lines.'
    # Sort by count in descending order
    ret_lines = ['hybrid_type' + sep + 'count']
    sorted_pairs = sorted(analysis_dict['hybrid_type_counts'].items(),
                          key=lambda item: item[1], reverse=True)
    other_threshhold = 0.90
    total_count = 0
    for key, count in sorted_pairs:
        total_count += count
    labels = []
    counts = []
    running_total = 0
    for key, count in sorted_pairs:
        running_total += count
        if running_total/total_count > other_threshhold:
            labels.append('Other')
            counts.append(total_count - running_total)
            break
        else:
            labels.append(key)
            counts.append(count)

    _plot_pie_chart(labels,counts)


# Public Methods : HybRecord Type Analysis Parsing
def format_type_analysis_all_seg_types(analysis_dict, sep=','):
    'Return the results of all_seg_types in a list of sep-delimited lines.'
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
                ('combined_types', format_type_analysis),
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
                ('combined', format_mirna_analysis),
               ]
    for analysis_name, analysis_method in analyses:
        analysis_file_name = file_name_base + '_' + analysis_name + file_suffix
        with open(analysis_file_name, 'w') as out_file:
            out_file.write('\n'.join(analysis_method(analysis_dict, sep)))


# Private Methods : Utility
def _check_matplotlib():
    # lazy-importing is being used to allow matplotlib to serve as an 
    #   optional dependency for the package.
    try:
        import matplotlib
    except ImportError:
        message = 'The python package matplotlib is required for plotting funciton.'
        message += '\nPlease install this package before using plotting features.'
        print(message)
        raise

# Private Methods : Pie Chart
def _plot_pie_chart(labels, sizes):
    _check_matplotlib()
    import matplotlib.pyplot as plot
    # Data to plot
    #explode = (0.1, 0, 0, 0)  # explode 1st slice
    
    # Plot
    plot.pie(sizes, labels=labels, # explode=explode, colors=colors,
             autopct='%1.1f%%', shadow=False, startangle=90)
    
    plot.axis('equal')
    plot.show()
    #patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    #plt.legend(patches, labels, loc="best")
    #plt.axis('equal')
    #plt.tight_layout()
    #plt.show()
    
