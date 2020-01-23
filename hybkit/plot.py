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

# Private Methods : Utility
def _check_matplotlib():
    # lazy-importing is being used to allow matplotlib to serve as an 
    #   optional dependency for the package.
    try:
        import matplotlib
    except ModuleNotFoundError:
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
    
