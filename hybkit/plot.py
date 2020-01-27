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
import matplotlib.pyplot as plot

# Public Constants
DEFAULT_PIE_MATPLOTLIB_SETTINGS = {
    autopct='%1.1f%%',
    shadow=False,
    startangle=90
    }
DEFAULT_PIE_DPI = 600
DEFAULT_PIE_FILE_TYPE = 'png'
DEFAULT_PIE_OTHER_THRESHHOLD = 0.1
DEFAULT_PIE_MIN_WEDGE_SIZE = 0.05

# Public Methods : HybRecord Type Analysis Plotting
def plot_hybrid_type_counts(analysis_dict, plot_file_name, 
                    other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                    min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                    plot_file_type=DEFAULT_PIE_FILE_TYPE,
                    dpi=DEFAULT_PIE_DPI,
                    matplotlib_settings=DEFAULT_PIE_MATPLOTLIB_SETTINGS):
    'Plot the results of hybrid_type_counts'
    # Sort hybrid_type_counts counter object in descending order
    for key, count in analysis_dict['hybrid_type_counts'].most_common() 
        labels.append(key)
        counts.append(count)

    _plot_pie_chart(label=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    other_threshhold=other_threshhold,
                    min_wedge_size=min_wedge_size,
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    matplotlib_settings=matplotlib_settings,
                   )


# Private Methods : Pie Chart
def _plot_pie_chart(labels, sizes, plot_file_name,
                    other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD
                    min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                    plot_file_type=DEFAULT_PIE_FILE_TYPE,
                    dpi=DEFAULT_PIE_DPI,
                    matplotlib_settings=DEFAULT_PIE_MATPLOTLIB_SETTINGS):
    total_size = sum(counts)
    fraction_sizes = [size/total_size for size in sizes]
    use_labels = []
    use_sizes = []
    for i in range(len(labels)):
        total_fraction = sum(use_sizes) / total_size
        if ((fraction_size[i] < min_wedge_size)
            or (total_fraction > (1 - other_threshhold) and i != (len(labels) - 1)):
            break
        else:
            use_labels.append(labels[i])
            use_sizes.append(sizes[i])

    if len(use_labels) != len(labels):
        other_size = total_size - sum(sizes)
        use_labels.append('other')
        use_sizes.append(other_size)

    # Reverse ordering to create expected clockwise big->small order.
    use_labels_ordered = reversed(use_labels)
    use_sizes_ordered = reversed(use_sizes)

    plot.pie(use_sizes_ordered, 
             labels=use_labels_ordered, 
             **matplotlib_settings)
    plot.axis('equal')
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plot.savefig(plot_file_name)  

    #patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    #plt.legend(patches, labels, loc="best")
    #plt.axis('equal')
    #plt.tight_layout()
    #plt.show()
    
