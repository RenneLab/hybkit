#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""Methods for plotting analyses of HybRecord and FoldRecord Objects."""

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import hybkit
import matplotlib
import matplotlib.pyplot as plot
import copy

# Public Constants
DEFAULT_PIE_MATPLOTLIB_SETTINGS = {
    'autopct':'%1.1f%%',
    'shadow':False,
    'startangle':90,
    'counterclock':False,
    }
DEFAULT_PIE_RC_PARAMS = {
    }
DEFAULT_DPI = 600
DEFAULT_FILE_TYPE = 'png'
DEFAULT_PIE_OTHER_THRESHHOLD = 0.1
DEFAULT_PIE_MIN_WEDGE_SIZE = 0.04
DEFAULT_LINE_RC_PARAMS = {
    'axes.titlesize':'x-large',
    'axes.labelsize':'large',
    'xtick.labelsize':'large',
    'ytick.labelsize':'large',
    }
DEFAULT_LINE_DATA_FORMAT = '-'
DEFAULT_LINE_MIN_FRACTION_SIZE = 0.01
DEFAULT_HYBRID_TYPE_COUNTS_TITLE = 'Hybrid Types'
DEFAULT_ALL_SEG_TYPE_COUNTS_TITLE = 'Total Segment Portions'
DEFAULT_MIRNA_COUNTS_TITLE = 'Hybrid Types'
DEFAULT_FIG_SIZE = matplotlib.rcParams['figure.figsize']  # 6.4, 4.8 in
TITLE_PAD = 15
FORMAT_NAME_MAP = {'5p_mirna_hybrids':"5'_miRNA_Hybrids",
                   '3p_mirna_hybrids':"3'_miRNA_Hybrids",
                   'mirna_dimer_hybrids':'miRNA_Duplexes',
                   'no_mirna_hybrids':'Non-miRNA_Hybrids'}



# Public Methods : HybRecord Type Analysis Plotting
def hybrid_type_count(plot_file_name, analysis_dict,
                      name=None,
                      title=DEFAULT_HYBRID_TYPE_COUNTS_TITLE,
                      other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                      min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                      plot_file_type=DEFAULT_FILE_TYPE,
                      dpi=DEFAULT_DPI,
                      matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    """
    Plot the counts of types resulting from the :ref:`type_analysis`.

    Args:
        plot_file_name (str): File name for output plots.
        analysis_dict (dict): Dict from type analysis (see :func:`type_dict`).
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge. 
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options: 
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function
    """

    # Sort hybrid_type_counts counter object in descending order
    labels = []
    counts = []
    for key, count in analysis_dict['hybrid_type_counts'].most_common(): 
        labels.append(key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=other_threshhold,
                    min_wedge_size=min_wedge_size,
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    matplotlib_settings=matplotlib_settings,
                   )


# Public Methods : HybRecord Type Analysis Plotting
def all_seg_type(plot_file_name, analysis_dict,
                 name=None,
                 title=DEFAULT_ALL_SEG_TYPE_COUNTS_TITLE,
                 other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                 min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                 plot_file_type=DEFAULT_FILE_TYPE,
                 dpi=DEFAULT_DPI,
                 matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    """
    Plot the counts of segment types resulting from the :ref:`type_analysis`.

    Args:
        plot_file_name (str): File name for output plot.
        analysis_dict (dict): Dict from type analysis 
            (see :func:`type_dict`).
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.
    """

    # Sort hybrid_type_counts counter object in descending order
    labels = []
    counts = []
    for key, count in analysis_dict['all_seg_types'].most_common(): 
        labels.append(key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=other_threshhold,
                    min_wedge_size=min_wedge_size,
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    matplotlib_settings=matplotlib_settings,
                   )


# Public Methods : HybRecord Type Analysis Plotting
def mirna_count(plot_file_name, analysis_dict,
                name=None,
                title=DEFAULT_MIRNA_COUNTS_TITLE,
                other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                plot_file_type=DEFAULT_FILE_TYPE,
                dpi=DEFAULT_DPI,
                matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    """
    Plot the results of the :ref:`mirna_count_analysis`.

    Args:
        plot_file_name (str): File name for output plot.
        analysis_dict (dict): Dict from mirna count analysis 
            (see :func:`~hybkit.mirna_count_dict`).
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass 
            to the matplotlib.plot() function.
        """

    # Sort hybrid_type_counts counter object in descending order
    labels = []
    counts = []
    mirna_counts_keys = hybkit.analysis.MIRNA_COUNT_ANALYSIS_KEYS
    use_dict = {}
    for key in mirna_counts_keys:
        if key != 'all_mirna_hybrids':
            use_dict[key] = analysis_dict[key]
    types_by_count = sorted(use_dict.items(), key=lambda item: item[1], reverse=True)
    for key, count in types_by_count:
        if key in FORMAT_NAME_MAP:
            labels.append(FORMAT_NAME_MAP[key])
        else:
            labels.append(key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=(-1),
                    min_wedge_size=(-1),
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    matplotlib_settings=matplotlib_settings,
                   )


# Public Methods : HybRecord miRNA Target Analysis Plotting
def mirna_target(plot_file_name, mirna_name, mirna_targets_dict, 
                 title=None,
                 name=None,
                 other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                 min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                 plot_file_type=DEFAULT_FILE_TYPE,
                 dpi=DEFAULT_DPI,
                 matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    """
    Plot the targets of a single mirna from the :ref:`mirna_target_analysis`.

    Args:
        plot_file_name (str): File name for output plot.
        mirna_name (str): Name of miRNA to plot for title.
        mirna_targets_dict (dict): Dict from mirna_target analysis 
            (see :func:`~hybkit.mirna_target_dict`).
        title (str, optional): Title / header for plot (replaces default title).
        name (str, optional): Data source name to prepend to title.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.
    """

    labels = []
    counts = []
    for ((target, target_seg_type), count) in mirna_targets_dict.most_common(): 
        labels.append(target)
        counts.append(count)

    if title is None:
        title = 'Targets of ' + str(mirna_name) 

    if name is not None:
        title = str(name) + ': ' + title

    matplotlib_settings.update({'textprops':{'size':'small'},
                               })


    # Set figsize:
    # Width set at constant 9(inches)
    # Height 3 for >= 45 character labels, ratio of 9:3 (more rectangular)
    # Height 3 + (0.16 * count-20) for: 20 < count < 45 labels
    # Height 6.5 for <= 20 character labels, ratio of 9:6.5 (more square)
    longest = max((len(label) for label in labels))
    add_count = max(longest-20, 0)  
    max_height = 6.5
    min_height = 3.0
    height = max((max_height - (0.16 * add_count)), (3.0))  # Ensure no less than min_height
    height = min(height, max_height)                        # Ensure no more than max_height
    figsize = (9, height)

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=other_threshhold,
                    min_wedge_size=(0.025),
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    figsize=figsize,
                    matplotlib_settings=matplotlib_settings,
                    )


# Public Methods : HybRecord miRNA Target Analysis Plotting
def mirna_target_type(plot_file_name, mirna_name, mirna_target_type_counts_dict, 
                      title=None,
                      name=None,
                      other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                      min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                      plot_file_type=DEFAULT_FILE_TYPE,
                      dpi=DEFAULT_DPI,
                      matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    """
    Plot the targets types of a single mirna from the :ref:`mirna_target_analysis`.

    Args:
        plot_file_name (str): File name for output plot.
        mirna_name (str): Name of miRNA to plot for title.
        mirna_targets_dict (dict): Dict from mirna_target analysis
            (see :func:`~hybkit.mirna_target_dict`).
        title (str, optional): Title / header for plot (replaces default title).
        name (str, optional): Data source name to prepend to title.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.                                                          """

    labels = []
    counts = []
    for key, count in mirna_target_type_counts_dict.most_common(): 
        labels.append(key)
        counts.append(count)

    if title is None:
        title = str(mirna_name) + ' Target Types'

    if name is not None:
        title = str(name) + ': ' + title

    matplotlib_settings.update({'textprops':{'size':'small'},
                               })

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=other_threshhold,
                    # min_wedge_size=(0.025),
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    matplotlib_settings=matplotlib_settings,
                    )

# Public Methods : HybRecord miRNA Target Analysis Plotting
def mirna_fold(plot_file_name, fold_analysis_dict,
               title=None,
               name=None,
               data_format=DEFAULT_LINE_DATA_FORMAT,
               min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
               plot_file_type=DEFAULT_FILE_TYPE,
               dpi=DEFAULT_DPI,
               matplotlib_settings=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)):
    """
    Plot the bound percentage of mirna by base from the :ref:`mirna_fold_analysis`.

    Args:
        plot_file_name (str): File name for output plot.
        fold_analysis_dict (dict): Dict from mirna_fold analysis
            (see :func:`~hybkit.mirna_fold_dict`).
        title (str, optional): Title / header for plot (replaces default title).
        name (str, optional): Data source name to prepend to title.
        data_format (str, optional): matplotlib line/data format.
        min_fractione_size (float, optional): Minimum fraction to include at tail
            end of plot. Setting to 0 includes all bases evaluated.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        matplotlib_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.                                                          """

    labels = []
    fractions = []
    max_fraction = 0.0
    for seq_i, fraction in fold_analysis_dict['base_fractions'].items():
        max_fraction = max(max_fraction, fraction)
        if seq_i > 22:
            if fraction/max_fraction < min_fraction_size:
                break
        labels.append(seq_i)
        fractions.append(fraction)

    if title is None:
        title = 'Bound miRNA by Position'

    if name is not None:
        title = str(name) + ': ' + title

    #matplotlib_settings.update({'textprops':{'size':'small'},
    #                           })


    _plot_line(labels=labels,
               sizes=fractions,
               plot_file_name=plot_file_name,
               title=title,
               data_format=data_format,
               min_fraction_size=min_fraction_size,
               plot_file_type=plot_file_type,
               dpi=dpi,
               matplotlib_settings=matplotlib_settings,
               )


# Private Methods : Pie Chart
def _plot_pie_chart(labels, sizes, plot_file_name,
                    title=None,
                    other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                    min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                    plot_file_type=DEFAULT_FILE_TYPE,
                    dpi=DEFAULT_DPI,
                    figsize=DEFAULT_FIG_SIZE,
                    matplotlib_settings=copy.deepcopy(DEFAULT_PIE_MATPLOTLIB_SETTINGS)):
    total_size = sum(sizes)
    if total_size < 0.00000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    fraction_sizes = [size/total_size for size in sizes]
    use_labels = []
    use_sizes = []
    for i in range(len(labels)):
        total_fraction = sum(use_sizes) / total_size
        if fraction_sizes[i] < min_wedge_size:
            break
        elif total_fraction > (1 - other_threshhold) and i != (len(labels) - 1):
            break
        else:
            use_labels.append(labels[i])
            use_sizes.append(sizes[i])

    if len(use_labels) != len(labels):
        other_size = total_size - sum(use_sizes)
        use_labels.append('other')
        use_sizes.append(other_size)

    fig = plot.gcf()
    fig.set_size_inches(figsize)
    patches, texts, autotexts = plot.pie(use_sizes,
                                         labels=use_labels,
                                         **matplotlib_settings)
    plot.axis('equal')
    if title is not None:
        plot.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plot.savefig(plot_file_name, dpi=dpi)  
    plot.clf()
    plot.close()

    #patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    #plt.legend(patches, labels, loc="best")
    #plt.axis('equal')
    #plt.tight_layout()
    #plt.show()
   
 
# Private Methods : Pie Chart
def _plot_line(labels, sizes, plot_file_name,
                    title=None,
                    min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
                    plot_file_type=DEFAULT_FILE_TYPE,
                    data_format=DEFAULT_LINE_DATA_FORMAT,
                    dpi=DEFAULT_DPI,
                    figsize=DEFAULT_FIG_SIZE,
                    matplotlib_settings=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)):
    max_size = max(sizes)
    if abs(max_size) < 0.00000000000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    fraction_sizes = [size/max_size for size in sizes]
    use_labels = []
    use_sizes = []
    for i in range(len(labels)):
        if fraction_sizes[i] > min_fraction_size:
            use_labels.append(labels[i])
            use_sizes.append(sizes[i] * 100)

    plot.rcParams.update(matplotlib_settings)

    fig = plot.gcf()
    fig.set_size_inches(figsize)
    lines = plot.plot(use_labels, use_sizes, data_format)
    plot.xlim(left=1)
    plot.xticks(range(1,len(labels)+1,2))
    plot.xlabel('miRNA Base Position Index')
    plot.ylabel('Percent Bound')
    if title is not None:
        plot.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plot.savefig(plot_file_name, dpi=dpi)  
    plot.clf()
    plot.close()
    matplotlib.rcdefaults()

    #patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    #plt.legend(patches, labels, loc="best")
    #plt.axis('equal')
    #plt.tight_layout()
    #plt.show()
    
