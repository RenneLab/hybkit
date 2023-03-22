#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Methods for plotting analyses of HybRecord and FoldRecord objects."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import Counter
import copy
import hybkit
import hybkit.analysis

# Import module-level dunder-names:
from hybkit.__about__ import (__author__, __contact__, __credits__, __date__, __deprecated__,
                              __email__, __license__, __maintainer__, __status__, __version__)

# Public Constants
DEFAULT_PIE_PLOT_SETTINGS = {
    'autopct': '%1.1f%%',
    'shadow': False,
    'startangle': 90,
    'counterclock': False,
}
DEFAULT_PIE_PLOT_TARGET_SETTINGS = copy.deepcopy(DEFAULT_PIE_PLOT_SETTINGS)
DEFAULT_PIE_PLOT_TARGET_SETTINGS.update(
    {'textprops': {'size': 'small'}}
)

DEFAULT_PIE_RC_PARAMS = {
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
}
DEFAULT_DPI = 1200
DEFAULT_FILE_TYPE = 'png'
DEFAULT_PIE_OTHER_THRESHHOLD = 0.1
DEFAULT_PIE_MIN_WEDGE_SIZE = 0.04
DEFAULT_LINE_RC_PARAMS = {
    'axes.titlesize': 'x-large',
    'axes.titleweight': 'bold',
    'axes.labelsize': 'large',
    'axes.labelweight': 'bold',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'large',
}
DEFAULT_LINE_DATA_FORMAT = '-'
DEFAULT_LINE_MIN_FRACTION_SIZE = 0.01
DEFAULT_BAR_RC_PARAMS = {
    'axes.titlesize': 'x-large',
    #'axes.titleweight': 'bold',
    #'axes.labelsize': 'large',
    'axes.labelweight': 'bold',
    #'xtick.labelsize': 'large',
    #'ytick.labelsize': 'large',
}
DEFAULT_ENERGY_MIN_COUNT = 0
DEFAULT_ENERGY_MIN_DENSITY = 0.0
DEFAULT_MATCH_MIN_COUNT = 0
DEFAULT_BAR_WIDTH = 0.9
DEFAULT_BAR_ALIGN = 'edge'
DEFAULT_ENERGY_XLABEL = 'Hybrid Gibbs Free Energy (kcal/mol)'
DEFAULT_ENERGY_YLABEL = 'Hybrid Count'
DEFAULT_MATCH_XLABEL = 'Predicted Base-Pairs'
DEFAULT_MATCH_YLABEL = 'Hybrid Count'
DEFAULT_HYBRID_TYPE_COUNTS_TITLE = 'Hybrid Types'
DEFAULT_TARGET_TITLE = 'Targets'
DEFAULT_TARGET_TYPE_TITLE = 'Target Types'
DEFAULT_MIRNA_COUNTS_TITLE = 'Hybrid Types'
DEFAULT_PATTERN_TITLE = 'Bound miRNA by Position'
DEFAULT_FOLD_ENERGY_COUNTS_TITLE = 'Hybrid Predicted Energies'
DEFAULT_FOLD_ENERGY_DENSITIES_TITLE = 'Normalized Hybrid Predicted Energies'
DEFAULT_FOLD_MATCH_TITLE = 'Hybrid Predicted Base Pairs'
DEFAULT_FIG_SIZE = matplotlib.rcParams['figure.figsize']  # 6.4, 4.8 in
TITLE_PAD = 15
FORMAT_NAME_MAP = {'mirnas_5p': "5'_miRNA_Hybrids",
                   'mirnas_3p': "3'_miRNA_Hybrids",
                   'mirna_dimers': 'miRNA_Duplexes',
                   'non_mirnas': 'Non-miRNA_Hybrids'
                   }


# Public Methods : HybRecord Type Analysis Plotting
def type_count(type_counter,
               plot_file_name,
               name=None,
               title=DEFAULT_HYBRID_TYPE_COUNTS_TITLE,
               other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
               min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
               plot_file_type=DEFAULT_FILE_TYPE,
               dpi=DEFAULT_DPI,
               rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
               plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_SETTINGS)):
    """
    Plot the counts of types resulting from the :class:`~hybkit.analysis.TypeAnalysis`.

    Args:
        type_counter (Counter): Counter containing type information from
            :class:`~hybkit.analysis.hybkit.analysis.TypeAnalysis`.
        plot_file_name (str): File name for output plots.
        name (str, optional): name to prepend to title.
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
        plot_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function
    """
    # Sort hybrid_type_counts counter object in descending order
    labels = []
    counts = []
    for key, count in type_counter.most_common():
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
                    rc_params=rc_params,
                    plot_settings=plot_settings,
                    )


# Public Methods : HybRecord Type Analysis Plotting
def mirna(analysis,
          plot_file_name,
          title=DEFAULT_MIRNA_COUNTS_TITLE,
          other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
          min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
          plot_file_type=DEFAULT_FILE_TYPE,
          dpi=DEFAULT_DPI,
          rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
          plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_SETTINGS)
          ):
    """
    Plot the results of the :class:`~hybkit.analysis.MirnaAnalysis`.

    Args:
        analysis (Analysis): Analysis that includes attributes of
            :class:`~hybkit.analysis.MirnaAnalysis`.
        plot_file_name (str): File name for output plot.
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
        plot_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.
    """
    # Sort hybrid_type_counts counter object in descending order
    labels = []
    counts = []
    use_items = Counter()
    for key in ['mirnas_5p', 'mirnas_3p', 'mirna_dimers', 'non_mirnas']:
        use_items[key] = getattr(analysis, key)

    for key, count in use_items.most_common():
        if key in FORMAT_NAME_MAP:
            labels.append(FORMAT_NAME_MAP[key])
        else:
            labels.append(key)
        counts.append(count)

    if analysis.name is not None:
        title = str(analysis.name) + ': ' + title

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=(-1),
                    min_wedge_size=(-1),
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    rc_params=rc_params,
                    plot_settings=plot_settings,
                    )


# Public Methods : HybRecord miRNA Target Analysis Plotting
def target(mirna_counts, plot_file_name,
           name=None,
           title=DEFAULT_TARGET_TITLE,
           other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
           min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
           plot_file_type=DEFAULT_FILE_TYPE,
           dpi=DEFAULT_DPI,
           rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
           plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_TARGET_SETTINGS)):
    """
    Plot the targets of a single mirna from the :class:`~hybkit.analysis.MirnaAnalysis`.

    Args:
        mirna_counts (~collection.Counter): Counter of targets of a mirna produced by
            :class:`~hybkit.analysis.TargetAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str, optional): Name of analysis
        title (str, optional): Title / header for plot.
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
        plot_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.
    """
    labels = []
    counts = []
    for ((target, target_seg_type), count) in mirna_counts.most_common():
        labels.append(target)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    # Set figsize:
    # Width set at constant 9(inches)
    # Height 3 for >= 45 character labels, ratio of 9:3 (more rectangular)
    # Height 3 + (0.16 * count-20) for: 20 < count < 45 labels
    # Height 6.5 for <= 20 character labels, ratio of 9:6.5 (more square)
    longest = max([0] + [len(label) for label in labels])
    add_count = max((longest - 20), 0)
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
                    rc_params=rc_params,
                    plot_settings=plot_settings,
                    )


# Public Methods : HybRecord miRNA Target Analysis Plotting
def target_type(target_type_count, plot_file_name,
                name=None,
                title=DEFAULT_TARGET_TYPE_TITLE,
                other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                plot_file_type=DEFAULT_FILE_TYPE,
                dpi=DEFAULT_DPI,
                rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
                plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_TARGET_SETTINGS)):
    """
    Plot the targets types of a single mirna from the :class:`~hybkit.analysis.TargetAnalysis`.

    Args:
        target_type_count (~collections.Counter): Counter with types of targets for miRNAs from
            :class:`~hybkit.analysis.TargetAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str): Name of analysis to plot for title.
        title (str, optional): Title / header for plot (replaces default title).
        other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
            Setting to 0.0 disables the "other" wedge based on a threshhold.
        min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
            wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
        plot_settings (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.plot() function.
    """
    labels = []
    counts = []
    for key, count in target_type_count.most_common():
        labels.append(key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    _plot_pie_chart(labels=labels,
                    sizes=counts,
                    plot_file_name=plot_file_name,
                    title=title,
                    other_threshhold=other_threshhold,
                    # min_wedge_size=(0.025),
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    rc_params=rc_params,
                    plot_settings=plot_settings,
                    )


# Public Methods : HybRecord miRNA Binding Pattern Plotting
def pattern(analysis,
            plot_file_name,
            name=None,
            title=DEFAULT_PATTERN_TITLE,
            data_format=DEFAULT_LINE_DATA_FORMAT,
            min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
            plot_file_type=DEFAULT_FILE_TYPE,
            dpi=DEFAULT_DPI,
            rc_params=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)
            ):
    """
    Plot the bound percentage of mirna by base from the :class:`~hybkit.analysis.PatternAnalysis`.

    Args:
        analysis (~hybkit.analysis.PatternAnalysis): Analysis that includes attributes of
            :class:`~hybkit.analysis.PatternAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str, optional): name to prepend to title.
        title (str, optional): Title / header for plot (replaces default title).
        data_format (str, optional): matplotlib line/data format.
        min_fraction_size (float, optional): Minimum fraction to include at tail
            end of plot. Setting to 0 includes all bases evaluated.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
    """
    labels = []
    fractions = []
    max_fraction = 0.0
    for seq_i, fraction in analysis.mirna_fold_frac.items():
        max_fraction = max(max_fraction, fraction)

    if max_fraction < 0.00000000000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return

    for seq_i, fraction in analysis.mirna_fold_frac.items():
        if seq_i > 22 and (fraction / max_fraction) < min_fraction_size:
            break
        labels.append(seq_i)
        fractions.append(fraction)

    if analysis.name is not None:
        title = str(analysis.name) + ': ' + title

    _plot_line(labels=labels,
               sizes=fractions,
               plot_file_name=plot_file_name,
               title=title,
               data_format=data_format,
               min_fraction_size=min_fraction_size,
               plot_file_type=plot_file_type,
               dpi=dpi,
               rc_params=rc_params,
               )


# Public Methods : HybRecord Fold Analysis Plotting
def fold_energy_counts(analysis,
                       plot_file_name,
                       name=None,
                       title=DEFAULT_FOLD_ENERGY_COUNTS_TITLE,
                       min_count=DEFAULT_ENERGY_MIN_COUNT,
                       plot_file_type=DEFAULT_FILE_TYPE,
                       dpi=DEFAULT_DPI,
                       rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
                       ):
    """
    Plot barogram of hybrid energy counts from a `~hybkit.analysis.FoldAnalysis`.

    Args:
        analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
            :class:`~hybkit.analysis.FoldAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str, optional): name to prepend to title.
        title (str, optional): Title / header for plot (replaces default title).
        min_count (float, optional): Minimum count to include at right
            side of plot of energy counts. Setting to 0 includes all values evaluated.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
    """
    # Check if empty
    max_count = 0
    for energy, count in analysis.energy_bins.items():
        max_count = max(max_count, count)

    if max_count < 1:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    
    # Plot energy bins
    energies, counts = [], []
    for energy, count in analysis.energy_bins.items():
        if float(energy) < -5 and count < min_count:
            break
        energies.append(float(energy))
        counts.append(count)

    if analysis.name is not None:
        title = str(analysis.name) + ': ' + title

    energy_width = float(energies[1]) + 0.02

    _plot_energy_bar(x_vals=energies,
                     sizes=counts,
                     plot_file_name=plot_file_name,
                     title=title,
                     width=energy_width,
                     align=DEFAULT_BAR_ALIGN,
                     plot_file_type=plot_file_type,
                     dpi=dpi,
                     rc_params=rc_params,
                     )


# Public Methods : HybRecord Fold Analysis Plotting
def fold_energy_densities(analysis,
                          plot_file_name,
                          name=None,
                          title=DEFAULT_FOLD_ENERGY_DENSITIES_TITLE,
                          min_density=DEFAULT_ENERGY_MIN_DENSITY,
                          plot_file_type=DEFAULT_FILE_TYPE,
                          dpi=DEFAULT_DPI,
                          rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
                          ):
    """
    Plot barogram of hybrid energy counts from a `~hybkit.analysis.FoldAnalysis`.

    Args:
        analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
            :class:`~hybkit.analysis.FoldAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str, optional): name to prepend to title.
        title (str, optional): Title / header for plot (replaces default title).
        min_density (float, optional): Minimum density to include at right
            side of plot of energy counts. Setting to 0 includes all values evaluated.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
    """
    # Check if empty
    max_density = 0
    for energy, density in analysis.energy_bin_densities.items():
        max_density = max(max_density, density)

    if max_density < 0.0000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    
    # Plot energy bins
    energies, densities = [], []
    for energy, density in analysis.energy_bin_densities.items():
        if float(energy) < -5 and density < min_density:
            break
        energies.append(float(energy))
        densities.append(density)

    if analysis.name is not None:
        title = str(analysis.name) + ': ' + title

    energy_width = float(energies[1]) + 0.02

    _plot_energy_bar(x_vals=energies,
                     sizes=densities,
                     plot_file_name=plot_file_name,
                     title=title,
                     width=energy_width,
                     align=DEFAULT_BAR_ALIGN,
                     plot_file_type=plot_file_type,
                     dpi=dpi,
                     rc_params=rc_params,
                     )


# Public Methods : HybRecord Fold Analysis Plotting
def fold_matches(analysis,
                 plot_file_name,
                 name=None,
                 title=DEFAULT_FOLD_MATCH_TITLE,
                 min_count=DEFAULT_MATCH_MIN_COUNT,
                 plot_file_type=DEFAULT_FILE_TYPE,
                 dpi=DEFAULT_DPI,
                 rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
                 ):
    """
    Plot predicted hybrid match counts from a `~hybkit.analysis.FoldAnalysis`.

    Args:
        analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
            :class:`~hybkit.analysis.FoldAnalysis`.
        plot_file_name (str): File name for output plot.
        name (str, optional): name to prepend to title.
        title (str, optional): Title / header for plot (replaces default title).
        min_count (float, optional): Minimum count to include at right
            side of plot of energy counts. Setting to 0 includes all values evaluated.
        plot_file_type (str, optional): File type for saving of plots. Options:
            {'png', 'ps', 'pdf', 'svg'}
        dpi (int, optional): DPI for saving of plots.
        rc_params (dict, optional): Dict of keys and values of settings to pass
            to the matplotlib.rcParams.update() method.
    """
    # Check if empty
    max_count = 0
    for base, count in analysis.match_counts.items():
        max_count = max(max_count, count)

    if max_count < 1:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    
    # Plot match count bins
    total_count = sum(analysis.match_counts.values())
    bases, counts, densities = [], [], []
    for base, count in analysis.match_counts.items():
        if base > 10 and count < min_count:
            break
        bases.append(base)
        counts.append(count)
        densities.append(count/total_count)
 
    density_title = 'Normalized ' + title
    if analysis.name is not None:
        title = str(analysis.name) + ': ' + title
        density_title = str(analysis.name) + ': ' + density_title

    counts_plot_file_name = plot_file_name + '_counts' 
    densities_plot_file_name = plot_file_name + '_densities' 

    _plot_match_bar(x_vals=bases,
                    sizes=counts,
                    plot_file_name=counts_plot_file_name,
                    title=title,
                    align='center',
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    rc_params=rc_params,
                    )

    _plot_match_bar(x_vals=bases,
                    sizes=densities,
                    plot_file_name=densities_plot_file_name,
                    title=density_title,
                    ylabel='Normalized ' + DEFAULT_MATCH_YLABEL,
                    align='center',
                    plot_file_type=plot_file_type,
                    dpi=dpi,
                    rc_params=rc_params,
                    )


# Private Methods : Pie Chart
def _plot_pie_chart(labels, sizes, plot_file_name,
                    title=None,
                    other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
                    min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
                    plot_file_type=DEFAULT_FILE_TYPE,
                    dpi=DEFAULT_DPI,
                    figsize=DEFAULT_FIG_SIZE,
                    rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
                    plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_SETTINGS)):
    total_size = sum(sizes)
    if total_size < 0.00000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    fraction_sizes = [(size / total_size) for size in sizes]
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

    plt.rcParams.update(rc_params)

    fig = plt.gcf()
    fig.set_size_inches(figsize)
    patches, texts, autotexts = plt.pie(use_sizes,
                                        labels=use_labels,
                                        **plot_settings
                                        )
    plt.axis('equal')
    if title is not None:
        plt.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plt.savefig(plot_file_name, dpi=dpi)
    plt.clf()
    plt.close()

    # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    # plt.legend(patches, labels, loc="best")
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()


# Private Methods : Plot Line Graph
def _plot_line(labels, sizes, plot_file_name,
               title=None,
               min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
               plot_file_type=DEFAULT_FILE_TYPE,
               data_format=DEFAULT_LINE_DATA_FORMAT,
               dpi=DEFAULT_DPI,
               figsize=DEFAULT_FIG_SIZE,
               rc_params=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)
               ):
    max_size = max(sizes)
    if abs(max_size) < 0.00000000000001:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return
    fraction_sizes = [(size / max_size) for size in sizes]
    use_labels = []
    use_sizes = []
    for i in range(len(labels)):
        if fraction_sizes[i] > min_fraction_size:
            use_labels.append(labels[i])
            use_sizes.append(sizes[i] * 100)

    plt.rcParams.update(rc_params)

    fig = plt.gcf()
    fig.set_size_inches(figsize)
    lines = plt.plot(use_labels, use_sizes, data_format)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    plt.xlim(left=1)
    plt.xticks(range(1, (len(labels) + 1), 2))
    plt.xlabel('miRNA Base Index')
    plt.yticks(range(10, 110, 10))
    plt.ylabel('Percent Bound')
    if title is not None:
        plt.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plt.savefig(plot_file_name, dpi=dpi)
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()

    # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    # plt.legend(patches, labels, loc="best")
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()

# Private Methods : Plot Energy Bar Graph
def _plot_energy_bar(x_vals, sizes, plot_file_name,
                     title=None,
                     width=DEFAULT_BAR_WIDTH,
                     align=DEFAULT_BAR_ALIGN,
                     xlabel=DEFAULT_ENERGY_XLABEL,
                     ylabel=DEFAULT_ENERGY_YLABEL,
                     plot_file_type=DEFAULT_FILE_TYPE,
                     dpi=DEFAULT_DPI,
                     figsize=DEFAULT_FIG_SIZE,
                     rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
                     ):

    plt.rcParams.update(rc_params)

    fig = plt.gcf()
    fig.set_size_inches(figsize)
    bars = plt.bar(x_vals, sizes, width=width, align=align)
    ax = plt.gca()
    ax.invert_xaxis()
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    plt.xlim(left=0)
    plt.xticks(range(0, (int(x_vals[-1])-1), (-2)), rotation=-30)
    plt.xlabel(xlabel)
    #plt.yticks(range(10, 110, 10))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plt.savefig(plot_file_name, dpi=dpi)
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()

    # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    # plt.legend(patches, labels, loc="best")
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()

# Private Methods : Plot Bases Bar Graph
def _plot_match_bar(x_vals, sizes, plot_file_name,
                    title=None,
                    width=DEFAULT_BAR_WIDTH,
                    align=DEFAULT_BAR_ALIGN,
                    xlabel=DEFAULT_MATCH_XLABEL,
                    ylabel=DEFAULT_MATCH_YLABEL,
                    plot_file_type=DEFAULT_FILE_TYPE,
                    dpi=DEFAULT_DPI,
                    figsize=DEFAULT_FIG_SIZE,
                    rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
                    ):

    plt.rcParams.update(rc_params)

    fig = plt.gcf()
    fig.set_size_inches(figsize)
    bars = plt.bar(x_vals, sizes, width=width, align=align)
    ax = plt.gca()
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    plt.xlim(left=1)
    plt.xlabel(xlabel)
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.ylabel(ylabel)
    if title is not None:
        plt.title(title, pad=TITLE_PAD)
    if not plot_file_name.endswith(plot_file_type):
        plot_file_name += '.' + plot_file_type
    plt.savefig(plot_file_name, dpi=dpi)
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()

    # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
    # plt.legend(patches, labels, loc="best")
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()
