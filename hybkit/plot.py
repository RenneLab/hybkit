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

# Import module-level dunder-names:
from hybkit.__about__ import (__author__, __contact__, __credits__, __date__, __deprecated__,
                              __email__, __license__, __maintainer__, __status__, __version__)

# ----- Begin Plot Constants -----
GENERAL_DEFAULTS = {
    'DPI': 1200,
    'FILE_TYPE': 'png',
    'FIG_SIZE': matplotlib.rcParams['figure.figsize'],  # 6.4, 4.8 in
    'TITLE_PAD': 15,
}


PIE_DEFAULTS = {
    'SETTINGS':  {
        'autopct': '%1.1f%%',
        'shadow': False,
        'startangle': 90,
        'counterclock': False,
    },
    'RC_PARAMS': {
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold',
    },
    'OTHER_THRESHHOLD': 0.1,
    'MIN_WEDGE_SIZE': 0.04,
}
# PIE['DEFAULT_TARGET_SETTINGS'] = copy.deepcopy(PIE['DEFAULT_SETTINGS'])
# PIE['DEFAULT_TARGET_SETTINGS'].update({'textprops': {'size': 'small'}})

LINE_DEFAULTS = {
    'RC_PARAMS': {
        'axes.titlesize': 'x-large',
        'axes.titleweight': 'bold',
        'axes.labelsize': 'large',
        'axes.labelweight': 'bold',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'large',
    },
    'LINE_DATA_FORMAT': '-',
    'LINE_MIN_FRACTION_SIZE': 0.01,
}

BAR_DEFAULTS = {
    'RC_PARAMS': {
        'axes.titlesize': 'x-large',
        # 'axes.titleweight': 'bold',
        # 'axes.labelsize': 'large',
        'axes.labelweight': 'bold',
        # 'xtick.labelsize': 'large',
        # 'ytick.labelsize': 'large',
    },
    'BAR_WIDTH': 0.9,
    'BAR_ALIGN': 'edge',
}

ENERGY_DEFAULTS = {
    'MIN_COUNT': 0,
    'MIN_DENSITY': 0.0,
    'XLABEL': 'Hybrid Gibbs Free Energy (kcal/mol)',
    'YLABEL': 'Hybrid Count',
    'COUNTS_TITLE': 'Hybrid Predicted Energies',
    'PROPS_TITLE': 'Normalized Predicted Energies',
}

MATCH_DEFAULTS = {
    'MATCH_MIN_COUNT': 0,
    'MATCH_XLABEL': 'Predicted Base-Pairs',
    'MATCH_YLABEL': 'Hybrid Count',
}

TARGET_DEFAULTS = {
    'TARGET_TITLE': 'Targets',
    'TARGET_TYPE_TITLE': 'Target Types',
}

FOLD_DEFAULTS = {
    'ENERGY_DENSITIES_TITLE': 'Normalized Hybrid Predicted Energies',
    'MATCH_TITLE': 'Hybrid Predicted Base Pairs',
}
# DEFAULT_HYBRID_TYPE_COUNTS_TITLE = 'Hybrid Types'
# DEFAULT_MIRNA_COUNTS_TITLE = 'Hybrid Types'
# DEFAULT_PATTERN_TITLE = 'Bound miRNA by Position'

FORMAT_NAME_MAP = {
    'mirnas_5p': "5'_miRNA_Hybrids",
    'mirnas_3p': "3'_miRNA_Hybrids",
    'mirna_dimers': 'miRNA_Duplexes',
    'non_mirnas': 'Non-miRNA_Hybrids'
}


# ----- Begin Plotting Methods -----
# Public Methods : Energy : energy_histogram
def energy_histogram(results, plot_file_name, name=None):
    """
    Plot histogram of hybrid energy counts from a `~hybkit.analysis.Analysis` Fold Analysis.

    Args:
        results (dict): Dictionary of energy counts from a `~hybkit.analysis.Analysis`
            Fold Analysis.
        plot_file_name (str): Name of output file.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
    """
    # Check if empty
    max_count = 0
    for energy, count in results.items():
        max_count = max(max_count, count)
    if not max_count:
        message = 'Warning: Attempted to create empty plot to name: %s' % plot_file_name
        print(message % plot_file_name)
        raise RuntimeError(message)

    # Plot energy bins
    energies, counts = [], []
    for energy, count in results.energy_bins.items():
        # if float(energy) > ENERGY_DEFAULTS['ENERGY_MIN_DENSITY']:
        #     break
        energies.append(float(energy))
        counts.append(count)

    title = ENERGY_DEFAULTS['COUNTS_TITLE']
    if name is not None:
        title = name + ': ' + title

    energy_width = float(energies[1]) + 0.02

    plot_params = {
        'x_vals': energies,
        'y_vals': counts,
        'plot_file_name': plot_file_name,
        'title': title,
        'xlabel': ENERGY_DEFAULTS['XLABEL'],
        'ylabel': ENERGY_DEFAULTS['YLABEL'],
        'width': energy_width,
        'align': BAR_DEFAULTS['BAR_ALIGN'],
        'rc_params': copy.deepcopy(BAR_DEFAULTS['RC_PARAMS']),
        'dpi': GENERAL_DEFAULTS['DPI'],
        'plot_file_type': GENERAL_DEFAULTS['FILE_TYPE'],
        'figsize': GENERAL_DEFAULTS['FIG_SIZE'],
        'title_pad': GENERAL_DEFAULTS['TITLE_PAD'],
    }
    _plot_energy_histogram(plot_params)


# Public Methods : Type : type_count
# Plot a pie plot from a hybkit Type Analysis
def type_count(results, plot_file_name, title, name=None):
    """
    Plot pie chart of hybrid type counts from a `~hybkit.analysis.Analysis` Type Analysis.

    Args:
        results (~collections.Counter): Counter Object of type counts from a
            `~hybkit.analysis.Analysis` Type Analysis.
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
    """

    # Sort hybrid_type_counts counter object in descending order
    labels, counts = [], []
    for key, count in results.most_common():
        labels.append(key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    plot_params = {
        'labels': labels,
        'sizes': counts,
        'plot_file_name': plot_file_name,
        'title': title,
        'other_threshhold': PIE_DEFAULTS['OTHER_THRESHHOLD'],
        'min_wedge_size': PIE_DEFAULTS['MIN_WEDGE_SIZE'],
        'rc_params': copy.deepcopy(PIE_DEFAULTS['RC_PARAMS']),
        'dpi': GENERAL_DEFAULTS['DPI'],
        'plot_file_type': GENERAL_DEFAULTS['FILE_TYPE'],
        'figsize': GENERAL_DEFAULTS['FIG_SIZE'],
        'title_pad': GENERAL_DEFAULTS['TITLE_PAD'],
        'settings': copy.deepcopy(PIE_DEFAULTS['SETTINGS']),
    }

    _plot_pie_chart(plot_params)


# ----- Begin Private Plotting Methods -----
# Private Methods : Energy : _plot_energy_histogram
def _plot_energy_histogram(plot_params):
    plt.rcParams.update(plot_params['rc_params'])

    fig = plt.gcf()
    fig.set_size_inches(plot_params['figsize'])
    bars = plt.bar(
        plot_params['x_vals'],
        plot_params['y_vals'],
        width=plot_params['width'],
        align=plot_params['align']
    )
    ax = plt.gca()
    ax.invert_xaxis()
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    plt.xlim(left=0)
    plt.xticks(range(0, (int(plot_params['x_vals'][-1]) - 1), (-2)), rotation=-30)
    plt.xlabel(plot_params['xlabel'])
    # plt.yticks(range(10, 110, 10))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    plt.ylabel(plot_params['ylabel'])
    plt.title(title, pad=plot_params['title_pad'])
    plt.savefig(plot_file_name, dpi=plot_params['dpi'])
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()

# Private Methods : Energy : _plot_types_pie_chart
def _plot_types_pie_chart(plot_params):
    total_size = sum(plot_params['sizes'])
    if total_size < 0.00000001:
        message = 'Warning: Attempted to create empty plot to name:'
        message += ' %s' % plot_params['plot_file_name']
        print(message)
        raise RuntimeError(message)
    fraction_sizes = [(size / total_size) for size in plot_params['sizes']]
    use_labels = []
    use_sizes = []
    for i in range(len(plot_params['labels'])):
        total_fraction = sum(use_sizes) / total_size
        if fraction_sizes[i] < plot_params['min_wedge_size']:
            break
        elif (total_fraction > (1 - plot_params['other_threshhold'])
              and i != (len(plot_params['labels']) - 1)):
            break
        else:
            use_labels.append(plot_params['labels'][i])
            use_sizes.append(plot_params['sizes'][i])

    if len(use_labels) != len(plot_params['labels']):
        other_size = total_size - sum(use_sizes)
        use_labels.append('other')
        use_sizes.append(other_size)

    plt.rcParams.update(plot_params['rc_params'])

    fig = plt.gcf()
    fig.set_size_inches(plot_params['figsize'])
    patches, texts, autotexts = plt.pie(
        use_sizes,
        labels=use_labels,
        **plot_params['settings']
    )
    plt.axis('equal')
    plt.title(plot_params['title'], pad=plot_params['title_pad'])
    plt.savefig(plot_params['plot_file_name'], dpi=plot_params['dpi'])
    plt.clf()
    plt.close()

# TODO: Implement fold_match_counts_histogram
# TODO: Implement fold_mirna_nt_counts_histogram

# ----- Begin Old Methods

# # Public Methods : HybRecord Type Analysis Plotting
# def mirna(analysis,
#           plot_file_name,
#           title=DEFAULT_MIRNA_COUNTS_TITLE,
#           other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
#           min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
#           plot_file_type=DEFAULT_FILE_TYPE,
#           dpi=DEFAULT_DPI,
#           rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
#           plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_SETTINGS)
#           ):
#     """
#     Plot the results of the :class:`~hybkit.analysis.MirnaAnalysis`.

#     Args:
#         analysis (Analysis): Analysis that includes attributes of
#             :class:`~hybkit.analysis.MirnaAnalysis`.
#         plot_file_name (str): File name for output plot.
#         title (str, optional): Title / header for plot.
#         other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
#             Setting to 0.0 disables the "other" wedge based on a threshhold.
#         min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
#             wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#         plot_settings (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.plot() function.
#     """
#     # Sort hybrid_type_counts counter object in descending order
#     labels = []
#     counts = []
#     use_items = Counter()
#     for key in ['mirnas_5p', 'mirnas_3p', 'mirna_dimers', 'non_mirnas']:
#         use_items[key] = getattr(analysis, key)

#     for key, count in use_items.most_common():
#         if key in FORMAT_NAME_MAP:
#             labels.append(FORMAT_NAME_MAP[key])
#         else:
#             labels.append(key)
#         counts.append(count)

#     if analysis.name is not None:
#         title = str(analysis.name) + ': ' + title

#     _plot_pie_chart(labels=labels,
#                     sizes=counts,
#                     plot_file_name=plot_file_name,
#                     title=title,
#                     other_threshhold=(-1),
#                     min_wedge_size=(-1),
#                     plot_file_type=plot_file_type,
#                     dpi=dpi,
#                     rc_params=rc_params,
#                     plot_settings=plot_settings,
#                     )


# # Public Methods : HybRecord miRNA Target Analysis Plotting
# def target(mirna_counts, plot_file_name,
#            name=None,
#            title=DEFAULT_TARGET_TITLE,
#            other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
#            min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
#            plot_file_type=DEFAULT_FILE_TYPE,
#            dpi=DEFAULT_DPI,
#            rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
#            plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_TARGET_SETTINGS)):
#     """
#     Plot the targets of a single mirna from the :class:`~hybkit.analysis.MirnaAnalysis`.

#     Args:
#         mirna_counts (~collection.Counter): Counter of targets of a mirna produced by
#             :class:`~hybkit.analysis.TargetAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str, optional): Name of analysis
#         title (str, optional): Title / header for plot.
#         other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
#             Setting to 0.0 disables the "other" wedge based on a threshhold.
#         min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
#             wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#         plot_settings (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.plot() function.
#     """
#     labels = []
#     counts = []
#     for ((target, target_seg_type), count) in mirna_counts.most_common():
#         labels.append(target)
#         counts.append(count)

#     if name is not None:
#         title = str(name) + ': ' + title

#     # Set figsize:
#     # Width set at constant 9(inches)
#     # Height 3 for >= 45 character labels, ratio of 9:3 (more rectangular)
#     # Height 3 + (0.16 * count-20) for: 20 < count < 45 labels
#     # Height 6.5 for <= 20 character labels, ratio of 9:6.5 (more square)
#     longest = max([0] + [len(label) for label in labels])
#     add_count = max((longest - 20), 0)
#     max_height = 6.5
#     min_height = 3.0
#     height = max((max_height - (0.16 * add_count)), (3.0))  # Ensure no less than min_height
#     height = min(height, max_height)                        # Ensure no more than max_height
#     figsize = (9, height)

#     _plot_pie_chart(labels=labels,
#                     sizes=counts,
#                     plot_file_name=plot_file_name,
#                     title=title,
#                     other_threshhold=other_threshhold,
#                     min_wedge_size=(0.025),
#                     plot_file_type=plot_file_type,
#                     dpi=dpi,
#                     figsize=figsize,
#                     rc_params=rc_params,
#                     plot_settings=plot_settings,
#                     )


# # Public Methods : HybRecord miRNA Target Analysis Plotting
# def target_type(target_type_count, plot_file_name,
#                 name=None,
#                 title=DEFAULT_TARGET_TYPE_TITLE,
#                 other_threshhold=DEFAULT_PIE_OTHER_THRESHHOLD,
#                 min_wedge_size=DEFAULT_PIE_MIN_WEDGE_SIZE,
#                 plot_file_type=DEFAULT_FILE_TYPE,
#                 dpi=DEFAULT_DPI,
#                 rc_params=copy.deepcopy(DEFAULT_PIE_RC_PARAMS),
#                 plot_settings=copy.deepcopy(DEFAULT_PIE_PLOT_TARGET_SETTINGS)):
#     """
#     Plot the targets types of a single mirna from the :class:`~hybkit.analysis.TargetAnalysis`.

#     Args:
#         target_type_count (~collections.Counter): Counter with types of targets for miRNAs from
#             :class:`~hybkit.analysis.TargetAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str): Name of analysis to plot for title.
#         title (str, optional): Title / header for plot (replaces default title).
#         other_threshhold (float, optional): Total fraction at which to begin the "other" wedge.
#             Setting to 0.0 disables the "other" wedge based on a threshhold.
#         min_wedge_size (float, optional): Minimum wedge fraction at which to add to the "other"
#             wedge. Setting to 0.0 disables the "other" wedge based on minimum wedge size.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#         plot_settings (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.plot() function.
#     """
#     labels = []
#     counts = []
#     for key, count in target_type_count.most_common():
#         labels.append(key)
#         counts.append(count)

#     if name is not None:
#         title = str(name) + ': ' + title

#     _plot_pie_chart(labels=labels,
#                     sizes=counts,
#                     plot_file_name=plot_file_name,
#                     title=title,
#                     other_threshhold=other_threshhold,
#                     # min_wedge_size=(0.025),
#                     plot_file_type=plot_file_type,
#                     dpi=dpi,
#                     rc_params=rc_params,
#                     plot_settings=plot_settings,
#                     )


# # Public Methods : HybRecord miRNA Binding Pattern Plotting
# def pattern(analysis,
#             plot_file_name,
#             name=None,
#             title=DEFAULT_PATTERN_TITLE,
#             data_format=DEFAULT_LINE_DATA_FORMAT,
#             min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
#             plot_file_type=DEFAULT_FILE_TYPE,
#             dpi=DEFAULT_DPI,
#             rc_params=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)
#             ):
#     """
#     Plot the bound percentage of mirna by base from the :class:`~hybkit.analysis.PatternAnalysis`.

#     Args:
#         analysis (~hybkit.analysis.PatternAnalysis): Analysis that includes attributes of
#             :class:`~hybkit.analysis.PatternAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str, optional): name to prepend to title.
#         title (str, optional): Title / header for plot (replaces default title).
#         data_format (str, optional): matplotlib line/data format.
#         min_fraction_size (float, optional): Minimum fraction to include at tail
#             end of plot. Setting to 0 includes all bases evaluated.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#     """
#     labels = []
#     fractions = []
#     max_fraction = 0.0
#     for seq_i, fraction in analysis.mirna_fold_frac.items():
#         max_fraction = max(max_fraction, fraction)

#     if max_fraction < 0.00000000000001:
#         print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
#         return

#     for seq_i, fraction in analysis.mirna_fold_frac.items():
#         if seq_i > 22 and (fraction / max_fraction) < min_fraction_size:
#             break
#         labels.append(seq_i)
#         fractions.append(fraction)

#     if analysis.name is not None:
#         title = str(analysis.name) + ': ' + title

#     _plot_line(labels=labels,
#                sizes=fractions,
#                plot_file_name=plot_file_name,
#                title=title,
#                data_format=data_format,
#                min_fraction_size=min_fraction_size,
#                plot_file_type=plot_file_type,
#                dpi=dpi,
#                rc_params=rc_params,
#                )


# # Public Methods : HybRecord Fold Analysis Plotting
# def fold_energy_counts(analysis,
#                        plot_file_name,
#                        name=None,
#                        title=DEFAULT_FOLD_ENERGY_COUNTS_TITLE,
#                        min_count=DEFAULT_ENERGY_MIN_COUNT,
#                        plot_file_type=DEFAULT_FILE_TYPE,
#                        dpi=DEFAULT_DPI,
#                        rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
#                        ):
#     """
#     Plot barogram of hybrid energy counts from a `~hybkit.analysis.FoldAnalysis`.

#     Args:
#         analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
#             :class:`~hybkit.analysis.FoldAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str, optional): name to prepend to title.
#         title (str, optional): Title / header for plot (replaces default title).
#         min_count (float, optional): Minimum count to include at right
#             side of plot of energy counts. Setting to 0 includes all values evaluated.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#     """
#     # Check if empty
#     max_count = 0
#     for energy, count in analysis.energy_bins.items():
#         max_count = max(max_count, count)

#     if max_count < 1:
#         print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
#         return

#     # Plot energy bins
#     energies, counts = [], []
#     for energy, count in analysis.energy_bins.items():
#         if float(energy) < -5 and count < min_count:
#             break
#         energies.append(float(energy))
#         counts.append(count)

#     if analysis.name is not None:
#         title = str(analysis.name) + ': ' + title

#     energy_width = float(energies[1]) + 0.02

#     _plot_energy_bar(x_vals=energies,
#                      sizes=counts,
#                      plot_file_name=plot_file_name,
#                      title=title,
#                      width=energy_width,
#                      align=DEFAULT_BAR_ALIGN,
#                      plot_file_type=plot_file_type,
#                      dpi=dpi,
#                      rc_params=rc_params,
#                      )


# # Public Methods : HybRecord Fold Analysis Plotting
# def fold_energy_densities(analysis,
#                           plot_file_name,
#                           name=None,
#                           title=DEFAULT_FOLD_ENERGY_DENSITIES_TITLE,
#                           min_density=DEFAULT_ENERGY_MIN_DENSITY,
#                           plot_file_type=DEFAULT_FILE_TYPE,
#                           dpi=DEFAULT_DPI,
#                           rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
#                           ):
#     """
#     Plot barogram of hybrid energy counts from a `~hybkit.analysis.FoldAnalysis`.

#     Args:
#         analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
#             :class:`~hybkit.analysis.FoldAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str, optional): name to prepend to title.
#         title (str, optional): Title / header for plot (replaces default title).
#         min_density (float, optional): Minimum density to include at right
#             side of plot of energy counts. Setting to 0 includes all values evaluated.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#     """
#     # Check if empty
#     max_density = 0
#     for energy, density in analysis.energy_bin_densities.items():
#         max_density = max(max_density, density)

#     if max_density < 0.0000001:
#         print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
#         return

#     # Plot energy bins
#     energies, densities = [], []
#     for energy, density in analysis.energy_bin_densities.items():
#         if float(energy) < -5 and density < min_density:
#             break
#         energies.append(float(energy))
#         densities.append(density)

#     if analysis.name is not None:
#         title = str(analysis.name) + ': ' + title

#     energy_width = float(energies[1]) + 0.02

#     _plot_energy_bar(x_vals=energies,
#                      sizes=densities,
#                      plot_file_name=plot_file_name,
#                      title=title,
#                      width=energy_width,
#                      align=DEFAULT_BAR_ALIGN,
#                      plot_file_type=plot_file_type,
#                      dpi=dpi,
#                      rc_params=rc_params,
#                      )


# # Public Methods : HybRecord Fold Analysis Plotting
# def fold_matches(analysis,
#                  plot_file_name,
#                  name=None,
#                  title=DEFAULT_FOLD_MATCH_TITLE,
#                  min_count=DEFAULT_MATCH_MIN_COUNT,
#                  plot_file_type=DEFAULT_FILE_TYPE,
#                  dpi=DEFAULT_DPI,
#                  rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
#                  ):
#     """
#     Plot predicted hybrid match counts from a `~hybkit.analysis.FoldAnalysis`.

#     Args:
#         analysis (~hybkit.analysis.FoldAnalysis): Analysis that includes attributes of
#             :class:`~hybkit.analysis.FoldAnalysis`.
#         plot_file_name (str): File name for output plot.
#         name (str, optional): name to prepend to title.
#         title (str, optional): Title / header for plot (replaces default title).
#         min_count (float, optional): Minimum count to include at right
#             side of plot of energy counts. Setting to 0 includes all values evaluated.
#         plot_file_type (str, optional): File type for saving of plots. Options:
#             {'png', 'ps', 'pdf', 'svg'}
#         dpi (int, optional): DPI for saving of plots.
#         rc_params (dict, optional): Dict of keys and values of settings to pass
#             to the matplotlib.rcParams.update() method.
#     """
#     # Check if empty
#     max_count = 0
#     for base, count in analysis.match_counts.items():
#         max_count = max(max_count, count)

#     if max_count < 1:
#         print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
#         return

#     # Plot match count bins
#     total_count = sum(analysis.match_counts.values())
#     bases, counts, densities = [], [], []
#     for base, count in analysis.match_counts.items():
#         if base > 10 and count < min_count:
#             break
#         bases.append(base)
#         counts.append(count)
#         densities.append(count/total_count)

#     density_title = 'Normalized ' + title
#     if analysis.name is not None:
#         title = str(analysis.name) + ': ' + title
#         density_title = str(analysis.name) + ': ' + density_title

#     counts_plot_file_name = plot_file_name + '_counts'
#     densities_plot_file_name = plot_file_name + '_densities'

#     _plot_match_bar(x_vals=bases,
#                     sizes=counts,
#                     plot_file_name=counts_plot_file_name,
#                     title=title,
#                     align='center',
#                     plot_file_type=plot_file_type,
#                     dpi=dpi,
#                     rc_params=rc_params,
#                     )

#     _plot_match_bar(x_vals=bases,
#                     sizes=densities,
#                     plot_file_name=densities_plot_file_name,
#                     title=density_title,
#                     ylabel='Normalized ' + DEFAULT_MATCH_YLABEL,
#                     align='center',
#                     plot_file_type=plot_file_type,
#                     dpi=dpi,
#                     rc_params=rc_params,
#                     )



# # Private Methods : Plot Line Graph
# def _plot_line(labels, sizes, plot_file_name,
#                title=None,
#                min_fraction_size=DEFAULT_LINE_MIN_FRACTION_SIZE,
#                plot_file_type=DEFAULT_FILE_TYPE,
#                data_format=DEFAULT_LINE_DATA_FORMAT,
#                dpi=DEFAULT_DPI,
#                figsize=DEFAULT_FIG_SIZE,
#                rc_params=copy.deepcopy(DEFAULT_LINE_RC_PARAMS)
#                ):
#     max_size = max(sizes)
#     if abs(max_size) < 0.00000000000001:
#         print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
#         return
#     fraction_sizes = [(size / max_size) for size in sizes]
#     use_labels = []
#     use_sizes = []
#     for i in range(len(labels)):
#         if fraction_sizes[i] > min_fraction_size:
#             use_labels.append(labels[i])
#             use_sizes.append(sizes[i] * 100)

#     plt.rcParams.update(rc_params)

#     fig = plt.gcf()
#     fig.set_size_inches(figsize)
#     lines = plt.plot(use_labels, use_sizes, data_format)
#     ax = plt.gca()
#     ax.yaxis.set_major_formatter(mtick.PercentFormatter())
#     plt.xlim(left=1)
#     plt.xticks(range(1, (len(labels) + 1), 2))
#     plt.xlabel('miRNA Base Index')
#     plt.yticks(range(10, 110, 10))
#     plt.ylabel('Percent Bound')
#     if title is not None:
#         plt.title(title, pad=TITLE_PAD)
#     if not plot_file_name.endswith(plot_file_type):
#         plot_file_name += '.' + plot_file_type
#     plt.savefig(plot_file_name, dpi=dpi)
#     plt.clf()
#     plt.close()
#     matplotlib.rcdefaults()

#     # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
#     # plt.legend(patches, labels, loc="best")
#     # plt.axis('equal')
#     # plt.tight_layout()
#     # plt.show()

# # Private Methods : Plot Energy Bar Graph
# def _plot_energy_bar(x_vals, sizes, plot_file_name,
#                      title=None,
#                      width=DEFAULT_BAR_WIDTH,
#                      align=DEFAULT_BAR_ALIGN,
#                      xlabel=DEFAULT_ENERGY_XLABEL,
#                      ylabel=DEFAULT_ENERGY_YLABEL,
#                      plot_file_type=DEFAULT_FILE_TYPE,
#                      dpi=DEFAULT_DPI,
#                      figsize=DEFAULT_FIG_SIZE,
#                      rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
#                      ):

#     plt.rcParams.update(rc_params)

#     fig = plt.gcf()
#     fig.set_size_inches(figsize)
#     bars = plt.bar(x_vals, sizes, width=width, align=align)
#     ax = plt.gca()
#     ax.invert_xaxis()
#     ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
#     plt.xlim(left=0)
#     plt.xticks(range(0, (int(x_vals[-1])-1), (-2)), rotation=-30)
#     plt.xlabel(xlabel)
#     #plt.yticks(range(10, 110, 10))
#     ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
#     plt.ylabel(ylabel)
#     if title is not None:
#         plt.title(title, pad=TITLE_PAD)
#     if not plot_file_name.endswith(plot_file_type):
#         plot_file_name += '.' + plot_file_type
#     plt.savefig(plot_file_name, dpi=dpi)
#     plt.clf()
#     plt.close()
#     matplotlib.rcdefaults()

#     # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
#     # plt.legend(patches, labels, loc="best")
#     # plt.axis('equal')
#     # plt.tight_layout()
#     # plt.show()

# # Private Methods : Plot Bases Bar Graph
# def _plot_match_bar(x_vals, sizes, plot_file_name,
#                     title=None,
#                     width=DEFAULT_BAR_WIDTH,
#                     align=DEFAULT_BAR_ALIGN,
#                     xlabel=DEFAULT_MATCH_XLABEL,
#                     ylabel=DEFAULT_MATCH_YLABEL,
#                     plot_file_type=DEFAULT_FILE_TYPE,
#                     dpi=DEFAULT_DPI,
#                     figsize=DEFAULT_FIG_SIZE,
#                     rc_params=copy.deepcopy(DEFAULT_BAR_RC_PARAMS)
#                     ):

#     plt.rcParams.update(rc_params)

#     fig = plt.gcf()
#     fig.set_size_inches(figsize)
#     bars = plt.bar(x_vals, sizes, width=width, align=align)
#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
#     plt.xlim(left=1)
#     plt.xlabel(xlabel)
#     ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
#     plt.ylabel(ylabel)
#     if title is not None:
#         plt.title(title, pad=TITLE_PAD)
#     if not plot_file_name.endswith(plot_file_type):
#         plot_file_name += '.' + plot_file_type
#     plt.savefig(plot_file_name, dpi=dpi)
#     plt.clf()
#     plt.close()
#     matplotlib.rcdefaults()

#     # patches, texts = plt.pie(sizes, colors=colors, shadow=True, startangle=90)
#     # plt.legend(patches, labels, loc="best")
#     # plt.axis('equal')
#     # plt.tight_layout()
#     # plt.show()
