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
RC_PARAMS_DEFAULTS = {
    'figure.dpi': 1200,
    'figure.figsize': matplotlib.rcParams['figure.figsize'],  # 6.4, 4.8 in
    'axes.titlepad': 15,
}

PIE_RC_PARAMS = copy.deepcopy(RC_PARAMS_DEFAULTS)
PIE_RC_PARAMS.update({
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    # 'figure.figsize': (10, 4.8),
}),
PIE_DEFAULTS = {
    'SETTINGS': {
        'autopct': '%1.1f%%',
        'shadow': False,
        'startangle': 90,
        'counterclock': False,
    },
    'RC_PARAMS': PIE_RC_PARAMS,
    'OTHER_THRESHHOLD': 0.1,
    'MIN_WEDGE_SIZE': 0.04,
}
# PIE['DEFAULT_TARGET_SETTINGS'] = copy.deepcopy(PIE['DEFAULT_SETTINGS'])
# PIE['DEFAULT_TARGET_SETTINGS'].update({'textprops': {'size': 'small'}})

LINE_RC_PARAMS = copy.deepcopy(RC_PARAMS_DEFAULTS)
LINE_RC_PARAMS.update({
    'axes.titlesize': 'x-large',
    'axes.titleweight': 'bold',
    'axes.labelsize': 'large',
    'axes.labelweight': 'bold',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'large',
})
LINE_DEFAULTS = {
    'RC_PARAMS': LINE_RC_PARAMS,
    'LINE_DATA_FORMAT': '-',
    'LINE_MIN_FRACTION_SIZE': 0.01,
}

BAR_RC_PARAMS = copy.deepcopy(RC_PARAMS_DEFAULTS)
BAR_RC_PARAMS.update({
    'axes.titlesize': 'x-large',
    'axes.titleweight': 'bold',
    # 'axes.labelsize': 'large',
    'axes.labelweight': 'bold',
    # 'xtick.labelsize': 'large',
    # 'ytick.labelsize': 'large',
}),
BAR_DEFAULTS = {
    'RC_PARAMS': BAR_RC_PARAMS,
    'BAR_WIDTH': 0.9,
    'BAR_ALIGN': 'edge',
}

ENERGY_HIST_RC_PARAMS = copy.deepcopy(BAR_RC_PARAMS)
ENERGY_HIST_RC_PARAMS.update({})
ENERGY_DEFAULTS = {
    'MIN_COUNT': 0,
    'MIN_DENSITY': 0.0,
    'XLABEL': 'Hybrid Gibbs Free Energy (kcal/mol)',
    'YLABEL': 'Hybrid Count',
    'RC_PARAMS': ENERGY_HIST_RC_PARAMS,
}

TYPE_PIE_RC_PARAMS = copy.deepcopy(PIE_RC_PARAMS)
TYPE_PIE_RC_PARAMS.update({})
# MATCH_DEFAULTS = {
#     'MATCH_MIN_COUNT': 0,
#     'MATCH_XLABEL': 'Predicted Base-Pairs',
#     'MATCH_YLABEL': 'Hybrid Count',
# }

TARGET_PIE_RC_PARAMS = copy.deepcopy(PIE_RC_PARAMS)
TARGET_PIE_RC_PARAMS.update({
    'figure.figsize': (10, 4.8),
})
# TARGET_DEFAULTS = {
#     'TARGET_TITLE': 'Targets',
#     'TARGET_TYPE_TITLE': 'Target Types',
# }

FOLD_MATCH_HIST_RC_PARAMS = copy.deepcopy(BAR_RC_PARAMS)
FOLD_MATCH_HIST_RC_PARAMS.update({})
# FOLD_DEFAULTS = {
#     'ENERGY_DENSITIES_TITLE': 'Normalized Hybrid Predicted Energies',
#     'MATCH_TITLE': 'Hybrid Predicted Base Pairs',
# }

FOLD_NT_COUNTS_HIST_RC_PARAMS = copy.deepcopy(BAR_RC_PARAMS)
FOLD_NT_COUNTS_HIST_RC_PARAMS.update({})

FORMAT_NAME_MAP = {
    'mirnas_5p': "5'_miRNA_Hybrids",
    'mirnas_3p': "3'_miRNA_Hybrids",
    'mirna_dimers': 'miRNA_Duplexes',
    'non_mirnas': 'Non-miRNA_Hybrids'
}


# ----- Begin Plotting Methods -----
# Public Methods : Energy : energy_histogram
def energy_histogram(results,
                     plot_file_name,
                     title,
                     name=None,
                     rc_params=copy.deepcopy(ENERGY_HIST_RC_PARAMS),
                     ):
    """
    Plot histogram of hybrid energy counts from a `~hybkit.analysis.Analysis` Fold Analysis.

    Args:
        results (dict): Dictionary of energy counts from a `~hybkit.analysis.Analysis`
            Fold Analysis (Key: 'binned_energy_vals').
        plot_file_name (str): Name of output file.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        rc_params (:obj:`dict`, optional): Dictionary of matplotlib rcParams. Defaults
            to :obj:`~hybkit.plot.ENERGY_HIST_RC_PARAMS`.
    """
    # Check if empty
    max_count = 0
    for energy, count in results.items():
        max_count = max(max_count, count)
    if not max_count:
        message = 'Warning: Attempted to create empty plot to name: %s' % plot_file_name
        print(message)
        raise RuntimeError(message)

    # Plot energy bins
    energies, counts = [], []
    for energy, count in results.items():
        # if float(energy) > ENERGY_DEFAULTS['ENERGY_MIN_DENSITY']:
        #     break
        energies.append(float(energy))
        counts.append(count)

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
        'rc_params': rc_params,
    }
    _plot_energy_histogram(plot_params)


# Public Methods : Type : type_count
# Plot a pie plot from a hybkit Type Analysis
def type_count(results,
               plot_file_name,
               title,
               name=None,
               join_entries=False,
               rc_params=copy.deepcopy(TYPE_PIE_RC_PARAMS),
               ):
    """
    Plot pie chart of hybrid type counts from a `~hybkit.analysis.Analysis` Type Analysis.

    Args:
        results (~collections.Counter): Counter Object of type counts from a
            `~hybkit.analysis.Analysis` Type Analysis.
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        join_entries (:obj:`bool`, optional): If True, join two-tuple pairs into a
            single string for plot labels.
        rc_params (:obj:`dict`, optional): Dictionary of matplotlib rcParams. Defaults
            to :obj:`~hybkit.plot.TYPE_PIE_RC_PARAMS`.
    """
    # Sort hybrid_type_counts counter object in descending order
    labels, counts = [], []
    for key, count in results.most_common():
        if join_entries:
            key = '--'.join(key)
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
        'rc_params': rc_params,
        'settings': copy.deepcopy(PIE_DEFAULTS['SETTINGS']),
    }

    _plot_types_pie_chart(plot_params)


def target_count(*args, **kwargs):
    """Hold Place for replaced docstring."""
    return type_count(*args, rc_params=TARGET_PIE_RC_PARAMS, **kwargs)


# Modify target_count doc
target_count.__doc__ = type_count.__doc__
for q, r in [('hybrid type counts', 'target counts'),
             ('TYPE_PIE_RC_PARAMS', 'TARGET_PIE_RC_PARAMS')]:
    target_count.__doc__ = target_count.__doc__.replace(q, r)


def fold_match_counts_histogram(results,
                                plot_file_name,
                                title,
                                name=None,
                                join_entries=False,
                                rc_params=copy.deepcopy(FOLD_MATCH_HIST_RC_PARAMS),
                                ):
    """
    Plot histogram of predicted miRNA/target match count from `~hybkit.analysis.Analysis` Analysis.

    Args:
        results (~collections.Counter): Counter Object of match counts from a
            `~hybkit.analysis.Analysis` Type Analysis.
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        rc_params (:obj:`dict`, optional): Dictionary of matplotlib rcParams. Defaults
            to :obj:`~hybkit.plot.FOLD_MATCH_HIST_RC_PARAMS`.
    """
    x_vals = []
    y_vals = []
    min_x = 10000
    max_x = -10000
    max_val = -10000
    for x_val, y_val in results.items():
        min_x = min(min_x, x_val)
        max_x = max(max_x, x_val)
        max_val = max(max_val, y_val)

    if max_val == -10000:
        print('Warning: Attempted to create empty plot to name: %s' % plot_file_name)
        return

    for x_val in range(min_x, (max_x + 1)):
        if x_val in results:
            y_val = results[x_val]
        else:
            y_val = 0
        # if nt_i > 22 and (y_val / max_fraction) < min_fraction_size:
        #     break
        x_vals.append(x_val)
        y_vals.append(y_val)

    if name is not None:
        title = str(name) + ': ' + title

    plot_params = {
        'x_vals': x_vals,
        'y_vals': y_vals,
        'plot_file_name': plot_file_name,
        'title': title,
        'xlabel': 'Predicted miRNA/Target Matches',
        'ylabel': 'Hybrid Count',
        # 'width': 1,
        # 'align': BAR_DEFAULTS['BAR_ALIGN'],
        'rc_params': rc_params,
    }

    _plot_int_hist(plot_params=plot_params)


def fold_mirna_nt_counts_histogram(*args, **kwargs):
    """Hold Place for replaced docstring."""
    return fold_match_counts_histogram(*args, rc_params=FOLD_NT_COUNTS_HIST_RC_PARAMS, **kwargs)


# Modify target_count doc
fold_mirna_nt_counts_histogram.__doc__ = fold_match_counts_histogram.__doc__
for q, r in [('predicted miRNA/target match count ', 'predicted binding at nt position'),
             ('FOLD_MATCH_HIST_RC_PARAMS', 'FOLD_NT_COUNTS_HIST_RC_PARAMS')]:
    fold_mirna_nt_counts_histogram.__doc__ = fold_mirna_nt_counts_histogram.__doc__.replace(q, r)


# ----- Begin Private Plotting Methods -----
# Private Methods : Energy : _plot_energy_histogram
def _plot_energy_histogram(plot_params):
    # Update plot parameters
    plt.rcParams.update(plot_params['rc_params'])

    # Create the figure with the specified size
    fig, ax = plt.subplots()

    # Create the bars for the histogram
    bars = ax.bar(
        plot_params['x_vals'],
        plot_params['y_vals'],
        width=plot_params['width'],
        align=plot_params['align']
    )

    # Invert the x-axis and set tick locators
    ax.invert_xaxis()
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

    # Set x-axis limits and ticks
    plt.xlim(left=0)
    plt.xticks(range(0, (int(plot_params['x_vals'][-1]) - 1), (-2)), rotation=-30)

    # Set x-axis and y-axis labels
    plt.xlabel(plot_params['xlabel'])
    plt.ylabel(plot_params['ylabel'])

    # Set y-axis minor tick locators
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    # Set the title with custom padding
    plt.title(plot_params['title'])

    # Add space for the x-axis title
    plt.subplots_adjust(bottom=0.15)

    # Adjust layout to ensure labels are not cut off
    plt.tight_layout()

    # Save the figure to a file with specified DPI
    plt.savefig(plot_params['plot_file_name'])

    # Clear the current figure and close it to release resources
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()


# Private Methods : Energy : _plot_types_pie_chart
def _plot_types_pie_chart(plot_params):
    # Calculate the total size and the size of each fraction
    total_size = sum(plot_params['sizes'])
    if total_size < 0.00000001:
        message = 'Warning: Attempted to create empty plot to name:'
        message += ' %s' % plot_params['plot_file_name']
        print(message)
        raise RuntimeError(message)

    fraction_sizes = [(size / total_size) for size in plot_params['sizes']]

    # Filter labels and sizes based on minimum wedge size and other threshold
    use_labels, use_sizes = [], []
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

    # Add the "other" category if necessary
    if len(use_labels) != len(plot_params['labels']):
        other_size = total_size - sum(use_sizes)
        use_labels.append('other')
        use_sizes.append(other_size)

    # Update plot parameters and create the figure with the specified size
    plt.rcParams.update(plot_params['rc_params'])

    # fig, ax = plt.subplots(figsize=plot_params['figsize'])
    fig, ax = plt.subplots()

    # Create the pie chart
    patches, texts, autotexts = ax.pie(
        use_sizes,
        labels=use_labels,
        **plot_params['settings']
    )

    # Set equal aspect ratio to ensure the pie chart is circular
    plt.axis('equal')

    # Set the title with custom padding
    plt.title(plot_params['title'])

    # Adjust layout to ensure labels are not cut off
    fig.draw_without_rendering()
    # plt.tight_layout(pad=1.5)
    plt.tight_layout(pad=1.5)

    # Save the figure to a file with specified DPI
    plt.savefig(plot_params['plot_file_name'])

    # Clear the current figure and close it to release resources
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()


# Private Methods : Energy : _plot_energy_histogram
def _plot_int_hist(plot_params):
    # Update plot parameters
    plt.rcParams.update(plot_params['rc_params'])

    # Create the figure with the specified size
    fig, ax = plt.subplots()

    # Create the bars for the histogram
    bars = ax.bar(
        plot_params['x_vals'],
        plot_params['y_vals'],
        # width=plot_params['width'],
        # align=plot_params['align']
    )

    # Invert the x-axis and set tick locators
    # ax.invert_xaxis()
    # ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))

    # Set x-axis limits and ticks
    # plt.xlim(left=0)
    # plt.xticks(range(0, (int(plot_params['x_vals'][-1]) - 1), (-2)), rotation=-30)

    # Set x-axis and y-axis labels
    plt.xlabel(plot_params['xlabel'])
    plt.ylabel(plot_params['ylabel'])

    # Set y-axis minor tick locators
    # ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    # Set the title with custom padding
    plt.title(plot_params['title'])

    # Add space for the x-axis title
    plt.subplots_adjust(bottom=0.15)

    # Adjust layout to ensure labels are not cut off
    plt.tight_layout()

    # Save the figure to a file with specified DPI
    plt.savefig(plot_params['plot_file_name'])

    # Clear the current figure and close it to release resources
    plt.clf()
    plt.close()
    matplotlib.rcdefaults()
