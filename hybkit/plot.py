#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Methods for plotting analyses of HybRecord and FoldRecord objects."""

import copy
import logging
from collections import Counter
from typing import Any, Callable, Dict, Generic, Optional, TypeVar

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

import hybkit
from hybkit.__about__ import (
    __author__,
    __contact__,
    __credits__,
    __date__,
    __deprecated__,
    __email__,
    __license__,
    __maintainer__,
    __status__,
    __version__,
)
from hybkit.errors import HybkitMiscError

# ----- File-Specific Linting Directives:
# ruff: noqa: F401 SLF001

# ----- Begin Typing Assist -----
# F = TypeVar('F', bound=Callable[..., Any])
# class copy_signature(Generic[F]):
#     def __init__(self, target: F) -> None: ...
#     def __call__(self, wrapped: Callable[..., Any]) -> F: ...

mpl.use('Agg')

# ----- Begin Plot Constants -----
#: Default Colors for colored plots:
#: Colors selected based on "Points of view: Color blindness" by Bang Wong, Nature Methods, 2011.
#: Colors in RGB nomenclature (1-255):
#: Black (0,0,0), Orange (230,159,0), Sky Blue (86,180,233), Bluish Green (0,158,115),
#: Yellow (240,228,66), Blue (0,114,178), Vermilion (213,94,0), Reddish Purple (204,121,167)
COLOR_DICT = ({
    'Blue': '#0072B2',
    'Vermilion': '#D55E00',
    'Bluish Green': '#009E73',
    'Reddish Purple': '#CC79A7',
    'Orange': '#E69F00',
    'Sky Blue': '#56B4E9',
    # 'Reddish Purple': '#CC79A7',
    # 'Orange': '#E69F00',
    'Yellow': '#F0E442',
    # 'Black': '#000000',
})

#: List of default colors for colored plots.
COLOR_LIST = list(COLOR_DICT.values())

# Base defaults for mpl rcParams.
RC_PARAMS_DEFAULTS = ({
    'figure.dpi': 1200,
    'figure.figsize': mpl.rcParams['figure.figsize'],  # 6.4, 4.8 in
    'axes.titlepad': 15,
})

# Default mpl rcParams for pie charts.
PIE_RC_PARAMS = dict(copy.deepcopy(RC_PARAMS_DEFAULTS))
PIE_RC_PARAMS.update({
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    # 'figure.figsize': (10, 4.8),
})
PIE_RC_PARAMS = (copy.deepcopy(PIE_RC_PARAMS))

# Default mpl rcParams for bar charts.
BAR_RC_PARAMS = dict(copy.deepcopy(RC_PARAMS_DEFAULTS))
BAR_RC_PARAMS.update({
    'axes.titlesize': 'large',
    # 'axes.titlesize': 'x-large',
    'axes.titleweight': 'bold',
    # 'axes.labelsize': 'large',
    'axes.labelweight': 'bold',
    # 'xtick.labelsize': 'large',
    # 'ytick.labelsize': 'large',
})
BAR_RC_PARAMS = (copy.deepcopy(BAR_RC_PARAMS))

#: Default mpl rcParams for energy analysis histograms.
ENERGY_HIST_RC_PARAMS = dict(copy.deepcopy(BAR_RC_PARAMS))
ENERGY_HIST_RC_PARAMS.update({})
ENERGY_HIST_RC_PARAMS = (copy.deepcopy(ENERGY_HIST_RC_PARAMS))

#: Default mpl rcParams for type analysis pie charts.
TYPE_PIE_SINGLE_RC_PARAMS = dict(copy.deepcopy(PIE_RC_PARAMS))
TYPE_PIE_SINGLE_RC_PARAMS.update({})
# MATCH_DEFAULTS = {
#     'MATCH_MIN_COUNT': 0,
#     'MATCH_XLABEL': 'Predicted Base-Pairs',
#     'MATCH_YLABEL': 'Hybrid Count',
# }
TYPE_PIE_SINGLE_RC_PARAMS = (copy.deepcopy(TYPE_PIE_SINGLE_RC_PARAMS))


#: Default mpl rcParams for type analysis pie charts.
TYPE_PIE_DUAL_RC_PARAMS = dict(copy.deepcopy(PIE_RC_PARAMS))
TYPE_PIE_DUAL_RC_PARAMS.update({
    'figure.figsize': (8, 4.8),
})
TYPE_PIE_DUAL_RC_PARAMS = (copy.deepcopy(TYPE_PIE_DUAL_RC_PARAMS))

#: Default mpl rcParams for target analysis pie charts.
TARGET_PIE_RC_PARAMS = dict(copy.deepcopy(PIE_RC_PARAMS))
TARGET_PIE_RC_PARAMS.update({
    'figure.figsize': (9.6, 4.8),
})
# TARGET_DEFAULTS = {
#     'TARGET_TITLE': 'Targets',
#     'TARGET_TYPE_TITLE': 'Target Types',
# }
TARGET_PIE_RC_PARAMS = (copy.deepcopy(TARGET_PIE_RC_PARAMS))

#: Default mpl rcParams for fold match analysis histograms.
FOLD_MATCH_HIST_RC_PARAMS = dict(copy.deepcopy(BAR_RC_PARAMS))
FOLD_MATCH_HIST_RC_PARAMS.update({})
# FOLD_DEFAULTS = {
#     'ENERGY_DENSITIES_TITLE': 'Normalized Hybrid Predicted Energies',
#     'MATCH_TITLE': 'Hybrid Predicted Base Pairs',
# }
FOLD_MATCH_HIST_RC_PARAMS = (copy.deepcopy(FOLD_MATCH_HIST_RC_PARAMS))

#: Default mpl rcParams for fold nt counts analysis histograms.
FOLD_NT_COUNTS_HIST_RC_PARAMS = dict(copy.deepcopy(BAR_RC_PARAMS))
FOLD_NT_COUNTS_HIST_RC_PARAMS.update({})
FOLD_NT_COUNTS_HIST_RC_PARAMS = (copy.deepcopy(FOLD_NT_COUNTS_HIST_RC_PARAMS))

# _FORMAT_NAME_MAP = {
#     'mirnas_5p': "5'_miRNA_Hybrids",
#     'mirnas_3p': "3'_miRNA_Hybrids",
#     'mirna_dimers': 'miRNA_Duplexes',
#     'non_mirnas': 'Non-miRNA_Hybrids'
# }


#: Default Pie Chart Plot Settings.
PIE_DEFAULTS = ({
    'SETTINGS': {
        'autopct': '%1.1f%%',
        'shadow': False,
        'startangle': 90,
        'counterclock': False,
    },
    'COLORS': COLOR_LIST,
    'OTHER_THRESHOLD': 0.05,
    'MIN_WEDGE_SIZE': 0.04,
})
# PIE['DEFAULT_TARGET_SETTINGS'] = copy.deepcopy(PIE['DEFAULT_SETTINGS'])
# PIE['DEFAULT_TARGET_SETTINGS'].update({'textprops': {'size': 'small'}})

#: Default Bar Chart Plot Settings.
BAR_DEFAULTS = ({
    'BAR_WIDTH': 0.9,
    'BAR_ALIGN': 'edge',
    'BAR_EDGE_COLOR': None,
})

#: Default Bar Chart of Integer Plot Settings.
BAR_INT_DEFAULTS = dict(copy.deepcopy(BAR_DEFAULTS))
BAR_INT_DEFAULTS.update({
    'BAR_ALIGN': 'center',
})
# BAR_INT_DEFAULTS = (BAR_INT_DEFAULTS)

#: Default Bar Chart Plot Settings for Energy Histograms.
ENERGY_DEFAULTS = ({
    'MIN_COUNT': 0,
    'MIN_DENSITY': 0.0,
    'XLABEL': 'Hybrid Gibbs Free Energy (kcal/mol)',
    'YLABEL': 'Hybrid Count',
})

# LINE_RC_PARAMS = copy.deepcopy(RC_PARAMS_DEFAULTS)
# LINE_RC_PARAMS.update({
#     'axes.titlesize': 'x-large',
#     'axes.titleweight': 'bold',
#     'axes.labelsize': 'large',
#     'axes.labelweight': 'bold',
#     'xtick.labelsize': 'medium',
#     'ytick.labelsize': 'large',
# })
# LINE_DEFAULTS = {
#     'RC_PARAMS': LINE_RC_PARAMS,
#     'LINE_DATA_FORMAT': '-',
#     'LINE_MIN_FRACTION_SIZE': 0.01,
# }


# ----- Begin Plotting Methods -----
# Public Methods : Energy : energy_histogram
def energy_histogram(
        results: Dict[str, Any],
        plot_file_name: str,
        title: str,
        name: Optional[str] = None,
        rc_params: Dict[str, Any] = ENERGY_HIST_RC_PARAMS,
        bar_params: Dict[str, Any] = BAR_DEFAULTS,
        ) -> None:
    """
    Plot histogram of hybrid energies from an :class:`~hybkit.analysis.Analysis` fold analysis.

    Args:
        results (dict): Dictionary of energy counts from an :class:`~hybkit.analysis.Analysis`
            fold analysis (Key: ``binned_energy_vals``).
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        rc_params (:obj:`dict`, optional): Dictionary of mpl rcParams. Defaults
            to :obj:`~hybkit.plot.ENERGY_HIST_RC_PARAMS`.
        bar_params (:obj:`dict`, optional): Dictionary of bar plot parameters. Defaults
            to :obj:`~hybkit.plot.BAR_DEFAULTS`.
    """
    # Check if empty
    max_count = 0
    for _energy, count in results.items():
        max_count = max(max_count, count)
    if not max_count:
        message = 'Warning: Attempted to create empty plot to name: %s' % plot_file_name
        logging.warning(message)

    # Plot energy bins
    energies, counts = [], []
    for energy, count in results.items():
        # if float(energy) > ENERGY_DEFAULTS['ENERGY_MIN_DENSITY']:
        #     break
        energies.append(float(energy))
        counts.append(count)

    if name is not None:
        title = name + ': ' + title

    # energy_width = float(energies[1]) + 0.02

    plot_params = {
        'x_vals': energies,
        'y_vals': counts,
        'plot_file_name': plot_file_name,
        'title': title,
        'xlabel': ENERGY_DEFAULTS['XLABEL'],
        'ylabel': ENERGY_DEFAULTS['YLABEL'],
        'width': bar_params['BAR_WIDTH'],
        # 'width': energy_width,
        'align': bar_params['BAR_ALIGN'],
        'edgecolor': bar_params['BAR_EDGE_COLOR'],
        'rc_params': rc_params,
    }
    _plot_energy_histogram(plot_params)


# Public Methods : Type : type_count
# Plot a pie plot from a hybkit Type Analysis
def type_count(
        results: Counter,
        plot_file_name: str,
        title: str,
        name: Optional[str] = None,
        join_entries: bool = False,
        rc_params: Dict[str, Any] = TYPE_PIE_SINGLE_RC_PARAMS,
        ) -> None:
    """
    Plot pie chart of hybrid type counts from an :class:`~hybkit.analysis.Analysis` type analysis.

    Args:
        results (~collections.Counter): Counter Object of type counts from an
            :class:`~hybkit.analysis.Analysis` type analysis.
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        join_entries (:obj:`bool`, optional): If True, join two-tuple pairs into a
            single string for plot labels.
        rc_params (:obj:`dict`, optional): Dictionary of mpl rcParams. Defaults
            to :obj:`~hybkit.plot.TYPE_PIE_RC_PARAMS`.
    """
    # Sort hybrid_type_counts counter object in descending order
    labels, counts = [], []
    for key, count in results.most_common():
        use_key = key
        if join_entries:
            use_key = '--'.join(key)
        labels.append(use_key)
        counts.append(count)

    if name is not None:
        title = str(name) + ': ' + title

    plot_params = {
        'labels': labels,
        'sizes': counts,
        'plot_file_name': plot_file_name,
        'title': title,
        'other_threshold': PIE_DEFAULTS['OTHER_THRESHOLD'],
        'min_wedge_size': PIE_DEFAULTS['MIN_WEDGE_SIZE'],
        'colors': PIE_DEFAULTS['COLORS'],
        'rc_params': rc_params,
        'settings': copy.deepcopy(PIE_DEFAULTS['SETTINGS']),
    }

    _plot_types_pie_chart(plot_params)


# Public Methods : Type : type_count
# Plot a pie plot for two types from a hybkit Type Analysis
def type_count_dual(
        results: Counter,
        plot_file_name: str,
        title: str,
        name: Optional[str] = None,
        join_entries: bool = False,
        rc_params: Dict[str, Any] = TYPE_PIE_DUAL_RC_PARAMS,
        ) -> None:
    """Hold Place for replaced docstring."""
    return type_count(results, plot_file_name, title, name, join_entries, rc_params)


# Modify type_count_dual
type_count_dual.__doc__ = type_count.__doc__


# Public Methods : Target : target_count
# Plot a pie plot for targets from a hybkit target analysis
# @copy_signature(type_count)
def target_count(*args, **kwargs) -> None:  # noqa: ANN003, ANN002
    """Hold Place for replaced docstring."""
    return type_count(*args, rc_params=TARGET_PIE_RC_PARAMS, **kwargs)


# Modify target_count doc
target_count.__doc__ = type_count.__doc__
for q, r in [('hybrid type counts', 'target counts'),
             ('TYPE_PIE_RC_PARAMS', 'TARGET_PIE_RC_PARAMS')]:
    target_count.__doc__ = target_count.__doc__.replace(q, r)


def fold_match_counts_histogram(
        results: Counter,
        plot_file_name: str,
        title: str,
        name: Optional[str] = None,
        is_prop: bool = False,
        rc_params: Dict[str, Any]= FOLD_MATCH_HIST_RC_PARAMS,
        bar_params: Dict[str, Any] = BAR_INT_DEFAULTS,
        ) -> None:
    """
    Plot histogram of predicted miRNA/target match count.

    Args:
        results (~collections.Counter): Counter Object of match counts from an
            :class:`~hybkit.analysis.Analysis` type analysis.
        plot_file_name (str): Name of output file.
        title (str): Title of plot.
        is_prop (:obj:`bool`, optional): If True, y axis is proportion.
        name (:obj:`str`, optional): Name of analysis to be included in plot title.
        rc_params (:obj:`dict`, optional): Dictionary of mpl rcParams. Defaults
            to :obj:`~hybkit.plot.FOLD_MATCH_HIST_RC_PARAMS`.
        bar_params (:obj:`dict`, optional): Dictionary of bar plot parameters. Defaults
            to :obj:`~hybkit.plot.BAR_INT_DEFAULTS`.
    """
    x_vals = []
    y_vals = []
    init_max_val = -10000
    min_x = 10000
    max_x = -10000
    max_val = init_max_val
    for x_val, y_val in results.items():
        min_x = min(min_x, x_val)
        max_x = max(max_x, x_val)
        max_val = max(max_val, y_val)

    if max_val == init_max_val:
        message = f'Warning: Attempted to create empty plot to name: {plot_file_name}'
        logging.warning(message)
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

    if is_prop:
        y_label = 'Proportion of Hybrids'
    else:
        y_label = 'Hybrid Count'

    plot_params = {
        'x_vals': x_vals,
        'y_vals': y_vals,
        'plot_file_name': plot_file_name,
        'title': title,
        'xlabel': 'Predicted miRNA/Target Matches',
        'ylabel': y_label,
        'width': bar_params['BAR_WIDTH'],
        'align': bar_params['BAR_ALIGN'],
        'rc_params': rc_params,
        'edgecolor': bar_params['BAR_EDGE_COLOR'],
    }

    _plot_int_hist(plot_params=plot_params)

# @copy_signature(fold_match_counts_histogram)
def fold_mirna_nt_counts_histogram(*args, **kwargs) -> None:  # noqa: ANN002, ANN003
    """Hold Place for replaced docstring."""
    if 'rc_params' not in kwargs:
        kwargs['rc_params'] = copy.deepcopy(FOLD_NT_COUNTS_HIST_RC_PARAMS)
    return fold_match_counts_histogram(*args, **kwargs)


# Modify target_count doc
fold_mirna_nt_counts_histogram.__doc__ = fold_match_counts_histogram.__doc__
for q, r in [('predicted miRNA/target match count ', 'predicted binding at nt position'),
             ('FOLD_MATCH_HIST_RC_PARAMS', 'FOLD_NT_COUNTS_HIST_RC_PARAMS')]:
    fold_mirna_nt_counts_histogram.__doc__ = fold_mirna_nt_counts_histogram.__doc__.replace(q, r)


# ----- Begin Private Plotting Methods -----
# Private Methods : Energy : _plot_energy_histogram
def _plot_energy_histogram(plot_params: Dict[str, Any]) -> None:
    # Update plot parameters
    plt.rcParams.update(plot_params['rc_params'])

    # Create the figure with the specified size
    fig, ax = plt.subplots()

    # Create the bars for the histogram
    _bars = ax.bar(
        plot_params['x_vals'],
        plot_params['y_vals'],
        width=plot_params['width'],
        align=plot_params['align'],
        edgecolor=plot_params['edgecolor'],
    )

    # for items in zip(plot_params['x_vals'], plot_params['y_vals']):
    #     pass

    # Invert the x-axis and set tick locators
    ax.invert_xaxis()
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))

    # Set x-axis limits and ticks
    left_xlim = max(0, int(plot_params['x_vals'][0]))
    last_x_tick = int(plot_params['x_vals'][-1]) - 1
    plt.xlim(left=left_xlim, right=last_x_tick)
    plt.xticks(range(0, last_x_tick, -2), rotation=-30)

    # Set x-axis and y-axis labels
    plt.xlabel(plot_params['xlabel'])
    plt.ylabel(plot_params['ylabel'])

    # Set y-axis minor tick locators
    ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

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
    mpl.rcdefaults()


# Private Methods : Energy : _plot_types_pie_chart
def _plot_types_pie_chart(plot_params: Dict[str, Any]) -> None:
    # Calculate the total size and the size of each fraction
    total_size = sum(plot_params['sizes'])
    if total_size < 0.00000001:  #noqa: PLR2004
        message = 'Warning: Attempted to create empty plot to name:'
        message += ' %s' % plot_params['plot_file_name']
        logging.warning(message)

    fraction_sizes = [(size / total_size) for size in plot_params['sizes']]

    # Filter labels and sizes based on minimum wedge size and other threshold
    use_labels, use_sizes = [], []
    for i in range(len(plot_params['labels'])):
        total_fraction = sum(use_sizes) / total_size
        if (fraction_sizes[i] < plot_params['min_wedge_size']
            or (total_fraction > (1 - plot_params['other_threshold'])
              and i != (len(plot_params['labels']) - 1))):
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

    # Set colors, and match length if more slices than colors.
    colors = plot_params['colors']
    if len(use_sizes) > len(colors):
        colors = colors * (len(use_sizes) // len(colors)) + colors[:len(use_sizes) % len(colors)]

    # Create the pie chart
    patches, texts, autotexts = ax.pie(
        use_sizes,
        colors=colors,
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
    mpl.rcdefaults()


# Private Methods : Energy : _plot_energy_histogram
def _plot_int_hist(
        plot_params: Dict[str, Any],
        truncate_to_first_int: bool = True
        ) -> None:
    # Update plot parameters
    plt.rcParams.update(plot_params['rc_params'])

    # Create the figure with the specified size
    fig, ax = plt.subplots()

    # Create the bars for the histogram
    _bars = ax.bar(
        plot_params['x_vals'],
        plot_params['y_vals'],
        width=plot_params['width'],
        align=plot_params['align'],
        edgecolor=plot_params['edgecolor'],

    )

    # Invert the x-axis and set tick locators
    # ax.invert_xaxis()
    # ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))

    # Set x-axis limits and ticks
    if truncate_to_first_int and int(plot_params['x_vals'][0]) > 0:
        left_xlim = int(plot_params['x_vals'][0]) - 1
        plt.xlim(left=left_xlim)
    plt.xticks(range(2, int(plot_params['x_vals'][-1]) + 1, 2))

    # Set x-axis and y-axis labels
    plt.xlabel(plot_params['xlabel'])
    plt.ylabel(plot_params['ylabel'])

    # Set y-axis minor tick locators
    # ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

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
    mpl.rcdefaults()
