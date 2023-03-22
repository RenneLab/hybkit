#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Functions for analysis of HybRecord and FoldRecord objects."""

import sys
import copy
from collections import Counter
import hybkit
import numpy as np

# Import module-level dunder-names:
from hybkit.__about__ import (__author__, __contact__, __credits__, __date__, __deprecated__,
                              __email__, __license__, __maintainer__, __status__, __version__)

# --- Base Analysis --- #
class BaseAnalysis(object):
    """Base class for specific analysis class types."""

    #: Class-level settings. See :attr:`hybkit.settings.Analysis_settings` for descriptions.
    settings = hybkit.settings.Analysis_settings

    # BaseAnalysis : Public Methods
    def __init__(self, *args, **kwargs):
        """Stub method to be replaced by subclasses."""
        message = 'BaseAnalysis is a base class and is not meant to be used directly.'
        print(message)
        raise NotImplementedError(message)

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _ensure_same_class(self, cmp_obj):
        if not type(self) == type(cmp_obj):
            message = 'Objects are of different classes and cannot be compared.\n'
            message += '    %s (%s)\n' % (str(self), type(self))
            message += '    %s (%s)\n' % (str(cmp_obj), type(cmp_obj))
            print(message)
            raise RuntimeError(message)

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _format_all_seg_types(self, out_delim):
        """Return the results of all_seg_types in a list of delimited lines."""
        sorted_pairs = self.all_seg_types.most_common()
        ret_lines = [''.join(['seg_type', out_delim, 'count'])]
        ret_lines += ['%s%s%i' % (key, out_delim, count) for (key, count) in sorted_pairs]
        return ret_lines

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _format_hybrid_type_counts(self, out_delim):
        """Return the results of hybrid_type_counts in a list of delimited lines."""
        ret_lines = ['hybrid_type' + out_delim + 'count']
        sorted_pairs = self.hybrid_types.most_common()
        ret_lines += ['%s%s%i' % (key, out_delim, count) for (key, count) in sorted_pairs]
        return ret_lines

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _init_type_analysis(self):
        self.hybrid_types = Counter()
        self.seg1_types = Counter()
        self.seg2_types = Counter()
        self.all_seg_types = Counter()

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _add_type_analysis(self, hyb_record):
        hyb_record._ensure_set('eval_types')
        count = hyb_record.get_count(self.settings['count_mode'])
        seg1_type, seg2_type = hyb_record.get_seg_types()

        if self.settings['mirna_sort'] and hyb_record.is_set('eval_mirna'):
            if hyb_record.has_prop('5p_mirna'):
                join_types = seg1_type, seg2_type
            elif hyb_record.has_prop('3p_mirna'):
                join_types = seg2_type, seg1_type
            else:
                join_types = sorted((seg1_type, seg2_type))
        else:
            join_types = sorted((seg1_type, seg2_type))
        hybrid_type = self.settings['type_sep'].join(join_types)

        self.hybrid_types[hybrid_type] += count
        self.seg1_types[seg1_type] += count
        self.seg2_types[seg2_type] += count

        for seg_type in seg1_type, seg2_type:
            self.all_seg_types[seg_type] += count

    # BaseAnalysis : Private Methods : TypeAnalysis
    def _update_type_analysis(self, add_analysis):
        self.hybrid_types.update(add_analysis.hybrid_types)
        self.seg1_types.update(add_analysis.seg2_types)
        self.seg2_types.update(add_analysis.seg1_types)
        self.all_seg_types.update(add_analysis.all_seg_types)

    # BaseAnalysis : Private Methods : MirnaAnalysis
    def _init_mirna_analysis(self):
        self.mirnas_5p = 0
        self.mirnas_3p = 0
        self.mirna_dimers = 0
        self.non_mirnas = 0
        self.has_mirna = 0
        self._mirna_counts = ['mirnas_5p', 'mirnas_3p', 'mirna_dimers', 'non_mirnas', 'has_mirna']

    # BaseAnalysis : Private Methods : MirnaAnalysis
    def _add_mirna_analysis(self, hyb_record):
        hyb_record._ensure_set('eval_mirna')
        count = hyb_record.get_count(self.settings['count_mode'])

        if hyb_record.has_prop('mirna_dimer'):
            self.mirna_dimers += count
            self.has_mirna += count
        elif hyb_record.has_prop('5p_mirna'):
            self.mirnas_5p += count
            self.has_mirna += count
        elif hyb_record.has_prop('3p_mirna'):
            self.mirnas_3p += count
            self.has_mirna += count
        elif hyb_record.has_prop('no_mirna'):
            self.non_mirnas += count
        else:
            raise RuntimeError()

    # BaseAnalysis : Private Methods : MirnaAnalysis
    def _update_mirna_analysis(self, add_analysis):
        self.mirnas_5p += add_analysis.mirnas_5p
        self.mirnas_3p += add_analysis.mirnas_3p
        self.mirna_dimers += add_analysis.mirna_dimers
        self.non_mirnas += add_analysis.non_mirnas
        self.has_mirna += add_analysis.has_mirna

    # BaseAnalysis : Private Methods : MirnaAnalysis
    def _format_mirna_counts(self, out_delim):
        ret_lines = ['miRNA_type' + out_delim + 'count']
        for key in self._mirna_counts:
            ret_lines.append('%s%s%i' % (key, out_delim, getattr(self, key)))
        return ret_lines

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _init_target_analysis(self):
        self.mirna_target_info = {}
        self.mirna_target_count = 0
        self.mirna_target_total_counts = None
        self.mirna_target_type_counts = None

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _add_target_analysis(self, hyb_record, count_mode, allow_mirna_dimers):
        hyb_record._ensure_set('eval_mirna')
        count = hyb_record.get_count(count_mode)

        mirna_name = None

        if hyb_record.has_prop('has_mirna'):
            if hyb_record.has_prop('mirna_dimer'):
                if allow_mirna_dimers:
                    mirna_name = hyb_record.seg1_props['ref_name']
                    mirna_type = hyb_record.get_seg1_type()
                    target_name = hyb_record.seg2_props['ref_name']
                    target_type = hyb_record.get_seg2_type()
            else:
                mirna_details = hyb_record.mirna_detail()
                mirna_name = mirna_details['mirna_ref']
                mirna_type = mirna_details['mirna_seg_type']
                target_name = mirna_details['target_ref']
                target_type = mirna_details['target_seg_type']

        if mirna_name is not None:
            mirna_id = (mirna_name, mirna_type)
            target_id = (target_name, target_type)
            if mirna_id not in self.mirna_target_info:
                self.mirna_target_info[mirna_id] = Counter()
            self.mirna_target_info[mirna_id][target_id] += count
            self.mirna_target_count += count

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _update_target_analysis(self, add_analysis):
        for mirna_id in add_analysis.mirna_target_info:
            if mirna_id not in self.mirna_target_info:
                self.mirna_target_info[mirna_id] = add_analysis.mirna_target_info[mirna_id]
            else:
                self.mirna_target_info[mirna_id].update(add_analysis.mirna_target_info[mirna_id])
        self.mirna_target_count += add_analysis.mirna_target_count

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _process_target_analysis_info(self):
        self.mirna_target_total_counts = Counter()
        self.mirna_target_type_counts = {}
        # Iterate over each miRNA
        for mirna_id in sorted(self.mirna_target_info.keys()):
            mirna_name, mirna_type = mirna_id
            # Sum target counts for each mirna
            total_count = sum(self.mirna_target_info[mirna_id].values())
            self.mirna_target_total_counts[mirna_name] = total_count
            # Count target types for each mirna
            mirna_i_target_type_counts = Counter()
            for target_id in self.mirna_target_info[mirna_id].keys():
                target_name, target_type = target_id
                target_count = self.mirna_target_info[mirna_id][target_id]
                mirna_i_target_type_counts[target_type] += target_count
            self.mirna_target_type_counts[mirna_name] = mirna_i_target_type_counts
            # Sort mirna targets by most common

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _format_target_analysis(self, out_delim, only_mirna_id=None):
        ret_lines = [out_delim.join(['mirna', 'mirna_type', 'target', 'target_type', 'count'])]
        if only_mirna_id is not None:
            use_mirna_ids = [only_mirna_id]
        else:
            use_mirna_ids = self.mirna_target_info.keys()

        for mirna_id in use_mirna_ids:
            mirna_name, mirna_seg_type = mirna_id
            for target_id, count in self.mirna_target_info[mirna_id].most_common():
                target_name, target_seg_type = target_id
                count_str = str(count)
                line_items = [mirna_name, mirna_seg_type, target_name, target_seg_type, count_str]
                ret_lines.append(out_delim.join(line_items))
            ret_lines.append('')
        return ret_lines

    # BaseAnalysis : Private Methods : TargetAnalysis
    def _format_target_analysis_totals(self, out_delim):
        ret_lines = [out_delim.join(['mirna', 'target_count'])]
        for mirna_id in self.mirna_target_info:
            mirna_name, mirna_seg_type = mirna_id
            mirna_total_count = str(self.mirna_target_total_counts[mirna_name])
            ret_lines.append(out_delim.join([mirna_name, mirna_total_count]))
        return ret_lines

    # BaseAnalysis : Private Methods : PatternAnalysis
    def _init_pattern_analysis(self):
        self.mirna_folds = 0
        self.mirna_fold_count = {i: 0 for i in range(1, 25)}
        self.mirna_fold_frac = {i: 0.0 for i in range(1, 25)}

    # BaseAnalysis : Private Methods : PatternAnalysis
    def _add_pattern_analysis(self, hyb_record, count_mode, allow_mirna_dimers=False):
        hyb_record._ensure_set('eval_mirna')
        hyb_record._ensure_set('fold_record')
        count = hyb_record.get_count(count_mode)

        mirna_fold = None

        if hyb_record.has_prop('has_mirna'):
            if hyb_record.has_prop('mirna_dimer'):
                if allow_mirna_dimers:
                    mirna_fold = hyb_record.fold_record._get_seg_fold(hyb_record.seg1_props,
                                                                      hyb_record)
            else:
                mirna_fold = hyb_record.mirna_detail('mirna_fold')

            if (('(' in mirna_fold and ')' in mirna_fold)
                    or ('(' not in mirna_fold and ')' not in mirna_fold)
                    or len(mirna_fold) < 5):
                message = 'WARNING: addto_mirna_fold: Record: %s, ' % str(hyb_record)
                message += 'Bad Fold: %s' % mirna_fold
                print(message)
                return

        if mirna_fold is not None:
            self.mirna_folds += count
            for i in range(len(mirna_fold)):
                if mirna_fold[i] in {'(', ')'}:
                    self.mirna_fold_count[(i + 1)] += count
            for i in self.mirna_fold_count:
                self.mirna_fold_frac[i] = (self.mirna_fold_count[i] / self.mirna_folds)

    # BaseAnalysis : Private Methods : PatternAnalysis
    def _update_pattern_analysis(self, add_analysis):
        self.mirna_folds += add_analysis.mirna_folds
        for i in add_analysis.mirna_fold_count.keys():
            if i not in self.mirna_fold_count:
                self.mirna_fold_count[i] = 0
            self.mirna_fold_count[i] += add_analysis.mirna_fold_count[i]
        for i in self.mirna_fold_count.keys():
            self.mirna_fold_frac[i] = (self.mirna_fold_count[i] / self.mirna_folds)

    # BaseAnalysis : Private Methods : PatternAnalysis
    def _format_pattern_analysis(self, out_delim):
        ret_lines = [out_delim.join(['data', 'count'])]
        ret_lines.append(out_delim.join(['mirna_folds', str(self.mirna_folds)]))
        ret_lines.append('')

        line_values = ['index']
        max_i = 24
        while self.mirna_fold_count[max_i] == 0:
            max_i -= 1
        line_values += [str(i) for i in range(1, (max_i + 1))]
        ret_lines.append(out_delim.join(line_values))

        line_values = ['i_count']
        line_values += [str(self.mirna_fold_count[i]) for i in range(1, (max_i + 1))]
        ret_lines.append(out_delim.join(line_values))

        line_values = ['i_fraction']
        line_values += ["%f.3" % self.mirna_fold_frac[i] for i in range(1, (max_i + 1))]
        ret_lines.append(out_delim.join(line_values))

        return ret_lines

    # -- Fold Analysis Methods --
    # BaseAnalysis : Private Methods : FoldAnalysis
    def _init_fold_analysis(self):
        self.energies = Counter()
        self.match_counts = Counter()
        self.energy_bins = None
        self.energy_bin_densities = None
        for i in range(26):
            self.match_counts[i] = 0
        self._energy_bin_strf = '%.1f'
        energy_bin_size = int(round(((float(self.settings['energy_bin_size']) + 0.0001) * 10), 0))
        energy_strs_num = int((-10 * float(self.settings['energy_min_bin'])) + 0.0001)
       
        # Initialize bins sorted by least -> most magnitude
        self._energy_bins = [(-1 * i/10) for i 
                             in range(0, ((energy_strs_num) + 1), energy_bin_size)] 
        self._energy_bin_strs = [(self._energy_bin_strf % f) for f in self._energy_bins]
        # Reverse bins into lowest -> highest order
        self._energy_bins.reverse()

    # BaseAnalysis : Private Methods : FoldAnalysis
    def _add_fold_analysis(self, hyb_record, count_mode, allow_mirna_dimers=False):
        hyb_record._ensure_set('fold_record')
        count = hyb_record.get_count(count_mode)

        record_energy = hyb_record.fold_record.energy
        if record_energy is None:
            message = 'ERROR: add_fold_analysis: Record: %s, ' % str(hyb_record)
            message += 'Bad Energy: %s' % record_energy
            print(message)
            raise RuntimeError(message)

        self.energies[record_energy] += count

        record_fold = hyb_record.fold_record.fold
        l_matches = record_fold.count('(')
        r_matches = record_fold.count(')')
        if l_matches != r_matches:
            message = 'ERROR: add_fold_analysis: Record: %s, ' % str(hyb_record)
            message += 'Malformed Fold, L/R Matches do not match: %s' % record_fold
            print(message)
            raise RuntimeError(message)
        self.match_counts[l_matches] += count       
        
    # BaseAnalysis : Private Methods : FoldAnalysis
    def _update_fold_analysis(self, add_analysis):
        self.energies += add_analysis.energies
        self.match_counts += add_analysis.match_counts
        self.energy_bins = None
        self.energy_bin_densities = None

    # BaseAnalysis : Private Methods : FoldAnalysis
    def _process_fold_analysis(self):
        # Create binned energy analysis with counts
        energies = [float(i) for i in self.energies.keys()]
        weights = [i for i in self.energies.values()]
        bin_vals, keys = np.histogram(energies, bins=self._energy_bins, weights=weights)
        # Reverse Histogram to put in low -> high magnitude order
        bin_vals = bin_vals[::-1]
        keys = keys[::-1]
        energy_bins = {}
        for i in range(len(bin_vals)):
            energy_bins[self._energy_bin_strf % keys[i]] = bin_vals[i]
        energy_bins[keys[-1]] = 0

        self.energy_bins = energy_bins

        # Create binned energy analysis with density
        bin_vals, keys = np.histogram(energies, bins=self._energy_bins, weights=weights, 
                         density=True)
        # Reverse Histogram to put in low -> high magnitude order
        bin_vals = bin_vals[::-1]
        keys = keys[::-1]
        energy_bin_densities = {}
        for i in range(len(bin_vals)):
            energy_bin_densities[self._energy_bin_strf % keys[i]] = bin_vals[i]
        energy_bin_densities[keys[-1]] = 0

        self.energy_bin_densities = energy_bin_densities
        
    # BaseAnalysis : Private Methods : FoldAnalysis
    def _format_fold_analysis_matches(self, out_delim):
        ret_lines = [out_delim.join(['data'] + [str(i) for i in self.match_counts.keys()])]
        out_line = ['hybrid_match_counts']
        out_line += [str(v) for v in self.match_counts.values()]
        ret_lines.append(out_delim.join(out_line))
        return ret_lines

    # BaseAnalysis : Private Methods : FoldAnalysis
    def _format_fold_analysis_energies(self, out_delim):
        ret_lines = [out_delim.join(['data'] + self._energy_bin_strs)]
        out_line = ['hybrid_energies']
        out_line += [str(v) for v in self.energy_bins.values()]
        ret_lines.append(out_delim.join(out_line))
        return ret_lines

    # BaseAnalysis : Private Methods : FoldAnalysis
    def _format_fold_analysis_energy_densities(self, out_delim):
        ret_lines = [out_delim.join(['data'] + self._energy_bin_strs)]
        out_line = ['hybrid_energy_densities']
        out_line += ['%.4f' % v for v in self.energy_bin_densities.values()]
        ret_lines.append(out_delim.join(out_line))
        return ret_lines
    
    
    # BaseAnalysis : Private Classmethods
    @classmethod
    def _sanitize_name(self, file_name):
        for char, replace in [('*', 'star'), (',', 'com')]:
            file_name = file_name.replace(char, replace)
        return file_name


# --- Type Analysis --- #
class TypeAnalysis(BaseAnalysis):
    """
    Analysis of segment types included in the analyzed hyb_records.

    Before using the analysis, the :ref:`seg1_type <seg1_type>` and
    :ref:`seg2_type <seg2_type>` flags must be set for the hyb_record,
    as is done by :func:`hybkit.HybRecord.eval_types`.
    A count is added to
    the analysis class object for each hybrid type (Ex: "miRNA-mRNA") with
    segments placed in sorted order for non-redundant type-combinations.
    The analysis additionally reports the number of individual segment
    types.

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        hybrid_types (~collections.Counter): Counter containing annotated ordered types of
            seg1 and seg2
        seg1_types (~collections.Counter): Counter containing annotated type of
            segment in position seg1
        seg2_types (~collections.Counter): Counter containing annotated type of
            segment in position seg2
        all_seg_types (~collections.Counter): Counter containing position-independent
            annotated types

    """

    # TypeAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_type_analysis()

    # TypeAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add information from a :class:`~hybkit.HybRecord` record to analysis.

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        self._add_type_analysis(hyb_record)

    # TypeAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another TypeAnalysis object to this one by combining counts.

        Args:
            add_analysis (TypeAnalysis): :class:`TypeAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_type_analysis(add_analysis)

    # TypeAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a type analysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        ret_lines = []
        ret_lines += self._format_hybrid_type_counts(out_delim)
        ret_lines.append('')
        ret_lines += self._format_all_seg_types(out_delim)
        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # TypeAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of a type analysis to a file.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        analyses = [('type_hybrids', self._format_hybrid_type_counts),
                    ('type_segs', self._format_all_seg_types),
                    ]

        for analysis_name, analysis_method in analyses:
            analysis_file_name = file_name_base + '_' + analysis_name + '.csv'
            with open(analysis_file_name, 'w') as out_file:
                out_lines = analysis_method(out_delim=out_delim)
                out_str = '\n'.join(out_lines) + '\n'
                out_file.write(out_str)

    # TypeAnalysis : Public Methods
    def plot(self, file_name_base):
        """
        Create plots of the results using :func:`hybkit.plot.type_count`.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        hybkit.plot.type_count(self.hybrid_types, file_name_base + '_type_hybrids',
                               title='Hybrid Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_segall',
                               title='Total Segment Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_seg1',
                               title='5p Segment Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_seg2',
                               title='3p Segment Types Types', name=self.name)


# --- miRNA Count Analysis --- #
class MirnaAnalysis(BaseAnalysis):
    """
    Analyze counts/types of miRNA in hyb records.

    The mirna_count analysis determines what type each record is
    with regard to mirna and counts them accordingly.
    This includes:

        | :obj:`mirnas_5p`: Hybrids with a |5p| miRNA.
        | :obj:`mirnas_3p`: Hybrids with a |3p| miRNA.
        | :obj:`mirna_dimers`: Hybrids with both a |5p| and |3p| miRNA.
        | :obj:`non_mirnas`: Hybrids with no miRNA.
        | :obj:`has_mirna`: Hybrids with a miRNA (one of first three categories).

    Before using the analysis, the :ref:`mirna_seg <mirna_seg>` flag
    must be set for each record as can be done by sequential use of the
    :func:`hybkit.HybRecord.eval_types` and :func:`hybkit.HybRecord.eval_mirna`
    methods.

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        mirnas_5p (int): Count of |5p| miRNAs detected
        mirnas_3p (int): Count of |3p| miRNAs detected
        mirna_dimers (int): Count of miRNA dimers (|5p| + |3p|) detected
        non_mirnas (int): Count of non-miRNA hybrids detected
        has_mirna (int): Hybrids with |5p|, |3p|, or both as miRNA
    """

    # MirnaAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_mirna_analysis()

    # MirnaAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add information from a :class:`~hybkit.HybRecord` record to analysis.

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        self._add_mirna_analysis(hyb_record)

    # MirnaAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another MirnaAnalysis object to this one by combining counts.

        Args:
            add_analysis (MirnaAnalysis): :class:`MirnaAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_mirna_analysis(add_analysis)

    # MirnaAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a miRNA analysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        ret_lines = self._format_mirna_counts(out_delim)

        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # MirnaAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of a type analysis to a file.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        analysis_file_name = file_name_base + '_mirna_counts.csv'
        with open(analysis_file_name, 'w') as out_file:
            out_lines = self.results(out_delim=out_delim)
            out_str = '\n'.join(out_lines) + '\n'
            out_file.write(out_str)

    # MirnaAnalysis : Public Methods
    def plot(self, file_name_base):
        """
        Create plots of the results using :func:`hybkit.plot.mirna`.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        hybkit.plot.mirna(self, file_name_base + '_mirna_counts')


# --- Full Summary Analysis --- #

class SummaryAnalysis(BaseAnalysis):
    """
    Analysis of segment types and miRNA counts in hyb records.

    This analysis includes the components of both the :class:`TypeAnalysis` and
    :class:`MirnaAnalysis` analyses, performed simultaneously.

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        hybrid_types (~collections.Counter): Counter containing annotated ordered types of
            seg1 and seg2
        seg1_types (~collections.Counter): Counter containing annotated type of
            segment in position seg1
        seg2_types (~collections.Counter): Counter containing annotated type of
            segment in position seg2
        all_seg_types (~collections.Counter): Counter containing position-independent
            annotated types
        mirnas_5p (int): Count of |5p| miRNAs detected
        mirnas_3p (int): Count of |3p| miRNAs detected
        mirna_dimers (int): Count of miRNA dimers (|5p| + |3p|) detected
        non_mirnas (int): Count of non-miRNA hybrids detected
        has_mirna (int): Hybrids with |5p|, |3p|, or both as miRNA
    """

    # SummaryAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_type_analysis()
        self._init_mirna_analysis()

    # SummaryAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add information from a :class:`hybkit.HybRecord` record to analysis.

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        self._add_type_analysis(hyb_record)
        self._add_mirna_analysis(hyb_record)

    # SummaryAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another SummaryAnalysis object to this one by combining counts.

        Args:
            add_analysis (SummaryAnalysis): :class:`SummaryAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_type_analysis(add_analysis)
        self._update_mirna_analysis(add_analysis)

    # SummaryAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a miRNA analysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        ret_lines = []
        ret_lines += self._format_hybrid_type_counts(out_delim)
        ret_lines.append('')
        ret_lines += self._format_all_seg_types(out_delim)
        ret_lines.append('')
        ret_lines += self._format_mirna_counts(out_delim)

        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # SummaryAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of a type analysis to a file.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        analyses = [('type_hybrids', self._format_hybrid_type_counts),
                    ('type_segs', self._format_all_seg_types),
                    ('mirna_counts', self._format_mirna_counts),
                    ]

        for analysis_name, analysis_method in analyses:
            analysis_file_name = file_name_base + '_' + analysis_name + '.csv'
            with open(analysis_file_name, 'w') as out_file:
                out_lines = analysis_method(out_delim=out_delim)
                out_str = '\n'.join(out_lines) + '\n'
                out_file.write(out_str)

    # SummaryAnalysis : Public Methods
    def plot(self, file_name_base):
        """
        Create plots of results via :func:`hybkit.plot.type_count` and :func:`hybkit.plot.mirna`.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        hybkit.plot.type_count(self.hybrid_types, file_name_base + '_type_hybrids',
                               title='Hybrid Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_segall',
                               title='Total Segment Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_seg1',
                               title='5p Segment Types', name=self.name)
        hybkit.plot.type_count(self.all_seg_types, file_name_base + '_type_seg2',
                               title='3p Segment Types Types', name=self.name)
        hybkit.plot.mirna(self, file_name_base + '_mirna_counts')


# --- Target Analysis --- #

class TargetAnalysis(BaseAnalysis):
    """
    Class for analysis of targets in miRNA-containing hybrids.

    The mirna_target analysis provides an analysis of what sequences are targeted
    by each respective miRNA within the hyb records. The attribute has keys
    of each miRNA, with each value being a dict of targeted sequences and their
    associated count of times targeted.

    Before using the analysis, the :ref:`seg1_type <seg1_type>`,
    :ref:`seg2_type <seg2_type>`, and :ref:`mirna_seg <mirna_seg>` flags
    must be set for each record as can be done by sequential use of the
    :func:`hybkit.HybRecord.eval_types` and :func:`hybkit.HybRecord.eval_mirna`
    methods.

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        targets_5p (int): Count of |5p| miRNAs detected
        targets_3p (int): Count of |3p| miRNAs detected
        target_dimers (int): Count of miRNA dimers (|5p| + |3p|) detected
        non_targets (int): Count of non-miRNA hybrids detected
        has_target (int): Hybrids with |5p|, |3p|, or both, as miRNA
    """

    # TargetAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_target_analysis()

    # TargetAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add the information from a :class:`~hybkit.HybRecord` to a mirna_target analysis.

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        count_mode = self.settings['count_mode']
        allow_mirna_dimers = self.settings['allow_mirna_dimers']
        self._add_target_analysis(hyb_record, count_mode, allow_mirna_dimers)

    # TargetAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another TargetAnalysis object to this one by combining counts.

        Args:
            add_analysis (TargetAnalysis): :class:`TargetAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_target_analysis(add_analysis)

    # TargetAnalysis : Public Methods
    def process_target_analysis(self):
        """
        Summarize results of the target analysis.

        This fills the analysis variables:

            | :obj:`mirna_target_total_counts` : :class:`Counter` with keys of miRNA names and
              values of miRNA total counts.
            | :obj:`mirna_target_type_counts` : Dict with keys of miRNA names and values of
              :class:`Counter` objects, with each counter object containin the count of
              segment types targeted by that miRNA.
        """
        self._process_target_analysis_info()

    # TargetAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a miRNA analysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        if self.mirna_target_total_counts is None:
            self.process_target_analysis()

        ret_lines = self._format_target_analysis(out_delim)
        ret_lines.append('')
        ret_lines += self._format_target_analysis_totals(out_delim)

        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # TargetAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of the target analysis to a file.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        if self.mirna_target_total_counts is None:
            self.process_target_analysis()

        analyses = [('targets_all', self._format_target_analysis),
                    ('targets_total', self._format_target_analysis_totals),
                    ]

        for analysis_name, analysis_method in analyses:
            analysis_file_name = file_name_base + '_' + analysis_name + '.csv'
            with open(analysis_file_name, 'w') as out_file:
                out_lines = analysis_method(out_delim=out_delim)
                out_str = '\n'.join(out_lines) + '\n'
                out_file.write(out_str)

    # TargetAnalysis : Public Methods
    def write_individual(self, file_name_base, out_delim=None):
        """
        Write the results of the target analysis to individual files per miRNA.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        for mirna_id in self.mirna_target_info:
            mirna, mirna_seg_type = mirna_id
            analysis_file_name = file_name_base + '_targets_' + mirna + '.csv'
            analysis_file_name = self._sanitize_name(analysis_file_name)
            out_lines = self._format_target_analysis(out_delim, only_mirna_id=mirna_id)
            with open(analysis_file_name, 'w') as out_file:
                out_file.write('\n'.join(out_lines) + '\n')

    # TargetAnalysis : Public Methnds
    def plot(self, file_name_base):
        """
        Plot targets of all miRNAs in analysis.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        if self.mirna_target_total_counts is None:
            self.process_target_analysis()

        combined_targets = Counter()
        combined_target_types = Counter()
        for mirna_id in self.mirna_target_info:
            mirna_name, _ = mirna_id
            combined_targets.update(self.mirna_target_info[mirna_id])
            combined_target_types.update(self.mirna_target_type_counts[mirna_name])
        hybkit.plot.target(combined_targets,
                           file_name_base + '_targets_all',
                           title='Combined miRNA Targets',
                           )
        hybkit.plot.target_type(combined_target_types,
                                file_name_base + '_targets_types_all',
                                title='Combined miRNA Target Types',
                                )

    # TargetAnalysis : Public Methnds
    def plot_individual(self, file_name_base):
        """
        Plot targets of each individual miRNA in analysis.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        if self.mirna_target_total_counts is None:
            self.process_target_analysis()

        for mirna_id in self.mirna_target_info:
            mirna_name, mirna_seg_type = mirna_id
            analysis_file_name = file_name_base + '_targets_' + mirna_name
            analysis_file_name = self._sanitize_name(analysis_file_name)
            analysis_target_file_name = file_name_base + '_targets_' + mirna_name + '_types'
            analysis_target_file_name = self._sanitize_name(analysis_target_file_name)

            hybkit.plot.target(self.mirna_target_info[mirna_id],
                               analysis_file_name,
                               name=mirna_name,
                               )
            hybkit.plot.target_type(self.mirna_target_type_counts[mirna_name],
                                    analysis_target_file_name,
                                    name=mirna_name,
                                    )


# --- miRNA Pattern Analysis --- #

class PatternAnalysis(BaseAnalysis):
    """
    Evaluate and quantify predicted miRNA binding patterns.

    This analysis evaluates the predicted binding of miRNA within hyb records
    that contain a miRNA and have an associated :class:`~hybkit.FoldRecord` object
    as the attribute :attr:`~hybkit.HybRecord.fold_record`. This includes an analysis and
    plotting of the predicted binding by position among the provided miRNA.

    Before using the analysis, the :ref:`mirna_seg <mirna_seg>` flag
    must be set for each record as can be done by sequential use of the
    :func:`hybkit.HybRecord.eval_types` and :func:`hybkit.HybRecord.eval_mirna`
    methods.

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        mirna_folds (int): Number of miRNA folds represented
        mirna_fold_counts (dict): Dict with keys of miRNA index and values of number of miRNAs
            folded at that index
        mirna_fold_frac (dict): Dict with keys of miRNA index and values of percent of miRNAs
            folded at that index
    """
    # PatternAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_pattern_analysis()

    # PatternAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add the information from a :class:`~hybkit.HybRecord` to a pattern analysis.

        If the record contains a single miRNA, the miRNA fold is identified.
        miRNA Dimers are skipped unless the :attr:`settings['all_mirna_dimers'] <settings>`
        setting is True.
        The count for this miRNA and its target is then added to :obj:`mirna_fold_counts` and
        :obj:`mirna_fold_frac` is recalculated.

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        count_mode = self.settings['count_mode']
        allow_mirna_dimers = self.settings['allow_mirna_dimers']
        self._add_pattern_analysis(hyb_record, count_mode, allow_mirna_dimers)

    # PatternAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another PatternAnalysis object to this one by combining counts.

        Args:
            add_analysis (PatternAnalysis): :class:`PatternAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_pattern_analysis(add_analysis)

    # PatternAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a miRNA analysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        ret_lines = self._format_pattern_analysis(out_delim)

        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # PatternAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of the pattern analysis to files.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        pattern_counts_name = file_name_base + '_pattern_counts.csv'
        with open(pattern_counts_name, 'w') as pattern_counts_file:
            pattern_counts_file.write(out_delim.join(['data', 'count']) + '\n')
            pattern_counts_file.write(out_delim.join(['mirna_folds', str(self.mirna_folds)]) + '\n')

        line_values = ['index']
        max_i = 24
        while self.mirna_fold_count[max_i] == 0:
            max_i -= 1
        line_values += [str(i) for i in range(1, (max_i + 1))]
        write_lines = [out_delim.join(line_values)]

        line_values = ['i_count']
        line_values += [str(self.mirna_fold_count[i]) for i in range(1, (max_i + 1))]
        write_lines.append(out_delim.join(line_values))

        line_values = ['i_fraction']
        line_values += ["%f.3" % self.mirna_fold_frac[i] for i in range(1, (max_i + 1))]
        write_lines.append(out_delim.join(line_values))

        pattern_bases_name = file_name_base + '_pattern_bases.csv'
        with open(pattern_bases_name, 'w') as pattern_bases_file:
            pattern_bases_file.write('\n'.join(write_lines) + '\n')

    # PatternAnalysis : Public Methods
    def plot(self, file_name_base):
        """
        Create plots of the results.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        hybkit.plot.pattern(self, file_name_base + '_pattern_bases')


# --- Hybrid Fold Analysis --- #

class FoldAnalysis(BaseAnalysis):
    """
    Evaluate hybrid fold energy and match count characteristics.

    This analysis evaluates the predicted Gibbs Free Energy of folding 
    and count of predicted Watson-Crick Base Pairing Interactions within
    hyb records that have an associated :class:`~hybkit.FoldRecord` object
    as the attribute :attr:`~hybkit.HybRecord.energy_record`. 

    Args:
        name (str, optional): Analysis name

    Attributes:
        name (str): Analysis name
        energies (collections.Counter): Python counter object containing strings representing
            hybrid energy counts.
        matches (collections.Counter): Python counter object containing counts of hybrids
            with a predicted number of match interactions.

    """
    # FoldAnalysis : Public Methods
    def __init__(self,
                 name=None,
                 ):
        self.name = name
        self._init_fold_analysis()
    
    
    # FoldAnalysis : Public Methods
    def add(self, hyb_record):
        """
        Add the information from a :class:`~hybkit.HybRecord` to the fold analysis.
       
        Fold values > 0.0 are assumed to be artifacts and are capped at 0.0. 

        Args:
            hyb_record (HybRecord): Record with information to add.
        """
        count_mode = self.settings['count_mode']
        self._add_fold_analysis(hyb_record, count_mode)

    # FoldAnalysis : Public Methods
    def update(self, add_analysis):
        """
        Add another FoldAnalysis object to this one by combining counts.

        Args:
            add_analysis (FoldAnalysis): :class:`FoldAnalysis` object to add.
        """
        self._ensure_same_class(add_analysis)
        self._update_fold_analysis(add_analysis)

    # FoldAnalysis : Public Methods
    def process_fold_analysis(self):
        """
        Calculate binned results of fold analysis.

        This fills the analysis variables:

            | :obj:`energy_bins` : `dict` with keys of a Gibbs Free Energy quantity 
              and values of counts of hybrids with that predicted energy
            | :obj:`energy_bin_densities` : `dict` with keys of a Gibbs Free Energy quantity 
              and values of counts of hybrids with that predicted energy normalized by total
              counts
        """
        self._process_fold_analysis()

    # FoldAnalysis : Public Methods
    def results(self, out_delim=None, newline=False):
        r"""
        Return the results of a FoldAnalysis in a list of delimited lines.

        Args:
            out_delim (str, optional): Delimiter for entries within lines, such as ',' or '\\t'.
                If not provided, defaults to :attr:`settings['out_delim'] <settings>`
            newline (bool, optional): Terminate lines with a newline_character. Default False

        Returns:
            list of string objects representing the results of the analysis.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        if self.energy_bins is None:
            self.process_fold_analysis()

        ret_lines = self._format_fold_analysis_matches(out_delim)
        ret_lines += [''] + self._format_fold_analysis_energies(out_delim)
        ret_lines += [''] + self._format_fold_analysis_energy_densities(out_delim)

        if newline:
            for i in range(len(ret_lines)):
                ret_lines[i] += '\n'
        return ret_lines

    # FoldAnalysis : Public Methods
    def write(self, file_name_base, out_delim=None):
        """
        Write the results of the FoldAnalysis to files.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
            out_delim (str, optional): Delimiter to write between fields.
                Defaults to :attr:`settings['out_delim'] <settings>` if not provided.
        """
        if out_delim is None:
            out_delim = self.settings['out_delim']

        if self.energy_bins is None:
            self.process_fold_analysis()

        matches_name = file_name_base + '_fold_match.csv'
        energies_name = file_name_base + '_fold_energy.csv'
        energy_densities_name = file_name_base + '_fold_energy_densities.csv'

        with open(matches_name, 'w') as matches_file:
            write_lines = self._format_fold_analysis_matches(out_delim)
            matches_file.write('\n'.join(write_lines) + '\n')

        with open(energies_name, 'w') as energies_file:
            write_lines = self._format_fold_analysis_energies(out_delim)
            energies_file.write('\n'.join(write_lines) + '\n')

        with open(energy_densities_name, 'w') as energy_densities_file:
            write_lines = self._format_fold_analysis_energy_densities(out_delim)
            energy_densities_file.write('\n'.join(write_lines) + '\n')


    # FoldAnalysis : Public Methods
    def plot(self, file_name_base):
        """
        Create plots of the results.

        Args:
            file_name_base (str): "Base" name for output files. Final file names will be generated
                based on analysis type and provided parameters.
        """
        hybkit.plot.fold_energy_counts(self, file_name_base + '_fold_energy_counts')
        hybkit.plot.fold_energy_densities(self, file_name_base + '_fold_energy_densities')
        hybkit.plot.fold_matches(self, file_name_base + '_fold_match_counts')


