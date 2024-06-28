#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""Functions for analysis of HybRecord and FoldRecord objects."""

import copy
from collections import Counter
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np

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
from hybkit.errors import HybkitArgError

# ----- File-Specific Linting Directives:
# ruff: noqa: F401 SLF001

# ----- Begin Typing Directives -----
AnalysisOptions = Literal['energy', 'type', 'mirna', 'target']
AnalysisArg = Union[AnalysisOptions, List[AnalysisOptions]]
QuantModeArg = Literal['single', 'reads', 'records']

# --- Hybkit Analysis --- #
class Analysis:
    """
    Class for analysis of hybkit HybRecord and FoldRecord objects.

    .. _Analyses:

    This class contains multiple conceptual analyses for HybRecord/FoldRecord Data:

        | :ref:`Energy <EnergyAnalysis>`: Analysis of values of predicted intra-hybrid folding
            energy
        | :ref:`Type <TypeAnalysis>`: Analysis of segment types
        | :ref:`miRNA <MirnaAnalysis>`: Analysis of miRNA segments distributions
        | :ref:`Target <TargetAnalysis>`: Analysis of mirna target segment names and types
        | :ref:`Fold <FoldAnalysis>`: Analysis of folding data included in the analyzed
            hyb_records.

    This class used by selecting the desired analysis types on object
    initialization. Analyses are performed either by using either the :meth:`add_record`
    or the :meth:`add_all_records` methods. The results of the analysis are
    then available through the :meth:`get_all_results`, :meth:`get_analysis_results`,
    :meth:`get_specific_result`, and :meth:`plot_analysis_results`
    methods, which can return (or plot) the results
    of all analyses or of a specific subset of analyses.

    Details for each respective analysis are provided here:

    .. _EnergyAnalysis:

    **Energy Analysis:**

        This analysis evaluates the energy of each :class:`~hybkit.HybRecord` object
        and provides a binned-histogram of all energy values represented.

        Output Results:
            | ``energy_analysis_count`` (:obj:`int`): Count of energy values evaluated
            | ``has_energy_val`` (:obj:`int`): Count of hyb_records with an energy value
            | ``no_energy_val`` (:obj:`int`): Count of hyb_records without an energy value
            | ``energy_min`` (:obj:`float`): Minimum energy value
            | ``energy_max`` (:obj:`float`): Maximum energy value
            | ``energy_mean`` (:obj:`float`): Mean energy value
            | ``energy_std`` (:obj:`float`): Standard deviation of energy values
            | ``binned_energy_vals`` (:obj:`~collections.Counter`): Counter with integer keys of
              energy values from ``energy_min`` to ``energy_max`` storing the count of any
              hyb_records with energy values that fall within that range
              (rounded to the next highest integer (e.g. -12.5 -> -12).


    .. _TypeAnalysis:

    **Type Analysis:**

        This analysis evaluates the counts of each type of segment included in the
        :class:`~hybkit.HybRecord` objects. The types of segments are determined by
        the :ref:`seg1_type <seg1_type>` and :ref:`seg2_type <seg2_type>` flags, which
        are set by the :func:`hybkit.HybRecord.eval_types` method.

        Requirements:

            | :ref:`seg1_type <seg1_type>` and
              :ref:`seg2_type <seg2_type>` flags must be set for each HybRecord,
              (can be done by :func:`hybkit.HybRecord.eval_types`).

        Output Results:
            | ``types_analysis_count`` (:obj:`int`): Count of hybrid types analyzed
            | ``hybrid_types`` (:obj:`~collections.Counter`): Counter containing annotated
              types of seg1 and seg (in original |5p| / |3p| order)
            | ``reordered_hybrid_types`` (:obj:`~collections.Counter`): Counter containing
              annotated types of seg1 and seg2. This is provided in "sorted" order, where
              types are sorted alphabetically (independent of |5p| / |3p| position).
            | ``mirna_hybrid_types`` (:obj:`~collections.Counter`): Counter containing annotated
              types of seg1 and seg2. This is provided in "miRNA-prime" order, where a
              miRNA type is always listed before other types, and then remaining
              types are sorted alphabetically (independent of |5p| / |3p| position).
            | ``seg1_types`` (:obj:`~collections.Counter`): Counter containing annotated type of
              segment in position seg1
            | ``seg2_types`` (:obj:`~collections.Counter`): Counter containing annotated type of
              segment in position seg2
            | ``all_seg_types`` (:obj:`~collections.Counter`): Counter containing
                position-independent
              annotated types

    .. _MirnaAnalysis:

    **miRNA Analysis:**

        Analysis of miRNA segments in hybrids.

        The mirna_analysis provides an analysis of what miRNA types are present in
        the hyb records. If a miRNA dimer is present in a hybrid, this is counted
        in ``mirna_dimers``. If a single miRNA is present in a hybrid, this is
        counted in ``mirnas_5p`` or ``mirnas_3p`` depending on the miRNA location.

        Requirements:
            | :ref:`mirna_seg <mirna_seg>` flag
              must be set for each HybRecord
              (can be done by :func:`hybkit.HybRecord.eval_mirna`).

        Output Results:
            | ``mirna_analysis_count`` (:obj:`int`): Count of miRNA hybrids analyzed
            | ``mirnas_5p`` (:obj:`int`): Count of |5p| miRNAs detected
            | ``mirnas_3p`` (:obj:`int`): Count of |3p| miRNAs detected
            | ``mirna_dimers`` (:obj:`int`): Count of miRNA dimers (|5p| + |3p|) detected
            | ``non_mirna`` (:obj:`int`): Count of non-miRNA hybrids detected
            | ``has_mirna`` (:obj:`int`): Hybrids with |5p|, |3p|, or both as miRNA

    .. _TargetAnalysis:

    **Target Analysis:**

        Analysis of targets in miRNA-containing hybrids.

        The target analysis provides an analysis of what annotated sequences
        and sequence types are targeted by any miRNA within the hyb records. If a
        miRNA is not present in a hybrid, the hybrid is not included in the analysis.
        If a miRNA dimer is present in a hybrid, the |5p| miRNA is used for the analysis,
        and the |3p| miRNA is considered the "target."

        Requirements:
            | :ref:`mirna_seg <mirna_seg>` flag
              must be set for each HybRecord
              (can be done by :func:`hybkit.HybRecord.eval_mirna`).

        Output Results:
            | ``target_analysis_count`` (:obj:`int`): Count of hybrids analyzed
            | ``target_evals`` (:obj:`int`): Count of target evaluations performed
            | ``target_names`` (:obj:`~collections.Counter`): Counter containing names of
              miRNA targets detected.
            | ``target_types`` (:obj:`~collections.Counter`): Counter containing types of
              miRNA targets detected.

    .. _FoldAnalysis:

    **Fold Analysis:**

        This analysis evaluates the predicted binding of miRNA within hyb records
        that contain a miRNA and have an associated :class:`~hybkit.FoldRecord` object
        as the attribute :attr:`~hybkit.HybRecord.fold_record`. This includes an analysis and
        plotting of the predicted binding by position among the provided miRNA.

        Requirements:
            | The :ref:`mirna_seg <mirna_seg>` flag
              must be set for each HybRecord
              (can be done by :func:`hybkit.HybRecord.eval_mirna`).
            | The :ref:`fold_record <HybRecord-Attributes>` attribute must be set for each
              HybRecord with a corresponding :class:`~hybkit.FoldRecord` object. This can be done
              using the :meth:`hybkit.HybRecord.set_fold_record()` method.

        Output Results:
            | ``fold_analysis_count`` (:obj:`int`): Count of miRNA fold predictions analyzed
            | ``folds_recorded`` (:obj:`int`): Count of fold predictions with a mirna fold
            | ``mirna_nt_fold_counts`` (:obj:`~collections.Counter`) : Counter with keys of
              miRNA position index and values of number of miRNAs with a predicted
              bound state at that index.
            | ``mirna_nt_fold_props`` (:obj:`~collections.Counter`) : Counter with keys of
              miRNA position index and values of proportion (0.0 - 1.0) of miRNAs
              with a predicted bound state at that index.
            | ``fold_match_counts`` (:obj:`~collections.Counter`) : Counter with keys of
              count of predicted matches between miRNA and target with
              values of count of miRNAs with that number of predicted matches.

    Args:
        analysis_types (:obj:`str` or :obj:`list` of :obj:`str`): Analysis types to perform
        name (:obj:`str`, optional): Name of the analysis
        quant_mode (:obj:`str`, optional): Mode to use for record quantification.
            Options are "single": One count per record; "reads": If "read_count" flag is set, count
            all reads in record (else count 1); "records": if the "record_count" flag is set, count
            all individual records within combined record (else count 1). If not provided,
            defaults to the value in :attr:`Analysis.settings['quant_mode'].`

    .. _Analysis-Attributes:

    Attributes:
        name (:obj:`str`): Name of the analysis
        analysis_types (:obj:`list` of :obj:`str`): List of analysis types to perform
        quant_mode (:obj:`str`): Mode to use for record quantification.
    """

    #: Class-level settings. See :attr:`hybkit.settings.Analysis_settings` for descriptions.
    settings = hybkit.settings.Analysis_settings

    #: Possible options for analyses
    # analysis_options = ['energy', 'type', 'mirna', 'target', 'fold']
    analysis_options = hybkit.settings.ANALYSIS_TYPE_OPTIONS

    # Class private variables:
    _result_keys = {
        'energy': (
            'energy_analysis_count', 'has_energy_val', 'no_energy_val',
            'energy_min', 'energy_max', 'energy_mean', 'energy_std',
            'binned_energy_vals',
        ),
        'type': (
            'types_analysis_count', 'hybrid_types', 'reordered_hybrid_types',
            'mirna_hybrid_types', 'seg1_types', 'seg2_types', 'all_seg_types',
        ),
        'mirna': (
            'mirna_analysis_count', 'mirnas_5p', 'mirnas_3p', 'mirna_dimers',
            'non_mirna', 'has_mirna',
        ),
        'target': (
            'target_analysis_count', 'target_evals', 'target_names', 'target_types',
        ),
        'fold': (
            'fold_analysis_count', 'mirna_nt_fold_counts', 'mirna_nt_fold_props',
            'fold_match_counts',
        ),
    }
    _all_result_keys_list_temp = []  # noqa: RUF012
    for key in _result_keys:
        _all_result_keys_list_temp += _result_keys[key]

    _all_result_keys_list = tuple(_all_result_keys_list_temp)
    _all_result_keys_set = frozenset(_all_result_keys_list_temp)
    del _all_result_keys_list_temp

    _quant_mode_options = frozenset(
        hybkit.settings.Analysis_settings_info['quant_mode'][4]['choices']
        )

    # ----- Begin Analysis Class -----
    # Start Analysis Public Methods
    # Analysis : Public Methods
    def __init__(
            self,
            analysis_types: AnalysisArg,
            name: Optional[str] = None,
            quant_mode: Optional[QuantModeArg] = None,
            ) -> None:
        """Describe in class docstring."""
        if analysis_types is None or not analysis_types:
            message = 'No analysis types provided. Analysis types must be provided'
            message += ' as a list of strings.\nOptions: %s' % ', '.join(self.analysis_options)
            raise HybkitArgError(message)
        if isinstance(analysis_types, str):
            analysis_types = [analysis_types]
        self.analysis_types = []
        for analysis_type in analysis_types:
            if (not isinstance(analysis_type, str)
                    or analysis_type.lower() not in self.analysis_options):
                message = (
                    f'Analysis type "{analysis_type!s}" not recognized.'
                    '\nChoices: {}'.format(', '.join(self.analysis_options))
                )
                raise HybkitArgError(message)
            else:
                self.analysis_types.append(analysis_type.lower())
        if name is not None:
            self.name = self._sanitize_name(name)
        else:
            self.name = name 

        if quant_mode is None:
            self.quant_mode = self.settings['quant_mode']
        elif quant_mode not in self._quant_mode_options:
            message = (
                f'Quantification mode "{quant_mode!s}" not recognized.'
                '\nChoices: {}'.format(', '.join(self._quant_mode_options))
            )
            raise HybkitArgError(message)
        else:
            self.quant_mode = quant_mode

        for analysis_type in self.analysis_types:
            getattr(self, '_init_' + analysis_type)()

    # Analysis : Public Methods : Add HybRecord
    def add_hyb_record(self, hyb_record: hybkit.HybRecord) -> None:
        """
        Add a HybRecord object to the analysis.

        Args:
            hyb_record (:class:`~hybkit.HybRecord`): HybRecord object to be added to the
                analysis.
        """
        for analysis_type in self.analysis_types:
            getattr(self, '_add_' + analysis_type)(hyb_record)

    # Analysis : Public Methods : Add HybRecords
    def add_hyb_records(
            self,
            hyb_records: List[hybkit.HybRecord],
            eval_types: bool = False,
            eval_mirna: bool = False
            ) -> None:
        """
        Add a list of HybRecord objects to the analysis.

        Args:
            hyb_records (:class:`~hybkit.HybFile` or :obj:`list` of :class:`~hybkit.HybRecord`):
                HybFile
                to iterate over, or iterable
                of HybRecord objects to be added to the analysis.
            eval_types (bool): If ``True``, evaluate the hybrid type of the HybRecord before adding
                it to the analysis using :meth:`hybkit.HybRecord.eval_types`.
            eval_mirna (bool): If ``True``, evaluate the miRNA segment of the HybRecord before
                adding it to the analysis using :meth:`hybkit.HybRecord.eval_mirna`.

        """
        for hyb_record in hyb_records:
            if eval_types:
                hyb_record.eval_types()
            if eval_mirna:
                hyb_record.eval_mirna()
            self.add_hyb_record(hyb_record)

    # Start Results Methods
    # Analysis : Public Methods : Results : get_all_results
    def get_all_results(self) -> dict:
        """
        Return a dictionary with all results for all active analyses.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Returns:
            dict: Dictionary with keys of analysis type and values of
                dictionaries with results for that analysis type.
        """
        results = {}
        for analysis_type in self.analysis_types:
            results[analysis_type] = getattr(self, '_get_' + analysis_type + '_results')()
        return results

    # Analysis : Public Methods : Results : get_analysis_result
    def get_analysis_results(self, analysis: AnalysisOptions) -> Dict:
        """
        Return a dictionary with all results for a specific analysis.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            analysis (str): Analysis type to return results for.

        Returns:
            dict: Dictionary with results for the specified analysis type.
                see :ref:Analyses for details.
        """
        self._ensure_analyses_active(analysis)
        return getattr(self, '_get_' + analysis + '_results')()

    # Analysis : Public Methods : Results : get_specific_result
    def get_specific_result(self, result_key: str) -> Any:  # noqa: ANN401
        """
        Get a specific result from the analysis.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            result_key (str): Result key to return from one of the enabled analyses.

        Returns:
            Result value for the specified result key.
        """
        if result_key not in self._all_result_keys_set:
            message = (
                f'Result key "{result_key!s}" not recognized.'
                '\nChoices: {}'.format(', '.join(self._all_result_keys_list))
            )
            raise HybkitArgError(message)
        for analysis_type in self.analysis_options:
            if result_key in self._result_keys[analysis_type]:
                if analysis_type not in self.analysis_types:
                    message = (
                        f'Result "{result_key}" cannot be gotten because analysis '
                        f'type "{analysis_type}" is not active'
                    )
                    raise HybkitArgError(message)
                return getattr(self, '_get_' + analysis_type + '_results')()[result_key]
        raise HybkitArgError('Result key "%s" not found.' % result_key)

    # Analysis : Public Methods : Results : get_analysis_delim_str
    def get_analysis_delim_str(
            self,
            analysis: AnalysisOptions = None,
            out_delim: Optional[str] = None,
            ) -> str:
        """
        Return a delimited string containing the results of the analysis.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            analysis (:obj:`str` or :obj:`list` of :obj:`str`): Analysis type for return results.
                If not provided, return the results for all active analyses.
            out_delim (str): Delimiter to use for output. If not provided, defaults to
                the value in :attr:`settings['out_delim'] <settings>`.
        """
        if analysis is None:
            use_analyses = self.analysis_types
        elif isinstance(analysis, str):
            use_analyses = [analysis]
        elif isinstance(analysis, (list, tuple)):
            use_analyses = analysis
        self._ensure_analyses_active(use_analyses)
        if out_delim is None:
            out_delim = self.settings['out_delim']

        ret_str = ''
        if self.name is not None:
            ret_str += out_delim.join(['analysis_name', self.name]) + '\n'
        for analysis_type in use_analyses:
            ret_str += self._get_analysis_results_delim_str(
                analysis_type, out_delim=out_delim
            )
        return ret_str

    # Start Write Methods
    # Analysis : Public Methods : Write Results : All Analyses
    def write_analysis_delim_str(
            self,
            out_file_name: Optional[str] = None,
            analysis: Optional[AnalysisArg] = None,
            out_delim: Optional[str] = None,
            ) -> None:
        """
        Write the results of the analysis to a delimited text file.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            out_file_name (str): Path to output file. If not provided, defaults to:
                ./<analysis_name>_<analysis>.csv if analysis/analyses provided, or
                ./<analysis_name>_multi_analysis.csv if no analysis/analyses provided.
            analysis (:obj:`str` or :obj:`list` of :obj:`str`): Analysis type for return results.
                If not provided, return the results for all active analyses.
            out_delim (str): Delimiter to use for output. If not provided, defaults to
                the value in :attr:`settings['out_delim'] <settings>`.
        """
        if analysis is None:
            use_analyses = self.analysis_types
            out_suffix = 'multi_analysis'
        elif isinstance(analysis, str):
            use_analyses = [analysis]
            out_suffix = 'analysis'
        elif isinstance(analysis, (list, tuple)):
            use_analyses = analysis
            out_suffix = '-'.join(use_analyses)
        self._ensure_analyses_active(use_analyses)

        if out_file_name is None:
            if self.name is not None:
                out_file_name = self.name + '_' + out_suffix + '.csv'
            else:
                out_file_name = 'analysis_' + out_suffix + '.csv'
        if out_delim is None:
            out_delim = self.settings['out_delim']

        out_delim_str = self.get_analysis_delim_str(analysis, out_delim=out_delim)

        with open(out_file_name, 'w') as out_file:
            out_file.write(out_delim_str)

    # Analysis : Public Methods : Write Results : write_analysis_results_special
    def write_analysis_results_special(
            self,
            out_basename: Optional[str] = None,
            analysis: Optional[AnalysisArg] = None,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        """
        Write the results of the analyses to specialized text files.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            out_basename (str): Path for basename of output file. Files will be renamed
                using the provided path as the base name. If not provided, defaults to:
                ./<analysis_name>_<analysis> if :attr:`name` is set, or
                ./Analysis_multi_<analysis>  if name not set.
            analysis (:obj:`str` or :obj:`list` of :obj:`str`): Analysis type to write results
                files for.
                If not provided, write results files for all active analyses.
            out_delim (str): Delimiter to use for output where applicable.
                If not provided, defaults to
                the value in :attr:`settings['out_delim'] <settings>`.
        """
        if analysis is None:
            use_analyses = self.analysis_types
        elif isinstance(analysis, str):
            use_analyses = [analysis]
        elif isinstance(analysis, (list, tuple)):
            use_analyses = analysis
        self._ensure_analyses_active(use_analyses)

        if out_basename is None:
            if self.name is not None:
                out_basename = self.name
            else:
                out_basename = 'analysis'

        if out_delim is None:
            out_delim = self.settings['out_delim']

        all_out_files = []
        for analysis in use_analyses:
            out_files = getattr(self, '_write_' + analysis + '_results_special')(
                basename=out_basename,
                out_delim=out_delim
            )
            all_out_files += out_files
        return all_out_files

    # Analysis : Public Methods : Plot Results : plot_analysis_results
    def plot_analysis_results(
            self,
            out_basename: Optional[str] = None,
            analysis: Optional[AnalysisArg] = None,
            ) -> List[str]:
        """
        Plot the results of the analyses.

        See :ref:`Analyses <Analyses>` for details on the results for each analysis type.

        Args:
            analysis (:obj:`str` or :obj:`list` of :obj:`str`): Analysis type to plot results for.
                If not provided, plot results for all active analyses.
            out_basename (str): Path to output file. If not provided, defaults to:
                ./<analysis_name> if :attr:`name` provided or
                ./analysis if no name provided.
        """
        if analysis is None:
            use_analyses = self.analysis_types
        elif isinstance(analysis, str):
            use_analyses = [analysis]
        elif isinstance(analysis, (list, tuple)):
            use_analyses = analysis
        self._ensure_analyses_active(use_analyses)

        if out_basename is None:
            if self.name is not None:
                out_basename = self.name
            else:
                out_basename = 'analysis'

        all_out_files = []
        for analysis in use_analyses:
            if hasattr(self, '_plot_' + analysis + '_results'):
                out_files = getattr(self, '_plot_' + analysis + '_results')(
                    basename=out_basename
                )
                all_out_files += out_files
        return all_out_files

    # Start Helper Methods
    # Analysis : Private Methods : Helper Methods
    def _get_quant(self, hyb_record: hybkit.HybRecord) -> int:
        if self.quant_mode == 'single':
            return 1
        elif self.quant_mode == 'reads':
            return hyb_record.get_read_count(require=True)
        elif self.quant_mode == 'records':
            return hyb_record.get_record_count(require=True)
        else:
            raise HybkitArgError('Quantification mode "%s" not recognized.' % self.quant_mode)

    # Start Init Methods
    # Analysis : Private Methods : Init Methods : Energy Analysis
    def _init_energy(self) -> None:
        self._energy_analysis_count = 0
        self._has_energy_val = 0
        self._no_energy_val = 0
        self._energy_vals = np.array([], dtype=float)
        self._binned_energy_vals = Counter()
        for i in range(0, -31, -1):
            self._binned_energy_vals[i] = 0

    # Analysis : Private Methods : Init Methods : Type Analysis
    def _init_type(self) -> None:
        self._types_analysis_count = 0
        self._hybrid_types = Counter()
        self._reordered_hybrid_types = Counter()
        self._mirna_hybrid_types = Counter()
        self._seg1_types = Counter()
        self._seg2_types = Counter()
        self._all_seg_types = Counter()

    # Analysis : Private Methods : Init Methods : miRNA Analysis
    def _init_mirna(self) -> None:
        self._mirna_analysis_count = 0
        self._mirnas_5p = 0
        self._mirnas_3p = 0
        self._mirna_dimers = 0
        self._non_mirna = 0
        self._has_mirna = 0

    # Analysis : Private Methods : Init Methods : Target Analysis
    def _init_target(self) -> None:
        self._target_analysis_count = 0
        self._target_evals = 0
        self._target_names = Counter()
        self._target_types = Counter()

    # Analysis : Private Methods : Init Methods : Fold Analysis
    def _init_fold(self) -> None:
        self._fold_analysis_count = 0
        self._folds_recorded = 0
        self._mirna_nt_fold_counts = Counter()
        self._fold_match_counts = Counter()

    # Start Add Methods
    # Analysis : Private Methods : Add Methods : Energy Analysis
    def _add_energy(self, hyb_record: hybkit.HybRecord) -> None:
        count = self._get_quant(hyb_record)
        self._energy_analysis_count += 1
        if hyb_record.energy is not None:
            self._has_energy_val += count
            self._energy_vals = np.append(self._energy_vals, float(hyb_record.energy))
            if hyb_record.energy is not None:
                energy_bin = int(np.ceil(float(hyb_record.energy)))
            self._binned_energy_vals[energy_bin] += count
        else:
            self._no_energy_val += count

    # Analysis : Private Methods : Add Methods : Type Analysis
    def _add_type(self, hyb_record: hybkit.HybRecord) -> None:
        hyb_record._ensure_set('eval_types')
        count = self._get_quant(hyb_record)
        self._types_analysis_count += 1
        seg1_type = hyb_record.get_seg1_type()
        seg2_type = hyb_record.get_seg2_type()
        hybrid_type = (seg1_type, seg2_type)
        reordered_hybrid_type = tuple(sorted([seg1_type, seg2_type]))
        if hyb_record.is_set('eval_mirna') and hyb_record.prop('has_mirna'):
            if hyb_record.prop('5p_mirna'):
                mirna_hybrid_type = (seg1_type, seg2_type)
            else:
                mirna_hybrid_type = (seg2_type, seg1_type)
        else:
            mirna_hybrid_type = reordered_hybrid_type
        self._hybrid_types[hybrid_type] += count
        self._reordered_hybrid_types[reordered_hybrid_type] += count
        self._mirna_hybrid_types[mirna_hybrid_type] += count
        self._seg1_types[seg1_type] += count
        self._seg2_types[seg2_type] += count
        self._all_seg_types[seg1_type] += count
        self._all_seg_types[seg2_type] += count

    # Analysis : Private Methods : Add Methods : miRNA Analysis
    def _add_mirna(self, hyb_record: hybkit.HybRecord) -> None:
        hyb_record._ensure_set('eval_mirna')
        count = self._get_quant(hyb_record)
        self._mirna_analysis_count += 1
        if hyb_record.prop('has_mirna'):
            self._has_mirna += count
            if hyb_record.prop('mirna_dimer'):
                self._mirna_dimers += count
            elif hyb_record.prop('5p_mirna'):
                self._mirnas_5p += count
            elif hyb_record.prop('3p_mirna'):
                self._mirnas_3p += count
        else:
            self._non_mirna += count

    # Analysis : Private Methods : Add Methods : Target Analysis
    def _add_target(self, hyb_record: hybkit.HybRecord) -> None:
        hyb_record._ensure_set('eval_mirna')
        self._target_analysis_count += 1
        count = self._get_quant(hyb_record)
        if hyb_record.prop('has_mirna'):
            self._target_evals += 1
            mirna_details = hyb_record.mirna_details(allow_mirna_dimers=True)
            self._target_names[mirna_details['target_ref']] += count
            self._target_types[mirna_details['target_seg_type']] += count

    # Analysis : Private Methods : Add Methods : Fold Analysis
    def _add_fold(self, hyb_record: hybkit.HybRecord) -> None:
        hyb_record._ensure_set('fold_record')
        hyb_record._ensure_set('eval_mirna')
        self._fold_analysis_count += 1
        count = self._get_quant(hyb_record)
        if hyb_record.prop('has_mirna'):
            self._folds_recorded += count
            mirna_fold = hyb_record.mirna_details('mirna_fold')
            match_count = 0
            for pos_i, nt_i in enumerate(range(len(mirna_fold)), start=1):
                if mirna_fold[nt_i] in {'(', ')'}:
                    self._mirna_nt_fold_counts[pos_i] += count
                    match_count += 1
        self._fold_match_counts[match_count] += count

    # Start Get Results Methods
    # Analysis : Private Methods : Get Methods : Energy Analysis
    def _get_energy_results(self) -> dict:
        energy_results = {}
        energy_results['energy_analysis_count'] = copy.deepcopy(self._energy_analysis_count)
        energy_results['has_energy_val'] = copy.deepcopy(self._has_energy_val)
        energy_results['no_energy_val'] = copy.deepcopy(self._no_energy_val)
        if energy_results['has_energy_val'] > 0:
            energy_results['energy_min'] = self._energy_vals.min()
            energy_results['energy_max'] = self._energy_vals.max()
            energy_results['energy_mean'] = self._energy_vals.mean()
            energy_results['energy_std'] = self._energy_vals.std()
        else:
            energy_results['energy_min'] = None
            energy_results['energy_max'] = None
            energy_results['energy_mean'] = None
            energy_results['energy_std'] = None
        ret_vals = Counter()
        if energy_results['energy_max'] is None:
            ret_range_max = 0
        else:
            ret_range_max = max(0, int(np.ceil(energy_results['energy_max'])))
        if energy_results['energy_min'] is None:
            ret_range_min = -1
        else:
            ret_range_min = int(np.floor(energy_results['energy_min'])) - 1
        for i in range(ret_range_max, ret_range_min, -1):
            if i in self._binned_energy_vals:
                ret_vals[i] = self._binned_energy_vals[i]
            else:
                ret_vals[i] = 0
        # Conditionally trim artifact of floor/ceil:
        if (ret_vals[ret_range_min + 1]) == 0:
            del ret_vals[ret_range_min + 1]
        energy_results['binned_energy_vals'] = ret_vals
        return energy_results

    # Analysis : Private Methods : Get Methods : Type Analysis
    def _get_type_results(self) -> dict:
        type_results = {}
        type_results['types_analysis_count'] = copy.deepcopy(self._types_analysis_count)
        type_results['hybrid_types'] = copy.deepcopy(self._hybrid_types)
        type_results['reordered_hybrid_types'] = copy.deepcopy(self._reordered_hybrid_types)
        type_results['mirna_hybrid_types'] = copy.deepcopy(self._mirna_hybrid_types)
        type_results['seg1_types'] = copy.deepcopy(self._seg1_types)
        type_results['seg2_types'] = copy.deepcopy(self._seg2_types)
        type_results['all_seg_types'] = copy.deepcopy(self._all_seg_types)
        return type_results

    # Analysis : Private Methods : Get Methods : miRNA Analysis
    def _get_mirna_results(self) -> dict:
        mirna_results = {}
        mirna_results['mirna_analysis_count'] = copy.deepcopy(self._mirna_analysis_count)
        mirna_results['has_mirna'] = copy.deepcopy(self._has_mirna)
        mirna_results['non_mirna'] = copy.deepcopy(self._non_mirna)
        mirna_results['mirnas_5p'] = copy.deepcopy(self._mirnas_5p)
        mirna_results['mirnas_3p'] = copy.deepcopy(self._mirnas_3p)
        mirna_results['mirna_dimers'] = copy.deepcopy(self._mirna_dimers)
        return mirna_results

    # Analysis : Private Methods : Get Methods : Target Analysis
    def _get_target_results(self) -> dict:
        target_results = {}
        target_results['target_analysis_count'] = copy.deepcopy(self._target_analysis_count)
        target_results['target_evals'] = copy.deepcopy(self._target_evals)
        target_results['target_names'] = copy.deepcopy(self._target_names)
        target_results['target_types'] = copy.deepcopy(self._target_types)
        return target_results

    # Analysis : Private Methods : Get Methods : Fold Analysis
    def _get_fold_results(self) -> dict:
        fold_results = {}
        fold_results['fold_analysis_count'] = copy.deepcopy(self._fold_analysis_count)
        fold_results['folds_recorded'] = copy.deepcopy(self._folds_recorded)
        fold_results['fold_match_counts'] = copy.deepcopy(self._fold_match_counts)
        fold_results['mirna_nt_fold_counts'] = copy.deepcopy(self._mirna_nt_fold_counts)
        fold_results['mirna_nt_fold_props'] = {
            k: (v / self._folds_recorded) for k, v in self._mirna_nt_fold_counts.items()
        }
        return fold_results

    # Start Get Results String Methods
    # Analysis : Private Methods : Result String Methods : All Analyses
    def _get_analysis_results_delim_str(
            self,
            analysis: AnalysisOptions,
            out_delim: Optional[str] = None,
            ) -> str:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        analysis_results = getattr(self, f'_get_{analysis}_results')()

        analysis_results_str = ''
        for key, val in analysis_results.items():
            if isinstance(val, dict):
                for subkey, subval in val.items():
                    if isinstance(subkey, tuple):
                        use_subkey = '--'.join(subkey)
                    else:
                        use_subkey = str(subkey)
                    if isinstance(subval, float):
                        use_subval = f'{subval:.3f}'
                    else:
                        use_subval = str(subval)
                    analysis_results_str += out_delim.join([key, use_subkey, use_subval]) + '\n'
            else:
                if isinstance(val, float):
                    use_val = f'{val:.3f}'
                else:
                    use_val = str(val)
                analysis_results_str += out_delim.join([key, use_val]) + '\n'
        return analysis_results_str

    # Start Write Methods
    # Analysis : Private Methods : Result Special Writing Methods : Energy Analysis
    def _write_energy_results_special(
            self,
            basename: str,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        energy_results = self._get_energy_results()

        out_file_names = []
        main_result_keys = [
            'energy_analysis_count', 'has_energy_val', 'no_energy_val',
            'energy_min', 'energy_max', 'energy_mean', 'energy_std',
        ]
        main_result_file_name = basename + '_energy_main_results.csv'
        with open(main_result_file_name, 'w') as main_result_file:
            for key in main_result_keys:
                val = energy_results[key]
                main_result_file.write(out_delim.join([key, str(val)]) + '\n')
        out_file_names.append(main_result_file_name)

        binned_energy_vals_file_name = basename + '_energy_binned_vals.csv'
        with open(binned_energy_vals_file_name, 'w') as binned_energy_vals_file:
            for val, count in energy_results['binned_energy_vals'].items():
                binned_energy_vals_file.write(out_delim.join([str(val), str(count)]) + '\n')
        out_file_names.append(binned_energy_vals_file_name)
        return out_file_names

    # Analysis : Private Methods : Result Special Writing Methods : Type Analysis
    def _write_type_results_special(
            self,
            basename: str,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        type_results = self._get_type_results()

        out_file_names = []
        name_keys = {
            'hybrid_types': 'type_hybrid_types',
            'reordered_hybrid_types': 'type_reordered_hybrid_types',
            'mirna_hybrid_types': 'type_mirna_hybrid_types',
            'seg1_types': 'type_seg1_types',
            'seg2_types': 'type_seg2_types',
            'all_seg_types': 'type_all_seg_types',
        }
        for name_key, name_val in name_keys.items():
            out_file_name = basename + '_' + name_val + '.csv'
            val_dict = type_results[name_key]
            if hasattr(val_dict, 'most_common'):
                iter_method = 'most_common'
            else:
                iter_method = 'items'

            with open(out_file_name, 'w') as out_file:
                for info_key, info_val in getattr(val_dict, iter_method)():
                    if isinstance(info_key, tuple):
                        use_info_key = '--'.join(info_key)
                    else:
                        use_info_key = str(info_key)
                    out_file.write(out_delim.join([use_info_key, str(info_val)]) + '\n')
            out_file_names.append(out_file_name)

        return out_file_names

    # Analysis : Private Methods : Result Special Writing Methods : miRNA Analysis
    def _write_mirna_results_special(
            self,
            basename: str,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        mirna_results = self._get_mirna_results()

        out_file_names = []
        main_keys = [
            'mirna_analysis_count', 'has_mirna', 'non_mirna',
            'mirnas_5p', 'mirnas_3p', 'mirna_dimers',
        ]
        main_outfile_name = basename + '_mirna_results.csv'
        with open(main_outfile_name, 'w') as out_file:
            for key in main_keys:
                val = mirna_results[key]
                out_file.write(out_delim.join([key, str(val)]) + '\n')
        out_file_names.append(main_outfile_name)

        return out_file_names

    # Analysis : Private Methods : Result Special Writing Methods : Target Analysis
    def _write_target_results_special(
            self,
            basename: str,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        target_results = self._get_target_results()

        out_file_names = []
        main_keys = [
            'target_analysis_count', 'target_evals'
        ]
        main_outfile_name = basename + '_target_results.csv'
        with open(main_outfile_name, 'w') as out_file:
            for key in main_keys:
                val = target_results[key]
                out_file.write(out_delim.join([key, str(val)]) + '\n')
        out_file_names.append(main_outfile_name)

        names_outfile_name = basename + '_target_names.csv'
        with open(names_outfile_name, 'w') as out_file:
            for name, count in target_results['target_names'].most_common():
                out_file.write(out_delim.join([name, str(count)]) + '\n')
        out_file_names.append(names_outfile_name)

        types_outfile_name = basename + '_target_types.csv'
        with open(types_outfile_name, 'w') as out_file:
            for name, count in target_results['target_types'].most_common():
                out_file.write(out_delim.join([name, str(count)]) + '\n')
        out_file_names.append(types_outfile_name)

        return out_file_names

    # Analysis : Private Methods : Result Special Writing Methods : Fold Analysis
    def _write_fold_results_special(
            self,
            basename: str,
            out_delim: Optional[str] = None,
            ) -> List[str]:
        if out_delim is None:
            out_delim = self.settings['out_delim']

        fold_results = self._get_fold_results()

        out_file_names = []
        name_keys = {
            'folds_recorded': 'fold_main_results',
            'fold_match_counts': 'fold_fold_match_counts',
            'mirna_nt_fold_counts': 'fold_mirna_nt_fold_counts',
            'mirna_nt_fold_props': 'fold_mirna_nt_fold_props',
        }
        for name_key, name_val in name_keys.items():
            out_file_name = basename + '_' + name_val + '.csv'
            if name_key == 'folds_recorded':
                val_dict = {name_key: fold_results[name_key]}
                iter_method = 'items'
            elif name_key in {'fold_match_counts', 'mirna_nt_fold_counts', 'mirna_nt_fold_props'}:
                val_dict = {k: fold_results[name_key][k]
                            for k in sorted([*fold_results[name_key].keys()])}
                iter_method = 'items'
            else:
                val_dict = fold_results[name_key]
                if hasattr(val_dict, 'most_common'):
                    iter_method = 'most_common'
                else:
                    iter_method = 'items'

            with open(out_file_name, 'w') as out_file:
                for info_key, info_val in getattr(val_dict, iter_method)():
                    out_file.write(out_delim.join([str(info_key), str(info_val)]) + '\n')
            out_file_names.append(out_file_name)

        return out_file_names

    # Start Plot Methods
    # Analysis : Private Methods : Result Plotting Methods : Energy Analysis
    def _plot_energy_results(self, basename: str) -> List[str]:
        energy_results = self._get_energy_results()
        out_files = []
        # Plot Histogram
        energy_histogram_file_name = basename + '_energy_histogram.png'
        hybkit.plot.energy_histogram(
            results=energy_results['binned_energy_vals'],
            plot_file_name=energy_histogram_file_name,
            title='Hybrid Gibbs Free Energy (kcal/mol)',
            name=self.name,
        )
        out_files.append(energy_histogram_file_name)
        return out_files

    # Analysis : Private Methods : Result Plotting Methods : Type Analysis
    def _plot_type_results(self, basename: str) -> List[str]:
        type_results = self._get_type_results()
        out_files = []
        type_results = {}
        type_results['types_analysis_count'] = copy.deepcopy(self._types_analysis_count)
        type_results['hybrid_types'] = copy.deepcopy(self._hybrid_types)
        type_results['reordered_hybrid_types'] = copy.deepcopy(self._reordered_hybrid_types)
        type_results['mirna_hybrid_types'] = copy.deepcopy(self._mirna_hybrid_types)
        type_results['seg1_types'] = copy.deepcopy(self._seg1_types)
        type_results['seg2_types'] = copy.deepcopy(self._seg2_types)
        type_results['all_seg_types'] = copy.deepcopy(self._all_seg_types)

        # Plot hybrid_types
        hybrid_types_file_name = basename + '_types_hybrid_types.png'
        hybkit.plot.type_count_dual(
            type_results['hybrid_types'],
            hybrid_types_file_name,
            title='Hybrid Types',
            name=self.name,
            join_entries=True,
        )
        out_files.append(hybrid_types_file_name)

        # Plot reordered_hybrid_types
        reordered_hybrid_types_file_name = basename + '_types_reordered_hybrid_types.png'
        hybkit.plot.type_count_dual(
            type_results['reordered_hybrid_types'],
            reordered_hybrid_types_file_name,
            title='Reordered Hybrid Types',
            name=self.name,
            join_entries=True,
        )
        out_files.append(reordered_hybrid_types_file_name)

        # Plot mirna_hybrid_types
        mirna_hybrid_types_file_name = basename + '_types_mirna_hybrids.png'
        hybkit.plot.type_count_dual(
            type_results['mirna_hybrid_types'],
            mirna_hybrid_types_file_name,
            title="Reordered miRNA' Hybrid Types",
            name=self.name,
            join_entries=True,
        )
        out_files.append(mirna_hybrid_types_file_name)

        # Plot All Seg Types
        all_seg_types_file_name = basename + '_types_all_seg.png'
        hybkit.plot.type_count(
            type_results['all_seg_types'],
            all_seg_types_file_name,
            title='All Segment Types',
            name=self.name,
        )
        out_files.append(all_seg_types_file_name)

        # Plot Seg1 Types
        seg1_types_file_name = basename + '_types_seg1.png'
        hybkit.plot.type_count(
            type_results['seg1_types'],
            seg1_types_file_name,
            title='5p Segment Types',
            name=self.name,
        )
        out_files.append(seg1_types_file_name)

        # Plot Seg2 Types
        seg2_types_file_name = basename + '_types_seg2.png'
        hybkit.plot.type_count(
            type_results['seg2_types'],
            seg2_types_file_name,
            title='3p Segment Types',
            name=self.name,
        )
        out_files.append(seg2_types_file_name)
        return out_files

    # Analysis : Private Methods : Result Plotting Methods : Target Analysis
    def _plot_target_results(self, basename: str) -> List[str]:
        target_results = self._get_target_results()
        out_files = []

        # Plot target_names
        target_names_file_name = basename + '_target_names.png'
        hybkit.plot.target_count(
            target_results['target_names'],
            target_names_file_name,
            title='miRNA Target Names',
            name=self.name,
        )
        out_files.append(target_names_file_name)

        # Plot target_types
        target_types_file_name = basename + '_target_types.png'
        hybkit.plot.type_count(
            target_results['target_types'],
            target_types_file_name,
            title='miRNA Target Types',
            name=self.name,
        )
        out_files.append(target_types_file_name)
        return out_files

    # Analysis : Private Methods : Result Plotting Methods : Fold Analysis
    def _plot_fold_results(self, basename: str) -> List[str]:
        fold_results = self._get_fold_results()
        out_files = []

        # Plot Match Counts Histogram
        match_counts_histogram_file_name = basename + '_fold_match_counts_histogram.png'
        hybkit.plot.fold_match_counts_histogram(
            results=fold_results['fold_match_counts'],
            plot_file_name=match_counts_histogram_file_name,
            title='Predicted Fold Match Count',
            name=self.name,
        )
        out_files.append(match_counts_histogram_file_name)

        # Plot miRNA nt Fold Counts Histogram
        mirna_nt_fold_counts_histogram_file_name = basename + '_fold_mirna_nt_counts_histogram.png'
        hybkit.plot.fold_mirna_nt_counts_histogram(
            results=fold_results['mirna_nt_fold_counts'],
            plot_file_name=mirna_nt_fold_counts_histogram_file_name,
            title='Predicted Matches per miRNA Nucleotide',
            name=self.name,
        )
        out_files.append(mirna_nt_fold_counts_histogram_file_name)

        # Plot miRNA nt Fold Props Histogram
        mirna_nt_fold_props_histogram_file_name = basename + '_fold_mirna_nt_props_histogram.png'
        hybkit.plot.fold_mirna_nt_counts_histogram(
            results=fold_results['mirna_nt_fold_props'],
            plot_file_name=mirna_nt_fold_props_histogram_file_name,
            title='Predicted Match Proportion per miRNA Nucleotide',
            is_prop=True,
            name=self.name,
        )
        out_files.append(mirna_nt_fold_props_histogram_file_name)
        return out_files

    # Analysis : Private Methods : Utility Methods
    def _ensure_analyses_active(self, analyses: AnalysisArg) -> None:
        if isinstance(analyses, str):
            all_test_analyses = [analyses]
        elif isinstance(analyses, list):
            all_test_analyses = analyses
        else:
            message = 'Analyses must be a string or list of strings.'
            raise TypeError(message)
        for test_analysis in all_test_analyses:
            if test_analysis not in self.analysis_types:
                message = (
                    f'Analysis type "{test_analysis!s}" not an active analysis.'
                    '\nActive choices: {}'.format(', '.join(self.analysis_types))
                )
                raise HybkitArgError(message)

    # Analysis : Private Classmethods
    @classmethod
    def _sanitize_name(cls, file_name: str) -> str:
        for char, replace in [('*', 'star'), (',', 'com')]:
            file_name = file_name.replace(char, replace)
        return file_name
