#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""This module contains settings information for hybkit classes and methods."""

import copy
from typing import Dict, List


# ----- Begin Settings Helper Functions -----
# Util : Settings Helper Functions
def _all_str_cases(in_str: str) -> List[str]:
    return [in_str.lower(), in_str.title(), in_str.upper()]


# Util : Settings Helper Functions
def _settings_info_to_settings(settings_info_dict: Dict) -> dict:
    ret_dict = {}
    for key in settings_info_dict:
        ret_dict[key] = copy.deepcopy(settings_info_dict[key][0])
    return ret_dict


# ----- Begin Settings Constants -----
# Util : Global Variables
#: Allowed suffixes for "Hyb" files.
HYB_SUFFIXES = _all_str_cases('.hyb')

#: Allowed suffixes for "Vienna" files.
VIENNA_SUFFIXES = _all_str_cases('.vienna')

#: Allowed suffixes for "Connection-Table" files.
CT_SUFFIXES = _all_str_cases('.ct')

#: Allowed suffixes for "Vienna" and "Connection-Table" files.
FOLD_SUFFIXES = VIENNA_SUFFIXES + CT_SUFFIXES

_USE_ABSPATH = False

_FILTER_OUT_SUFFIX = '_filtered'
_EVAL_OUT_SUFFIX = '_evaluated'

#: Default miRNA types for use in :func:`mirna_analysis`.
MIRNA_TYPES = ['miRNA', 'microRNA']

#: Analysis types for use by :class:`hybkit.Analysis`.
ANALYSIS_TYPE_OPTIONS_HYB = ['energy', 'type', 'mirna', 'target']
ANALYSIS_TYPE_OPTIONS_FOLD = ['fold']
ANALYSIS_TYPE_OPTIONS = [*ANALYSIS_TYPE_OPTIONS_HYB, *ANALYSIS_TYPE_OPTIONS_FOLD]

#: Minimum number of fields in hyb line:
MIN_RECORD_FIELDS = 15

#: Maximum number of fields in hyb line:
MAX_RECORD_FIELDS = 16

# #: Default minimum Gibbs Free Energy for bins in :class:`EnergyAnalysis`
# #: (range: ENERGY_MIN_BIN <= 0).
# ENERGY_MIN_BIN = '-45.0'

# #: Default Gibbs Free Energy bin size for use by :class:`EnergyAnalysis` ('0.1' == no binning).
# ENERGY_BIN_SIZE = '0.4'

# ----- Begin settings_info variables -----
# Start settings_info : HybRecord
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.HybRecord` class.
#: Copied into :data:`HybRecord_settings` for use at runtime.
HybRecord_settings_info = {
    'mirna_types': [
        copy.deepcopy(MIRNA_TYPES),
        '"seg_type" fields identifying a miRNA',
        'str',
        None,
        {'nargs': '+'}
    ],
    'custom_flags': [
        [],
        """
        Custom flags to allow in addition to those specified in the hybkit
        specification.
        """,
        'str',
        None,
        {'nargs': '+'}
    ],
    'hyb_placeholder': [
        '.',
        """
        placeholder character/string for missing data in hyb files.
        """,
        'str',
        None,
        {}
    ],
    'reorder_flags': [
        True,
        """
        Re-order flags to the hybkit-specification order when
        writing hyb records.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'allow_undefined_flags': [
        False,
        """
        Allow use of flags not defined in the hybkit-specification order when
        reading and writing hyb records. As the preferred alternative to
        using this setting,
        the --custom_flags argument can be be used to supply custom allowed flags.
        """,
        'custom_bool_from_str',
        None,
        {'nargs': '?', 'const': True}
    ],
    'allow_unknown_seg_types': [
        False,
        """
        Allow unknown segment types when assigning segment types.
        """,
        'custom_bool_from_str',
        None,
        {'nargs': '?', 'const': True}
    ],
    # 'check_complete_seg_types': [
    #     False,
    #     """
    #     Check every segment possibility when assigning segment types, rather than
    #     breaking after the first match is found. If True, finding segment types
    #     is slower but better at catching errors.
    #     """,
    #     'custom_bool_from_str',
    #     None,
    #     {'nargs': '?', 'const': True}
    # ],
}

# Start settings_info : HybFile
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.HybFile` class.
#: Copied into :data:`HybFile_settings` for use at runtime.
HybFile_settings_info = {
    'hybformat_id': [
        False,
        """
        The Hyb Software Package places further information in the "id" field of the
        hybrid record that can be used to infer the number of contained read counts.
        When set to True, the identifiers will be parsed as: "<read_id>_<read_count>"
        """,
        'custom_bool_from_str',
        None,
        {'nargs': '?', 'const': True}
    ],
    'hybformat_ref': [
        False,
        """
        The Hyb Software Package uses a reference database with identifiers
        that contain sequence type and other sequence information.
        When set to True, all hyb file identifiers will be parsed as:
        "<gene_id>_<transcript_id>_<gene_name>_<seg_type>"
        """,
        'custom_bool_from_str',
        None,
        {'nargs': '?', 'const': True}
    ],
}

# Start settings_info : FoldRecord
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.FoldRecord` class.
#: Copied into :data:`FoldRecord_settings` for use at runtime.
FoldRecord_settings_info = {
    'allowed_mismatches': [
        0,
        """
        For DynamicFoldRecords, allowed number of mismatches with a HybRecord.
        """,
        'int',
        None,
        {}
    ],
    'fold_placeholder': [
        '.',
        """
        Placeholder character/string for missing data for reading/writing fold records.
        """,
        'str',
        None,
        {}
    ],
    'seq_type': [
        'static',
        """
        Type of fold record object to use. Options:
        "static": FoldRecord, requires an exact sequence
        match to be paired with a HybRecord; "dynamic": DynamicFoldRecord, requires a sequence
        match to the "dynamic" annotated regions of a HybRecord, and may be shorter/longer
        than the original sequence.
        """,
        'str',
        '-y',
        {'choices': ['static', 'dynamic']}
    ],
    'error_mode': [
        'raise',
        """
        Mode for handling errors during reading of HybFiles
        (overridden by HybFoldIter.settings['iter_error_mode'] when using HybFoldIter).
        Options: "raise": Raise an error when encountered and exit program ;
        "warn_return": Print a warning and return the error_value ;
        "return": Return the error value with no program output.
        record is encountered.
        """,
        'str',
        None,
        {'choices': ['raise', 'warn_return', 'return']}
    ],
}

# Start settings_info : FoldFile
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.FoldFile` class.
#: Copied into :data:`FoldFile_settings` for use at runtime.
FoldFile_settings_info = {
}

# Start settings_info : HybFoldIter
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.HybFoldIter` class.
#: Copied into :data:`HybFoldIter_settings` for use at runtime.
HybFoldIter_settings_info = {
    'error_checks': [
        ['hybrecord_indel', 'foldrecord_nofold', 'max_mismatch', 'energy_mismatch'],
        """
        Error checks for simultaneous HybFile and FoldFile parsing.
        Options: "hybrecord_indel": Error for HybRecord objects where one/both sequences have
        insertions/deletions in alignment, which prevents matching of sequences;
        "foldrecord_nofold": Error when failure in reading a fold_record object;
        "max_mismatch": Error when mismatch between hybrecord and foldrecord sequences is
        greater than FoldRecord "allowed_mismatches" setting; "energy_mismatch": Error when
        a mismatch exists between HybRecord and FoldRecord energy values.
        """,
        'str',
        None,
        {'choices': ['hybrecord_indel', 'foldrecord_nofold', 'max_mismatch', 'energy_mismatch']}
    ],
    'iter_error_mode': [
        'warn_skip',
        """
        Mode for handling errors found during error checks.
        Overrides HybRecord "error_mode" setting when using HybFoldIter.
        Options:
        "raise": Raise an error when encountered;
        "warn_return": Print a warning and return the value;
        "warn_skip": Print a warning and continue to the next iteration;
        "skip": Continue to the next iteration without any output;
        "return": return the value without any error output;
        """,
        'str',
        None,
        {'choices': ['raise', 'warn_return', 'warn_skip', 'skip', 'return']}
    ],
    'max_sequential_skips': [
        100,
        """
        Maximum number of record(-pairs) to skip in a row. Limited as several sequential skips
        usually indicates an issue with record formatting or a desynchronization between files.
        """,
        'int',
        None,
        {}
    ],
}

# Start settings_info : Analysis
# setting_info format contains structure::
# {
#     setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
# }
#: Information for settings of :class:`~hybkit.Analysis` class.
#: Copied into :data:`Analysis_settings` for use at runtime.
Analysis_settings_info = {
    'quant_mode': [
        'single',
        """
        Method for counting records. Options:
        "single": Count each record as a single entry;
        "reads": Use the number of reads per hyb record as the count (may contain PCR duplicates);
        "records": Count the number of records represented by each
        hyb record entry (1 for "unmerged" records, >= 1 for "merged" records)
        """,
        'str',
        None,
        {'choices': ['single', 'reads', 'records']}
    ],
    'out_delim': [
        ',',
        """
        Delimiter-string to place between fields in analysis output.
        """,
        'str',
        None,
        {}
    ],
    # 'mirna_sort': [
    #     True,
    #     """
    #     During TypeAnalysis, sort miRNAs first for "miRNA"-"Other" segtype pairs.
    #     If False, sort alphabetically.
    #     """,
    #     'custom_bool_from_str',
    #     None,
    #     {}
    # ],
    # 'allow_mirna_dimers': [
    #     False,
    #     """
    #     Include miRNA / miRNA dimers in TargetAnalysis.
    #     If False, exclude these from analysis results.
    #     """,
    #     'custom_bool_from_str',
    #     None,
    #     {'nargs': '?', 'const': True}
    # ],
    # 'type_sep': [
    #     '-',
    #     """
    #     Separator-string to place between types in analysis output.
    #     """,
    #     'str',
    #     None,
    #     {}
    # ],
    # 'energy_min_bin': [
    #     ENERGY_MIN_BIN,
    #     """
    #     Minimum Gibbs Free Energy value for binned analysis.
    #     """,
    #     'str',
    #     None,
    #     {}
    # ],
    # 'energy_bin_size': [
    #     ENERGY_BIN_SIZE,
    #     """
    #     Size of increment to bin energy values for binned energy analysis
    #     (allowed >= 0.1). A value of '0.1' represents no value binning.
    #     """,
    #     'str',
    #     None,
    #     {}
    # ],

}

# Settings (active) : HybRecord
#: Settings for :class:`~hybkit.HybRecord`,
#: created from :data:`HybRecord_settings_info`
HybRecord_settings = _settings_info_to_settings(HybRecord_settings_info)

# Settings (active) : HybFile
#: Settings for :class:`~hybkit.HybFile`,
#: created from :data:`HybFile_settings_info`
HybFile_settings = _settings_info_to_settings(HybFile_settings_info)

# Settings (active) : FoldRecord
#: Settings for :class:`~hybkit.FoldRecord`,
#: created from :data:`FoldRecord_settings_info`
FoldRecord_settings = _settings_info_to_settings(FoldRecord_settings_info)

# Settings (active) : FoldFile
#: Settings for :class:`~hybkit.FoldFile`,
#: created from :data:`FoldFile_settings_info`
FoldFile_settings = _settings_info_to_settings(FoldFile_settings_info)

# Settings (active) : HybFoldIter
#: Settings for :class:`~hybkit.HybFoldIter`,
#: created from :data:`HybFoldIter_settings_info`
HybFoldIter_settings = _settings_info_to_settings(HybFoldIter_settings_info)

# Settings (active) : Analysis
#: Settings for :class:`~hybkit.analysis.BaseAnalysis`,
#: created from :data:`Analysis_settings_info`
Analysis_settings = _settings_info_to_settings(Analysis_settings_info)
