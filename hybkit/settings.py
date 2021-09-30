#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""This module contains settings information for hybkit classes and methods."""

import copy

# Util : Settings Helper Functions
def _all_str_cases(in_str):
    return [in_str.lower(), in_str.title(), in_str.upper()]

# Util : Settings Helper Functions
def _settings_info_to_settings(settings_info_dict):
    ret_dict = {}
    for key in settings_info_dict:
        ret_dict[key] = copy.deepcopy(settings_info_dict[key][0]) 
    return ret_dict

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
_ANALYSIS_OUT_SUFFIX = '_analyzed'

#: Default miRNA types for use in :func:`mirna_analysis`.
MIRNA_TYPES = ['miRNA', 'microRNA']

#: Default coding sequence types for use in the :func:`target_region_analysis`.
CODING_TYPES = ['mRNA']

# settings_info : HybRecord
#: Information for settings of HybRecord class.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
HybRecord_settings_info = {
    'mirna_types': [ 
        copy.deepcopy(MIRNA_TYPES),
        '"seg_type" fields identifying a miRNA',
        'str',
        None,
        {'nargs':'+'}
    ],
    'coding_types': [ 
        copy.deepcopy(CODING_TYPES),
        '"seg_type" fields identifying a coding sequence',
        'str',
        None,
        {'nargs':'+'}
    ],
    'custom_flags': [
        [],
        """
        Custom flags to allow in addition to those specified in the hybkit 
        specification.
        """,
        'str',
        None,
        {'nargs':'+'}
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
        Re-order flags to the hybkit-specificiation order when 
        writing hyb records.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'allow_undefined_flags': [
        False,
        """
        Allow use of flags not definied in the hybkit-specificiation order when 
        reading and writing hyb records. As the preferred alternative to 
        using this setting,
        the --custom_flags arguement can be be used to supply custom allowed flags.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'allow_unknown_seg_types': [
        False,
        """
        Allow unknown segment types when assigning segment types.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'check_complete_seg_types': [
        False,
        """
        Check every segment possibility when assigning segment types, rather than
        breaking after the first match is found. If True, finding segment types
        is slower but better at catching errors.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'allow_unknown_target_regions': [
        False,
        """
        Allow unknown mRNA regions when performing target region analysis.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'warn_unknown_target_regions': [
        False,
        """
        Print a warning message for unknown regions when performing
        target region analysis.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
}

# settings_info : HybFile
#: Information for settings of HybFile class.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
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
        {}
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
        {}
    ],
}

# settings_info : FoldRecord
#: Information for settings of FoldRecord class.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
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
}

# settings_info : FoldFile
#: Information for settings of FoldFile class.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
FoldFile_settings_info = {
    'fold_record_type': [
        'strict',
        """
        Type of fold record object to use. Options: "strict": FoldRecord, requires an exact sequence
        match to be paired with a HybRecord; "dynamic": DynamicFoldRecord, requires a sequence 
        match to the "dynamic" annotated regions of a HybRecord, and may be shorter/longer
        than the original sequence.
        """,
        'str',
        None,
        {'choices':['strict','dynamic']}
    ],
    'error_mode': [ 
        'raise',
        """
        Mode for handling errors during reading of HybFiles 
        (overridden by HybFoldIter.settings['error_mode'] when using HybFoldIter).
        Options: "raise": Raise an error when encountered and exit program;
        "warn_return": Print a warning and return the error_value ;
        "return": Return the error value with no program output.
        record is encountered.
        """,
        'str',
        None,
        {'choices':{'raise', 'warn_return', 'return'}}
    ],
}

# settings_info : HybFoldIter
#: Information for settings of HybFoldIter class.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
HybFoldIter_settings_info = {
    'error_checks': [
        ['hybrecord_indel', 'foldrecord_nofold', 'max_mismatch'],
        """
        Error checks for simultaneous HybFile and FoldFile parsing.
        Options: "hybrecord_indel": Error for HybRecord objects where one/both sequences have 
        insertions/deletions in alignment, which prevents matching of sequences; 
        "foldrecord_nofold": Error when failure in reading a fold_record object; 
        "max_mismatch": Error when mismatch between hybrecord and foldrecord sequences is 
        """,
        'str',
        None,
        {'choices':{'hybrecord_indel', 'foldrecord_nofold', 'max_mismatch'}}
    ],
    'error_mode': [
        'warn_skip',
        """
        Mode for handling errors found during error checks. Options: 
        "raise": Raise an error when encountered and exit program;
        "warn_return": Print a warning and return the value ;
        "warn_skip": Print a warning and continue to the next iteration;
        "skip": Continue to the next iteration without any output.
        """,
        'custom_bool_from_str',
        None,
        {'chioces':{'raise', 'warn_return', 'warn_skip', 'skip'}}
    ],
    'max_sequential_skips': [
        20,
        """
        Maximum number of record(-pairs) to skip in a row. Limited as several sequential skips
        usually indicates an issue with record formatting or a desynchrnoization between files.
        """,
        'int',
        None,
        {}
    ],
}
# settings_info : Analysis
#: Information for settings of Analysis classes.
#: setting_info format contains structure:: 
#:
#:     {
#:         setting_name : [ default_value, description, type_str, short_flag, argparse_fields ]
#:     }                   
Analysis_settings_info = {
    'count_mode': [
        'record',
        """
        Method for counting records. "read": use the number of reads per hyb record as the count 
        (may contain PCR duplicates); "record" count the number of records represented by each 
        (hyb record entry 1 for "unmerged" records, >= 1 for "merged" records)
        """,
        'str',
        None,
        {'choices':['read','record']}
    ],
    'mirna_sort': [
        True,
        """
        During TypeAnalysis, sort miRNAs first for "miRNA"-"Other" segtype pairs. 
        If False, sort alphabetically.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'allow_mirna_dimers': [
        False,
        """
        Include miRNA / miRNA dimers in TargetAnalysis.
        If False, exclude these from analysis results.
        """,
        'custom_bool_from_str',
        None,
        {}
    ],
    'type_sep': [
        '-',
        """
        Separator-string to place between types in analysis output.
        """,
        'str',
        None,
        {}
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

}

# Settings (active) : HybRecord
#: HybRecord Active Settings
HybRecord_settings = _settings_info_to_settings(HybRecord_settings_info)

# Settings (active) : HybFile
#: HybFile Active Settings
HybFile_settings = _settings_info_to_settings(HybFile_settings_info)

# Settings (active) : FoldRecord
#: FoldRecord Active Settings
FoldRecord_settings = _settings_info_to_settings(FoldRecord_settings_info)

# Settings (active) : FoldFile
#: FoldFile Active Settings
FoldFile_settings = _settings_info_to_settings(FoldFile_settings_info)

# Settings (active) : Analysis
#: Analysis Active Settings
Analysis_settings = _settings_info_to_settings(Analysis_settings_info)

# Settings (active) : HybFoldIter
#: HybFoldIter Active Settings
HybFoldIter_settings = _settings_info_to_settings(HybFoldIter_settings_info)