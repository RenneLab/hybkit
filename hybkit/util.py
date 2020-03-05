#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""This module contains helper functions for hybkit's command line scripts."""

import argparse
import os

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import hybkit


# Util : Global Variables
HYB_SUFFIXES = ['.hyb', '.Hyb', '.HYB']
VIENNA_SUFFIXES = ['.vienna', '.Vienna', '.VIENNA']
VIENNAD_SUFFIXES = ['.viennad', '.Viennad', '.VIENNAD']
CT_SUFFIXES = ['.ct', '.Ct', '.CT']
FOLD_SUFFIXES = VIENNA_SUFFIXES + VIENNAD_SUFFIXES + CT_SUFFIXES

USE_ABSPATH = False


# Util : Path Helper Functions
def dir_exists(dir_name):
    """
    Check if a directory exists at the provided path, and return a normalized path.

    Args:
        dir_name (str): Name of directory to check for existence.

    Returns:
        A normalized version of the path passed to dir_name.
    """
    if '~' in dir_name:
        dir_name = os.path.expanduser(dir_name)
    if not os.path.isdir(dir_name):
        message = 'Provided Directory: %s Does not exist.\n' % dir_name
        raise argparse.ArgumentTypeError(message)
    dir_name = os.path.abspath(dir_name)
    if USE_ABSPATH:
        dir_name = os.path.abspath(dir_name)
    return dir_name

# Util : Path Helper Functions 
def file_exists(file_name, required_suffixes=[]):
    """
    Check if a file exists at the provided path, and return a normalized path.

    Args:
        file_name (str): Name of file to check for existence.
        required_suffixes (list, optional): List of strings containing file-name suffixes.
            If provided, a file passed to file-exists must end with one of the strings 
            provided. Otherwise an error will be raised.

    Returns:
        A normalized version of the path passed to file_name.
    """
    # Check that directory exists.
    dir_exists(os.path.dirname(file_name))

    # Check that file exists.
    if not os.path.isfile(file_name):
        message = 'Provided Input File: %s does not exist.' % file_name
        raise argparse.ArgumentTypeError(message)
     
    # If required_suffixes provided, ensure the file has a required suffix.
    if required_suffixes:
        if not any(file_name.endswith(suffix) for suffix in required_suffixes):
            message = ('Provided Input File: %s' % file_name
                       + ' does not have an allowed suffix.'
                       +' {%s} ' % ', '.join(required_suffixes)
                      )
            print(message)
            raise argparse.ArgumentTypeError(message)

    # Normalize the file path
    file_name = os.path.normpath(file_name)

    # If global option set, then find the absolute path.
    if USE_ABSPATH:
        file_name = os.path.abspath(file_name)

    return os.path.abspath(file_name)

# Util : Path Helper Functions 
def hyb_exists(file_name):
    """
    Check if a .hyb file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes: %s.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = HYB_SUFFIXES
    return file_exists(file_name, use_suffixes)

hyb_exists.__doc__ = hyb_exists.__doc__ % HYB_SUFFIXES

# Util : Path Helper Functions 
def vienna_exists(file_name):
    """
    Check if a .vienna file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes: %s.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = VIENNA_SUFFIXES
    return file_exists(file_name, use_suffixes)

hyb_exists.__doc__ = hyb_exists.__doc__ % VIENNA_SUFFIXES

# Util : Path Helper Functions 
def viennad_exists(file_name):
    """
    Check if a .viennad file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes: %s.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = VIENNAD_SUFFIXES
    return file_exists(file_name, use_suffixes)

hyb_exists.__doc__ = hyb_exists.__doc__ % VIENNAD_SUFFIXES

# Util : Path Helper Functions 
def ct_exists(file_name):
    """
    Check if a .ct file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes: %s.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = CT_SUFFIXES
    return file_exists(file_name, use_suffixes)

hyb_exists.__doc__ = hyb_exists.__doc__ % CT_SUFFIXES

# Util : Path Helper Functions 
def fold_exists(file_name):
    """
    Check if a fold-representing file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes: %s.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = FOLD_SUFFIXES
    return file_exists(file_name, use_suffixes)

hyb_exists.__doc__ = hyb_exists.__doc__ % FOLD_SUFFIXES

# Util : Path Helper Functions 
def out_path_exists(file_name):
    """
    Check if the directory of the specified output path exists, and return a normalized path.

    Args:
        file_name (str): Name of path to an output file to check.

    Returns:
        A normalized version of the path passed to file_name.
    """
    # Check that directory exists.
    dir_exists(os.path.dirname(file_name))

    # Normalize the file path
    file_name = os.path.normpath(file_name)

    # If global option set, then find the absolute path.
    if USE_ABSPATH:
        file_name = os.path.abspath(file_name)

    return os.path.abspath(file_name)


# Argument Parser : Input/Output Options
in_hyb_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to hyb-format file with a ".hyb" suffix for use in the analysis.
                 """
in_hyb_parser.add_argument('-i', '--in_hyb', type=hyb_exists,
                           metavar='PATH_TO/MY_FILE.HYB',
                           required=True,
                           #nargs='1', 
                           help=_this_arg_help)

# Argument Parser : Input/Output Options
in_hybs_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to one or more hyb-format files with a ".hyb" suffix for use 
                 in the analysis.
                 """
in_hybs_parser.add_argument('-i', '--in_hyb', type=hyb_exists,
                            metavar='PATH_TO/MY_FILE.HYB',
                            required=True,
                            nargs='+', 
                            help=_this_arg_help)

# Argument Parser : Input/Output Options
out_hyb_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Optional path to hyb-format file for output (should include a ".hyb" suffix).
                 If not provided, the output for input file "PATH_TO/MY_FILE.HYB"
                 will be used as a template for the output "PATH_TO/MY_FILE_OUT.hyb".
                 """
out_hyb_parser.add_argument('-o', '--out_hyb', type=out_path_exists,
                           metavar='PATH_TO/OUT_FILE.HYB',
                           #required=True,
                           #nargs='1', 
                           help=_this_arg_help)

# Argument Parser : Input/Output Options
req_out_hyb_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to hyb-format file for output (should include a ".hyb" suffix).
                 """
req_out_hyb_parser.add_argument('-o', '--out_hyb', type=out_path_exists,
                                metavar='PATH_TO/OUT_FILE.HYB',
                                required=True,
                                #nargs='1', 
                                help=_this_arg_help)

# Argument Parser : Input/Output Options
out_dir_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to directory to output analysis files.
                 """
out_dir_parser.add_argument('-d', '--out_dir', type=dir_exists,
                           metavar='PATH_TO/MY_FILE.HYB',
                           #required=True,
                           #nargs='1', 
                           help=_this_arg_help)

# Argument Parser : Input/Output Options
# Stub
gen_opts_parser = argparse.ArgumentParser(add_help=False)

# Argument Parser : HybRecord  
# Create parser for HybRecord options
hybrecord_parser = argparse.ArgumentParser(add_help=False)
hr_group = hybrecord_parser.add_argument_group('Hyb Record Settings')
hr_defaults = hybkit.HybRecord.DEFAULTS

# Argument Parser : HybRecord
_this_arg_key  = 'custom_flags'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Custom flags to allow in addition to those specified in the hybkit 
                 specification.
                 """
hr_group.add_argument(_this_arg_flag, type=str,
                      nargs='+',
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'mirna_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Segment types to identifiy as microRNAs (miRNAs) during miRNA analysis.
                 """
hr_group.add_argument(_this_arg_flag, type=str,
                      nargs='+',
                      default=hybkit.HybRecord.MIRNA_TYPES,
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'coding_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Segment types to identify as coding/messenger RNA (mRNA) during
                 target region analysis.
                 """
hr_group.add_argument(_this_arg_flag, type=str,
                      nargs='+',
                      default=hybkit.HybRecord.CODING_TYPES,
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'hyb_placeholder'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Placeholder character/string for missing data in hyb files.
                 """
hr_group.add_argument(_this_arg_flag, type=str,
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord 
_this_arg_key  = 'reorder_flags'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Re-order flags to the hybkit-specificiation order when 
                 writing hyb records.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'allow_undefined_flags'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Allow use of flags not definied in the hybkit-specificiation order when 
                 reading and writing hyb records. As the preferred alternative to 
                 using this setting,
                 the --custom_flags arguement can be be used to supply custom allowed flags.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'allow_unknown_seg_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Allow unknown segment types when assigning segment types.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'check_complete_seg_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Check every segment possibility when assigning segment types, rather than
                 breaking after the first match is found. If True, finding segment types
                 is slower but better at catching errors.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'allow_unknown_regions'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Allow unknown mRNA regions when performing target region analysis.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybRecord
_this_arg_key  = 'warn_unknown_regions'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Print a warning message for unknown regions when performing
                 target region analysis.
                 """
hr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)


# Argument Parser : HybFile
# Create parser for HybFile options
hybfile_parser = argparse.ArgumentParser(add_help=False)
hf_group = hybfile_parser.add_argument_group('Hyb File Settings')
hf_defaults = hybkit.HybFile.DEFAULTS

# Argument Parser : HybFile 
_this_arg_key  = 'hybformat_id'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 The Hyb Software Package places further information in the "id" field of the
                 hybrid record that can be used to infer the number of contained read counts.
                 When set to True, the identifiers will be parsed as: "<read_id>_<read_count>" 
                 """
hf_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hf_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : HybFile 
_this_arg_key  = 'hybformat_ref'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 The Hyb Software Package uses a reference database with identifiers
                 that contain sequence type and other sequence information.
                 When set to True, all hyb file identifiers will be parsed as: 
                 "<gene_id>_<transcript_id>_<gene_name>_<seg_type>"
                 """
hf_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=hf_defaults[_this_arg_key],
                      help=_this_arg_help)


# Argument Parser : FoldRecord
# Create parser for FoldRecord options
foldrecord_parser = argparse.ArgumentParser(add_help=False)
fr_group = foldrecord_parser.add_argument_group('Fold Record Settings')
fr_defaults = hybkit.FoldRecord.DEFAULTS

# Argument Parser : FoldRecord 
_this_arg_key  = 'skip_bad_fold_records'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 When reading a fold record, bad fold records should be skipped.
                 When set to False, an error will be raised for improperly formatted fold records.
                 When True, improperly formatted fold records will be skipped.
                 """
fr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=fr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : FoldRecord 
_this_arg_key  = 'warn_bad_fold_records'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 When reading a fold record, print a warning when an improperly-formatted fold
                 record is encountered.
                 """
fr_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=fr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : FoldRecord 
_this_arg_key  = 'fold_placeholder'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Placeholder character/string for missing data for reading/writing fold records.
                 """
fr_group.add_argument(_this_arg_flag, type=str,
                      default=fr_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : FoldFile
# Create parser for FoldFile options
foldfile_parser = argparse.ArgumentParser(add_help=False)
ff_group = foldfile_parser.add_argument_group('Fold File Settings')
ff_defaults = hybkit.FoldFile.DEFAULTS

# Argument Parser : FoldFile 
_this_arg_key  = 'hybformat_file'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 The Hyb Software Package contains further information in the "name" field of the
                 viennad record that can be used to infer further information 
                 about the fold divisions.
                 """
ff_group.add_argument(_this_arg_flag, type=bool,
                      choices=[True, False],
                      default=ff_defaults[_this_arg_key],
                      help=_this_arg_help)


# Allow execution of module for testing purposes.
if __name__ == '__main__':
    all_parsers = [#in_hyb_parser, 
                   in_hybs_parser,
                   out_hyb_parser,
                   out_dir_parser,
                   hybrecord_parser, hybfile_parser, foldrecord_parser, foldfile_parser]
    test_parser = argparse.ArgumentParser(parents=all_parsers,
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    test_parser.print_help()

