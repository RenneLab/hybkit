#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""This module contains helper functions for hybkit's command line scripts."""

import os
import sys
import argparse
import textwrap

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

FILTER_OUT_SUFFIX = '_filtered'
ANALYSIS_OUT_SUFFIX = '_analyzed'

# Util : Argparse Helper Functions
def bool_from_string(value):
    if isinstance(value, bool):
       return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    # Snippet Credit goes to @Maxim at https://stackoverflow.com/questions/15008758/
    #   /parsing-boolean-values-with-argparse

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
    if '$' in dir_name:
        dir_name = os.path.expandvars(dir_name)
    #if not dir_name:
    #    dir_name = os.getcwd()
    dir_name = dir_name.strip()
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

vienna_exists.__doc__ = vienna_exists.__doc__ % VIENNA_SUFFIXES

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

viennad_exists.__doc__ = viennad_exists.__doc__ % VIENNAD_SUFFIXES

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

ct_exists.__doc__ = ct_exists.__doc__ % CT_SUFFIXES

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

fold_exists.__doc__ = fold_exists.__doc__ % FOLD_SUFFIXES

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

# Util : Path Helper Functions 
def make_out_file_name(in_file_name, name_suffix='out', in_suffix='', out_suffix='', 
                       out_dir='', seg_sep='_'):
    """
    Given an input file name, generate an output file name.

    Args:
        in_file_name (str): Name of input file as template.
        file_suffix (str): File type suffix.
        add_suffix (str): Name suffix to add to indicate file is output.
        out_dir (str): Directory path in which to place output file. 
        seg_sep (str): Separator string between file name segements.

    Returns:
        An output file path based on the input file template.
    """
    # Normalize the file path
    in_file_basename = os.path.basename(in_file_name)
    if in_suffix and in_file_basename.lower().endswith(in_suffix.lower()):
        suffix_len = len(in_suffix)
        in_file_basename = in_file_basename[:(-1 * suffix_len)]

    full_name_suffix = name_suffix
    if out_suffix and not full_name_suffix.lower().endswith(out_suffix.lower()):
        full_name_suffix += out_suffix
    if seg_sep and not full_name_suffix.startswith(seg_sep):
        full_name_suffix = seg_sep + full_name_suffix
    
    out_file_basename = in_file_basename + full_name_suffix

    out_file_full = os.path.join(out_dir, out_file_basename)

    return out_file_full

    # If global option set, then find the absolute path.
    if USE_ABSPATH:
        file_name = os.path.abspath(file_name)

    return os.path.abspath(file_name)


# Util : Path Helper Functions 
def validate_args(args, parser=None):
    """
    Check supplied arguments to make sure there are no hidden contradictions.

    Current checks:
        | If explicit output file names supplied, be sure that they match the number of 
          input files provided.

    Args:
        args (namespace): The arguments produced by argparse.
    """
   
    message = '\nArgument validation error: '
    if parser is not None:
        suffix = '\n\n' + parser.format_usage() + '\n' 
    else:
        suffix = 'Please use the -h or --help options for input requirements.'

    if hasattr(args, 'in_hyb') and hasattr(args, 'out_hyb') and args.out_hyb is not None:
        len_in_hyb = len(args.in_hyb)
        len_out_hyb = len(args.out_hyb)
        if len_in_hyb != len_out_hyb:
            message += 'The number of input files and number of output files provided '
            message += 'do not match. ( %i and %i )' % (len_in_hyb, len_out_hyb)
            message += '\n\nInput Files:\n    '
            message += '\n    '.join([f for f in args.in_hyb])
            message += '\n\nOutput Files:\n    '
            message += '\n    '.join([f for f in args.out_hyb])
            print(message + suffix)
            sys.exit(1)


# Argument Parser : Input/Output Options
in_hyb_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to hyb-format file with a ".hyb" suffix for use in the analysis.
                 """
in_hyb_parser.add_argument('-i', '--in_hyb', type=hyb_exists,
                           metavar='PATH_TO/MY_FILE.HYB',
                           required=True,
                           # nargs='1', 
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
                 Optional path to a hyb-format file for 
                 output (should include a ".hyb" suffix).
                 If not provided, the output for input file "PATH_TO/MY_FILE.HYB"
                 will be used as a template for the output "OUT_DIR/MY_FILE_OUT.HYB".
                 """
out_hyb_parser.add_argument('-o', '--out_hyb', type=out_path_exists,
                            metavar='PATH_TO/OUT_FILE.HYB',
                            # required=True,
                            # nargs='+', 
                            help=_this_arg_help)

# Argument Parser : Input/Output Options
out_hybs_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Optional path to one or more hyb-format file for 
                 output (should include a ".hyb" suffix).
                 If not provided, the output for input file "PATH_TO/MY_FILE.HYB"
                 will be used as a template for the output "OUT_DIR/MY_FILE_OUT.HYB".
                 """
out_hybs_parser.add_argument('-o', '--out_hyb', type=out_path_exists,
                             metavar='PATH_TO/OUT_FILE.HYB',
                             # required=True,
                             nargs='+', 
                             help=_this_arg_help)

# Argument Parser : Input/Output Options
req_out_hyb_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to hyb-format file for output (should include a ".hyb" suffix).
                 """
req_out_hyb_parser.add_argument('-o', '--out_hyb', type=out_path_exists,
                                metavar='PATH_TO/OUT_FILE.HYB',
                                required=True,
                                # nargs='1', 
                                help=_this_arg_help)

# Argument Parser : Input/Output Options
out_dir_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Path to directory for output of analysis files. 
                 Defaults to the current working directory.
                 """
out_dir_parser.add_argument('-d', '--out_dir', type=dir_exists,
                           # required=True,
                           # nargs='1',
                           default='$PWD', 
                           help=_this_arg_help)


# Argument Parser : Input/Output Options
out_suffix_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Suffix to add to the name of output hyb files.
                 (Note: If the provided suffix does not end with ".hyb", then ".hyb" will be added)
                 """
out_suffix_parser.add_argument('-u', '--out_suffix',
                               # required=True,
                               # nargs='1',
                               help=_this_arg_help)


out_opts_parser = argparse.ArgumentParser(add_help=False, 
                                          parents=[
                                                   out_hybs_parser,
                                                   out_dir_parser,
                                                   out_suffix_parser,
                                                  ],)

# Argument Parser : General Options
gen_opts_parser = argparse.ArgumentParser(add_help=False)
verbosity_group = gen_opts_parser.add_mutually_exclusive_group()

# Argument Parser : General Options
_this_arg_help = """
                 Print verbose output during run.
                 """
verbosity_group.add_argument('-v', '--verbose', action='store_true',
                             # nargs='+',
                             help=_this_arg_help)

# Argument Parser : General Options
_this_arg_help = """
                 Print no output during run.
                 """
verbosity_group.add_argument('-s', '--silent', action='store_true',
                             # nargs='+',
                             help=_this_arg_help)


# Argument Parser : HybRecord  
# Create parser for HybRecord options
hybrecord_parser = argparse.ArgumentParser(add_help=False)
hr_group = hybrecord_parser.add_argument_group('General Hyb Record Settings')
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
hr_group.add_argument(_this_arg_flag, type=bool_from_string,
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
hr_group.add_argument(_this_arg_flag, type=bool_from_string,
                      choices=[True, False],
                      default=hr_defaults[_this_arg_key],
                      help=_this_arg_help)




# Argument Parser : HybFile
# Create parser for HybFile options
hybfile_parser = argparse.ArgumentParser(add_help=False)
hf_group = hybfile_parser.add_argument_group('General Hyb File Settings')
hf_defaults = hybkit.HybFile.DEFAULTS

# Argument Parser : HybFile 
_this_arg_key  = 'hybformat_id'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 The Hyb Software Package places further information in the "id" field of the
                 hybrid record that can be used to infer the number of contained read counts.
                 When set to True, the identifiers will be parsed as: "<read_id>_<read_count>" 
                 """
hf_group.add_argument(_this_arg_flag, type=bool_from_string,
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
hf_group.add_argument(_this_arg_flag, type=bool_from_string,
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
fr_group.add_argument(_this_arg_flag, type=bool_from_string,
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
fr_group.add_argument(_this_arg_flag, type=bool_from_string,
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
ff_group.add_argument(_this_arg_flag, type=bool_from_string,
                      choices=[True, False],
                      default=ff_defaults[_this_arg_key],
                      help=_this_arg_help)

# Argument Parser : hyb_analysis
hyb_analysis_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Types of analyses to perform on input hyb file.
                 (Note: analyses can be combined, such as "--analysis_types segtype mirna")
                 """
hyb_analysis_parser.add_argument('-t', '--analysis_types',
                                 # required=True,
                                 nargs='+',
                                 default=['segtype'],
                                 choices=['segtype', 'mirna', 'target_region'],
                                 help=_this_arg_help)


# Argument Parser : hyb_analysis : segtype
segtype_opts_group = hyb_analysis_parser.add_argument_group('segtype Analysis Options')
# Argument Parser : hyb_analysis : segtype
_this_arg_help = """
                 Segment-type finding method to use for segtype analysis.
                 For a description of the different methods, see the HybRecord documentation
                 for the find_seg_types method.
                 """
segtype_opts_group.add_argument('--segtype_method',
                                # required=True,
                                # nargs='?',
                                default='hyb',
                                choices=hybkit.HybRecord.find_type_methods.keys(),
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : segtype
_this_arg_help = """
                 Segment-type finding paramaters to use for segtype analysis with some segtype
                 finding methods: {string_match, id_map}.
                 For a description of the different methods, see the HybRecord documentation
                 for the find_seg_types method.
                 """
segtype_opts_group.add_argument('--segtype_params', type=file_exists,
                                metavar='PATH_TO/PARAMATERS_FILE',
                                # required=True,
                                # nargs='?',
                                # default='hyb',
                                # choices=hybkit.HybRecord.find_type_methods,
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : segtype
_this_arg_key  = 'allow_unknown_seg_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Allow unknown segment types when assigning segment types.
                 """
segtype_opts_group.add_argument(_this_arg_flag, type=bool_from_string,
                                choices=[True, False],
                                default=hr_defaults[_this_arg_key],
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : segtype
_this_arg_key  = 'check_complete_seg_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Check every segment possibility when assigning segment types, rather than
                 breaking after the first match is found. If True, finding segment types
                 is slower but better at catching errors.
                 """
segtype_opts_group.add_argument(_this_arg_flag, type=bool_from_string,
                                choices=[True, False],
                                default=hr_defaults[_this_arg_key],
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : mirna
mirna_opts_group = hyb_analysis_parser.add_argument_group('mirna Analysis Options')

# Argument Parser : hyb_analysis : mirna
_this_arg_key  = 'mirna_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Segment types to identifiy as microRNAs (miRNAs) during miRNA analysis.
                 """
mirna_opts_group.add_argument(_this_arg_flag, type=str,
                              nargs='+',
                              default=hybkit.HybRecord.MIRNA_TYPES,
                              help=_this_arg_help)

# Arguemnt Parser : hyb_analysis : target_region
regions_opts_group = hyb_analysis_parser.add_argument_group('target_region Analysis Options')

# Argument Parser : HybRecord
_this_arg_key  = 'coding_types'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Segment types to identify as coding/messenger RNA (mRNA) during
                 target region analysis.
                 """
regions_opts_group.add_argument(_this_arg_flag, type=str,
                                nargs='+',
                                default=hybkit.HybRecord.CODING_TYPES,
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : target_region
_this_arg_key  = 'allow_unknown_regions'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Allow unknown mRNA regions when performing target region analysis.
                 """
regions_opts_group.add_argument(_this_arg_flag, type=bool_from_string,
                                choices=[True, False],
                                default=hr_defaults[_this_arg_key],
                                help=_this_arg_help)

# Argument Parser : hyb_analysis : target_region
_this_arg_key  = 'warn_unknown_regions'
_this_arg_flag = '--' + _this_arg_key
_this_arg_help = """
                 Print a warning message for unknown regions when performing
                 target region analysis.
                 """
regions_opts_group.add_argument(_this_arg_flag, type=bool_from_string,
                                choices=[True, False],
                                default=hr_defaults[_this_arg_key],
                                help=_this_arg_help)


# Argument Parser : hyb_filter
hyb_filter_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = """
                 Modes for evaluating multiple filters. 
                 The "all" mode requires all provided filters to be true for inclusion.
                 The "any" mode requires only one provided filter to be true for inclusion.
                 (Note: matching any exclusion filter is grounds for exclusion of record.)
                 """
hyb_filter_parser.add_argument('-m', '--filter_mode',
                               # required=True,
                               # nargs='+',
                               default='all',
                               choices={'all', 'any'},
                               help=_this_arg_help)

for i in range(1,4):
    _this_arg_help = """
                     Filter criteria #%i. 
                     Records matching the criteria will be included in output.
                     Includes a filter type, Ex: "seg_name_contains",
                     and an argument, Ex: "ENST00000340384".
                     (Note: not all filter types require a second argument, 
                     for Example: "has_mirna_seg")
                     """ % i
    if i == 1:
        flag_suffix = ''
    else:
        flag_suffix = '_' + str(i)

    hyb_filter_parser.add_argument('--filter' + flag_suffix,
                                   # required=True,
                                   nargs='+',
                                   # default='all',
                                   # choices={'all', 'any'},
                                   help=_this_arg_help)

for i in range(1,4):
    _this_arg_help = """
                     Exclusion filter criteria #%i. 
                     Records matching the criteria will be excluded from output.
                     Includes a filter type, Ex: "seg_name_contains",
                     and an argument, Ex: "ENST00000340384".
                     (Note: not all filter types require a second argument, 
                     for Example: "has_mirna_seg")
                     """ % i
    if i == 1:
        flag_suffix = ''
    else:
        flag_suffix = '_' + str(i)

    hyb_filter_parser.add_argument('--exclude' + flag_suffix,
                                   # required=True,
                                   nargs='+',
                                   # default='all',
                                   # choices={'all', 'any'},
                                   help=_this_arg_help)

# Argument Parser : Standardized Documentation Settings
output_description = textwrap.dedent("""
Output File Naming:
    Output files can be named in two fashions: via automatic name generation,
    or by providing specific out file names.

    Automatic Name Generation:
        For output name generation, the default respective naming scheme is used::
        
            hyb_script -i PATH_TO/MY_FILE_1.HYB [...]
                -->  OUT_DIR/MY_FILE_1_ADDSUFFIX.HYB
    
        This output file path can be modified with the arguments {--out_dir, --out_suffix}
        described below.
    
        The output directory defaults to the current working directory ``($PWD)``, and 
        can be modified with the ``--out_dir <dir>`` argument. 
        Note: The provided directory must exist, or an error will be raised.
        For Example::
    
            hyb_script -i PATH_TO/MY_FILE_1.HYB [...] --out_dir MY_OUT_DIR
                -->  MY_OUT_DIR/MY_FILE_1_ADDSUFFIX.HYB
    
        The suffix used for output files is based on the primary actions of the script.
        It can be specified using ``--out_suffix <suffix>``. This can optionally include
        the ".hyb" final suffix.
        for Example::
    
            hyb_script -i PATH_TO/MY_FILE_1.HYB [...] --out_suffix MY_SUFFIX
                -->  OUT_DIR/MY_FILE_1_MY_SUFFIX.HYB 
            #OR
            hyb_script -i PATH_TO/MY_FILE_1.HYB [...] --out_suffix MY_SUFFIX.HYB
                -->  OUT_DIR/MY_FILE_1_MY_SUFFIX.HYB
    
    Specific Output Names:
        Alternatively, specific file names can be provided via the -o/--out_hyb argument,
        ensuring that the same number of input and output files are provided. This argument
        takes precedence over all automatic output file naming options 
        (--out_dir, --out_suffix), which are ignored if -o/--out_hyb is provided.
        For Example::
    
            hyb_script [...] --out_hyb MY_OUT_DIR/OUT_FILE_1.HYB MY_OUT_DIR/OUT_FILE_2.HYB
                -->  MY_OUT_DIR/OUT_FILE_1.hyb
                -->  MY_OUT_DIR/OUT_FILE_2.hyb
        
        Note: The directory provided with output file paths (MY_OUT_DIR above) must exist, 
        otherwise an error will be raised.
    """)

# Argument Parser : hyb_type_analysis
hyb_type_analysis_parser = argparse.ArgumentParser(add_help=False)
# _this_arg_help = """
#                  Modes for evaluating multiple filters. 
#                  """
# hyb_type_analysis_parser.add_argument('-m', '--filter_mode',
#                                       # required=True,
#                                       # nargs='+',
#                                       default='all',
#                                       choices={'all', 'any'},
#                                       help=_this_arg_help)



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

