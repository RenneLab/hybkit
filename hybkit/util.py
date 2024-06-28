#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""This module contains helper functions for hybkit's command line scripts."""

import argparse
import copy
import os
import sys
import textwrap
from typing import Any, List, Optional, Union

from hybkit import settings, type_finder
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

# ----- Linting Directives:
# ruff: noqa: F401 E402 SLF001

# ----- Begin Argparse Helper Functions -----
# Util : Argparse Helper Functions
def get_argparse_doc(docstring: str) -> str:
    """
    Get the argparse description from a docstring.

    Args:
        docstring (str): A docstring.

    Returns:
        A string containing the argparse description.
    """
    initial_replace_pairs = [('|', ''), ('==', ''), ('``', '"'), ('`', ''),
                             (':ref:', '')]
    final_replace_pairs = [('::', ':')]
    argparse_doc = docstring
    for query_str, replace_str in initial_replace_pairs:
        argparse_doc = argparse_doc.replace(query_str, replace_str)
    retain_lines = []
    for line in argparse_doc.splitlines():
        if line.strip() == '::':
            continue
        retain_lines.append(line.rstrip())
    argparse_doc = '\n'.join(retain_lines) + '\n'
    for query_str, replace_str in final_replace_pairs:
        argparse_doc = argparse_doc.replace(query_str, replace_str)
    return argparse_doc


# Util : Argparse Helper Functions
def _bool_from_string(value: Optional[Union[str, bool]]) -> bool:
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected. Provided: %s' % value)

_custom_types = {
    'custom_bool_from_str': _bool_from_string,
    'str': str,
    'int': int,
}
_custom_type_choices = {
    'custom_bool_from_str': [True, False],
    'str': None,
    'int': None,
}


# Util : Path Helper Functions
def dir_exists(dir_name: str) -> str:
    """
    Check if a directory exists at the provided path (else raise), and return a normalized path.

    Args:
        dir_name (str): Name of directory to check for existence.

    Returns:
        A normalized version of the path passed to dir_name.
    """
    if '~' in dir_name:
        dir_name = os.path.expanduser(dir_name)
    if '$' in dir_name:
        dir_name = os.path.expandvars(dir_name)
    # if not dir_name:
    #    dir_name = os.getcwd()
    dir_name = dir_name.strip()
    if not os.path.isdir(dir_name):
        message = 'Provided Directory: %s Does not exist.\n' % dir_name
        raise argparse.ArgumentTypeError(message)
    dir_name = os.path.abspath(dir_name)
    if settings._USE_ABSPATH:
        dir_name = os.path.abspath(dir_name)
    return dir_name


# Util : Path Helper Functions
def file_exists(file_name: str, required_suffixes: Optional[List[str]]=None) -> str:
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
    # If there is a directory, check that directory exists.
    if str(os.path.dirname(file_name)).strip():
        dir_exists(os.path.dirname(file_name))

    # Check that file exists.
    if not os.path.isfile(file_name):
        message = 'Provided Input File: %s does not exist.' % file_name
        raise argparse.ArgumentTypeError(message)

    # If required_suffixes provided, ensure the file has a required suffix.
    if (required_suffixes is not None
        and not any(file_name.endswith(suffix) for suffix in required_suffixes)):
            message = ('Provided Input File: %s' % file_name
                       + ' does not have an allowed suffix.'
                       + ' {%s} ' % ', '.join(required_suffixes)
                       )
            raise argparse.ArgumentTypeError(message)

    # Normalize the file path
    file_name = os.path.normpath(file_name)

    # If global option set, then find the absolute path.
    if settings._USE_ABSPATH:
        file_name = os.path.abspath(file_name)

    return os.path.abspath(file_name)


# Util : Path Helper Functions
def hyb_exists(file_name: str) -> str:
    """
    Check if a .hyb file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes in
    :attr:`hybkit.settings.HYB_SUFFIXES`.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = settings.HYB_SUFFIXES
    return file_exists(file_name, use_suffixes)


# Util : Path Helper Functions
def vienna_exists(file_name: str) -> str:
    """
    Check if a .vienna file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes in
    :attr:`hybkit.settings.VIENNA_SUFFIXES`.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = settings.VIENNA_SUFFIXES
    return file_exists(file_name, use_suffixes)


# Util : Path Helper Functions
def ct_exists(file_name: str) -> str:
    """
    Check if a .ct file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes in
    :attr:`hybkit.settings.CT_SUFFIXES`.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = settings.CT_SUFFIXES
    return file_exists(file_name, use_suffixes)


# Util : Path Helper Functions
def fold_exists(file_name: str) -> str:
    """
    Check if a fold-representing file exists at the provided path, and return a normalized path.

    Wrapper for :func:`file_exists` that includes the required suffixes in
    :attr:`hybkit.settings.FOLD_SUFFIXES`.

    Args:
        file_name (str): Name of file to check for existence.

    Returns:
        A normalized version of the path passed to file_name.
    """
    use_suffixes = settings.FOLD_SUFFIXES
    return file_exists(file_name, use_suffixes)


# Util : Path Helper Functions
def out_path_exists(file_name: str) -> str:
    """
    Check if the directory of the specified output path exists, and return a normalized path.

    Args:
        file_name (str): Name of path to an output file to check.

    Returns:
        A normalized version of the path passed to file_name.
    """
    # Check that directory exists.
    if str(os.path.dirname(file_name)).strip():
        dir_exists(os.path.dirname(file_name))

    # Normalize the file path
    file_name = os.path.normpath(file_name)

    # If global option set, then find the absolute path.
    if settings._USE_ABSPATH:
        file_name = os.path.abspath(file_name)

    return os.path.abspath(file_name)


# Util : Path Helper Functions
def make_out_file_name(
        in_file_name: str,
        name_suffix: str = 'out',
        in_suffix: str = '',
        out_suffix: str = '',
        out_dir: str = '',
        seg_sep: str = '_'
        ) -> str:
    """
    Given an input file name, generate an output file name.

    Args:
        in_file_name (str): Name of input file as template.
        name_suffix (str): Suffix to add to name before file type.
        in_suffix (str): File type suffix on in_file_name (to remove).
        out_suffix (str): File type suffix to add to final output file.
        out_dir (str): Directory path in which to place output file.
        seg_sep (str): Separator string between file name segments.

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
    if seg_sep and full_name_suffix and not full_name_suffix.startswith(seg_sep):
        full_name_suffix = seg_sep + full_name_suffix

    out_file_basename = in_file_basename + full_name_suffix

    out_file_full = os.path.join(out_dir, out_file_basename)

    # If global option set, then find the absolute path.
    if settings._USE_ABSPATH:
        out_file_full = os.path.abspath(out_file_full)

    return out_file_full


# Util : Path Helper Functions
def validate_args(
        args: argparse.Namespace,
        parser: Optional[argparse.ArgumentParser] = None
        ) -> bool:
    """
    Check supplied arguments to make sure there are no hidden contradictions.

    Current checks:
        * If explicit output file names supplied, be sure that they match the number of
          input files provided.
        * If fold files provided, make sure that they match the number of input hyb
          files provided.

    Args:
        args (argparse.Namespace): The arguments produced by argparse.
        parser (argparse.ArgumentParser, optional): Argparse parser object to use for
            verbose outputting of help message.
    """
    message = '\nArgument validation error: '
    suffix = '\n'
    if parser is not None:
        suffix = '\n\n' + parser.format_usage() + '\n'
    else:
        suffix = 'Please use the -h or --help options for input requirements.'

    ret_val = True
    if hasattr(args, 'in_hyb') and hasattr(args, 'out_hyb') and args.out_hyb is not None:
        len_in_hyb = len(args.in_hyb)
        len_out_hyb = len(args.out_hyb)
        if len_in_hyb != len_out_hyb:
            message += 'The number of input hyb files and output hyb files provided '
            message += 'do not match. ( %i and %i )' % (len_in_hyb, len_out_hyb)
            message += '\n\nInput Files:\n    '
            message += '\n    '.join(list(args.in_hyb))
            message += '\n\nOutput Files:\n    '
            message += '\n    '.join(list(args.out_hyb))
            print(message + suffix)
            ret_val = False

    if (hasattr(args, 'in_fold') and args.in_fold is not None
            and hasattr(args, 'out_fold') and args.out_fold is not None):
        len_in_fold = len(args.in_fold)
        len_out_fold = len(args.out_fold)
        if len_in_fold != len_out_fold:
            message += 'The number of input fold files and output fold files provided '
            message += 'do not match. ( %i and %i )' % (len_in_fold, len_out_fold)
            message += '\n\nInput Files:\n    '
            message += '\n    '.join(list(args.in_fold))
            message += '\n\nOutput Files:\n    '
            message += '\n    '.join(list(args.out_fold))
            print(message + suffix)
            ret_val = False

    if hasattr(args, 'in_hyb') and hasattr(args, 'in_fold') and args.in_fold is not None:
        len_in_hyb = len(args.in_hyb)
        len_in_fold = len(args.in_fold)
        if len_in_hyb != len_in_fold:
            message += 'The number of input hyb files and input fold files provided '
            message += 'do not match. ( %i and %i )' % (len_in_hyb, len_in_fold)
            message += '\n\nInput Files:\n    '
            message += '\n    '.join(list(args.in_hyb))
            message += '\n\nOutput Files:\n    '
            message += '\n    '.join(list(args.in_fold))
            print(message + suffix)
            ret_val = False

    return ret_val


# Util : Path Helper Functions
def validate_args_exit(
        args: argparse.Namespace,
        parser: Optional[argparse.ArgumentParser] = None
        ) -> None:
    """
    Check supplied arguments using :func:`validate_args`, and exit if a conflict exists.

    Args:
        args (argparse.Namespace): The arguments produced by argparse.
        parser (argparse.ArgumentParser, optional): Argparse parser object to use for
            verbose outputting of help message.
    """
    if not validate_args(args, parser):
        sys.exit(1)


# ----- Begin Argparse Helper Classes -----
# Util : Argparse Helper Formatter Class
class _HybkitFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
    ):
    pass


# ----- Begin Argparse Parsers -----
# Start I/O Files
# Argument Parser : Input/Output: I/O Files
in_hybs_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    REQUIRED path to one or more hyb-format files with a ".hyb" suffix for use
    in the evaluation.
    """
)
in_hybs_parser.add_argument(
    '-i', '--in_hyb', type=hyb_exists,
    metavar='PATH_TO/MY_FILE.HYB',
    required=True,
    nargs='+',
    help=_this_arg_help
)

# Argument Parser : Input/Output : I/O Files
in_folds_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    REQUIRED path to one or more RNA secondary-structure files with a
    ".vienna" or ".ct" suffix for use in the evaluation.
    """
)
in_folds_parser.add_argument(
    '-f', '--in_fold', type=fold_exists,
    metavar='PATH_TO/MY_FILE.VIENNA',
    # required=True,
    nargs='*',
    help=_this_arg_help
)

# Argument Parser : Input/Output : I/O Files
out_hybs_parser = argparse.ArgumentParser(add_help=False)
# In out_opts_... combined parsers
_this_arg_help = (
    """
    Optional path to one or more hyb-format file for
    output (should include a ".hyb" suffix).
    If not provided, the output for input file "PATH_TO/MY_FILE.HYB"
    will be used as a template for the output "OUT_DIR/MY_FILE_OUT.HYB".
    """
)
out_hybs_parser.add_argument(
    '-o', '--out_hyb', type=out_path_exists,
    metavar='PATH_TO/OUT_FILE.HYB',
    # required=True,
    nargs='+',
    help=_this_arg_help
)

# Argument Parser : Input/Output : I/O Files
out_folds_parser = argparse.ArgumentParser(add_help=False)
# In out_opts_... combined parsers
_this_arg_help = (
    """
    Optional path to one or more ".vienna" or ".ct"-format files for
    output (should include appropriate ".vienna"/".ct" suffix).
    If not provided, the output for input file "PATH_TO/MY_FILE.VIENNA"
    will be used as a template for the output "OUT_DIR/MY_FILE_OUT.VIENNA".
    """
)
out_folds_parser.add_argument(
    '-l', '--out_fold', type=out_path_exists,
    metavar='PATH_TO/OUT_FILE.VIENNA',
    # required=True,
    nargs='+',
    help=_this_arg_help
)

# Argument Parser : Input/Output : I/O Files
out_basenames_parser = argparse.ArgumentParser(add_help=False)
# In out_opts_... combined parsers
_this_arg_help = (
    """
    Optional path to one or more basename prefixes to use for
    output. The appropriate suffix will be added
    based on the specific name.
    If not provided, the output for input file "PATH_TO/MY_FILE.HYB"
    will be used as a template for the basename "OUT_DIR/MY_FILE".
    """
)
out_basenames_parser.add_argument(
    '-o', '--out_basename', type=out_path_exists,
    metavar='PATH_TO/OUT_BASENAME',
    # required=True,
    nargs='+',
    help=_this_arg_help
)

# Start I/O File Options
# Argument Parser : Input/Output : I/O File Options
# In out_opts_... combined parsers
out_dir_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Path to directory for output of files.
    Defaults to the current working directory.
    """
)
out_dir_parser.add_argument(
    '-d', '--out_dir', type=dir_exists,
    # required=True,
    # nargs='1',
    default='$PWD',
    help=_this_arg_help
)


# Argument Parser : Input/Output : I/O File Options
# In out_opts_... combined parsers
out_suffix_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Suffix to add to the name of output files, before any
    file- or analysis-specific suffixes. The file-type appropriate suffix
    will be added automatically.
    """
)
out_suffix_parser.add_argument(
    '-u', '--out_suffix',
    # required=True,
    # nargs='1',
    help=_this_arg_help
)

# Start I/O Combined
# Argument Parser : Input/Output : I/O Combined
# cmb_hyb_fold_inputs_parser = argparse.ArgumentParser(
#     add_help=False,
#     parents=[
#         in_hybs_parser,
#         in_folds_parser,
#     ],
# )

# Argument Parser : Input/Output : I/O Combined
cmb_hyb_fold_io_parser = argparse.ArgumentParser(
    add_help=False,
    parents=[
        in_hybs_parser,
        in_folds_parser,
        out_hybs_parser,
        out_folds_parser,
    ],
)

# Argument Parser : Input/Output : I/O Combined
cmb_out_opts_parser = argparse.ArgumentParser(
    add_help=False,
    parents=[
        out_dir_parser,
        out_suffix_parser,
    ],
)

# Argument Parser : Input/Output : I/O Combined
# out_opts_parser_hyb = argparse.ArgumentParser(
#     add_help=False,
#     parents=[
#         out_hybs_parser,
#         out_dir_parser,
#         out_suffix_parser,
#     ],
# )

# Argument Parser : Input/Output : I/O Combined
# out_opts_parser_basename = argparse.ArgumentParser(
#     add_help=False,
#     parents=[
#         out_basenames_parser,
#         out_dir_parser,
#         out_suffix_parser,
#     ],
# )

# Argument Parser : Input/Output : I/O Combined
# out_fold_opts_parser = argparse.ArgumentParser(
#     add_help=False,
#     parents=[
#         out_folds_parser,
#         out_dir_parser,
#         out_suffix_parser,
#     ],
# )

# Argument Parser : Input/Output : I/O Combined
# cmb_out_analysis_parser = argparse.ArgumentParser(
#     add_help=False,
#     parents=[
#         out_basenames_parser,
#         out_dir_parser,
#         out_suffix_parser,
#     ],
# )

# Start General Options
# Argument Parser : General Options
gen_opts_parser = argparse.ArgumentParser(add_help=False)

# Argument Parser : General Options
_this_arg_help = (
    """
    Print version and exit.
    """
)
_arg_version_str = f'{__version__}  (hybkit API: {__version__})'
gen_opts_parser.add_argument(
    '--version', action='version', version='    %(prog)s ' + _arg_version_str,
    help=_this_arg_help,
)

verbosity_group = gen_opts_parser.add_mutually_exclusive_group()
# Argument Parser : General Options : Verbosity
_this_arg_help = (
    """
    Print verbose output during run.
    """
)
verbosity_group.add_argument(
    '-v', '--verbose', action='store_true',
    # nargs='+',
    help=_this_arg_help
)

# Argument Parser : General Options : Verbosity
_this_arg_help = (
    """
    Print no output during run.
    """
)
verbosity_group.add_argument(
    '-s', '--silent', action='store_true',
    # nargs='+',
    help=_this_arg_help
)

# Start Record Manip Options
# Argument Parser : Record Manipulation Options
record_manip_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Set "dataset" flag to value of the input file name.
    """
)
record_manip_parser.add_argument(
    '--set_dataset',
    action='store_true',
    # required=True,
    # nargs='1',
    help=_this_arg_help
)

# Start Class Settings Parsers
# Argument Parser : Class Settings Parser
_class_settings_groups = {}
# Create parser for HybRecord options
hybrecord_parser = argparse.ArgumentParser(add_help=False)
hr_group = hybrecord_parser.add_argument_group('Hyb Record Settings')
_class_settings_groups['HybRecord'] = hr_group

# Create parser for HybFile options
hybfile_parser = argparse.ArgumentParser(add_help=False)
hf_group = hybfile_parser.add_argument_group('Hyb File Settings')
_class_settings_groups['HybFile'] = hf_group

# Create parser for FoldRecord options
foldrecord_parser = argparse.ArgumentParser(add_help=False)
fr_group = foldrecord_parser.add_argument_group('Fold Record Settings')
_class_settings_groups['FoldRecord'] = fr_group

# Create parser for FoldFile options
foldfile_parser = argparse.ArgumentParser(add_help=False)
ff_group = foldfile_parser.add_argument_group('Fold File Settings')
_class_settings_groups['FoldFile'] = ff_group

# Create parser for HybFoldIter options
hybfolditer_parser = argparse.ArgumentParser(add_help=False)
hfi_group = hybfolditer_parser.add_argument_group('Hyb-Fold Iterator Settings')
_class_settings_groups['HybFoldIter'] = hfi_group

# Create parser for Analysis options
analysis_parser = argparse.ArgumentParser(add_help=False)
a_group = analysis_parser.add_argument_group('Analysis Settings')
_class_settings_groups['Analysis'] = a_group

for _cls_name, _cls_group in _class_settings_groups.items():
    _cls_settings = getattr(settings, _cls_name + '_settings_info')
    for _key, _details in _cls_settings.items():
        _flag = '--' + _key
        _default, _description, _type_str, _short_flag, _extra_kwargs = _details
        _use_args = []
        if _short_flag is not None:
            _use_args.append(_short_flag)
        _use_args.append(_flag)
        _use_kwargs = copy.deepcopy(_extra_kwargs)
        _use_kwargs['help'] = _description
        _use_kwargs['type'] = _custom_types[_type_str]
        if _custom_type_choices[_type_str] is not None:
            _use_kwargs['choices'] = _custom_type_choices[_type_str]
        _use_kwargs['type'] = _custom_types[_type_str]
        _use_kwargs['default'] = _default
        _cls_group.add_argument(*_use_args, **_use_kwargs)

# Argument Parser : Record Manipulation Options
cmb_hyb_fold_class_settings_parser = argparse.ArgumentParser(
    add_help=False,
    parents=[
        hybrecord_parser,
        hybfile_parser,
        foldrecord_parser,
        foldfile_parser,
        hybfolditer_parser,
    ],
)


#  ----- Begin Task-specific Parsers -----
# Start eval
# Argument Parser : hyb_eval
hyb_eval_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Types of evaluations to perform on input hyb file.
    (Note: evaluations can be combined, such as "--eval_types type mirna")
    """
)
hyb_eval_parser.add_argument(
    '-t', '--eval_types',
    # required=True,
    nargs='+',
    default=['type'],
    choices=['type', 'mirna'],
    help=_this_arg_help
)

# Start type_opts
# Argument Parser : type_opts
type_opts_group = hyb_eval_parser.add_argument_group('type Analysis Options')
# Argument Parser : type_opts : type
_this_arg_help = (
    """
    Segment-type finding method to use for type evaluation.
    For a description of the different methods, see the HybRecord documentation
    for the eval_types method.
    """
)
type_opts_group.add_argument(
    '--type_method',
    # required=True,
    # nargs='?',
    default=type_finder.TypeFinder.default_method,
    choices=type_finder.TypeFinder.methods.keys(),
    help=_this_arg_help
)

# Argument Parser : type_opts : type
_this_arg_help = (
    """
    Segment-type finding parameters file to use for type evaluation with some type
    finding methods: {string_match, id_map}.
    For a description of the different methods, see the HybRecord documentation
    for the find_seg_types method.
    """
)
type_opts_group.add_argument(
    '--type_params_file', type=file_exists,
    metavar='PATH_TO/PARAMETERS_FILE',
    # required=True,
    # nargs='?',
    # default='hyb',
    # choices=HybRecord.find_type_methods,
    help=_this_arg_help
)

# Start filter
# Argument Parser : hyb_filter
hyb_filter_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Modes for evaluating multiple filters.
    The "all" mode requires all provided filters to be true for inclusion.
    The "any" mode requires only one provided filter to be true for inclusion.
    (Note: matching any exclusion filter is grounds for exclusion of record.)
    """
)
# Argument Parser : hyb_filter : filter_mode
hyb_filter_parser.add_argument(
    '-m', '--filter_mode',
    # required=True,
    # nargs='+',
    default='all',
    choices={'all', 'any'},
    help=_this_arg_help
)

_this_arg_help = (
    """
    Skip sequential duplicate read IDs before filtering.
    """
)
# Argument Parser : hyb_filter : skip_dup_id_before
hyb_filter_parser.add_argument(
    '--skip_dup_id_before',
    # required=True,
    # nargs='+',
    action='store_true',
    help=_this_arg_help
)

_this_arg_help = (
    """
    Skip sequential duplicate read IDs after filtering.
    """
)
# Argument Parser : hyb_filter : skip_dup_id_after
hyb_filter_parser.add_argument(
    '--skip_dup_id_after',
    # required=True,
    # nargs='+',
    action='store_true',
    help=_this_arg_help
)

# Argument Parser : hyb_filter : filter_criteria
for i in range(1, 4):
    _this_arg_help = (
        """
        Filter criteria #%i.
        Records matching the criteria will be included in output.
        Includes a filter type, Ex: "seg_name_contains",
        and an argument, Ex: "ENST00000340384".
        (Note: not all filter types require a second argument,
        for Example: "has_mirna_seg")
        """ % i
    )
    if i == 1:
        flag_suffix = ''
    else:
        flag_suffix = '_' + str(i)

    hyb_filter_parser.add_argument(
        '--filter' + flag_suffix,
        # required=True,
        nargs='+',
        # default='all',
        # choices={'all', 'any'},
        help=_this_arg_help
    )

# Argument Parser : hyb_filter : filter_criteria
for i in range(1, 4):
    _this_arg_help = (
        """
        Exclusion filter criteria #%i.
        Records matching the criteria will be excluded from output.
        Includes a filter type, Ex: "seg_name_contains",
        and an argument, Ex: "ENST00000340384".
        (Note: not all filter types require a second argument,
        for Example: "has_mirna_seg")
        """ % i
    )
    if i == 1:
        flag_suffix = ''
    else:
        flag_suffix = '_' + str(i)

    hyb_filter_parser.add_argument(
        '--exclude' + flag_suffix,
        # required=True,
        nargs='+',
        # default='all',
        # choices={'all', 'any'},
        help=_this_arg_help
    )

# Argument Parser : hyb_fold_analyze
hyb_analyze_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Analysis to perform on input hyb and fold files.
    """
)
hyb_analyze_parser.add_argument(
    '-a', '--analysis_types',
    # required=True,
    nargs='+',
    default=['fold'],
    action='store',
    choices=copy.deepcopy(settings.ANALYSIS_TYPE_OPTIONS),
    help=_this_arg_help
)

# Start  all_analyze
# Argument Parser : all_analyze : analysis_name
all_analyze_parser = argparse.ArgumentParser(add_help=False)
_this_arg_help = (
    """
    Name / title of analysis data.
    """
)
all_analyze_parser.add_argument(
    '-n', '--analysis_name',
    # required=True,
    # nargs='1',
    # default=None,
    # choices=[True, False],
    help=_this_arg_help
)

# Argument Parser : all_analyze : make_plots
_this_arg_help = (
    """
    Create plots of analysis output.
    """
)
all_analyze_parser.add_argument(
    '-p', '--make_plots',
    # required=True,
    type=_bool_from_string,
    default=True,
    choices=[True, False],
    help=_this_arg_help
)


# Start Documentation Settings
# Argument Parser : Standardized Documentation Settings
# output_description = textwrap.dedent(
output_description = textwrap.dedent(
    """
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

    """
)


# ------ Begin Settings Manipulation Functions ------
def set_setting(
        setting: str,
        set_value: Any,  # noqa: ANN401
        verbose: bool = False
        ) -> str:
    """
    Take a namespace object as from an argparse parser and update settings.

    Each setting in the following settings dictionaries are checked and set where applicable:

        ===================================== ============================================
        :class:`~hybkit.HybRecord` Settings   :attr:`hybkit.settings.HybRecord_settings`
        :class:`~hybkit.HybFile` Settings     :attr:`hybkit.settings.HybFile_settings`
        :class:`~hybkit.FoldRecord` Settings  :attr:`hybkit.settings.FoldRecord_settings`
        :class:`~hybkit.FoldFile` Settings    :attr:`hybkit.settings.FoldFile_settings`
        :class:`~hybkit.HybFoldIter` Settings :attr:`hybkit.settings.HybFoldIter_settings`
        :class:`~hybkit.Analysis` Settings    :attr:`hybkit.settings.Analysis_settings`
        ===================================== ============================================

    Args:
        setting (str): Name of setting to change
        set_value (str): New value for setting
        verbose (:obj:`bool`, optional): If True, print when changing setting.
    """
    out_report = ''
    setting_found = False
    for class_name in ['HybRecord', 'HybFile', 'FoldRecord',
                       'FoldFile', 'HybFoldIter', 'Analysis']:
        cls_settings_info = getattr(settings, class_name + '_settings_info')
        cls_settings = getattr(settings, class_name + '_settings')
        do_check_list = isinstance(set_value, list)
        if setting in cls_settings:
            setting_found = True
            old_setting = cls_settings[setting]
            if 'choices' in cls_settings_info[setting][4]:
                choices = cls_settings_info[setting][4]['choices']
            else:
                choices = None
            if choices is not None:
                if do_check_list:
                    for check_value in set_value:
                        if check_value not in choices:
                            message = f'Invalid value for {setting} setting: {check_value}'
                            message += '\nChoices are: %s' % str(choices)
                            raise HybkitArgError(message)
                elif set_value not in choices:
                    message = f'Invalid value for {setting} setting: {set_value}'
                    message += '\nChoices are: %s' % str(choices)
                    raise HybkitArgError(message)
            if old_setting is not None and set_value != old_setting:
                out_report += 'Setting %s Setting: ' % class_name
                out_report += f'"{setting}" to "{set_value!s}"\n'
                cls_settings[setting] = set_value
    if not setting_found:
        message = 'Setting "%s" not found' % setting
        raise HybkitArgError(message)
    if verbose and out_report.strip():
        print(out_report)
    if out_report:
        out_report = '\n' + out_report
    return out_report


def set_settings_from_namespace(
        nspace: argparse.Namespace,
        verbose: bool = False,
        ) -> None:
    """
    Take a namespace object as from an argparse parser and update settings.

    See :func:`set_setting` for details

    Args:
        nspace (argparse.Namespace): Namespace containing settings
        verbose (:obj:`bool`, optional): If True, print when changing setting.
    """
    out_report = '\n'
    for class_name in ['HybRecord', 'HybFile', 'FoldRecord',
                       'FoldFile', 'HybFoldIter', 'Analysis']:
        _cls_settings_info = getattr(settings, class_name + '_settings_info')
        cls_settings = getattr(settings, class_name + '_settings')
        for setting_key in cls_settings:
            if hasattr(nspace, setting_key):
                argparse_setting = getattr(nspace, setting_key)
                out_report += set_setting(setting_key, argparse_setting, verbose=False)

    if verbose and out_report.strip():
        print(out_report)

# Allow execution of module for testing purposes.
# if __name__ == '__main__':
#    all_parsers = [#in_hyb_parser,
#    test_parser = argparse.ArgumentParser(parents=all_parsers,
#                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    test_parser.print_help()
