#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""

# import sys
# import copy
# from contextlib import nullcontext as does_not_raise
import argparse
import os

import pytest

import hybkit
from auto_tests.test_helper_data import (
    test_ct_file_name,
    test_hyb_file_name,
    test_out_dir,
    test_vienna_file_name,
)
from auto_tests.test_helper_functions import get_expected_result_context

# ----- Linting Directives:
# ruff: noqa: SLF001 ARG001 FBT003

# ----- Start Test Util -----
def test_util_misc():
    """Test miscellaneous functions in hybkit.util."""
    _original_abspath = hybkit.settings._USE_ABSPATH
    hybkit.settings._USE_ABSPATH = True

    # Test get_argparse_doc
    source_str = ':ref:`|==``::\n::\n'
    assert hybkit.util.get_argparse_doc(source_str) == '":\n'

    # Test _bool_from_string
    assert hybkit.util._bool_from_string(True)
    assert hybkit.util._bool_from_string('yes')
    assert not hybkit.util._bool_from_string(False)
    assert not hybkit.util._bool_from_string('no')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util._bool_from_string('invalid')

    # Test hybkit.util.dir_exists
    assert hybkit.util.dir_exists('~')
    assert hybkit.util.dir_exists('${PWD}')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.dir_exists('nonexistent_dir')

    # Test hybkit.util.file_exists
    assert hybkit.util.file_exists(__file__)
    with pytest.raises(argparse.ArgumentTypeError):
        assert not hybkit.util.file_exists('nonexistent_file')

    # Test hybkit.util.hyb_exists
    assert hybkit.util.hyb_exists(test_hyb_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.hyb_exists(test_vienna_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.hyb_exists('nonexistent_file')

    # Test hybkit.util.vienna_exists
    assert hybkit.util.vienna_exists(test_vienna_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.vienna_exists(test_hyb_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.vienna_exists('nonexistent_file')

    # Test hybkit.util.ct_exists
    assert hybkit.util.ct_exists(test_ct_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.ct_exists(test_hyb_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.ct_exists('nonexistent_file')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.ct_exists(test_vienna_file_name)

    # Test hybkit.util.fold_exists
    assert hybkit.util.fold_exists(test_vienna_file_name)
    assert hybkit.util.fold_exists(test_ct_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.fold_exists(test_hyb_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.fold_exists('nonexistent_file')

    # Test hybkit.util.out_path_exists
    assert hybkit.util.out_path_exists(test_hyb_file_name)
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.out_path_exists('/nonexistent/path')

    test_path_1 = hybkit.util.make_out_file_name(
        test_hyb_file_name,
        name_suffix='out',
        in_suffix='.hyb',
        out_suffix='.new',
        out_dir='',
        seg_sep='_'
    )
    # In name: 'test_hybrid.hyb' -> 'test_hybrid_out.new'
    assert test_path_1 == os.path.join(os.path.abspath('.'), 'test_hybrid_out.new')

    test_path_2 = hybkit.util.make_out_file_name(
        test_hyb_file_name,
        name_suffix='out',
        in_suffix='.hyb',
        out_suffix='.new',
        out_dir=test_out_dir,
        seg_sep='.'
    )
    # In name: 'test_hybrid.hyb' -> 'test_hybrid_out.new'
    expected_name = os.path.join(test_out_dir, 'test_hybrid.out.new')
    assert test_path_2 == expected_name


# ----- Test hybkit.util.set_settings -----
settings_info_dicts = [
    (hybkit.settings.HybRecord_settings, hybkit.settings.HybRecord_settings_info),
    (hybkit.settings.HybFile_settings, hybkit.settings.HybFile_settings_info),
    (hybkit.settings.FoldRecord_settings, hybkit.settings.FoldRecord_settings_info),
    (hybkit.settings.FoldFile_settings, hybkit.settings.FoldFile_settings_info),
    (hybkit.settings.HybFoldIter_settings, hybkit.settings.HybFoldIter_settings_info),
    (hybkit.settings.Analysis_settings, hybkit.settings.Analysis_settings_info),
]
test_parameters = []
for source, settings_info_dict in settings_info_dicts:
    for setting in settings_info_dict:
        default_value, description, type_str = settings_info_dict[setting][:3]
        short_flag, argparse_fields = settings_info_dict[setting][3:]
        if 'choices' in argparse_fields:
            good_choice = next(iter(argparse_fields['choices']))
            bad_choice_1, bad_choice_2 = 'invalid', ['invalid']
        else:
            good_choice = default_value
            bad_choice_1, bad_choice_2 = None, None
        test_parameters.append((source, setting, 'Pass', good_choice, settings_info_dict[setting]))
        if bad_choice_1 is not None:
            test_parameters.append((source, setting, 'HybkitArgError',
                                    bad_choice_1, settings_info_dict[setting]))
            test_parameters.append((source, setting, 'HybkitArgError',
                                    bad_choice_2, settings_info_dict[setting]))
use_parameters = ('source', 'setting', 'expectation', 'set_val', 'setting_props')

@pytest.mark.parametrize(use_parameters, [*test_parameters])
def test_util_set_settings(source, setting, expectation, set_val, setting_props):
    """Test hybkit.util.set_settings_from_namespace."""
    expect_context = get_expected_result_context(expectation)
    use_namespace = argparse.Namespace()
    setattr(use_namespace, setting, set_val)
    first_val = source[setting]
    with expect_context:
        hybkit.util.set_setting(setting, set_val, verbose=True)
        assert source[setting] == set_val
        hybkit.util.set_setting(setting, first_val, verbose=True)
    with expect_context:
        hybkit.util.set_settings_from_namespace(use_namespace, verbose=True)
        assert source[setting] == set_val
        hybkit.util.set_setting(setting, first_val, verbose=True)


# ----- Test hybkit.util parser generation
def test_util_validate_args():
    """Test hybkit.util parser generation."""
    original_abspath = hybkit.settings._USE_ABSPATH
    hybkit.settings._USE_ABSPATH = True

    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.out_hybs_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', test_hyb_file_name, test_hyb_file_name, '-o', test_hyb_file_name]
    )
    assert not hybkit.util.validate_args(args)
    assert not hybkit.util.validate_args(args, script_parser)
    with pytest.raises(SystemExit):
        hybkit.util.validate_args_exit(args, script_parser)

    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.in_folds_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', test_hyb_file_name, test_hyb_file_name, '-f', test_vienna_file_name]
    )
    assert not hybkit.util.validate_args(args)
    assert not hybkit.util.validate_args(args, script_parser)
    with pytest.raises(SystemExit):
        hybkit.util.validate_args_exit(args, script_parser)

    hybkit.settings._USE_ABSPATH = original_abspath

# if __name__ == '__main__':
#     test_type_finder()
