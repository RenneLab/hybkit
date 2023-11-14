#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit ViennaFile Class.
"""

import os

import pytest

import hybkit
from auto_tests.test_helper_data import ART_HYB_VIENNA_PROPS_1, ART_HYB_VIENNA_PROPS_2, test_out_dir
from auto_tests.test_helper_functions import get_expected_result_context
from hybkit.errors import HybkitMiscError

# ----- Linting Directives:
# ruff: noqa: SLF001 ARG001

hybkit.util.set_setting('error_mode', 'raise')
hybkit.util.set_setting('iter_error_mode', 'raise')


# ----- Begin ViennaFile Tests -----
test_parameters = [
    ('bad_type', 'HybkitArgError', {'seq_type': 'badtype'}),
    ('bad_error_mode', 'HybkitArgError', {'error_mode': 'badmode'}),
]


@pytest.mark.parametrize(('test_name', 'expectation', 'test_kwargs'), [*test_parameters])
# ----- Test Misc Properties of ViennaFiles -----
def test_viennafile_constructor_misc(test_name, expectation, test_kwargs, tmp_path):
    """Test miscellaneous properties of ViennaFiles."""
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    expect_context = get_expected_result_context(expectation)
    with expect_context:
        _vienna_file = hybkit.ViennaFile.open(
            vienna_autotest_file_name, 'w',
            **test_kwargs
        )


# ----- Test IO of Vienna Strings / FoldRecords -----
test_parameters = []
for prop_set in [ART_HYB_VIENNA_PROPS_1, ART_HYB_VIENNA_PROPS_2]:
    if prop_set['overlapping']:
        test_name = 'Overlapping'
    else:
        test_name = 'Static'
    test_parameters.append(
        (test_name, 'Pass', prop_set)
    )


@pytest.mark.parametrize(('test_name', 'expectation', 'test_props'), [*test_parameters])
def test_viennafile_io(test_name, expectation, test_props, tmp_path):
    """Test reading and writing of ViennaFile objects."""
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    expect_context = get_expected_result_context(expectation)
    vienna_str = test_props['vienna_str']
    with open(vienna_autotest_file_name, 'w') as vienna_autotest_file:
        vienna_autotest_file.write(vienna_str)

    assert hybkit.util.vienna_exists(vienna_autotest_file_name)

    with expect_context:
        with hybkit.ViennaFile.open(vienna_autotest_file_name, 'r') as vienna_autotest_file:
            for vienna_record in vienna_autotest_file:
                gen_vienna_str = vienna_record.to_vienna_string()
                assert gen_vienna_str == vienna_str

    if expectation.lower() == 'pass':
        with hybkit.ViennaFile.open(vienna_autotest_file_name, 'r') as vienna_autotest_file:
            first_record = vienna_autotest_file.read_record()
        with hybkit.ViennaFile.open(vienna_autotest_file_name, 'r') as vienna_autotest_file:
            all_records = vienna_autotest_file.read_records()
        assert first_record == all_records[0]
        with hybkit.ViennaFile(vienna_autotest_file_name, 'w') as vienna_autotest_file:
            vienna_autotest_file.write_records([vienna_record, vienna_record])
            vienna_autotest_file.write_record(vienna_record)
            vienna_autotest_file.write_fh(vienna_str)


# ----- Test ViennaFile Misc -----
def test_viennafile_misc(tmp_path):
    """Test miscellaneous properties of ViennaFiles."""
    vienna_autotest_file_name = os.path.join(tmp_path, 'vienna_autotest_file.vienna')
    vienna_file = hybkit.ViennaFile.open(vienna_autotest_file_name, 'w')
    with pytest.raises(HybkitMiscError):
        vienna_file._ensure_foldrecord(None)
