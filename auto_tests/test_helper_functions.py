#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Helper functions for automatic testing of hybkit code.
"""

from contextlib import nullcontext as does_not_raise

import pytest
from test_helper_data import ERROR_TYPE_STRINGS

# import test_helper_data

# import hybkit

# Get expected result string for exception testing.
def get_expected_result_string(is_allowed=False, err_string='Raise'):
    """Return string identifying expected pass/error result."""
    if is_allowed:
        return 'Pass'
    else:
        return err_string


# Get expected result context for exception testing.
def get_expected_result_context(expect_str, error_types=(TypeError, RuntimeError)):
    """Return context for testing allowed types."""
    if expect_str.lower() == 'pass':
        return does_not_raise()
    elif expect_str.lower() == 'raise':
        if isinstance(error_types, list):
            error_types = tuple(error_types)
        elif not isinstance(error_types, tuple):
            error_types = (error_types,)
        return pytest.raises(error_types)
    elif expect_str.lower() in ERROR_TYPE_STRINGS:
        return pytest.raises(ERROR_TYPE_STRINGS[expect_str.lower()])
    else:
        message = 'Expected result string must be "Pass", "Raise", or one of: '
        message += f'{ERROR_TYPE_STRINGS.keys()}'
        raise ValueError(message)
