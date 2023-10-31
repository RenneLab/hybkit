#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Helper functions for automatic testing of hybkit code.
"""

from contextlib import nullcontext as does_not_raise

import pytest

# import test_helper_data

# import hybkit

# Get expected result string for exception testing.
def get_expected_result_string(is_allowed=False):
    """Return string identifying expected pass/error result."""
    if is_allowed:
        return 'Pass'
    else:
        return 'Raise'


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
    else:
        message = 'Expected result string must be "Pass" or "Raise".'
        raise ValueError(message)
