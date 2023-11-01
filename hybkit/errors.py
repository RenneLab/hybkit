#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

r"""
Module storing hybkit error classes.
"""

# ----- Begin Error Classes -----
class HybkitError(Exception):
    """
    Base class for Hybkit errors.

    Attributes:
        message (str): Human-readable string describing the error.
    """


class HybkitArgError(HybkitError):
    """
    Error raised when an invalid argument is provided to a Hybkit function.

    Subclass of :class:`HybkitError`.

    Attributes:
        message (str): Human-readable string describing the error.
    """


class HybkitConstructorError(HybkitError):
    """
    Error raised when a read error occurs.

    Subclass of :class:`HybkitError`.

    Attributes:
        message (str): Human-readable string describing the error.
    """


class HybkitIterError(HybkitError):
    """
    Error raised when an error is encountered during Hybkit iteration.

    Subclass of :class:`HybkitError`.

    Attributes:
        message (str): Human-readable string describing the error.
    """


class HybkitMiscError(HybkitError):
    """
    Error raised when an error is encountered during Hybkit usage.

    Subclass of :class:`HybkitError`.

    Attributes:
        message (str): Human-readable string describing the error.
    """
