#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
hybkit TypeFinder Class.

This module contains the TypeFinder class to work with :class:`HybRecord` to
parse sequence identifiers to identify sequence type.
"""

import os
import types
from typing import Any, Callable, Dict, List, NoReturn, Optional, Union

from hybkit.errors import HybkitArgError, HybkitMiscError

# ----- Begin Typing Variables ----- #
SegProps = Dict[str, Union[float, int, str]]


# ----- Begin TypeFinder Class ----- #
class TypeFinder:
    """
    Class for parsing identifiers to identify sequence 'type'.

    Designed to be used by the :class:`hybkit.HybRecord`

    .. _TypeFinder-Attributes:

    Attributes:
        params (dict): Stored parameters for string parsing, where applicable.
    """

    # TypeFinder : Public Attributes
    #: Placeholder for storing active method, set with :meth:`set_method`
    #: (see :meth:`set_method` for details).
    find_with_params = None

    # TypeFinder : Public Attributes
    #: Placeholder for parameters for active method, set with :meth:`set_method`
    #: (see :meth:`set_method` for details).
    params = None

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #:   Default method assigned using :meth:`check_set_method`
    default_method = 'hybformat'

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #:   Dict of provided methods available to assign segment types
    #:
    #:     ============== ==================================
    #:     'hybformat'    :meth:`method_hybformat`
    #:     'string_match' :meth:`method_string_match`
    #:     'id_map'       :meth:`method_id_map`
    #:     ============== ==================================
    methods = {
        'hybformat': 'method_hybformat',
        'string_match': 'method_string_match',
        'id_map': 'method_id_map'
    }

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #: Dict of param generation methods for type finding methods
    #:
    #:     ============== ===================================================
    #:     'hybformat'    'N/A'
    #:     'string_match' :meth:`make_string_match_params`
    #:     'id_map'       :meth:`make_id_map_params`
    #:     ============== ===================================================
    param_methods = {
        'hybformat': None,
        'string_match': 'make_string_match_params',
        'id_map': 'make_id_map_params'
    }

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #: Dict of whether parameter generation methods need an input file
    #:
    #:     ============== ===================================================
    #:     'hybformat'    :obj:`False`
    #:     'string_match' :obj:`True`
    #:     'id_map'       :obj:`True`
    #:     ============== ===================================================
    param_methods_needs_file = {
        'hybformat': False,
        'string_match': True,
        'id_map': True,
    }

    # TypeFinder : Public Methods : Initialization
    # STUB, class is designed to be used with class-level functions.
    def __init__(self) -> NoReturn:
        """Class not used with initialization."""
        message = 'TypeFinder class not intended to be initialized for use.'
        raise RuntimeError(message)

    # TypeFinder : Public Classmethods : method
    @classmethod
    def set_method(
            cls,
            method: str,
            params: Optional[Dict[str, Any]] = None
            ) -> None:
        """
        Select method to use when finding types.

        Available methods are listed in :attr:`methods`.

        Args:
            method (str): Method option from :attr:`methods` to set
                for use as :meth:`find`.
            params (dict, optional): Dict object of parameters to use by set method.
        """
        if method not in cls.methods:
            message = 'Selected method: %s is not defined.\n' % method
            message += 'Allowed Options:' + ', '.join(cls.methods.keys())
            raise HybkitArgError(message)
        cls.find_with_params = getattr(cls, cls.methods[method])
        if params is not None:
            use_params = params
        else:
            use_params = {}
        cls.params = use_params

    # TypeFinder : Public Classmethods : method
    @classmethod
    def method_is_set(cls) -> bool:
        """
        Return whether a TypeFinder method has been set.

        Methods should be set with :meth:`set_method`.

        Returns:
            bool: True if a method has been set, False otherwise.
        """
        return cls.find_with_params is not None

    # TypeFinder : Public Classmethods : method
    @classmethod
    def check_set_method(cls) -> None:
        """If no TypeFinder method set, set as :attr:`default_method`."""
        if not cls.method_is_set():
            cls.set_method(cls.default_method)

    # TypeFinder : Public Classmethods : method
    @classmethod
    def find(cls, seg_props: SegProps) -> Optional[str]:
        """
        Find type of segment using :meth:`TypeFinder.find_custom_method`.

        If a TypeFinder method has been set with :meth:`set_method`.
        use the current parameters set in
        :attr:`params` to find the type of the provided segment by calling::

            seg_type = :meth:`TypeFinder.find_custom_method`(seg_props, :attr`TypeFinder.params`)

        Args:
            seg_props (dict): :obj:`seg_props` from :class:`hybkit.HybRecord`

        Returns:
            str: Type of the provided segment, or None if a type cannot be identified.
        """
        if cls.find_with_params is None:
            message = 'TypeFinder method has not been set.'
            raise RuntimeError(message)
        return cls.find_with_params(seg_props, cls.params)

    # TypeFinder : Public Classmethods : method
    @classmethod
    def set_custom_method(
            cls,
            method: Callable,
            params: Optional[dict] = None
            ) -> None:
        """
        Set the method for use to find seg types.

        This method is for providing a custom function. To use the included functions,
        use :meth:`set_method`.
        Custom functions provided must have the signature::

            seg_type = custom_method(self, seg_props, params)

        This function should return the string of the assigned segment type if found, or a
        None object if the type cannot be found.
        It can also take a dictionary in the "params" argument that specifies
        additional or dynamic search properties, as desired.

        Args:
            method (method): Method to set for use.
            params (dict, optional): dict of custom parameters to set for use.
        """
        cls.find_with_params = types.MethodType(method, cls)
        if params is not None:
            cls.params = params
        else:
            cls.params = {}

    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def method_hybformat(
            seg_props: SegProps,
            params: Optional[dict] = None,
            ) -> Optional[str]:
        """
        Return the type of the provided segment, or None if segment cannot be identified.

        This method works with sequence / alignment mapping identifiers
        in the format of the reference database provided by the Hyb Software Package,
        specifically identifiers of the format::

            <gene_id>_<transcript_id>_<gene_name>_<seg_type>

        This method returns the last component of the identifier,
        split by "_", as the identified sequence type.
        (returns :obj:`None` if the segment identifier does not contain "_").

        Example:
            ::

                "MIMAT0000076_MirBase_miR-21_microRNA"  --->  "microRNA".

        Args:
            seg_props (dict): :obj:`seg_props` from :class:`hybkit.HybRecord`
            params (dict, optional): Unused in this method.
        """
        if params is not None and params:
            message = 'method_hybformat does not use params, but params were provided:\n'
            message += str(params)
            raise HybkitArgError(message)
        if '_' not in seg_props['ref_name']:
            return None
        else:
            split_id = seg_props['ref_name'].split('_')
            return split_id[-1]

    # TypeFinder : Public Staticmethods : methods
    @staticmethod
    def method_string_match(
            seg_props: SegProps,
            params: Optional[dict] = None,
            ) -> Optional[str]:
        """
        Return the type of the provided segment, or None if unidentified.

        This method attempts to find a string matching a specific pattern within the identifier
        of the aligned segment. Search options include "startswith", "contains", "endswith", and
        "matches", and returns the first type matching the criteria.
        The required params dict should contain a key for each desired
        search type, with a list of 2-tuples for each search-string with assigned-type.

        Example:
            ::

                params = {'endswith': [('_miR', 'microRNA'),
                                       ('_trans', 'mRNA')   ]}

        This dict can be generated with the associated :meth:`make_string_match_params`
        method and an associated csv legend file with format::

            #comment line
            #search_type,search_string,seg_type
            endswith,_miR,microRNA
            endswith,_trans,mRNA

        Args:
            seg_props (dict): :class:`~hybkit.HybRecord` segment properties dict
                to evaluate.
            params (dict, optional): Dict with search parameters as described above.
        """
        if params is None or not params:
            message = 'method_string_match requires params, but none were provided.'
            raise HybkitArgError(message)
        seg_name = seg_props['ref_name']
        found_type = None
        if 'startswith' in params and not found_type:
            for search_string, search_type in params['startswith']:
                if seg_name.startswith(search_string):
                    found_type = search_type
                    break
        if 'contains' in params and not found_type:
            for search_string, search_type in params['contains']:
                if search_string in seg_name:
                    found_type = search_type
                    break
        if 'endswith' in params and not found_type:
            for search_string, search_type in params['endswith']:
                if seg_name.endswith(search_string):
                    found_type = search_type
                    break
        if 'matches' in params and not found_type:
            for search_string, search_type in params['matches']:
                if search_string == seg_name:
                    found_type = search_type
                    break

        return found_type

    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def make_string_match_params(legend_file: str) -> dict:
        """
        Read csv and return a dict of search parameters for :meth:`method_string_match`.

        The my_legend.csv file should have the format::

            #comment line
            #search_type,search_string,seg_type
            endswith,_miR,microRNA
            endswith,_trans,mRNA

        Search_type options include "startswith", "contains", "endswith", and "matches"
        The produced dict object contains a key for each search type, with a list of
        2-tuples for each search-string and associated segment-type.

        For example::

            {'endswith': [('_miR', 'microRNA'),
                          ('_trans', 'mRNA')   ]}

        """
        allowed_search_types = {'startswith', 'contains', 'endswith', 'matches'}
        return_dict = {}
        if not isinstance(legend_file, (str)):
            message = 'legend_file must be a string, not %s' % type(legend_file)
            raise HybkitArgError(message)
        if not os.path.isfile(legend_file):
            message = 'File: %s for make_string_match_params() method ' % legend_file
            message += 'not found.'
            raise FileNotFoundError(message)

        with open(legend_file, 'r') as legend_file_obj:
            for line in legend_file_obj:
                # Skip Blank Lines
                if not line.split():
                    continue
                # Skip Commented Lines
                if line.lstrip().startswith('#'):
                    continue
                use_line = line.strip()
                split_line = use_line.split(',')
                if len(split_line) != 3: # noqa: PLR2004
                    message = f'Error reading legend line: \n{use_line!s}\n{split_line!s}'
                    message += '\nThree comma-separated entries expected.'
                    raise RuntimeError(message)
                search_type, search_string, seg_type = split_line
                if search_type not in allowed_search_types:
                    message = f'Read Search type: "{search_type}"\n'
                    message += 'Not in allowed types: {}'.format(', '.join(allowed_search_types))
                    message += f'\nFor legend line: \n{use_line!s}\n'
                    raise RuntimeError(message)

                if search_type not in return_dict:
                    return_dict[search_type] = []

                return_dict[search_type].append((search_string, seg_type))

        return return_dict

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    @staticmethod
    def method_id_map(
            seg_props: SegProps,
            params: dict = None,
            ) -> Optional[str]:
        """
        Return the type of the provided segment or None if it cannot be identified.

        This method checks to see if the identifier of the segment is present in the params dict.
        params should be formatted as a dict with keys as
        sequence identifier names, and the corresponding type as the respective values.

        Example:
            ::

                params = {'MIMAT0000076_MirBase_miR-21_microRNA': 'microRNA',
                          'ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA': 'rRNA'}

        This dict can be generated with the associated :meth:`make_id_map_params` method.

        Args:
            seg_props (dict): :class:`~hybkit.HybRecord` segment properties dict
                to evaluate.
            params (dict): Dict of mapping of sequence identifiers to sequence types.

        Returns:
           str: Identified sequence type, or None if it cannot be found.

        """
        seg_name = seg_props['ref_name']
        if params is None:
            params = {}
        if seg_name in params:
            return params[seg_name]
        else:
            return None

    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def make_id_map_params(mapped_id_files: List[str]) -> dict:
        """
        Read file(s) into a mapping of sequence identifiers.

        This method reads one or more files into a dict for use with the
        :meth:`method_id_map` method.
        The method requires passing a file path (or list/tuple of file paths)
        of mapped_id_files.
        Files listed in the mapped_id_files argument should have the format::

            #comment line
            #seg_id,seg_type
            segA_unique_id,segA_type
            segB_unique_id,segB_type

        Args:
            mapped_id_files (:obj:`str`, :obj:`list`, or :obj:`tuple`): Iterable
                object containing strings of paths
                to files containing id/type mapping information.

        """
        return_dict = {}
        if not isinstance(mapped_id_files, (str, list, tuple)):
            message = 'arguments passed to mapped_id_files and type_file_pairs must be '
            message += 'provided as a list or tuple.\n  Provided: "%s"' % str(mapped_id_files)
            raise TypeError(message)
        if isinstance(mapped_id_files, str):
            mapped_id_files = [mapped_id_files]

        for mapped_id_file in mapped_id_files:
            # Check if file not exists and raise error
            if not os.path.isfile(mapped_id_file):
                message = 'File: %s for make_id_map_params() method not found.' % mapped_id_file
                raise FileNotFoundError(message)
            with open(mapped_id_file, 'r') as mapped_id_file_obj:
                for raw_line in mapped_id_file_obj:
                    # Skip Blank Lines
                    line = raw_line.strip()
                    if not line.split():
                        continue
                    # Skip Commented Lines
                    if line.startswith('#'):
                        continue
                    split_line = line.split(',')
                    if len(split_line) != 2:  # noqa: PLR2004
                        message = 'Error reading mapped-id line: '
                        message += f'\n{line!s}\n{split_line!s}'
                        message += '\nTwo comma-separated entries expected.'
                        raise HybkitMiscError(message)
                    seq_id, seg_type = split_line

                    if seq_id in return_dict and seg_type != return_dict[seq_id]:
                        message = 'Conflicting types assigned for sequence id: %s\n' % seq_id
                        message += f'  {return_dict[seq_id]}  |  {seg_type}'
                        raise HybkitMiscError(message)
                    else:
                        return_dict[seq_id] = seg_type

        return return_dict

    # TypeFinder : Private classmethods : find_seg_type
    @classmethod
    def _reset(cls) -> None:
        """
        Reset the class to its initial state.

        This method is used for testing purposes.
        """
        cls.find_with_params = None
        cls.params = None
