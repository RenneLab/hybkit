#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
This module contains the TypeFinder class to work with :class:`HybRecord` to 
parse sequence identifiers to idenfity sequence type.
"""

import os
import io
import types
import csv
import copy
    

class TypeFinder(object):
    """
    Class for storing and using parsing methods to assign a sequence "type" by parsing a 
    sequence identifier.

    Designed to be used by the :class:`hybkit.HybRecord`

    .. _TypeFinder-Attributes:

    Attributes:
        params (dict): Stored parameters for string parsing, where applicable.
    """

    # TypeFinder : Public Attributes
    #: Placeholder for storing active method, set with :meth:`set_method`.
    find = None

    # TypeFinder : Public Methods : Initialization
    # STUB, class is designed to be used with class-level functions.
    def __init__(self):
        print('WARN: TypeFinder class not intended to be initialized for use.')

    # TypeFinder : Public Classmethods : method
    @classmethod
    def set_method(cls, method, params={}):
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
        cls.find = cls.methods[method]
        cls.params = params

    # TypeFinder : Public Classmethods : method
    @classmethod
    def set_custom_method(cls, method, params={}):
        """
        Set the method for use to find seg types.
        
        This method is for providing a custom function. To use the included functions, 
        use :meth:`set_method`.
        Custom functions provided must have the signature::

            seg_type = custom_method(self, seg_props, params, check_complete)
        
        This function should return the string of the assigned segment type if found, or a
        None object if the type cannot be found.
        It can also take a dictionary in the "params" argument that specifies
        additional or dynamic search properties, as desired.
        The if check_complete is true, the function should search for all possibilities for
        a given sequence, instead of stopping after the first is found.

        Args:
            method (method): Method to set for use.
            params (dict, optional): dict of custom parameters to set for use.
        """

        cls.find = types.MethodType(method, cls)
        cls.params = params


    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def method_hybformat(seg_props, params={}, check_complete=False):
        """Return the type of the provided segment, or None if segment cannot be identified.

        This method works with sequence / alignment mapping identifiers
        in the format of the reference database provided by the Hyb Software Package,
        specifically identifiers of the format:: 

            <gene_id>_<transcript_id>_<gene_name>_<seg_type>

        This method returns the last component of the identifier, 
        split by "_", as the identfied sequence type.

        Example:
            ::

                "MIMAT0000076_MirBase_miR-21_microRNA"  --->  "microRNA".

        Args:
            seg_props (dict): :obj:`seg_props` from :class:`hybkit.HybRecord`
            params (dict, optional): Unused in this method.
            check_complete (bool, optional): Unused in this method.

        """
        split_id = seg_props['ref_name'].split('_')
        return split_id[-1]


    # TypeFinder : Public Staticmethods : methods
    @staticmethod
    def method_string_match(seg_props, params={}, check_complete=False):
        """Return the type of the provided segment, or None if unidentified.
        
        This method attempts to find a string matching a specific pattern within the identifier
        of the aligned segment. Search options include "startswith", "contains", "endswith", and
        "matches". The required params dict should contain a key for each desired
        search type, with a list of 2-tuples for each search-string with assigned-type.

        Example:
            ::

                params = {'endswith': [('_miR', 'microRNA'),
                                       ('_trans', 'mRNA')   ]}

        This dict can be generated with the associated :meth:`make_string_match_params`
        method and an associated csv legend file with format::

            #commentline
            #search_type,search_string,seg_type
            endswith,_miR,microRNA
            endswith,_trans,mRNA

        Args:
            params (dict, optional): Dict with search paramaters as described above.
            check_complete (bool, optional): If true, the method will continue checking search
                options after an option has been found, to ensure that no options conflict
                (more sure method). If False, it will stop after the first match is found 
                (faster method). (Default: False)
        """
        seg_name = seg_props['ref_name']
        found_types = set()
        check_done = False
        if not check_done and 'startswith' in params:
            for search_string, search_type in params['startswith']:
                if seg_name.startswith(search_string):
                    found_types.add(search_type)
                    if not check_complete:
                        check_done = True
                        break
        if not check_done and 'contains' in params:
            for search_string, search_type in params['contains']:
                if search_string in seg_name:
                    found_types.add(search_type)
                    if not check_complete:
                        check_done = True
                        break
        if not check_done and 'endswith' in params:
            for search_string, search_type in params['endswith']:
                if seg_name.endswith(search_string):
                    found_types.add(search_type)
                    if not check_complete:
                        check_done = True
                        break
        if not check_done and 'matches' in params:
            for search_string, search_type in params['matches']:
                if search_string == seg_name:
                    found_types.add(search_type)
                    if not check_complete:
                        check_done = True
                        break

        if not found_types:
            return None
        elif len(found_types) == 1:
            return next(iter(found_types))
        elif len(found_types) > 1:
            message = 'Multiple sequence types found for item: %s' % seg_name
            message += '  ' + ', '.join(sorted(list(found_types)))
            print(message)
            raise Exception(message)

    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def make_string_match_params(
            legend_file
            ):
        """Read csv and return a dict of search parameters for :meth:`method_string_match`.

        The my_legend.csv file should have the format::

            #commentline
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

        ALLOWED_SEARCH_TYPES = {'startswith', 'contains', 'endswith', 'matches'}
        return_dict = {}
        with open(legend_file, 'r') as legend_file_obj:
            for line in legend_file_obj:
                # Skip Blank Lines
                if not line.split():
                    continue
                # Skip Commented Lines
                if line.lstrip().startswith('#'):
                    continue
                line = line.rstrip()
                split_line = line.split(',')
                if len(split_line) != 3:
                    message = 'Error reading legend line: \n%s\n%s' % (str(line), str(split_line))
                    message += '\nThree comma-separated entries expected.'
                    print(message)
                    raise Exception(message)
                search_type = split_line[0]
                search_string = split_line[1]
                seg_type = split_line[2]
                if search_type not in ALLOWED_SEARCH_TYPES:
                    message = 'Read Search type: "%s"\n' % search_type
                    message += 'Not in allowed types: %s' % ', '.join(ALLOWED_SEARCH_TYPES)
                    message += '\nFor legend line: \n%s\n' % (str(line))
                    print(message)
                    raise Exception(message)

                if search_type not in return_dict:
                    return_dict[search_type] = []

                return_dict[search_type].append((search_string, seg_type))

        return return_dict

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    @staticmethod
    def method_id_map(seg_props, params={}, check_complete=False):
        """Return the type of the provided segment or None if it cannot be identified.

        This method checks to see if the identifer of the segment is present in the params dict.
        params should be formatted as a dict with keys as 
        sequence identifier names, and the corresponding type as the respective values.

        Example:
            ::

                params = {'MIMAT0000076_MirBase_miR-21_microRNA': 'microRNA',
                          'ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA': 'rRNA'}

        This dict can be generated with the associated :meth:`make_id_map_params` method.

        Args:
            params (dict): Dict of mapping of sequence identifiers to sequence types.
            check_complete (bool, optional): Unused in this method.
     
        Returns:
           str: Identified sequence type, or None if it cannot be found.

        """
        seg_name = seg_props['ref_name']
        if seg_name in params:
            return params[seg_name]
        else:
            return None        

    # TypeFinder : Public Staticmethods : find_seg_type
    @staticmethod
    def make_id_map_params(mapped_id_files=None, type_file_pairs=None):
        """
        Read file(s) into a mapping of sequence identifiers.

        This method reads one or more files into a dict for use with the 
        :meth:`method_id_map` method.
        The method requires passing either a list/tuple of one or more files to mapped_id_files,
        or a list/tuple of one or more pairs of file lists and file types 
        passed to type_file_pairs.
        Files listed in the mapped_id_files argument should have the format:: 

            #commentline
            #seg_id,seg_type
            seg1_unique_id,seg1_type
            seg2_unique_id,seg2_type

        Entries in the list/tuple passed to type_file_pairs should have the format: 
        (seg1_type, file1_name)

        Example:
            ::

                [(seg1_type, file1_name), (seg2_type, file2_name),]

        The first entry in each (non-commented, non-blank) file line will be read and
        added to the mapping dictionary mapped to the provided seg_type.

        Args:
            mapped_id_files (list or tuple, optional) Iterable object containing strings of paths
                to files containing id/type mapping information.
            type_file_pairs (list or tuple, optional) Iterable object containing 2-tuple pairs
                containing id/type mapping information.
     
        """

        return_dict = {}
        if not any((arg is not None for arg in (mapped_id_files, type_file_pairs))):
            message = 'make_seg_type_id_map function requires either a mapped_id_files '
            message += 'or type_file_pairs argument.'
            print(message)
            raise Exception(message)
        for argument in mapped_id_files, type_file_pairs:
            if ((argument is not None) and not any((isinstance(argument, allowed_type) 
                                                     for allowed_type in (list, tuple)))):
                message = 'arguments passed to mapped_id_files and type_file_pairs must be '
                message += 'provided as a list or tuple.\n  Current passed aruement: '
                message += str(argument)
                print(message)
                raise Exception(message)

        if mapped_id_files is not None:
            for mapped_id_file in mapped_id_files:
                with open(mapped_id_file, 'r') as mapped_id_file_obj:
                    for line in mapped_id_file_obj:
                        # Skip Blank Lines
                        if not line.split():
                            continue
                        # Skip Commented Lines
                        if line.lstrip().startswith('#'):
                            continue
                        line = line.rstrip()
                        split_line = line.split(',')
                        if len(split_line) != 2:
                            message = 'Error reading mapped-id line: '
                            message += '\n%s\n%s' % (str(line), str(split_line))
                            message += '\nTwo comma-separated entries expected.'
                            print(message)
                            raise Exception(message)
                        seq_id = split_line[0]
                        seg_type = split_line[1]
        
                        if seq_id in return_dict and seg_type != return_dict[seq_id]:
                            message = 'Conflicting types assigned for sequence id: %s\n' % seq_id
                            message += '  %s  |  %s' % (return_dict[seq_id], seg_type)
                            print(message)
                            raise Exception(message)
                        else:
                            return_dict[seq_id] = seg_type
    
        if type_file_pairs is not None: 
            for seg_type, id_file in type_file_pairs:
                with open(id_file, 'r') as id_file_obj:
                    for line in id_file_obj:
                        # Skip Blank Lines
                        if not line.split():
                            continue
                        # Skip Commented Lines
                        if line.lstrip().startswith('#'):
                            continue
                        seq_id = line.strip().split()[0]
        
                        if seq_id in return_dict and seg_type != return_dict[seq_id]:
                            message = 'Conflicting types assigned for sequence id: %s\n' % seq_id
                            message += '  %s  |  %s' % (return_dict[seq_id], seg_type)
                            print(message)
                            raise Exception(message)
                        else:
                            return_dict[seq_id] = seg_type
    
        return return_dict

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #:   Dict of provided methods available to assign segment types
    #:   
    #:     ============== ==================================
    #:     'hybformat'    :meth:`method_hybformat`
    #:     'string_match' :meth:`method_string_match`
    #:     'id_map'       :meth:`method_id_map`
    #:     ============== ==================================
    methods = {'hybformat': method_hybformat,
               'string_match': method_string_match,
               'id_map': method_id_map}

    # TypeFinder : Public Methods : Flag_Info : find_seg_type
    #:   Dict of param generation methods for type finding methods
    #:   
    #:     ============== ===================================================
    #:     'hybformat'    :obj:`None`
    #:     'string_match' :meth:`make_string_match_params`
    #:     'id_map'       :meth:`make_id_map_params`
    #:     ============== ===================================================
    parameter_methods = {
        'hybformat': None,
        'string_match': make_string_match_params.__func__,
        'id_map': make_id_map_params.__func__,
        }

