#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
This module contains the TypeFinder class to work with :class:`HybRecord` to 
parse sequence identifiers to idenfity sequence type.
"""

import os
import csv

class RegionFinder(object):
    """
    Class for storing and using parsing methods to assign a sequence "type" by parsing a 
    sequence identifier.

    Designed to be used by the :class:`hybkit.HybRecord`

    .. _RegionFinder-Attributes:

    Attributes:
        region_info (dict): Sequence coding region information.
           
            Must have form {'seqid':info} where info has keys: 
            "cdna_coding_start", "cdna_coding_end" with string / interger locus values.
            Example::

                    region_info = {'ENST00000372098': {'cdna_coding_start':'45340255',
                                                       'cdna_coding_end':'45340388'}}

    """

    #: Dict of information required for use by :func:`find`.
    #: Populate dict using :func:`make_set_region_info`. 
    region_info = {}

    # RegionFinder : Public Methods : Initialization
    def __init__(self):
        print('WARN: RegionFinder class not intended to be initialized for use.')

    # RegionFinder : Public Classmethods : Prep Region Info
    @classmethod
    def find(cls, ref_name, ref_start, ref_end):
        """
        For miRNA/coding-target pairs, find the region of the coding transcript targeted.

        # Method is temporarily disabled pending conceptual review.

        The analysis requires :attr:`region_finder.RegionFinder.region_info` 
        to be populated by using 
        the :func:`set_region_info` method, with information prepared using 
        the :func:`make_region_info` method.

        Example:
            ::

                region_info = {'ENST00000372098': {'cdna_coding_start':'45340255',
                                                   'cdna_coding_end':'45340388'}}

        Args:
            ref_name (str): Sequence identifier for the sequence.
            ref_start (str or int): Base position of the start of the sequence.
            ref_end (str or int): Base position of end of the sequence.
        """
        raise NotImplementedError()
        if cls.region_info == {}:
            message = 'RegionFinder.region_info attribute has not been populated.\n'
            message += '  Please fill this variable with "set_region_info()" before use.'
            print(message)
            raise Exception(message)

        if ref_name not in region_info:
            ref_region = None
        else:
            ref_details = region_info[ref_name]
            required_keys = ['cdna_coding_start', 'cdna_coding_end']
            if not all(key in ref_details for key in required_keys):
                message = 'Problem with finding target for reference: %s\n' % ref_name
                message += 'Provided target region information: \n%s\n' % str(ref_details)
                message += 'Does not have required keys:\n    %s' % ', '.join(required_keys)
                print(message)
                raise Exception(message)

            # Setup Variables
            ref_middle = ref_start + ((ref_end - ref_start) // 2)
            ref_cds_start = int(target_regions['cdna_coding_start'])
            ref_cds_end = int(target_regions['cdna_coding_end'])

            if ref_cds_start >= ref_cds_end:
                message = 'Reverse oreintation reference: %s\n' % self.target_props['ref']
                message += '    "Start:", %i,  "End:" %i' % (ref_cds_start, ref_cds_end)
                print(message)
                raise Exception(message)

            if ref_middle < ref_cds_start:
                ref_region = '5pUTR'
            elif ref_middle > ref_cds_end:
                ref_region = '3pUTR'
            else:
                ref_region = 'coding'

            #print('cds:', ref_cds_start, ref_cds_end, 'seg:', ref_start, ref_end)
            #print(target_reg)
            #input()

            self.mirna_details['target_reg'] = target_reg 
            self.set_flag('target_reg', target_reg)           
 
        return ref_region


    # RegionFinder : Public Classmethods : Prep Region Info
    @classmethod
    def make_region_info(cls, region_csv_name, sep=','):
        """
        Return dict with information on coding transcript utr regions from an input csv.

        The input csv must contain a header line, and must contain the columns::

            identifier,cdna_coding_start,cdna_coding_end

        Example:
            Example return dict object::

                region_info = {'ENST00000372098': {'cdna_coding_start':'45340255',
                                                   'cdna_coding_end':'45340388'}}

        The return dict can then be passed to :func:`set_region_info` or supplied directly to
        the :func:`find` method.

        Args:
            region_csv_name (str): String of path to csv file to read information from.
            sep (str, optional): Separator for columns of input delimited file. (Default: ',')
        """

        index_keys = ['cdna_coding_start', 'cdna_coding_end']
        required_keys = ['identifier'] + index_keys

        if not os.path.isfile(region_csv_name):
            message = 'Problem with creation of target region information dict.\n'
            message += 'Provided input target region csv file:\n    %s\n' % region_csv_name
            message += 'cannot be found.'
            print(message)
            raise Exception(message)

        ret_dict = {}

        with open(region_csv_name, 'r', newline='') as region_csv:
            reader = csv.DictReader(region_csv, delimiter=sep) 
            for i, line_items in enumerate(reader, start=1):
                # If line does not contain required format or information, raise error.
                if (not all(key in line_items for key in required_keys)
                    or any (line_items[key] is None for key in required_keys)):
                    message = 'Problem with creation of target region information dict.\n'
                    message += 'Provided input target region csv file:\n    %s\n' % region_csv_name
                    message += 'Does not contain the required information at line: '
                    message += '%i\nRequired Keys: %s\n' % (i, ', '.join(required_keys))
                    print(message)
                    raise Exception(message)

                seq_id = line_items['identifier']
            
                # If information has already been read for line, raise error.
                if seq_id in ret_dict:
                    message = 'Problem with creation of target region information dict.\n'
                    message += 'Provided input target region csv file:\n    '
                    message +=  region_csv_name + '\n'
                    message += 'Contains duplicate information for entry:'
                    message += '%s\nProvided at line %i\n' % (seq_id, i)
                    message += 'Previous: %s\n' % str(ret_dict[seq_id])
                    message += 'Current: %s\n' % str(line_items)
                    print(message)
                    raise Exception(message)
                
                data_pairs = [(key, line_items[key]) for key in index_keys]
                for key, item in data_pairs:
                    if not item.isnumeric():
                        message = 'Data for index entry: %s contains non-numeric values.' % seq_id
                        message += 'Data: %s' % str(data_pairs)
                        print(message)
                        raise Exception(message)
                ret_dict[seq_id] = {key:int(line_items[key]) for key in index_keys}                

        return ret_dict

    # RegionFinder : Public Classmethods : Prep Region Info
    @classmethod
    def set_region_info(cls, region_info):
        """Set :attr:`region_info_dict` with information on coding transcript UTR regions. 

        This dict must have transcript identifiers as keys, with values of dicts with 
        containing (at least): cdna_coding_start, cdna_coding_end

        Example:
           ::

                region_info = {'ENST00000372098': {'cdna_coding_start':'45340255',
                                                   'cdna_coding_end':'45340388'}}

        Args:
            region_info (dict): Dict of region information to set as 
                :attr:`region_info`.
        """
        cls.region_info = region_info
    
    # RegionFinder : Public Classmethods: find
    @classmethod
    def make_set_region_info(cls, region_csv_name, sep=','):
        """
        Calling :func:`make_region_info` then :func:`set_region_info`.

        Args:
            region_csv_name (str): String of path to csv file to read information from.
            sep (str, optional): Separator for columns of input delimited file. (Default: ',')
        """
        region_info = cls.make_region_info(region_csv_name, sep=sep)
        cls.set_region_info(region_info)

