#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Classes and Functions for manipulating data in the ".hyb" genomic sequence format.
Public classes and methods are imported by hybkit/__init__.py so they are accessible
as hybkit.HybRecord() ... etc.
'''

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __depreciated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import os
import io
import types

class HybRecord(object):
    '''
    Class for storing information about a chimeric or hybrid sequence read from a .hyb format file.
    A minimum amount of data necessary for a HybRecord object is the genomic sequence and its
    corresponding identifier.
    '''

    # Class-level constants
    # Columns 1-3 defining parameters of the overall hybrid, defined by the Hyb format
    HYBRID_COLUMNS = [
                      'id',      # str, Hybrid read identifier, ex: "1257_12"
                      'seq',     # str, Hybrid nucleotide sequence, ex: "ATCGGCTAATCGGTCA..."
                      'energy',  # str(of float), Intra-hybrid folding energy, ex: "-11.33"
                     ]
    # Columns 4-9 and 10-15, reespectively, defining parameters of each respective segment mapping,
    # defined by the Hyb format
    SEGMENT_COLUMNS = [
                       'ref',         # str, Mapping Reference Identity: ex:
                                      #   "MIMAT0000076_MirBase_miR-21_microRNA"
                       'read_start',  # int, Start position of mapping within hybrid, ex: "0"
                       'read_end',    # int, End position of mapping within hybrid, ex: "21"
                       'ref_start',   # int, Start position of mapping within reference, ex: "1723"
                       'ref_end',     # int, End position of mapping within reference, ex: "1744"
                       'score',       # str(of int/float) Alignment Score, can be BLAST e-score
                                      #   or mapping alignment score, depending on analysis
                                      #   implementation type.
                      ]
    # Arbitrary details included in column 16 of the hyb format in the form:
    #   “feature1=value1; feature2=value2;..."
    #   Flags utilized in the Hyb software package
    HYB_FLAGS = [
                 'count_total',            # int, total represented hybrids
                 'count_last_clustering',  # int, total represented hybrids at last clustring
                 'two_way_merged',         # "0" or "1", boolean representation of whether
                                           #   entries with mirrored 5' and 3' hybrids were merged
                 'seq_IDs_in_cluster',     # str, comma-separated list of all ids of hybrids
                                           #   merged into this hybrid entry.
                ]
    # Additional flag specifications utilized by hybkit
    HYBKIT_FLAGS = [
                    'orient',     # str, orientation of strand. Options:
                                  #   "F" (Forward), "IF" (Inferred Forward),
                                  #   "R" (Reverse), "IR" (Inferred Reverse),
                                  #   "U" (Unknown), or "C" (Conflicting)
                    'seg1_type',  # str, assigned "type" of segment 1, ex: "miRNA" or "mRNA"
                    'seg2_type',  # str, assigned "type" of segment 2, ex: "miRNA" or "mRNA"
                    'miRNA_seg',  # int, indicates which (if any) mapping is a miRNA
                                  #   options are "N" (none), "3p" (seg1), or "5p" (seg2),
                                  #   or "B" (both)
                    'extended',   # "0" or "1", boolean representation of whether
                                  #   sequences were artificailly extended as is
                                  #   performed by the Hyb software package.
                    'vienna',     # Vienna/dot-bracket format representation of binding
                   ]

    DEFAULT_FLAGS = HYB_FLAGS

    # Class-level variables
    custom_flags = []
    all_flags = HYB_FLAGS + HYBKIT_FLAGS
    find_type_method = None
    #find_type_method = find_seg_type_hyb
    find_type_params = {}
    reorder_flags = True           # Reorder flags to default when outputting
    fill_flags = False             # Fill flags-dict with default flags
    allow_undefined_flags = False  # Allow undefined flags

    # Placeholder symbol for empty entries. Default is "." in the Hyb software package.
    placeholder = '.'

    def __init__(self, hyb_id, seq, energy=placeholder, 
                 seg1_info={}, seg2_info={}, flags={}, 
                 find_seg_types=False):
        self.id = hyb_id
        self.seq = seq
        self.energy = energy

        self.seg1_info = self._make_seg_info_dict(seg1_info)
        self.seg2_info = self._make_seg_info_dict(seg2_info)
        self.flags = self._make_flags_dict(flags)

        if find_seg_types:
            self.find_seg_types()

    def seg1_id(self):
        'Return a copy of the id for segment 1 (5p).'
        return self.seg1_info['ref']

    def seg1_info(self):
        'Return a copy of the info dict object for segment 1 (5p).'
        return self.seg1_info.copy()

    def seg2_id(self):
        'Return a copy of the id for segment 2 (3p).'
        return self.seg2_info['ref']

    def seg2_info(self):
        'Return a copy of the info dict object for segment 2 (3p).'
        return self.seg2_info.copy()

    def seg_ids(self):
        'Return a tuple of the ids of segment 1 (5p) and segment 2 (3p).'
        return (self.seg1_info['ref'], self.seg2_info['ref'])

    def initialize_hyb_flags(self):
        'Add Hyb-default flags with initial options to entry flags.'
        hyb_flags = self.default_hyb_flags()
        for flag in hyb_flags:
            if flag not in self.flags:
                self.flags[flag] = hyb_flags[flag]

    def initialize_hybkit_flags(self):
        'Add hybkit-default flags with initial options to entry flags.'
        hybkit_flags = self.default_hybkit_flags()
        for flag in hybkit_flags:
            if flag not in self.flags:
                self.flags[flag] = hybkit_flags[flag]

    def seg1_type(self):
        'If the "seg1_type" flag is defined, return it. Otherwise return "None"'
        return self._get_flag_or_none('seg1_type')

    def seg2_type(self):
        'If the "seg2_type" flag is defined, return it. Otherwise return "None"'
        return self._get_flag_or_none('seg2_type')

    def seg_types(self):
        '''
        If the "seg1_type" and "seg2_type" flags are defined, return a tuple with both values.
        Otherwise return "None"
        '''
        return self._get_flag_or_none('seg1_type')

    def set_seg1_type(self, seg1_type):
        'Set "seg1_type" flag in flags.'
        self._set_flag('seg1_type', seg1_type)

    def set_seg2_type(self, seg2_type):
        'Set "seg2_type" flag in flags.'
        self._set_flag('seg2_type', seg2_type)

    def set_seg_types(self, seg1_type, seg2_type):
        'Set "seg1_type" and "seg2_type" flags in flags.'
        self._set_flag('seg1_type', seg1_type)
        self._set_flag('seg2_type', seg2_type)

    def find_seg_types(self, allow_unknown=False):
        '''
        Find the types of each segment using the method currently set for the class. 
        The default supplied method is HybRecord.find_seg_type_hyb, and works with alignemnt 
        mapping identifiers in the format of the reference database provided by the Hyb 
        Software Package. Custom methods with more complex behavior can be supplied via the 
        "set_find_method" function. 
        If allow_unknown is False, an error will be raised if a segment type cannot be identified.
        If allow_unknown is True, unidentified segments will be designated as "unknown". 
        ''' 
        types = []
        for seg_info in [self.seg1_info, self.seg2_info]:
            seg_type = self.find_type_method(seg_info, self.find_type_params)
            if seg_type is None:
                if allow_unknown:
                    types.append('unknown')
                else:
                    message = 'Cannot identify segment type for segment:\n'
                    message += self._format_seg_info(seg_info, prefix=' '*2) + '\n'
                    print(message)
                    raise Exception(message)
            else:
                types.append(seg_type)
        self.set_seg_types(types[0], types[1])

    def find_seg_type_hyb(self, seg_info, find_type_params={}):
        '''
        Return the type of the provided segment, or return None if the segment cannot be
        identified. This method works with sequence / alignment mapping identifiers
        in the format of the reference database provided by the Hyb Software Package,
        specifically identifiers of the format: "AAAA_BBBB_CCCC_DDDD" This function returns
        the fourth component of the identifier, split by "_", as the identfied sequence.
        For Example, "MIMAT0000076_MirBase_miR-21_microRNA" is identified as "microRNA".
        This function does not utilize the find_type_params arg.
        '''
        split_id = seg_info['ref'].split('_')
        if len(split_id) != 4:
            return None
        elif not bool(split_id[3]):
            return None
        else:
            return split_id[-1]

    #find_type_method = find_seg_type_hyb

    def to_line(self, newline=False):
        'Return a Hyb-format string representation of the Hyb record.'
        line_items = []
        for item_key in self.HYBRID_COLUMNS:
            line_items.append(getattr(self, item_key, '.'))
        for seg_dict in [self.seg1_info, self.seg2_info]:
            for item_key in self.SEGMENT_COLUMNS:
                if item_key in seg_dict and seg_dict[item_key] is not None:
                    line_items.append(seg_dict[item_key])
                else:
                    line_items.append('.')

        flag_string = self._make_flag_string() 

        if flag_string:
            line_items.append(flag_string)

        ret_string = '\t'.join((str(x) for x in line_items))
        if newline:
            ret_string += '\n'
        return ret_string

    def _format_seg_info(self, seg_info, prefix='', suffix='', indent_str=''):
        #Returns a formatted string of the sgement info information
        ret_string = prefix
        ret_string += indent_str + 'Map Reference:  %s\n' % seg_info['ref']
        ret_string += indent_str + 'Read Start Pos: %s\n' % seg_info['read_start']
        ret_string += indent_str + 'Read End Pos:   %s\n' % seg_info['read_end']
        ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_info['ref_start']
        ret_string += indent_str + 'Map End Pos:    %s\n' % seg_info['ref_end']
        ret_string += indent_str + 'Map Score:      %s\n' % seg_info['score']
        ret_string += suffix
        return ret_string

    def _get_flag_or_none(self, flag_key):
        if flag_key in self.flags:
            return self.flags[flag_key]
        else:
            return None
         
    def _set_flag(self, flag_key, flag_val, allow_undefined_flags=None):
        # Set the value of self.flags: flag_key to value flag_val.
        # allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        # If argument is provided to the method, it overrides default behavior.
        # Otherwise, the method falls back to the object-defaults.
        if allow_undefined_flags is None:
           allow_undefined_flags = self.allow_undefined_flags

        if not allow_undefined_flags and flag_key not in self.all_flags:
            message = 'Flag "%s" is not defined. Please check flag key' % flag_key
            message += ' or run with: "allow_undefined_flags=True"'
            print(message)
            raise Exception(message)

        self.flags[flag_key] = flag_val       


    def _make_flag_string(self, use_all_flags=False):
        flag_string = ''
        for flag in self._get_flag_keys():
           flag_string += ('%s=%s;' % (flag, str(self.flags[flag])))
        return flag_string

    def _get_flag_keys(self, reorder_flags=None):
        # reorder_flags option returns flags in deafult ordering scheme.
        #  If reorder-flags argument provided, it overrides default behavior.
        #  Otherwise, the method falls back to the object-default.
        return_list = []
        if reorder_flags is None:
            reorder_flags = self.reorder_flags
        if reorder_flags:
            return_list = self._get_ordered_flag_keys()
        else:
            return_list = self.flags.keys()
        return return_list

    def _get_ordered_flag_keys(self):
        return_list = []
        for flag in self.all_flags:
            if flag in self.flags:
                return_list.append(flag)
        for flag in self.flags:
            if flag not in self.all_flags:
                return_list.append(flag)
        return return_list

    def _make_seg_info_dict(self, seg_info_obj={}):
        # Create a dictionary with mapping entries, ensuring each read data point is either a
        #   placeholder or is of the correct data type.
        return_dict = {}
        segment_column_types = {'ref': str,
                                'read_start': int,
                                'read_end': int,
                                'ref_start': int,
                                'ref_end': int,
                                'score': str,
                                }
        for column in self.SEGMENT_COLUMNS:
            if column in seg_info_obj:
                if seg_info_obj[column] in ['.']:
                    return_dict[column] = None
                else:
                    column_type = segment_column_types[column]
                    try:
                        return_dict[column] = column_type(seg_info_obj[column])
                    except TypeError:
                        message = 'Entry "%s" for column: %s' % (seg_info_obj[column], column)
                        message += 'could not be converted to type: %s' % str(column_type)
                        print(message)
                        raise
            else:
                return_dict[column] = None
        return return_dict

    def _make_flags_dict(self, flag_obj={}, fill_flags=None, allow_undefined_flags=None):
        #  fill_flags_dict option populates the entry with the default flags and corresponding
        #   values where appropriate.
        #  allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        #  If either argument is provided to the method, it overrides default behavior.
        #  Otherwise, the method falls back to the object-defaults.
        if fill_flags is None:
            fill_flags = self.fill_flags
        if allow_undefined_flags is None:
            allow_undefined_flags = self.allow_undefined_flags

        if not isinstance(flag_obj, dict):
            message = '"flag_obj" argument must be a dict obj. Defined keys are:'
            message += self.all_flags.join(', ')
            print(message)
            raise Exception(message)

        if fill_flags:
            message = 'This feature has not yet been implemented, to ensure that false data is not'
            message += ' unintentionally added to Hyb entries.'
            print(message)
            raise NotImplementedError(message)

        if not allow_undefined_flags:
            for flag in flag_obj:
                if flag not in self.all_flags:
                    message = 'Flag "%s" is not defined. Please check flag key' % flag
                    message += ' or run with: "allow_undefined_flags=True"'
                    print(message)
                    raise Exception(message)
        return flag_obj

    # Methods for hashing and comparison. The object ".id" value is presumed to be unique, and
    # can therefore be used for comparison and hashing.
    def __eq__(self, other):
        return self.id == other.id
    def __neq__(self, other):
        return self.id != other.id
    def __hash__(self):
        return hash(self.id)
    def __bool__(self):
        return True
    # Container length properties will refer to the length of the hybrid sequence.
    def __len__(self):
        return len(self.seq)

    @classmethod
    def set_find_type_method(cls, find_method, find_params={}):
        '''
        Set the class-level custom method for segment assignemnt to callable function 
        in find_method, that has the form: "def my_method(self, seg_info, find_params)".
        This method should return the string of the assigned segment type if found, or a
        None object if the type cannot be found.
        It can also take a dictionary in the "find_params" argument that specifies 
        additional or dynamic search properties, as desired.
        ''' 
        cls.find_type_method = types.MethodType(find_method, cls)
        cls.find_type_params = find_params

    @classmethod
    def set_custom_flags(cls, custom_flags):
        '''
        Set the class-level HybRecord.custom_flags variable, (and update the 
        HybRecord.all_flags variable) to allow custom flags in your Hyb file without
        causing an exception.
        '''
        cls.custom_flags = custom_flags
        cls.all_flags = HYB_FLAGS + HYBKIT_FLAGS + cls.custom_flags

    @classmethod
    def _read_flags(cls, flag_string, allow_undefined_flags=True):
        # allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        # undefined flags allowed in this method by default, to allow the object-level setting to
        # take precedence
        flag_string = flag_string.rstrip()
        flag_string = flag_string.rstrip(';')
        flag_pairs = [flag_pair.split('=') for flag_pair in flag_string.split(';')]
        flags = {}
        for flag_key, flag_value in flag_pairs:
            if not allow_undefined_flags and flag_key not in cls.all_flags:
                message = 'Problem: Unidefined Flag: %s' % flag_key
                print(message)
                raise Exception(message)
            flags[flag_key] = flag_value
        return flags

    @classmethod
    def default_hyb_flags(cls):
        flags = {}
        flags['count_total'] = 0
        flags['count_last_clustering'] = 0
        flags['two_way_merged'] = 0
        flags['seq_IDs_in_cluster'] = ''
        return flags

    @classmethod
    def default_hybkit_flags(cls):
        flags = {
                 'orient': 'U',
                 'seg1_type': 'unk',
                 'seg2_type': 'unk',
                 'miRNA_seg': 'U',
                 'extended': 0,
                }
        return flags

    @classmethod
    def from_line(cls, line):
        line_items = line.strip().split('\t')
        # print(line_items)
        hyb_id = line_items[0]
        seq = line_items[1]
        energy = line_items[2]
        seg1_info = {}
        seg1_info['ref'] = line_items[3]
        seg1_info['read_start'] = line_items[4]
        seg1_info['read_end'] = line_items[5]
        seg1_info['ref_start'] = line_items[6]
        seg1_info['ref_end'] = line_items[7]
        seg1_info['score'] = line_items[8]
        seg2_info = {}
        seg2_info['ref'] = line_items[9]
        seg2_info['read_start'] = line_items[10]
        seg2_info['read_end'] = line_items[11]
        seg2_info['ref_start'] = line_items[12]
        seg2_info['ref_end'] = line_items[13]
        seg2_info['score'] = line_items[14]
        flags = {}
        if len(line_items) > 15:
            flags = cls._read_flags(line_items[15])

        return_obj = cls(hyb_id, seq, energy, seg1_info, seg2_info, flags)
        return return_obj


class HybFile(io.FileIO):
    '''
    File-object subclass that provides abiltity to return file lines as HybRecord entries.
    io.FileIO (apparently) only supports reading and writing bytestrings, so any utilized
    base class methods (like .read() and .write() utilize bytestrings, not unicode strings.
    '''
    
    def __next__(self): 
        'Call io.FileIO __next__ method and return output as HybRecord object.'
        return HybRecord.from_line(str(super().__next__(), 'utf-8'))

    def read_record(self):
        'Return next line of hyb file as HybRecord object.'
        return next(self)

    def read_records(self):
        'Return list of all records in hyb file as HybRecord objects.'
        records = []
        for record in self:
            records.append(record)
        return records

    def write_record(self, write_record):
        '''
        Write a HybRecord object to file as a Hyb-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        '''
        self._ensure_HybRecord(write_record)
        record_string = write_record.to_line(newline=True)
        record_bytestring = bytearray(record_string, 'utf-8')
        self.write(record_bytestring)

    def write_records(self, write_records):
        '''
        Write a sequence of HybRecord objects as hyb-format lines to the Hyb file.
        Unlike the file.writelines() method, this method will add a newline to the 
        end of each written record line.
        '''
        for write_record in write_records:
            self._ensure_HybRecord(write_record)
            self.write(write_record.to_line(newline=True))
       
    def _ensure_HybRecord(self, record):
        if not isinstance(record, HybRecord):
            message = 'Item: "%s" is not a HybRecord object.' % record
            print(message)
            raise Exception(message) 
        
