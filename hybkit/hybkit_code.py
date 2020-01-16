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
    reorder_flags = True                # Reorder flags to default when outputting
    fill_flags = False                  # Fill flags-dict with default flags
    allow_undefined_flags = False       # Allow undefined flags

    #Ideally the following paramaters should be set to False, and True, respctively, however the 
    #  output of the Hyb program often provides vienna/viennad fold-records that do not 
    #  match the sequences in the corresponding hyb entry, as an artifact of collapsing 
    #  multiple hyb records into the same entry.
    allow_fold_record_mismatch = True   # Allow mismatch with self.fold_record sequence
    warn_fold_record_mismatch = False   # Warn if mismatch with self.fold_record sequence

    # Placeholder symbol for empty entries. Default is "." in the Hyb software package.
    placeholder = '.'

    def __init__(self, hyb_id, seq, energy=placeholder,
                 seg1_info={}, seg2_info={}, flags={},
                 find_seg_types=False,
                 fold_record=None):
        self.id = hyb_id
        self.seq = seq
        self.energy = energy

        self.seg1_info = self._make_seg_info_dict(seg1_info)
        self.seg2_info = self._make_seg_info_dict(seg2_info)
        self.flags = self._make_flags_dict(flags)

        if find_seg_types:
            self.find_seg_types()

        if fold_record is not None:
            self.set_fold_record(fold_record)
        else:
            self.fold_record = None

    def seg1_id(self):
        'Return a copy of the id for segment 1 (5p), or None if not defined.'
        if 'ref' in self.seg1_info:
            return self.seg1_info['ref']
        else:
            return None

    def seg2_id(self):
        'Return a copy of the id for segment 2 (3p), or None if not defined.'
        if 'ref' in self.seg2_info:
            return self.seg2_info['ref']
        else:
            return None

    def seg1_info(self):
        'Return a copy of the info dict object for segment 1 (5p).'
        return self.seg1_info.copy()

    def seg2_info(self):
        'Return a copy of the info dict object for segment 2 (3p).'
        return self.seg2_info.copy()

    def seg_ids(self):
        'Return a tuple of the ids of segment 1 (5p) and segment 2 (3p), or None if not defined.'
        return (self.seg1_id(), self.seg2_id())

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

    def check_fold_record_match(self, fold_record):
        '''
        Return True if the sequence (".seq") attribute of a FoldRecord instance matches the
        sequence (".seq") attribute of this instance.
        '''
        return (self.seq == fold_record.seq)

    def set_fold_record(self, fold_record, 
                        allow_fold_record_mismatch=None,
                        warn_fold_record_mismatch=None):
        '''
        Check to ensure that fold_record argument is an instance of FoldRecord, and that 
        it has a matching sequence to this HybRecord, then set it as self.fold_record.
        allow_fold_record_mismatch allows mismatches between the HybRecord sequence 
        and the FoldRecord sequence.
        warn_fold_record_mismatch prints a warning when there is a mismatch between the 
        HybRecord sequence and hte FoldRecord sequence.
        If either argument is provided to the method, it overrides default behavior.
        Otherwise, the method falls back to the class default setting.
        '''
        if allow_fold_record_mismatch is None:
            allow_fold_record_mismatch = self.allow_fold_record_mismatch
        if warn_fold_record_mismatch is None:
            warn_fold_record_mismatch = self.warn_fold_record_mismatch

        if not isinstance(fold_record, FoldRecord):
            message = 'Supplied argument to fold_record: %s' % str(fold_record)
            message += '\n   is not a FoldRecord object.'
            print(message)
            raise Exception(message)
        elif not self.check_fold_record_match(fold_record):
            if not allow_fold_record_mismatch:
                message = 'For HybRecord: %s\n' % str(self)
                message += 'Supplied FoldRecord %s\n' % str(fold_record)
                message += '   does not have a matching sequence.\n'
                message += '  %s\n  %s' % (self.seq, fold_record.seq)
                print(message)
                raise Exception(message)

            elif warn_fold_record_mismatch:
                shortest = min(len(self), len(fold_record))
                seq1 = self.seq[:]
                seq2 = fold_record.seq[:]
                while (shortest > 0) and seq1[0] == seq2[0]:
                    seq1 = seq1[1:]
                    seq2 = seq2[1:]
                    shortest = min(len(seq1), len(seq2))
                while (shortest > 0) and seq1[-1] == seq2[-1]:
                    seq1 = seq1[1:]
                    seq2 = seq2[1:]
                    shortest = min(len(seq1), len(seq2))
                num_mismatch = max(len(seq1), len(seq2))
                message = 'WARNING: Hybrid: %s' % self.id
                message += ' has fold_record sequence with a %i base mismatch.' % num_mismatch
                #message += '\n%s\n%s' % (self.seq, fold_record.seq)
                #message += '\nMismatch Region: %s | %s' % (seq1, seq2) 
                print(message)
 
        self.fold_record = fold_record

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

    # Magic methods for hashing and comparison. The object ".id" value and ".seq" value
    # is presumed to be unique, and can therefore be used for comparison and hashing.
    # These comparisons can also be used with FoldRecord objects.
    def __eq__(self, other):
        return (self.id == other.id and self.seq == other.seq)
    def __neq__(self, other):
        return (self.id != other.id or self.seq != other.seq)
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
            self.write_record(write_record)

    def _ensure_HybRecord(self, record):
        if not isinstance(record, HybRecord):
            message = 'Item: "%s" is not a HybRecord object.' % record
            print(message)
            raise Exception(message)



class FoldRecord(object):
    '''
    Class for storing secondary structure (folding) information for a nucleotide sequence.
    This class supports:
        The Vienna file format: http://unafold.rna.albany.edu/doc/formats.php#VIENNA
        ex: 34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-...
            TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
            .....((((((.((((((......)))))).))))))   (-11.1)

        (TODO)
        The Viennad file format utilizied in the Hyb Software package
        ex:
            34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-34-...
            TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
            TAGCTTATCAGACTGATGT------------------   miR-21_microRNA 1       19
            -------------------TAGCTTATCAGACTGATG   miR-21_microRNA 1       18
            .....((((((.((((((......)))))).))))))   (-11.1)
            [space-line]

    A minimum amount of data necessary for a FoldRecord object is a sequence identifier,
    a genomic sequence, and its fold representaiton.
    '''

    # Class-level constants

    # Placeholder symbol for empty entries. Default is "." in the Hyb software package.
    placeholder = '.'

    # Class-level variables
    def __init__(self, id, seq, fold, energy,
                 seg1_fold_info={},
                 seg2_fold_info={}):
        self.id = id
        self.seq = seq
        self.fold = fold
        self.energy = energy

        self.set_seg1_fold_info(seg1_fold_info)
        self.set_seg2_fold_info(seg2_fold_info)

    def seg1_id(self):
        'Return a copy of the id for segment 1 (5p), or None if not defined.'
        if 'ref' in self.seg1_fold_info:
            return self.seg1_fold_info['ref']
        else:
            return None

    def seg2_id(self):
        'Return a copy of the id for segment 2 (3p), or None if not defined.'
        if 'ref' in self.seg2_fold_info:
            return self.seg2_fold_info['ref']
        else:
            return None

    def seg_ids(self):
        'Return a tuple of the ids of segment 1 (5p) and segment 2 (3p), or None if not defined.'
        return (self.seg1_id(), self.seg2_id())

    def seg1_info(self):
        'Return a copy of the info dict object for segment 1 (5p).'
        return self.seg1_fold_info.copy()

    def seg2_info(self):
        'Return a copy of the info dict object for segment 2 (3p).'
        return self.seg2_fold_info.copy()

    def seg1_detail(self, detail):
        '''
        Return a copy of the detail for seg1 provided in by the key in "detail" parameter,
        or if it does not exist return None.
        '''
        return self._get_segN_detail(1, detail)

    def seg2_detail(self, detail):
        '''
        Return a copy of the detail for seg2 provided in by the key in "detail" parameter,
        or if it does not exist return None.
        '''
        return self._get_segN_detail(2, detail)

    def to_vienna_lines(self, newline=False):
        'Return a list of lines for the record in vienna format.'
        ret_lines = []
        suffix = ''
        if newline:
            suffix = '\n'
        ret_lines.append(self.id + suffix)   # Add line 1, id
        ret_lines.append(self.seq + suffix)  # Add line 2, sequence

        #Create formatted energy string which uses no decimal places for integer numbers
        if abs(self.energy - round(self.energy)) > 0.00001:
            energy_str = ("%.5f" % self.energy).rstrip('0')
        else:
            energy_str = "%i" % int(round(self.energy))

        line_3 = '%s\t(%s)' % (self.fold, energy_str)
        ret_lines.append(line_3 + suffix)    # Add line 3, fold representaiton and energy
        return ret_lines

    def to_vienna_string(self, newline=False):
        'return a 3-line string for the record in vienna format.'
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_vienna_lines(newline=False)) + suffix)

    def to_viennad_lines(self, newline=False):
        'Return a list of lines for the record in viennad format.'
        ret_lines = []
        if not newline:
            suffix = ''
        elif newline:
            suffix = '\n'
        ret_lines.append(self.id + suffix)   # Add line 1, id
        ret_lines.append(self.seq + suffix)  # Add line 2, sequence
        line_3_items = [self.seg1_detail('highlight'),
                        self.seg1_detail('ref'),
                        self.seg1_detail('start'),
                        self.seg1_detail('end')]
        line_3 = '\t'.join([str(item) if item is not None
                            else self.placeholder
                            for item in line_3_items ])
        ret_lines.append(line_3 + suffix)
        line_4_items = [self.seg2_detail('highlight'),
                        self.seg2_detail('ref'),
                        self.seg2_detail('start'),
                        self.seg2_detail('end')]
        line_4 = '\t'.join([str(item) if item is not None
                            else self.placeholder
                            for item in line_4_items])
        ret_lines.append(line_4 + suffix)

        #Create formatted energy string which uses no decimal places for integer numbers
        if abs(self.energy - round(self.energy)) > 0.00001:
            energy_str = ("%.5f" % self.energy).rstrip('0')
        else:
            energy_str = "%i" % int(round(self.energy))

        line_5 = '%s\t(%s)' % (self.fold, energy_str)
        ret_lines.append(line_5 + suffix)  # Add line 5, fold representaiton and energy
        ret_lines.append(suffix)           # Add line 6, blank/newline  
        return ret_lines

    def to_viennad_string(self, newline=False):
        'return a 6-line string for the record in viennad format.'
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_viennad_lines(newline=False)) + suffix)

    def set_seg1_fold_info(self, seg_info_obj):
        'Set folding information for segment 1'
        self._set_segN_fold_info(1, seg_info_obj)

    def set_seg2_fold_info(self, seg_info_obj):
        'Set folding information for segment 2'
        self._set_segN_fold_info(2, seg_info_obj)

    def set_segs_fold_info_from_hybrecord(self, hybrecord):
        raise NotImplementedError

    def find_seg1_fold_info(self, seg1_name, seg1_start_n, seg1_end_n):
        raise NotImplementedError

    def find_seg2_fold_info(self, seg2_name, seg2_start_n, seg2_end_n):
        raise NotImplementedError

    def check_hyb_record_match(self, hyb_record):
        '''
        Return True if the sequence (".seq") attribute of a HybRecord instance matches the
        sequence (".seq") attribute of this instance.
        '''
        return (self.seq == hyb_record.seq)

    def _find_segN_fold_details(self, seg_start_n, seg_end_n):
        raise NotImplementedError

    def _set_segN_fold_info(self, seg, seg_info_obj={}):
        # Create a dictionary with segment fold information.
        return_dict = {}
        segment_param_types = {'ref': str,
                               'start': int,
                               'end': int,
                               'highlight': str,
                               'seg_fold': str ,
                               }
        for param in segment_param_types.keys():
            if param in seg_info_obj:
                param_type = segment_param_types[param]
                try:
                    return_dict[param] = param_type(seg_info_obj[param])
                except TypeError:
                    message = 'Entry "%s" ' % seg_info_obj[param]
                    message += 'for segment parameter: %s' % param
                    message += 'could not be converted to type: %s' % str(param_type)
                    print(message)
                    raise
            else:
                return_dict[param] = None
        if seg == 1:
            self.seg1_fold_info = return_dict
        elif seg == 2:
            self.seg2_fold_info = return_dict
        else:
            raise Exception('Unidentified value "%s" for seg argument' % str(seg))

    def _get_segN_detail(self, seg, detail):
        if seg == 1:
            seg_dict = self.seg1_fold_info
        elif seg == 2:
            seg_dict = self.seg2_fold_info
        else:
            raise Exception('Unidentified value "%s" for seg argument' % str(seg))
        if detail in seg_dict:
            return seg_dict[detail]
        else:
            return None

    def _format_seg_info(self, seg_info, prefix='', suffix='', indent_str=''):
        raise NotImplementedError
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

    # Magic methods for hashing and comparison. The object ".id" value and ".seq" value
    # is presumed to be unique, and can therefore be used for comparison and hashing.
    # These comparisons can also be used with HybRecord objects.
    def __eq__(self, other):
        return (self.id == other.id and self.seq == other.seq)
    def __neq__(self, other):
        return (self.id != other.id or self.seq != other.seq)
    def __hash__(self):
        return hash(self.id)
    def __bool__(self):
        return True
    # Container length properties will refer to the length of the hybrid sequence.
    def __len__(self):
        return len(self.seq)

    @classmethod
    def from_vienna_lines(cls, record_lines):
        '''
        Create a FoldRecord entry from a list of 3 strings corresponding to lines in the
        Vienna format.
        '''

        if len(record_lines) != 3:
            message = 'Provided Vienna Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... are not in required 3-line format.'
            print(message)
            raise Exception(message)

        rec_id = record_lines[0].strip()
        seq = record_lines[1].strip()
        line_3 = record_lines[2].strip()
        line_3_split = line_3.split('\t')
        if len(line_3_split) != 2:
            message = 'Provided Vienna Record Line 3:\n'
            message += line_3.rstrip() + '\n'
            message += str(line_3_split) + '\n'
            message += '\n  ... does not have required ".(.).<tab>(-1.23)" format.'
            print(message)
            raise Exception(message)
        fold = line_3_split[0]
        energy = float(line_3_split[1].strip('()'))

        return_obj = cls(rec_id, seq, fold, energy)
        return return_obj

    @classmethod
    def from_vienna_string(cls, record_lines):
        '''
        Create a FoldRecord entry from a string containing 3 lines corresponding to lines in the
        Vienna format.
        '''
        lines = record_lines.strip().split('\n')
        return cls.from_vienna_lines(lines)

    @classmethod
    def from_viennad_lines(cls, record_lines):
        '''
        Create a FoldRecord entry from a list of 5 or 6 strings corresponding to lines in the
        Viennad format.
        '''

        if len(record_lines) not in [5, 6]:
            message = 'Provided Viennad Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... are not in required 5-line or 6-line format.'
            print(message)
            raise Exception(message)

        rec_id = record_lines[0].strip()
        seq = record_lines[1].strip()
        line_3 = record_lines[2].strip()
        line_3_split = line_3.split('\t')
        if len(line_3_split) != 4:
            message = 'Provided Vienna Record Line 3:\n'
            message += line_3.rstrip() + '\n'
            message += str(line_3_split) + '\n'
            message += '\n  ... does not have required "ACTG---<tab>'
            message += 'name<tab>start<tab>end" format.'
            print(message)
            raise Exception(message)
        seg1_fold_info = {
                          'highlight' : line_3_split[0],
                          'ref' : line_3_split[1],
                          'start' : line_3_split[2],
                          'end' : line_3_split[3],
                         }
        line_4 = record_lines[3].strip()
        line_4_split = line_4.split('\t')
        if len(line_4_split) != 4:
            message = 'Provided Vienna Record Line 4:\n'
            message += line_4.rstrip() + '\n'
            message += str(line_4_split) + '\n'
            message += '\n  ... does not have required "ACTG---<tab>'
            message += 'name<tab>start<tab>end" format.'
            print(message)
            raise Exception(message)
        seg2_fold_info = {
                          'highlight' : line_4_split[0],
                          'ref' : line_4_split[1],
                          'start' : line_4_split[2],
                          'end' : line_4_split[3],
                         }

        line_5 = record_lines[4].strip()
        line_5_split = line_5.split('\t')
        fold = line_5_split[0]
        energy = float(line_5_split[1].strip('()'))
        return_obj = cls(rec_id, seq, fold, energy,
                         seg1_fold_info, seg2_fold_info)
        return return_obj

    @classmethod
    def from_viennad_string(cls, record_lines):
        '''
        Create a FoldRecord entry from a string containing 5 or 6 lines corresponding 
        to lines in the Viennad format.
        '''
        lines = record_lines.strip().split('\n')
        return cls.from_vienna_lines(lines)


class ViennaFile(io.FileIO):
    '''
    File-object subclass that provides abiltity to return sets of three file lines as
    FoldRecord entries.
    io.FileIO (apparently) only supports reading and writing bytestrings, so any utilized
    base class methods (like .read() and .write() utilize bytestrings, not unicode strings.
    '''

    def __next__(self):
        '''
        Call io.FileIO __next__ method for next three lines and return
        output as HybRecord object.'
        '''
        line_1 = str(super().__next__(), 'utf-8')
        line_2 = str(super().__next__(), 'utf-8')
        line_3 = str(super().__next__(), 'utf-8')
        return FoldRecord.from_vienna_lines((line_1, line_2, line_3))

    def read_record(self):
        'Return next three line of vienna file as FoldRecord object.'
        return next(self)

    def read_records(self):
        'Return list of vienna records in vienna file as FoldRecord objects.'
        records = []
        for record in self:
            records.append(record)
        return records

    def write_record(self, write_record):
        '''
        Write a FoldRecord object to file as a vienna-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        '''
        self._ensure_FoldRecord(write_record)
        record_string = write_record.to_vienna_string(newline=True)
        record_bytestring = bytearray(record_string, 'utf-8')
        self.write(record_bytestring)

    def write_records(self, write_records):
        '''
        Write a sequence of FoldRecord objects as vienna-format lines to the vienna file.
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        '''
        for write_record in write_records:
            self.write_record(write_record)

    def _ensure_FoldRecord(self, record):
        if not isinstance(record, FoldRecord):
            message = 'Item: "%s" is not a FoldRecord object.' % record
            print(message)
            raise Exception(message)

class HybViennaIter(object):
    '''
    This class provides an iterator to iterate through a HybFile and ViennaFile simultaneously,
    returning a tuple of hyb_record, fold_record instances on each iteration.
    '''
    def __init__(self, hybfile_handle, viennafile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennafile_handle = viennafile_handle
        self.counter = 0

    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_vienna_record = next(self.viennafile_handle)
        return (next_hyb_record, next_vienna_record)


class HybViennaCmbIter(object):
    '''
    This class provides an iterator to iterate through a HybFile and ViennaFile simultaneously.
    It is presumed that each respective hyb entry corresponds to an aligned vienna entry.
    Each ViennaRecord will be added to the corresponding HybRecord. 
    Only the HybRecord entry will then be returned, containing the associated ViennaRecord entry.
    '''
    def __init__(self, hybfile_handle, viennafile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennafile_handle = viennafile_handle
        self.counter = 0

    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_vienna_record = next(self.viennafile_handle)
        try:
            next_hyb_record.set_fold_record(next_vienna_record)
        except:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return next_hyb_record

class ViennadFile(io.FileIO):
    '''
    File-object subclass that provides abiltity to return sets of six viennad file lines as
    FoldRecord entries.
    io.FileIO (apparently) only supports reading and writing bytestrings, so any utilized
    base class methods (like .read() and .write() utilize bytestrings, not unicode strings.
    '''

    def __next__(self):
        '''
        Call io.FileIO __next__ method for next three six lines and return
        output as HybRecord object.'
        '''
        line_1 = str(super().__next__(), 'utf-8')
        line_2 = str(super().__next__(), 'utf-8')
        line_3 = str(super().__next__(), 'utf-8')
        line_4 = str(super().__next__(), 'utf-8')
        line_5 = str(super().__next__(), 'utf-8')
        line_6 = str(super().__next__(), 'utf-8')
        return FoldRecord.from_viennad_lines((line_1, line_2, line_3, line_4, line_5, line_6))

    def read_record(self):
        'Return next six lines of viennad file as FoldRecord object.'
        return next(self)

    def read_records(self):
        'Return list of viennad records in viennad file as FoldRecord objects.'
        records = []
        for record in self:
            records.append(record)
        return records

    def write_record(self, write_record):
        '''
        Write a FoldRecord object to file as a viennad-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        '''
        self._ensure_FoldRecord(write_record)
        record_string = write_record.to_viennad_string(newline=True)
        record_bytestring = bytearray(record_string, 'utf-8')
        self.write(record_bytestring)

    def write_records(self, write_records):
        '''
        Write a sequence of FoldRecord objects as viennad-format lines to the viennad file.
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        '''
        for write_record in write_records:
            self.write_record(write_record)

    def _ensure_FoldRecord(self, record):
        if not isinstance(record, FoldRecord):
            message = 'Item: "%s" is not a FoldRecord object.' % record
            print(message)
            raise Exception(message)

class HybViennadIter(object):
    '''
    This class provides an iterator to iterate through a HybFile and ViennadFile simultaneously,
    returning a tuple of hyb_record, fold_record instances on each iteration.
    '''
    def __init__(self, hybfile_handle, viennadfile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennadfile_handle = viennadfile_handle
        self.counter = 0

    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_viennad_record = next(self.viennadfile_handle)
        return (next_hyb_record, next_viennad_record)


class HybViennadCmbIter(object):
    '''
    This class provides an iterator to iterate through a HybFile and ViennadFile simultaneously.
    It is presumed that each respective hyb entry corresponds to an aligned viennad entry.
    Each FoldRecord will be added to the corresponding HybRecord. 
    Only the HybRecord entry will then be returned, containing the associated FoldRecord entry.
    '''
    def __init__(self, hybfile_handle, viennadfile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennadfile_handle = viennadfile_handle
        self.counter = 0

    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_viennad_record = next(self.viennadfile_handle)
        try:
            next_hyb_record.set_fold_record(next_viennad_record)
        except:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return next_hyb_record
