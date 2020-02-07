#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Classes and Methods for manipulating data in the ".hyb" genomic sequence format.

Todo:
    Expand Definition of class.
    Check hybrecord/foldrecord match settings.
    Add Hybrecord.to_csv_line()
    Add Hybrecord.to_csv_header()
    Clean from_viennd_lines() method.
    Add analysis.target_region analyses.
    Overhaul docstrings and make beautiful documentation.
    Create hyb-format database download script.
    Add user-friendly individual scripts.
    Implement all sample analyses with bash workflows and individual scripts.
    Implement all sample analyses with nextflow workflows and individual scripts.
    Make decision and clean "extra" scripts.
"""

import os
import io
import types
import csv
import copy
from collections import OrderedDict
import hybkit
import hybkit.analysis
import hybkit.plot

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

class HybRecord(object):
    """Class for storing and analayzing information about a chimeric or hybrid 
    sequence read from a .hyb format file.  

    A minimum amount of data necessary 
    for a HybRecord object is the genomic sequence and its
    corresponding identifier. The .hyb file format is specified in
        "Travis, Anthony J., et al. "Hyb: a bioinformatics pipeline for the analysis of CLASH 
        (crosslinking, ligation and sequencing of hybrids) data." 
        Methods 65.3 (2014): 263-273."
    and consists of 15-16 columns separated by tabs, briefly described as: 
    "id, sequence, binding_energy, seg1_name, seg1_read_start, seg1_read_end, 
     seg1_ref_start, seg1_ref_end, seg1_map_score, seg2_read_start, seg2_read_end,
     seg2_ref_start, seg2_ref_end, seg2_map_score, [flag1=val1; flag2=val2;flag3=val3...]"
    """

    # HybRecord : Class-Level Constants
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
                       'score',       # str(int or float) Alignment Score, can be BLAST e-score
                                      #   or mapping alignment score, depending on analysis
                                      #   implementation type.
                      ]
    # Arbitrary details included in column 16 of the hyb format in the form:
    #   “feature1=value1; feature2=value2;..."
    #   Flags utilized in the Hyb software package
    HYB_FLAGS = [
                 'count_total',            # str(int), total represented hybrids
                 'count_last_clustering',  # str(int), total represented hybrids at last clustring
                 'two_way_merged',         # "0" or "1", boolean representation of whether
                                           #   entries with mirrored 5' and 3' hybrids were merged
                 'seq_IDs_in_cluster',     # str, comma-separated list of all ids of hybrids
                                           #   merged into this hybrid entry.
                ]
    # Additional flag specifications utilized by hybkit
    HYBKIT_FLAGS = [
                    'read_count',   # str(int), number of sequence reads represented by record
                                    #   if merged record, represents total for all merged entries
                    'orient',       # str, orientation of strand. Options:
                                    #   "F" (Forward), "IF" (Inferred Forward),
                                    #   "R" (Reverse), "IR" (Inferred Reverse),
                                    #   "U" (Unknown), or "C" (Conflicting)
                    'seg1_type',    # str, assigned "type" of segment 1, ex: "miRNA" or "mRNA"
                    'seg2_type',    # str, assigned "type" of segment 2, ex: "miRNA" or "mRNA"
                    'seg1_det',     # str, arbitrary detail about segment 1
                    'seg2_det',     # str, arbitrary detail about segment 2
                    'miRNA_seg',    # str, indicates which (if any) mapping is a miRNA
                                    #   options are "N" (none), "3p" (seg1), "5p" (seg2),
                                    #   "B" (both), or "U" (unknown)
                    'target_reg',   # str, assigned region of the miRNA target.
                                    #   options are "5pUTR", "coding", "3pUTR",
                                    #   "N" (none), or "U" (unknown)
                    'ext',          # int, "0" or "1", boolean representation of whether
                                    #   record sequences were bioinformatically extended as is
                                    #   performed by the Hyb software package.
                    'source',       # str, label for sequence source, when combining records 
                                    #   from different sources.
                   ]

    # miRNA Types for use in ".mirna_analysis" function
    MIRNA_TYPES = {'miRNA', 'microRNA'}

    # Coding types for use in the target_region_analysis function
    CODING_TYPES = {'mRNA'}


    # HybRecord : Class-Level Default Settings:
    # These settings can be changed in a script to change the analysis-specific behavior
    #   of the class.
    DEFAULTS = {}
    DEFAULTS['reorder_flags'] = True           # Reorder flags to default order when writing
    DEFAULTS['allow_undefined_flags'] = False  # Allow undefined flags

    # Ideally the following paramaters should be set to False, and True, respectively. However
    #   the output of the Hyb program often provides vienna/viennad fold-records that do not
    #   match the sequences in the corresponding hyb entry, perhaps as an artifact of collapsing
    #   multiple hyb records into the same entry.
    DEFAULTS['allow_fold_record_mismatch'] = True   # Allow mismatch in self.fold_record sequence
    DEFAULTS['warn_fold_record_mismatch'] = False   # Warn if mismatch, self.fold_record sequence

    DEFAULTS['allow_unknown_regions'] = False
    DEFAULTS['warn_unknown_regions'] = True

    # Placeholder symbol for empty entries. Default is "." in the Hyb software package.
    DEFAULTS['placeholder'] = '.'

    # HybRecord : Class-Level Variables
    custom_flags = []
    all_flags = HYB_FLAGS + HYBKIT_FLAGS
    find_type_method = None             # Will be populated after method definition
    find_type_params = {}               # May be populated during use.
    target_region_info = {}             # May be populated during use.

    # set class initial settings to defaults.
    settings = copy.deepcopy(DEFAULTS)

    # HybRecord : Public Methods : Initialization
    def __init__(self, hyb_id, seq, energy=None,
                 seg1_info={}, seg2_info={}, flags={},
                 read_count=None,
                 find_seg_types=False,
                 fold_record=None):
        self.id = hyb_id
        self.seq = seq

        if energy is not None:
            self.energy = energy
        else:
            energy = self.settings['placeholder']

        self.seg1_info = self._make_seg_info_dict(seg1_info)
        self.seg2_info = self._make_seg_info_dict(seg2_info)
        self.flags = self._make_flags_dict(flags)

        if read_count is not None:
            if read_count not in flags:
                self.set_flag('read_count', read_count)

        if find_seg_types:
            self.find_seg_types()

        self.fold_seq_match = None  # Placeholder to be filled during set_fold_record()
        if fold_record is not None:
            self.set_fold_record(fold_record)
        else:
            self.fold_record = None

        self.mirna_details = None  # Placeholder variable for mirna_analysis
        self.mirna_info = None     # Placeholder variable for mirna_analysis
        self.target_info = None    # Placeholder variable for mirna_analysis

        self._post_init_tasks()    # Method stub to enable convenient subclassing

    # HybRecord : Public Methods : Segment_Info
    def seg1_id(self):
        """Return a copy of the id for segment 1 (5p), or None if not defined."""
        if 'ref' in self.seg1_info:
            return self.seg1_info['ref']
        else:
            return None

    # HybRecord : Public Methods : Segment_Info
    def seg2_id(self):
        """Return a copy of the id for segment 2 (3p), or None if not defined."""
        if 'ref' in self.seg2_info:
            return self.seg2_info['ref']
        else:
            return None

    # HybRecord : Public Methods : Segment_Info
    def seg_ids(self):
        """Return a tuple of the ids of segment 1 (5p) and segment 2 (3p), or None if not defined."""
        return (self.seg1_id(), self.seg2_id())

    # HybRecord : Public Methods : Segment_Info
    def seg1_info(self):
        """Return a copy of the info dict object for segment 1 (5p)."""
        return self.seg1_info.copy()

    # HybRecord : Public Methods : Segment_Info
    def seg2_info(self):
        """Return a copy of the info dict object for segment 2 (3p)."""
        return self.seg2_info.copy()

    # HybRecord : Public Methods : flags
    def set_flag(self, flag_key, flag_val, allow_undefined_flags=None):
        """Set the value of self.flags: flag_key to value flag_val.
        allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        If argument is provided to the method, it overrides default behavior.
        Otherwise, the method falls back to the object-defaults.
        """

        if allow_undefined_flags is None:
            allow_undefined_flags = self.settings['allow_undefined_flags']

        if not allow_undefined_flags and flag_key not in self.all_flags:
            message = 'Flag "%s" is not defined. Please check flag key' % flag_key
            message += ' or run with: "allow_undefined_flags=True"'
            print(message)
            raise Exception(message)

        self.flags[flag_key] = flag_val

    # HybRecord : Public Methods : Flag_Info : seg_type
    def seg1_type(self, require=False):
        """If the "seg1_type" flag is defined, return it. 
        If require is provided as True, raise an error if the seg1_type is not defined.
        Otherwise return "None"
        """
        if require:
            return self._get_flag('seg1_type')
        else:
            return self._get_flag_or_none('seg1_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def seg2_type(self, require=False):
        """If the "seg2_type" flag is defined, return it.
        If require is provided as True, raise an error if the seg2_type flag is not defined.
        Otherwise return "None"
        """
        if require:
            return self._get_flag('seg2_type')
        else:
            return self._get_flag_or_none('seg2_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def seg_types(self, require=False):
        """If the "seg1_type" and "seg2_type" flags are defined, return a tuple with both values.
        If require is provided as True, raise an error if either of the seg1_type or seg2_type 
        flags are not defined. Otherwise return "None" for each undefined flag.
        """
        if require:
            return (self._get_flag('seg1_type'), self._get_flag('seg2_type'))
        else:
            return (self._get_flag_or_none('seg1_type'), self._get_flag_or_none('seg2_type'))

    # HybRecord : Public Methods : Flag_Info : seg_type
    def seg_types_sorted(self, require=False):
        """If the "seg1_type" and "seg2_type" flags are defined, return a tuple of both values
        sorted alphabetically. 
        If require is provided as True, raise an error if either of the seg1_type or seg2_type
        flags are not defined. Otherwise return "None" for each undefined flag. 
        """
        return sorted(self.seg_types(require=require))

    # HybRecord : Public Methods : Flag_Info : seg_type
    def set_seg1_type(self, seg1_type):
        """Set "seg1_type" flag in flags."""
        self.set_flag('seg1_type', seg1_type)

    # HybRecord : Public Methods : Flag_Info : seg_type
    def set_seg2_type(self, seg2_type):
        """Set "seg2_type" flag in flags."""
        self.set_flag('seg2_type', seg2_type)

    # HybRecord : Public Methods : Flag_Info : seg_type
    def set_seg_types(self, seg1_type, seg2_type):
        """Set "seg1_type" and "seg2_type" flags in flags."""
        self.set_flag('seg1_type', seg1_type)
        self.set_flag('seg2_type', seg2_type)

    # HybRecord : Public Methods : Flag_Info : read_count
    def read_count(self, require=False):
        """If the "read_count" flag is defined, return it in integer form.
        If require is provided as True, raise an error if the read_count flag is not defined.
        Otherwise return "None"
        """
        if require:
            return int(self._get_flag('read_count'))
        else:
            ret_val = self._get_flag_or_none('read_count')
            if ret_val is None:
                return None
            else:
                return int(ret_val)

    # HybRecord : Public Methods : Flag_Info : read_count
    def set_read_count(self, read_count):
        """Set "read_count" flag in flags."""
        self.set_flag('read_count', str(read_count))

    # HybRecord : Public Methods : Flag_Info : record_count
    def record_count(self, require=False):
        """If the "count_total" flag is defined, return it as an int.
        If require is provided as True, raise an error if the count_total is not defined.
        Otherwise return 1, indicating that the record contains only itself.
        """
        if require:
            return int(self._get_flag('count_total'))
        else:
            ret_val = self._get_flag_or_none('count_total')
            if ret_val is None:
                ret_val = 1
            return ret_val

    # HybRecord : Public Methods : Flag_Info : record_count
    count_total = record_count

    def count(self, count_mode):
        if count_mode in {'read', 'read_count'}:
            return self.read_count(require=True)
        elif count_mode in {'record', 'record_count', 'total', 'count_total'}:
            return self.record_count(require=False)
        else:
            message = 'Unrecognized Count Mode: "%s"\n' % count_mode
            message += 'Allowed options are: %s' % ', '.join(['read', 'record'])
            print(message)
            raise Exception(message)

    # HybRecord : Public Methods : Flag_Info : find_seg_type
    def find_seg_types(self, allow_unknown=False):
        """Find the types of each segment using the method currently set for the class.
        The default supplied method is HybRecord.find_seg_type_hyb, and works with alignemnt
        mapping identifiers in the format of the reference database provided by the Hyb
        Software Package. Custom methods with more complex behavior can be supplied via the
        "set_find_method" method.
        If allow_unknown is False, an error will be raised if a segment type cannot be identified.
        If allow_unknown is True, unidentified segments will be designated as "unknown".
        """
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

    # HybRecord : Public Methods : Flag_Info : find_seg_type
    def find_seg_type_hyb(self, seg_info, find_type_params={}):
        """Return the type of the provided segment, or return None if the segment cannot be
        identified. This method works with sequence / alignment mapping identifiers
        in the format of the reference database provided by the Hyb Software Package,
        specifically identifiers of the format: "AAAA_BBBB_CCCC_DDDD" This method returns
        the fourth component of the identifier, split by "_", as the identfied sequence.
        For Example, "MIMAT0000076_MirBase_miR-21_microRNA" is identified as "microRNA".
        This method does not utilize the find_type_params arg.
        """
        split_id = seg_info['ref'].split('_')
        if len(split_id) != 4:
            return None
        elif not bool(split_id[3]):
            return None
        else:
            return split_id[-1]

    # Set find_seg_type_method as the defulat seg_type finding method.
    find_type_method = find_seg_type_hyb

    # HybRecord : Public Methods : Flag_Info : find_seg_type
    def find_seg_type_string_match(self, seg_info, find_type_params={}, check_all=False):
        """Return the type of the provided segment, or return None if the segment cannot be
        identified.
        This method attempts to find a string matching a specific pattern within the identifier
        of the aligned segment. Search options include "prefix", "contains", "suffix", and
        "matches". The required find_type_params dict should contain a key for each desired
        search type, with a list of 2-tuples for each search-string with assigned-type.
        For example:
        find_type_params = {'suffix': [('_miR', 'microRNA'),
                                       ('_trans', 'mRNA')   ]}
        This dict can be generated with the associated make_string_match_parameters()
        method and an associated csv legend file with format:
            #commentline
            #search_type,search_string,seg_type
            suffix,_miR,microRNA
            suffix,_trans,mRNA
        If check_all is provided as true, the method will continue checking search oftens after
        an option has been found, to ensure that no options conflict.
        """
        seg_name = seg_info['ref']
        found_types = []
        check_done = False
        if not check_done and 'prefix' in find_type_params:
            for search_string, search_type in find_type_params['prefix']:
                if seg_name.startswith(search_string):
                    found_types.append(search_type)
                    if not check_all:
                        check_done = True
                        break
        if not check_done and 'contains' in find_type_params:
            for search_string, search_type in find_type_params['contains']:
                if search_string in seg_name:
                    found_types.append(search_type)
                    if not check_all:
                        check_done = True
                        break
        if not check_done and 'suffix' in find_type_params:
            for search_string, search_type in find_type_params['suffix']:
                if seg_name.endswith(search_string):
                    found_types.append(search_type)
                    if not check_all:
                        check_done = True
                        break
        if not check_done and 'matches' in find_type_params:
            for search_string, search_type in find_type_params['matches']:
                if search_string == seg_name:
                    found_types.append(search_type)
                    if not check_all:
                        check_done = True
                        break

        if not found_types:
            return None
        elif len(found_types) == 1:
            return found_types[0]
        elif len(found_types) > 1:
            # If multiple types found, check if they are the same.
            if all(((found_type == found_types[0]) for found_type in found_types)):
                return found_types[0]
            else:
                message = 'Multiple sequence types found for item: %s' % seg_name
                message += '  ' + ', '.join(found_types)
                print(message)
                raise Exception(message)


    # HybRecord : Public Methods : Flag_Info : find_seg_type
    def find_seg_type_from_id_map(self, seg_info, find_type_params={}):
        """Return the type of the provided segment, or return None if the segment cannot be
        identified.
        This method checks to see if the identifer of the segment is present in a list provided 
        in find_type_params. find_type_params should be formatted as a dict with keys as 
        sequence identifier names, and the corresponding type as the respective values.
        For example:
        find_type_params = {'MIMAT0000076_MirBase_miR-21_microRNA': 'microRNA',
                            'ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA': 'rRNA'}
        This dict can be generated with the associated make_seg_type_id_map()
        method.
        """
        seg_name = seg_info['ref']
        if seg_name in find_type_params:
            return find_type_params[seg_name]
        else:
            return None        

    # HybRecord : Public Methods : Flag_Info : find_seg_type
    # Populate the list of builtin methods available to assign segment types
    find_type_methods = {'hyb': find_seg_type_hyb,
                         'string_match': find_seg_type_string_match,
                         'id_map': find_seg_type_from_id_map}

    # HybRecord : Public Methods : fold_record
    def check_fold_seq_match(self, fold_record):
        """Return True if the sequence (".seq") attribute of a FoldRecord instance matches the
        sequence (".seq") attribute of this instance.
        """
        return (self.seq == fold_record.seq)

    # HybRecord : Public Methods : fold_record
    def set_fold_record(self, fold_record,
                        allow_fold_record_mismatch=None,
                        warn_fold_record_mismatch=None):
        """Check to ensure that fold_record argument is an instance of FoldRecord, and that
        it has a matching sequence to this HybRecord, then set it as self.fold_record.
        allow_fold_record_mismatch allows mismatches between the HybRecord sequence
        and the FoldRecord sequence.
        warn_fold_record_mismatch prints a warning when there is a mismatch between the
        HybRecord sequence and hte FoldRecord sequence.
        If either argument is provided to the method, it overrides default behavior.
        Otherwise, the method falls back to the class default setting.
        """
        if allow_fold_record_mismatch is None:
            allow_fold_record_mismatch = self.settings['allow_fold_record_mismatch']
        if warn_fold_record_mismatch is None:
            warn_fold_record_mismatch = self.settings['warn_fold_record_mismatch']

        if fold_record is None or (isinstance(fold_record, tuple) and fold_record[0] is None):
            if not allow_fold_record_mismatch:
                message = 'Trying to assign None-object as FoldRecord is disallowed with '
                message += 'option "allow_fold_record_mismatch" as true.'
                print(message)
                raise Exception(message)
            elif warn_fold_record_mismatch:
                print('WARNING: Assigning None-Type as fold_record.')        
            self.fold_seq_match = False
            self.fold_record = None
            self.fold_record_str = fold_record[1]
            return

        if fold_record is not None and not isinstance(fold_record, FoldRecord):
            message = 'Supplied argument to fold_record: %s' % str(fold_record)
            message += '\n   is not a FoldRecord object.'
            print(message)
            raise Exception(message)
        if self.check_fold_seq_match(fold_record):
            self.fold_seq_match = True
        else:
            self.fold_seq_match = False
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
                # message += '\n%s\n%s' % (self.seq, fold_record.seq)
                # message += '\nMismatch Region: %s | %s' % (seq1, seq2)
                print(message)

        self.fold_record = fold_record

    # HybRecord : Public Methods : mir_analysis
    def mirna_analysis(self, mirna_types=None):
        """Analyze and store miRNA properties from other properties in the hyb record.

        Perform an analysis of miRNA properties within the sequence record, and store
        the results in the miRNA_seg flag, and in the miRNA_analysis dict.
        This analysis requries the seg1_type and seg2_type flags to be populated, which
        can be performed by the ".find_seg_types()" method.
        If mirna_types is provided, miRNA types will be assigned using this list of types.
        Otherwise it will be assigned using the entries in HybRecord.MIRNA_TYPES .
        # If allow_mirna_dimers is provided as True, entries where seg1_type and seg2_type are
        # both identified as a miRNA will be included in miRNA analyses, where they would
        # otherwise be excluded to prevent ambiguity.
        """

        if mirna_types is None:
            mirna_types = self.MIRNA_TYPES

        # If miRNA_seg flag is not defined, find and set flag.
        if not self.has_property('has_mirna_seg'):
            seg_types = self.seg_types()
            if any(seg_type is None for seg_type in seg_types):
                message = 'Segment types not assigned for mirna analysis of record: %s' % str(self)
                message += 'Segment types must be previously defined for miRNA analysis.'
                message += 'Please use the "find_seg_types" method to assign segment types.'
                print(message)
                raise Exception(message)

            seg1_is_mirna = seg_types[0] in mirna_types
            seg2_is_mirna = seg_types[1] in mirna_types

            if seg1_is_mirna and seg2_is_mirna:
                mirna_flag = 'B'
            elif seg1_is_mirna:
                mirna_flag = '5p'
            elif seg2_is_mirna:
                mirna_flag = '3p'
            else:
                mirna_flag = 'N'
            self.set_flag('miRNA_seg', mirna_flag)

        # If miRNA_seg flag already defined, recall value.
        else:
            mirna_flag = self._get_flag_or_none('miRNA_seg')

        # Analyze miRNA details
        self.mirna_details = {}

        if mirna_flag == 'B':
            #if allow_dimers:
            #    self.mirna_details['mirna_hybrid'] = True
            #else:
            #    self.mirna_details['mirna_hybrid'] = False
            self.mirna_details['mirna_hybrid'] = True
            self.mirna_details['mirna_seg_type'] = self.seg1_type(require=True)
            self.mirna_details['target_seg_type'] = self.seg2_type(require=True)
            self.mirna_info = self.seg1_info
            self.target_info = self.seg2_info
        elif mirna_flag == '5p':
            self.mirna_details['mirna_hybrid'] = True
            self.mirna_details['mirna_seg_type'] = self.seg1_type(require=True)
            self.mirna_details['target_seg_type'] = self.seg2_type(require=True)
            self.mirna_info = self.seg1_info
            self.target_info = self.seg2_info
        elif mirna_flag == '3p':
            self.mirna_details['mirna_hybrid'] = True
            self.mirna_details['mirna_seg_type'] = self.seg2_type(require=True)
            self.mirna_details['target_seg_type'] = self.seg1_type(require=True)
            self.mirna_info = self.seg2_info
            self.target_info = self.seg1_info
        elif mirna_flag == 'N':
            self.mirna_details['mrina_hybird'] = False
            self.mirna_details['mirna_seg_type'] = None
            self.mirna_details['target_seg_type'] = None
        elif mirna_flag == 'U':
            self.mirna_details['mirna_hybrid'] = False
            self.mirna_details['mirna_seg_type'] = None
            self.mirna_details['target_seg_type'] = None
        else:
            message = 'Problem with mirna_analysis for hybrecord: %s ' % str(self)
            message += 'Undefined value: %s found for flag: miRNA_seg' % mirna_flag
            print(message)
            raise Exception(message)

        self.mirna_details['mirna_fold'] = None
        self.mirna_details['target_fold'] = None
        if self.fold_record is not None and mirna_flag in {'B', '5p', '3p'}:
            if mirna_flag in {'B', '5p'}:
                mirna_details = self.fold_record.seg1_info()
                target_details = self.fold_record.seg2_info()
            elif mirna_flag == '3p':
                target_details = self.fold_record.seg1_info()
                mirna_details = self.fold_record.seg2_info()
            self.mirna_details['mirna_fold'] = mirna_details['seg_fold']
            self.mirna_details['target_fold'] = target_details['seg_fold']

    def target_region_analysis(self, region_info=None, coding_types=None,
                               allow_unknown_regions=None, warn_unknown_regions=None):
        """If the record contains an identified mirna and coding target, 
        find the region in which the targeted sequence resides and store the results in the 
        "target_reg" flag and miRNA_analysis dict.
        This analysis requries the seg1_type and seg2_type flags to be populated, which
        can be performed by the ".find_seg_types()" method. It also requires ".mirna_analysis()"
        to have been performed. If the miRNA_seg flag is one of "3p" or "5p", the target type
        will be checked for membership in the coding types arguement if provided, and 
        default HybRecord.CODING_TYPES if not. if both are true, then the analysis 
        will be performed. 
        The region_info dict should contain keys of transcript identifiers with values as a 
        dict containing containing transcript region information.
        Required keys are cdna_coding_start, cdna_coding_end
        Optional keys are 5_utr_start, 5_utr_end, 3_utr_start, 3_utr_end
        Example:
           {'ENST00000372098': {'5_utr_start':'45340255', '5_utr_end':'45340388',
                                'cdna_coding_start':  'cdna_coding_end':,
                                '3_utr_start':'45329242', '3_utr_end':'45329305'}}
        This dict can be set at the class level using 
        ".make_region_info()', with ".set_region_info()", or with
        ".make_set_region_info()", or can be supplied directly to this method.         
        """

        if coding_types is None:
            coding_types = self.CODING_TYPES

        if region_info is None:
            region_info = self.target_region_info

        if allow_unknown_regions is None:
            allow_unknown_regions = self.settings['allow_unknown_regions']

        if warn_unknown_regions is None:
            warn_unknown_regions = self.settings['warn_unknown_regions']

        # Ensure mirna_analysis has been performed.
        self._ensure_mirna_analysis()
         
        # Get miRNA flag.
        mirna_flag = self._get_flag_or_none('miRNA_seg')
        if (mirna_flag in {'5p', '3p'} 
            and self.mirna_details['target_seg_type'] in coding_types):

            self.mirna_details['target_coding'] = True
            #target = self.target_info['ref']
            target = self.target_info['ref']
            warn = True  # TODO Make Fancy
            if target not in region_info:
                if allow_unknown_regions:
                    if warn_unknown_regions:
                        message = 'WARNING: Unknown transcript identifier: %s' % target
                        print(message)
                    self.mirna_details['target_reg'] = 'unknown'
                    self.set_flag('target_reg', 'Unk')               
                    return
                else:
                    message = 'Problem with target_region_analysis for hybrecord: %s \n' % str(self)
                    message += 'Coding Transcript identifier: "%s" ' % target
                    message += 'cannot be found.\n'
                    message += 'Please check the region_info dict provided to method/class.'
                    print(message)
                    raise Exception(message)

            target_regions = region_info[target]
            required_keys = ['cdna_coding_start', 'cdna_coding_end']
            if not all(key in target_regions for key in required_keys):
                message = 'Problem with target_region_analysis for hybrecord: %s \n' % str(self)
                message += 'Provided region information: \n%s\n' % str(target_regions)
                message += 'Does not have required keys:\n    %s' % ', '.join(required_keys)
                print(message)
                raise Exception(message)

            # Setup Variables
            # five_utr = ('5pUTR', target_regions['5_utr_start'], target_regions['5_utr_end'])
            # three_utr = ('3pUTR', target_regions['3_utr_start'], target_regions['3_utr_end'])
            ref_cds_start = int(target_regions['cdna_coding_start'])
            ref_cds_end = int(target_regions['cdna_coding_end'])
            seg_start = int(self.target_info['ref_start'])
            seg_end = int(self.target_info['ref_end'])

            # Assignment based on lower and higher position number, which can vary based on 
            #   sense orientation of transcript, so determine and assign order:
            # if five_utr[1] < three_utr[1]:
            #     first_utr, second_utr = five_utr, three_utr
            # else:
            #     first_utr, second_utr = three_utr, five_utr
            if ref_cds_start < ref_cds_end:
                rev_orientation = False
                low_cds_bound = ref_cds_start
                high_cds_bound = ref_cds_end
                low_utr = '5pUTR'
                high_utr = '3pUTR'
            else:
                rev_orientation = True
                low_cds_bound = ref_cds_end
                high_cds_bound = ref_cds_start
                low_utr = '3pUTR'
                high_utr = '5pUTR'
                print('Reverse Oreintation Detected!',self.target_info['ref'])
                print('"Start:", %i,  "End:" %i' % (ref_cds_start, ref_cds_end))
                raise Exception

            # Assign Categories
            # if seg_start < first_utr[2]:
            #     target_reg = first_utr[0]
            # elif seg_end > second_utr[1]:
            #     target_reg = second_utr[0]
            # else:
            #     target_reg = 'coding'

            if seg_end < low_cds_bound:
                target_reg = low_utr
            elif seg_end > high_cds_bound:
                target_reg = high_utr
            else:
                target_reg = 'coding'

            #print('cds:', ref_cds_start, ref_cds_end, 'seg:', seg_start, seg_end)
            #print(target_reg)
            #input()

            self.mirna_details['target_reg'] = target_reg 
            self.set_flag('target_reg', target_reg)           
 
        else:
            self.mirna_details['target_coding'] = False
            self.mirna_details['target_reg'] = None
            self.set_flag('target_reg', 'N')

        

    # TODO move to "constants section"
    # Set object of string-comparison properties for the ".has_property()" method.
    _STR_PROPERTIES = {
        'id', 'id_prefix', 'id_suffix', 'id_contains',
        'seg', 'seg_prefix', 'seg_suffix', 'seg_contains',
        'seg1', 'seg1_prefix', 'seg1_suffix', 'seg1_contains',
        'seg2', 'seg2_prefix', 'seg2_suffix', 'seg2_contains',
        'seq', 'seq_prefix', 'seq_suffix', 'seq_contains',
        'seg_type', 'seg_type_prefix', 'seg_type_suffix', 'seg_type_contains',
        'seg1_type', 'seg1_type_prefix', 'seg1_type_suffix', 'seg1_type_contains',
        'seg2_type', 'seg2_type_prefix', 'seg2_type_suffix', 'seg2_type_contains',
    }
    # Set object of non-sstring-comparison properties for the ".has_property()" method.
    _HAS_PROPERTIES = {
        'has_seg1_type', 'has_seg2_type', 'has_seg_types',
        'has_mirna_details', 'has_fold_record', 'has_mirna_seg', 'has_mirna_fold',
    }
    # Set object of miRNA-analysis properties for the ".has_property()" method.
    _MIRNA_PROPERTIES = {
        'has_mirna', 'has_mirna_dimer', 'has_mirna_not_dimer',
        '3p_mirna', '5p_mirna',
        '3p_target', '5p_target',
    }
    # Set object of miRNA-analysis properties for the ".has_property()" method.
    _TARGET_PROPERTIES = {
        'has_target',
        'target_3p_utr', 'target_coding', 'target_5p_utr',
    }
    

    # Set object of all allowed properties for the ".has_property()" method.
    PROPERTIES = _STR_PROPERTIES | _HAS_PROPERTIES | _MIRNA_PROPERTIES | _TARGET_PROPERTIES

    # HybRecord : Public Methods : has_property
    def has_property(self, prop_type, prop_compare=None, allow_unknown=False):
        """Check if HybRecord has property of prop_type defined in list of allowed properties
        stored in record.PROPERTIES. If query property has a comparator, provide this
        in prop_compare.
        If allow_unknown is False, an error will be raised if the requested property is undefined.
        """

        if prop_type not in self.PROPERTIES:
            message = 'Requested Property: %s is not defined. ' % prop_type
            message += 'Available proprties are:\n' + ', '.join(self.PROPERTIES)
            print(message)
            raise Exception(message)

        # Check if a substring compares to a desired property string.
        if prop_type in self._STR_PROPERTIES:
            if not prop_compare:
                message = 'Property: %s  requires a comparison string. ' % prop_type
                message += 'Please provide an argument to prop_compare.'

            prop_type_split = prop_type.split('_')
            if len(prop_type_split) == 1 or prop_type_split[-1] == 'type':
                check_attr = prop_type
                check_type = 'matches'
            else:
                check_attr = '_'.join(prop_type_split[:-1])
                check_type = prop_type_split[-1]

            check_info = None
            multi_check = False
            if check_attr in {'id', 'seq'}:
                check_info = getattr(self, check_attr)
            elif check_attr == 'seg1':
                check_info = self.seg1_info.get('ref')
            elif check_attr == 'seg2':
                check_info = self.seg2_info.get('ref')
            elif check_attr == 'seg':
                check_info = {self.seg1_info.get('ref'),
                              self.seg2_info.get('ref')}
                multi_check = True
            elif check_attr == 'seg1_type':
                check_info = self.seg1_type()
            elif check_attr == 'seg2_type':
                check_info = self.seg2_type()
            elif check_attr == 'seg_types':
                check_info = self.seg_types()
                multi_check = True

            # Wrap single check_info value in a list, if not already.
            if not multi_check:
                check_info = {check_info}

            if not allow_unknown:
                for info in check_info:
                    if info is None:
                        message = 'HybRecord Instance: %s does not have a ' % (str(self))
                        message += 'value for requested property: %s' % check_attr
                        message += '\nTo allow unknown values, require_known as False.'
                        print(message)
                        raise Exception(message)

            if any((val is None) for val in check_info):
                ret_val = None
            elif check_type == 'prefix':
                ret_val = any((val.startswith(prop_compare)) for val in check_info)
            elif check_type == 'contains':
                ret_val = any((prop_compare in val) for val in check_info)
            elif check_type == 'suffix':
                ret_val = any((val.endswith(prop_compare)) for val in check_info)
            elif check_type == 'matches':
                 ret_val = bool(prop_compare in check_info)
                # ret_val = any((prop_copmpare == val) for val in check_info)

        # Check whether a property exists and has a non-None / non-blank value.
        elif prop_type in self._HAS_PROPERTIES:
            if prop_type == 'has_seg1_type':
                datum = self._get_flag_or_none(self, 'seg1_type')
            elif prop_type == 'has_seg2_type':
                datum = self._get_flag_or_none(self, 'seg2_type')
            elif prop_type == 'has_seg_types':
                datum = all((bool(self._get_flag_or_none('seg1_type')),
                             bool(self._get_flag_or_none('seg2_type'))))
            elif prop_type == 'has_mirna_details':
                datum = self.mirna_details
            elif prop_type == 'has_fold_record':
                datum = self.fold_record
            elif prop_type == 'has_mirna_seg':
                datum = self._get_flag_or_none('miRNA_seg')
            elif prop_type == 'has_mirna_fold':
                datum = None
                if self.mirna_details is not None:
                    if 'mirna_fold' in self.mirna_details:
                        datum = self.mirna_details['mirna_fold']
            ret_val = bool(datum)

        # Check mirna-specific properties (requires mirna-analysis)
        elif prop_type in self._MIRNA_PROPERTIES:
            self._ensure_mirna_analysis()

            if prop_type == 'has_mirna':
                ret_val = self.flags['miRNA_seg'] in {'5p', '3p', 'B'}
            elif prop_type == 'has_mirna_dimer':
                ret_val = self.flags['miRNA_seg'] == 'B'
            elif prop_type == 'has_mirna_not_dimer':
                ret_val = self.flags['miRNA_seg'] in {'5p', '3p'}
            elif prop_type == '5p_mirna':
                ret_val = self.flags['miRNA_seg'] in {'5p', 'B'}
            elif prop_type == '3p_mirna':
                ret_val = self.flags['miRNA_seg'] in {'3p', 'B'}
            elif prop_type == '5p_target':
                ret_val = self.flags['miRNA_seg'] == '3p'
            elif prop_type == '3p_target':
                ret_val = self.flags['miRNA_seg'] == '5p'

        elif prop_type in self._TARGET_PROPERTIES:
            self._ensure_mirna_analysis()
            if prop_type == 'has_target':
                ret_val = bool(self.mirna_details['target_coding'])
            elif prop_type == 'target_5p_utr':
                ret_val = (self._get_flag('target_reg') == '5pUTR')
            elif prop_type == 'target_coding':
                ret_val = (self._get_flag('target_reg') == 'coding')
            elif prop_type == 'target_3p_utr':
                ret_val = (self._get_flag('target_reg') == '3pUTR')
        return ret_val         
                 

    # HybRecord : Public Methods : Record Parsing
    def to_line(self, newline=False):
        """Return a Hyb-format string representation of the Hyb record."""
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

    # HybRecord : Public MagicMethods : Comparison
    def __eq__(self, other):
        """Return True if ".id" and ".seq" attributes match."""
        return (self.id == other.id and self.seq == other.seq)

    # HybRecord : Public MagicMethods : Comparison
    def __neq__(self, other):
        """Return False if either ".id" or ".seq" attributes mismatch."""
        return (self.id != other.id or self.seq != other.seq)

    # HybRecord : Public MagicMethods : Evaluation
    def __hash__(self):
        """Return a hash of the record ".id" attribute"""
        return hash(self.id)

    # HybRecord : Public MagicMethods : Evaluation
    def __bool__(self):
        """Return True wherever the class is defined."""
        return True

    # HybRecord : Public MagicMethods : Evaluation
    def __len__(self):
        """Return the length of the genomic sequence"""
        return len(self.seq)

    # HybRecord : Public MagicMethods : Printing
    def __str__(self):
        """Print the identifier of the record."""
        return '<HybRecord ID: %s>' % self.id

    # HybRecord : Public Classmethods : find_type_method
    @classmethod
    def set_find_type_method(cls, find_method, find_params={}):
        """Set the class-level custom method for segment assignemnt to callable method
        in find_method, that has the form: "def my_method(self, seg_info, find_params)".
        This method should return the string of the assigned segment type if found, or a
        None object if the type cannot be found.
        It can also take a dictionary in the "find_params" argument that specifies
        additional or dynamic search properties, as desired.
        """
        cls.find_type_method = types.MethodType(find_method, cls)
        cls.find_type_params = find_params

    # HybRecord : Public Classmethods : find_type_method
    @classmethod
    def select_find_type_method(cls, find_method_name, find_params={}):
        """
        Set the class-level custom method for segment assignemnt to callable method
        in find_method, that has the form: "def my_method(self, seg_info, find_params)".
        Select this method from the methods available in HybRecord.find_type_methods.
        """
        if find_method_name not in cls.find_type_methods:
            message = 'Selected find_seg_type_method: %s is not defined.\n' % find_method_name
            message += 'Allowed Options:' + ', '.join(cls.find_type_methods.keys())
        cls.find_type_method = cls.find_type_methods[find_method_name]
        cls.find_type_params = find_params

    @classmethod
    def make_region_info(cls, region_csv_name, sep=','):
        """
        Make a dict containing information on a coding transcript utr regions from an
        input csv file. The input csv must contain a header line, and must have the columns:
        Required keys are identifier, cdna_coding_start, cdna_coding_end
        Optional keys are 5_utr_start, 5_utr_end, 3_utr_start, 3_utr_end
        Example return dict object:
           {'ENST00000372098': {'5_utr_start':'45340255', '5_utr_end':'45340388',
                                'cdna_coding_start':  'cdna_coding_end':,
                                '3_utr_start':'45329242', '3_utr_end':'45329305'}}
        This dict can then be passed to ".set_region_info()" or supplied directly to
        the ".target_region_analysis()" method.
        """

        data_keys = ['cdna_coding_start', 'cdna_coding_end']
        required_keys = ['identifier'] + data_keys

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
                allow_duplicate_info = False  # TODO Make Fancy
                warn = True   # TODO Make Fancy
                permissive = False  #T TODO Make Fancy
                if seq_id in ret_dict:
                    if allow_duplicate_info:
                        if warn:
                            message = 'WARNING: duplicate entry for transcript: %s detected.' % seq_id
                            print(message)
                            continue
                    else:
                        message = 'Problem with creation of target region information dict.\n'
                        message += 'Provided input target region csv file:\n    %s\n' % region_csv_name
                        message += 'Contains duplicate information for entry:'
                        message += '%s\nProvided at line %i\n' % (seq_id, i)
                        message += 'Previous: %s\n' % str(ret_dict[seq_id])
                        message += 'Current: %s\n' % str(line_items)
                        print(message)
                        raise Exception(message)
                
                data_pairs = [(key, line_items[key]) for key in data_keys]
                bad_entry = False
                for key, item in data_pairs:
                    if not item.isnumeric():
                        if permissive:
                            if warn:
                                message = 'WARNING: Skipping entry for %s ' % seq_id
                                message += 'with non-numeric values.'
                                print(message)
                            bad_entry = True
                            break
                        else:
                            message = 'Data for entry: %s contains non-numeric values.' % seq_id
                            message += 'Data: %s' % str(data_pairs)
                            print(message)
                            raise Exception(message)
                if bad_entry:
                    continue
                ret_dict[seq_id] = {key:int(line_items[key]) for key in data_keys}                

        return ret_dict

    @classmethod
    def set_region_info(cls, region_info_dict):
        """
        Set the class-level reference dict object with information on coding transcript utr 
        regions. This dict must have transcript identifiers as keys, with values of dicts with
        defined keys of: 
        Required keys are cdna_coding_start, cdna_coding_end
        Optional keys are 5_utr_start, 5_utr_end, 3_utr_start, 3_utr_end
        Example dict format:
           {'ENST00000372098': {'5_utr_start':'45340255', '5_utr_end':'45340388',
                                'cdna_coding_start':  'cdna_coding_end':,
                                '3_utr_start':'45329242', '3_utr_end':'45329305'}}
        """

        cls.target_region_info = region_info_dict
    
    @classmethod
    def make_set_region_info(cls, region_csv_name, sep=','):
        """
        Convenience wrapper, sequentially callse ".make_region_info()" followed by 
        ".set_region_info()" to set the HybRecord.target_region_info class dict object 
        required for target_region_analysis. 
        """
        region_info = cls.make_region_info(region_csv_name, sep=sep)
        cls.set_region_info(region_info)

    # HybRecord : Public Classmethods : flags
    @classmethod
    def set_custom_flags(cls, custom_flags):
        """
        Set the class-level HybRecord.custom_flags variable, (and update the
        HybRecord.all_flags variable) to allow custom flags in your Hyb file without
        causing an exception.
        """
        cls.custom_flags = custom_flags
        cls.all_flags = HYB_FLAGS + HYBKIT_FLAGS + cls.custom_flags

    # HybRecord : Public Classmethods : Record Construction
    @classmethod
    def from_line(cls, line, hybformat_id=False, hybformat_ref=False):
        """
        Takes as input a line in .hyb format and returns a HybRecord object containing the 
        line's data.
        The Hyb Software Package contains further information in the "id" field of the
        line that can be used to infer read counts represented by the hyb record.
        If hybformat_id is provided as True, this extra information will be read.
        The Hyb Software Package also utilizes a database by default that contains 
        further information in the names of each respective reference sequence. 
        If hybformat_ref is provided as True, this extra information will be read.
        """
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

        if hybformat_id:
            read_id, read_count = cls._parse_hybformat_id(hyb_id)
            if 'read_count' not in flags:
                flags['read_count'] = read_count
                 
        if hybformat_ref:
            for i, seg_info in enumerate([seg1_info, seg2_info], start=1):
                ref = seg_info['ref']
                seg_type_key = 'seg%i_type' % i
                gene_id, transcript_id, gene_name, seg_type = cls._parse_hybformat_ref(ref)
                if seg_type_key in flags and flags[seg_type_key] != seg_type:
                    message = 'Problem reading in hybformat ref for reference: %s\n' % ref
                    message += 'Inferred type: %s\n' % seg_type
                    message += 'Does not equal current type flag: %s' % flags[seg_string]
                    print(message)
                    raise Exception(message)
                elif seg_type_key not in flags:
                    flags[seg_type_key] = seg_type
              
        return_obj = cls(hyb_id, seq, energy, seg1_info, seg2_info, flags)
        return return_obj

    # HybRecord : Public Staticmethods : find_seg_type
    @staticmethod
    def make_string_match_parameters(
            legend_file=os.path.join(hybkit.__about__.code_dir, 'find_type_string_match.csv')
            ):
        """
        Read csv file provided in legend_file, and return a dict of search parameters
        for use with the find_seg_type_string_match method.
        The my_legend.csv file should have the format:
            #commentline
            #search_type,search_string,seg_type
            suffix,_miR,microRNA
            suffix,_trans,mRNA
        search_type options include "prefix", "contains", "suffix", and "matches"
        The produced dict object contains a key for each search type, with a list of
        2-tuples for each search-string and associated segment-type. For example:
          {'suffix': [('_miR', 'microRNA'),
                      ('_trans', 'mRNA')   ]}
        """

        ALLOWED_SEARCH_TYPES = {'prefix', 'contains', 'suffix', 'matches'}
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
                    message += 'Not in allowed types: %s' % ', '.join(ALLOWED_SEARCH_TYPES.keys())
                    message += '\nFor legend line: \n%s\n' % (str(line))
                    print(message)
                    raise Exception(message)

                if search_type not in return_dict:
                    return_dict[search_type] = []

                return_dict[search_type].append((search_string, seg_type))

        return return_dict


    # HybRecord : Public Staticmethods : find_seg_type
    @staticmethod
    def make_seg_type_id_map(mapped_id_files=None, type_file_pairs=None):
        """
        Read file(s) provided and return a mapping of sequence identifiers to types 
        for use with the find_seg_type_from_id_map method.
        The method requires passing either a list/tuple of one or more files to mapped_id_files,
        or a list/tuple of one or more pairs of file lists and file types 
        passed to type_file_pairs.
        Files listed in the mapped_id_files argument should have the format: 
            #commentline
            #seg_id,seg_type
            seg1_unique_id,seg1_type
            seg2_unique_id,seg2_type
        Entries in the list/tuple passed to type_file_pairs should have the format:
            (seg1_type, file1_name)
            Ex: [(seg1_type, file1_name), (seg2_type, file2_name),]
            The first entry in each (non-commented, non-blank) file line will be read and
            added to the mapping dictionary mapped to the provided seg_type.
        """

        return_dict = {}
        if not any((arg is not None for arg in (mapped_id_files, type_file_pairs))):
            message = 'make_seg_type_id_map function requires either a mapped_id_files '
            message += 'or type_file_pairs arguement.'
            print(message)
            raise Exception(message)
        for arguement in mapped_id_files, type_file_pairs:
            if ((arguement is not None) and not any((isinstance(arguement, allowed_type) 
                                                     for allowed_type in (list, tuple)))):
                message = 'arguements passed to mapped_id_files and type_file_pairs must be '
                message += 'provided as a list or tuple.\n  Current passed aruement: '
                message += str(arguement)
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

    # HybRecord : Private Methods : Initialization
    def _post_init_tasks(self):
        # Stub for subclassing
        pass

    # HybRecord : Private Methods : Record Parsing
    def _format_seg_info(self, seg_info, prefix='', suffix='', indent_str=''):
        # Returns a formatted string of the sgement info information
        ret_string = prefix
        ret_string += indent_str + 'Map Reference:  %s\n' % seg_info['ref']
        ret_string += indent_str + 'Read Start Pos: %s\n' % seg_info['read_start']
        ret_string += indent_str + 'Read End Pos:   %s\n' % seg_info['read_end']
        ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_info['ref_start']
        ret_string += indent_str + 'Map End Pos:    %s\n' % seg_info['ref_end']
        ret_string += indent_str + 'Map Score:      %s\n' % seg_info['score']
        ret_string += suffix
        return ret_string

    # HybRecord : Private Methods : flags
    def _get_flag(self, flag_key):
        if flag_key in self.flags:
            return self.flags[flag_key]
        else:
            message = 'Expected Flag Key: %s, but it is not present in record.' % flag_key
            print(message)
            raise Exception(message)

    # HybRecord : Private Methods : flags
    def _get_flag_or_none(self, flag_key):
        if flag_key in self.flags:
            return self.flags[flag_key]
        else:
            return None

    # HybRecord : Private Methods : flags
    def _make_flag_string(self, use_all_flags=False):
        flag_string = ''
        for flag in self._get_flag_keys():
            flag_string += ('%s=%s;' % (flag, str(self.flags[flag])))
        return flag_string

    # HybRecord : Private Methods : flags
    def _get_flag_keys(self, reorder_flags=None):
        # reorder_flags option returns flags in deafult ordering scheme.
        #  If reorder-flags argument provided, it overrides default behavior.
        #  Otherwise, the method falls back to the object-default.
        return_list = []
        if reorder_flags is None:
            reorder_flags = self.settings['reorder_flags']
        if reorder_flags:
            return_list = self._get_ordered_flag_keys()
        else:
            return_list = self.flags.keys()
        return return_list

    # HybRecord : Private Methods : flags
    def _get_ordered_flag_keys(self):
        return_list = []
        for flag in self.all_flags:
            if flag in self.flags:
                return_list.append(flag)
        for flag in self.flags:
            if flag not in self.all_flags:
                return_list.append(flag)
        return return_list

    # HybRecord : Private Methods : flags
    def _make_flags_dict(self, flag_obj, allow_undefined_flags=None):
        #  allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        #  If either argument is provided to the method, it overrides default behavior.
        #  Otherwise, the method falls back to the object-defaults.
        if allow_undefined_flags is None:
            allow_undefined_flags = self.settings['allow_undefined_flags']

        if not isinstance(flag_obj, dict):
            message = '"flag_obj" argument must be a dict obj. Defined keys are:'
            message += self.all_flags.join(', ')
            print(message)
            raise Exception(message)

        if not allow_undefined_flags:
            for flag in flag_obj:
                if flag not in self.all_flags:
                    message = 'Flag "%s" is not defined. Please check flag key' % flag
                    message += ' or run with: "allow_undefined_flags=True"'
                    print(message)
                    raise Exception(message)
        return flag_obj

    # HybRecord : Private Methods : seg_info
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

    # HybRecord : Private Methods : mirna_details
    def _ensure_mirna_analysis(self):
        if self.mirna_details is None:
            message = 'Problem with HybRecord instance: %s\n' % str(self)
            message += 'Method requries running ".mirna_analysis()" before use.'
            print(message)
            raise Exception(message)

    # HybRecord : Private Classmethods : hybformat record parsing
    @classmethod
    def _parse_hybformat_id(cls, record_id):
        # Parse id in format: "48_50002" into read_id, read_count
        split_id = record_id.split('_')
        if not len(split_id) == 2:
            message = 'Failed attempt to parse record id: %s in hyb format.\n' % record_id
            message += 'Hyb-Program format record ids have form: <read_id>_<read_count>'
            print(message)
            raise Exception(message)
        return (split_id[0], split_id[1])

    # HybRecord : Private Classmethods : hybformat record parsing
    @classmethod
    def _parse_hybformat_ref(cls, seg_ref):
        # Parse reference sequence identifier in format: 
        # "ENSG00000146425_ENST00000367089_DYNLT1_mRNA" 
        # into <gene_id>_<transcript_id>_<gene_name>_<seg_type> information.
        split_ref = seg_ref.split('_')
        if not len(split_ref) == 4:
            message = 'Failed attempt to parse segment reference id: "%s"'  % seg_ref
            message += ' in hyb format.\n'
            message += 'Hyb-Program format record ids have form:\n'
            message += '    <gene_id>_<transcript_id>_<gene_name>_<seg_type>'
            print(message)
            raise Exception(message)
        return (split_ref[0], split_ref[1], split_ref[2], split_ref[3])

    # HybRecord : Private Classmethods : flags
    @classmethod
    def _read_flags(cls, flag_string, allow_undefined_flags=None):
        # allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        # undefined flags allowed in this method by default, to allow the object-level setting to
        # take precedence
        if allow_undefined_flags is None:
            allow_undefined_flags = cls.settings['allow_undefined_flags']

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


class HybFile(object):
    """
    File-Object wrapper that provides abiltity to return file lines as HybRecord entries.
    """

    # HybFile : Class-Level Constants:
    DEFAULTS = {}
    # The Hyb Software Package places further information in the "id" field of the
    #   hybrid record that can be used to infer the number of contained read counts.
    #   Set this value to True with "hybkit.HybFile.hybformat_id = True" to read this
    #   extra information.
    DEFAULTS['hybformat_id'] = False
    # The Hyb Software Package by default uses a reference database with identifiers
    #   that contain sequence type and other information.
    #   Set this value to True with "hybkit.HybFile.hybformat_ref = True" to read this
    #   extra information.
    DEFAULTS['hybformat_ref'] = False

    # HybFile : Class-Level variables:
    settings = copy.deepcopy(DEFAULTS)

    # HybFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Wrapper for open() function that stores resulting file."""
        self.fh = open(*args, **kwargs)

    # HybFile : Public Methods : Initialization / Closing
    def __enter__(self, *args, **kwargs):
        """Open "with" syntax."""
        return self

    # HybFile : Public Methods : Initialization / Closing
    def __exit__(self, type, value, traceback):
        """Close "with" syntax"""
        self.close()

    # HybFile : Public Methods : Initialization / Closing
    def __iter__(self):
        """Return an iterator."""
        return self

    # HybFile : Public Methods : Reading
    def __next__(self):
        """Return next line as HybRecord object."""
        return HybRecord.from_line(self.fh.__next__(),
                                   hybformat_id=self.settings['hybformat_id'],
                                   hybformat_ref=self.settings['hybformat_ref'])

    # HybFile : Public Methods : Reading
    def close(self):
        """Close the file."""
        self.fh.close()

    # HybFile : Public Methods : Reading
    def read_record(self):
        """Return next line of hyb file as HybRecord object."""
        return next(self)

    # HybFile : Public Methods : Reading
    def read_records(self):
        """Return list of all records in hyb file as HybRecord objects."""
        records = []
        for record in self:
            records.append(record)
        return records

    # HybFile : Public Methods : Writing
    def write_record(self, write_record):
        """
        Write a HybRecord object to file as a Hyb-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        """
        self._ensure_HybRecord(write_record)
        record_string = write_record.to_line(newline=True)
        self.fh.write(record_string)

    # HybFile : Public Methods : Writing
    def write_records(self, write_records):
        """
        Write a sequence of HybRecord objects as hyb-format lines to the Hyb file.
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        """
        for write_record in write_records:
            self.fh.write_record(write_record)

    # HybFile : Public Classmethods : Initialization
    @classmethod
    def open(cls, *args, **kwargs):
        """Return a new HybFile object."""
        return cls(*args, **kwargs)

    # HybFile : Private Methods
    def _ensure_HybRecord(self, record):
        if not isinstance(record, HybRecord):
            message = 'Item: "%s" is not a HybRecord object.' % record
            print(message)
            raise Exception(message)


class FoldRecord(object):
    """
    Class for storing secondary structure (folding) information for a nucleotide sequence.
    This class supports:
        The Vienna file format: http://unafold.rna.albany.edu/doc/formats.php#VIENNA
        ex: 34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-...
            TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
            .....((((((.((((((......)))))).))))))   (-11.1)

        The Viennad file format utilizied in the Hyb Software package
        ex:
            34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-34-...
            TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
            TAGCTTATCAGACTGATGT------------------   miR-21_microRNA 1       19
            -------------------TAGCTTATCAGACTGATG   miR-21_microRNA 1       18
            .....((((((.((((((......)))))).))))))   (-11.1)
            [space-line]

        The Ct file format utilized by the UNAFold Software Package.
        ex:
            41	dG = -8	dH = -93.9	seq1_name-seq2_name
            1	A	0	2	0	1	0	0
            2	G	1	3	0	2	0	0
            ...
            40	G	39	41	11	17	39	41
            41	T	40	0	10	18	40	0

    A minimum amount of data necessary for a FoldRecord object is a sequence identifier,
    a genomic sequence, and its fold representaiton.
    """

    # FoldRecord : Class-Level Constants
    # Default settings:
    DEFAULTS = {}
    DEFAULTS['skip_bad'] = False
    DEFAULTS['warn_bad'] = True

    # Placeholder symbol for empty entries. Default is "." in the Hyb software package.
    DEFAULTS['placeholder'] = '.'

    # FoldRecord : Class-Level Settings
    settings = copy.deepcopy(DEFAULTS)

    # FoldRecord : Public Methods : Initialization
    def __init__(self, id, seq, fold, energy,
                 seg1_fold_info={},
                 seg2_fold_info={}):
        self.id = id          # Sequence Name (often seg1name-seg2name)
        self.seq = seq        # Genomic Sequence
        self.fold = fold      # Fold Representation, consisting of '(', '.', and ')' characters
        self.energy = energy  # Predicted energy of folding

        self.set_seg1_fold_info(seg1_fold_info)
        self.set_seg2_fold_info(seg2_fold_info)

    # FoldRecord : Public Methods : seg_info
    def seg1_id(self):
        """Return a copy of the id for segment 1 (5p), or None if not defined."""
        if 'ref' in self.seg1_fold_info:
            return self.seg1_fold_info['ref']
        else:
            return None

    # FoldRecord : Public Methods : seg_info
    def seg2_id(self):
        """Return a copy of the id for segment 2 (3p), or None if not defined."""
        if 'ref' in self.seg2_fold_info:
            return self.seg2_fold_info['ref']
        else:
            return None

    # FoldRecord : Public Methods : seg_info
    def seg_ids(self):
        """Return a tuple of the ids of segment 1 (5p) and segment 2 (3p), or None if not defined."""
        return (self.seg1_id(), self.seg2_id())

    # FoldRecord : Public Methods : seg_info
    def seg1_info(self):
        """Return a copy of the info dict object for segment 1 (5p)."""
        return self.seg1_fold_info.copy()

    # FoldRecord : Public Methods : seg_info
    def seg2_info(self):
        """Return a copy of the info dict object for segment 2 (3p)."""
        return self.seg2_fold_info.copy()

    # FoldRecord : Public Methods : seg_info
    def seg1_detail(self, detail):
        """
        Return a copy of the detail for seg1 provided in by the key in "detail" parameter,
        or if it does not exist return None.
        """
        return self._get_segN_detail(1, detail)

    # FoldRecord : Public Methods : seg_info
    def seg2_detail(self, detail):
        """
        Return a copy of the detail for seg2 provided in by the key in "detail" parameter,
        or if it does not exist return None.
        """
        return self._get_segN_detail(2, detail)

    # FoldRecord : Public Methods : seg_info
    def set_seg1_fold_info(self, seg_info_obj):
        """Set folding information for segment 1"""
        self._set_segN_fold_info(1, seg_info_obj)

    # FoldRecord : Public Methods : seg_info
    def set_seg2_fold_info(self, seg_info_obj):
        """Set folding information for segment 2"""
        self._set_segN_fold_info(2, seg_info_obj)

    # FoldRecord : Public Methods : seg_info
    def set_segs_fold_info_from_hybrecord(self, hybrecord):
        raise NotImplementedError

    # FoldRecord : Public Methods : Parsing : Vienna
    def to_vienna_lines(self, newline=False):
        """Return a list of lines for the record in vienna format."""
        ret_lines = []
        suffix = ''
        if newline:
            suffix = '\n'
        ret_lines.append(self.id + suffix)   # Add line 1, id
        ret_lines.append(self.seq + suffix)  # Add line 2, sequence

        # Create formatted energy string which uses no decimal places for integer numbers
        if abs(self.energy - round(self.energy)) > 0.00001:
            energy_str = ("%.5f" % self.energy).rstrip('0')
        else:
            energy_str = "%i" % int(round(self.energy))

        line_3 = '%s\t(%s)' % (self.fold, energy_str)
        ret_lines.append(line_3 + suffix)    # Add line 3, fold representaiton and energy
        return ret_lines

    # FoldRecord : Public Methods : Parsing : Vienna
    def to_vienna_string(self, newline=False):
        """return a 3-line string for the record in vienna format."""
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_vienna_lines(newline=False)) + suffix)

    # FoldRecord : Public Methods : Parsing : Viennad
    def to_viennad_lines(self, newline=False):
        """Return a list of lines for the record in viennad format."""
        ret_lines = []
        if not newline:
            suffix = ''
        elif newline:
            suffix = '\n'
        ret_lines.append(self.id + suffix)   # Add line 1, id
        ret_lines.append(self.seq + suffix)  # Add line 2, sequence
        line_3_items = [self.seg1_detail('highlight'),
                        self.seg1_detail('ref_short'),
                        self.seg1_detail('ref_start'),
                        self.seg1_detail('ref_end')]
        line_3 = '\t'.join([str(item) if item is not None
                            else self.settings['placeholder']
                            for item in line_3_items])
        ret_lines.append(line_3 + suffix)
        line_4_items = [self.seg2_detail('highlight'),
                        self.seg2_detail('ref_short'),
                        self.seg2_detail('ref_start'),
                        self.seg2_detail('ref_end')]
        line_4 = '\t'.join([str(item) if item is not None
                            else self.settings['placeholder']
                            for item in line_4_items])
        ret_lines.append(line_4 + suffix)

        # Create formatted energy string which uses no decimal places for integer numbers
        if abs(self.energy - round(self.energy)) > 0.00001:
            energy_str = ("%.5f" % self.energy).rstrip('0')
        else:
            energy_str = "%i" % int(round(self.energy))

        line_5 = '%s\t(%s)' % (self.fold, energy_str)
        ret_lines.append(line_5 + suffix)  # Add line 5, fold representaiton and energy
        ret_lines.append(suffix)           # Add line 6, blank/newline
        return ret_lines

    # FoldRecord : Public Methods : Parsing : Viennad
    def to_viennad_string(self, newline=False):
        """return a 6-line string for the record in viennad format."""
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_viennad_lines(newline=False)) + suffix)

    # FoldRecord : Public Methods : hyb_record
    def check_hyb_record_match(self, hyb_record):
        """
        Return True if the sequence (".seq") attribute of a HybRecord instance matches the
        sequence (".seq") attribute of this instance.
        """
        return (self.seq == hyb_record.seq)

    # FoldRecord : Public MagicMethods : Comparison
    def __eq__(self, other):
        return (self.id == other.id and self.seq == other.seq)

    # FoldRecord : Public MagicMethods : Comparison
    def __neq__(self, other):
        return (self.id != other.id or self.seq != other.seq)

    # FoldRecord : Public MagicMethods : Evaluation
    def __hash__(self):
        return hash(self.id)

    # FoldRecord : Public MagicMethods : Evaluation
    def __bool__(self):
        return True

    # FoldRecord : Public MagicMethods : Evaluation
    def __len__(self):
        return len(self.seq)

    # FoldRecord : Public MagicMethods : Printing
    def __str__(self):
        """Print the identifier of the record."""
        return '<FoldRecord ID: %s>' % self.id

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_lines(cls, record_lines, hybformat_file=False):
        """
        Create a FoldRecord entry from a list of 3 strings corresponding to lines in the
        Vienna format.
        The Hyb Software Package contains further information in the "name" field of the
        vienna record that can be used to infer further information about the fold divisions.
        If hybformat_file is provided as True, this extra information will be read.
        extra information.
        """

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

        seg1_fold_info = {}
        seg2_fold_info = {}

        if hybformat_file:
            seg1_fold_info, seg2_fold_info = self._parse_hybformat_name(rec_id, seq, fold)

        return_obj = cls(rec_id, seq, fold, energy,
                         seg1_fold_info=seg1_fold_info,
                         seg2_fold_info=seg2_fold_info)
        return return_obj

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_string(cls, record_string):
        """
        Create a FoldRecord entry from a string containing 3 lines corresponding to lines in the
        Vienna format.
        """
        lines = record_string.strip().split('\n')
        return cls.from_vienna_lines(lines)

    # FoldRecord : Public Classmethods : Construction : Viennad
    @classmethod
    def from_viennad_lines(cls, record_lines, hybformat_file=False, skip_bad=None, warn_bad=None):
        """
        Create a FoldRecord entry from a list of 5 or 6 strings corresponding to lines in the
        Viennad format.
        """
        if skip_bad is None:
            skip_bad = cls.settings['skip_bad']
        if warn_bad is None:
            warn_bad = cls.settings['warn_bad']

        if len(record_lines) not in {5, 6}:
            if skip_bad:
                if warn_bad:
                    message = 'WARNING: Improperly Viennad: Wrong-Line-Number'
                    print(message)
                return (None, ''.join(record_lines))
            else:
                message = 'Provided Viennad Record Lines:\n'
                message += '\n'.join([line.rstrip() for line in record_lines])
                message += '\n  ... are not in required 5-line or 6-line format.'
                print(message)
                raise Exception(message)

        rec_id = record_lines[0].strip()
        seq = record_lines[1].strip()
        line_3 = record_lines[2].strip()
        line_3_split = line_3.split('\t')
        if hybformat_file and len(line_3_split) != 4:
            message = 'Provided Vienna Record Line 3:\n'
            message += line_3.rstrip() + '\n'
            message += str(line_3_split) + '\n'
            message += '\n  ... does not have required "ACTG---<tab>'
            message += 'name<tab>refstart<tab>refend" format.'
            print(message)
            raise Exception(message)
        seg1_fold_info = {
                          'highlight': line_3_split[0],
                          'ref': line_3_split[1],
                          'ref_short': line_3_split[1],
                          'ref_start': None,
                          'ref_end': None,
                         }
        if len(line_3_split) == 4 and all(x.isnumeric() for x in line_3_split[2:]):
            seg1_fold_info['ref_start'] = line_3_split[2]
            seg1_fold_info['ref_end'] = line_3_split[3]

        line_4 = record_lines[3].strip()
        line_4_split = line_4.split('\t')
        if hybformat_file and len(line_4_split) != 4:
            message = 'Provided Vienna Record Line 4:\n'
            message += line_4.rstrip() + '\n'
            message += str(line_4_split) + '\n'
            message += '\n  ... does not have required "ACTG---<tab>'
            message += 'name<tab>refstart<tab>refend" format.'
            print(message)
            raise Exception(message)

        seg2_fold_info = {
                          'highlight': line_4_split[0],
                          'ref': line_4_split[1],
                          'ref_short': line_4_split[1],
                          'ref_start': None,
                          'ref_end': None,
                         }
        if len(line_4_split) == 4 and all(x.isnumeric() for x in line_4_split[2:]):
            seg2_fold_info['ref_start'] = line_4_split[2]
            seg2_fold_info['ref_end'] = line_4_split[3]

        if hybformat_file:
            seg1_hybformat_info, seg2_hybformat_info = cls._parse_hybformat_name(full_name)
            seg1_fold_info['ref_short'] = seg1_hybformat_info['ref_short']
            seg2_fold_info['ref_short'] = seg2_hybformat_info['ref_short']
         
        line_5 = record_lines[4].strip()
        line_5_split = line_5.split('\t')
        if len(line_5_split) != 2:
            message = 'WARNING: Improper viennad record: Line-5 Formatting Error.\n'
            message += 'Full Record:\n'
            message += ''.join(record_lines) + '\n'
            if skip_bad:
                if warn_bad:
                    print(message)
                return (None, ''.join(record_lines))
            else:
                print(message)
                raise Exception(message)
        fold = line_5_split[0]
        energy = float(line_5_split[1].strip('()'))

        for seg_info in (seg1_fold_info, seg2_fold_info):
            seg_info['seg_fold'] = ''
            for i, highlight_char in enumerate(seg_info['highlight']):
                if highlight_char != '-':
                    seg_info['seg_fold'] += fold[i]

        return_obj = cls(rec_id, seq, fold, energy,
                         seg1_fold_info, seg2_fold_info)
        return return_obj

    # FoldRecord : Public Classmethods : Construction : Viennad
    @classmethod
    def from_viennad_string(cls, record_string):
        """
        Create a FoldRecord entry from a string containing 5 or 6 lines corresponding
        to lines in the Viennad format.
        """
        lines = record_string.strip().split('\n')
        return cls.from_vienna_lines(lines)

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_lines(cls, record_lines, hybformat_file=False):
        """
        Create a FoldRecord entry from a list of an arbitrary number of strings
        corresponding to lines in the ".ct" file format.
        The Hyb Software Package contains further information in the "name" field of the
        ct record that can be used to infer further information about the fold divisions.
        If hybformat_file is provided as True, this extra information will be read.
        extra information.
        """
        header_line = record_lines[0].strip()
        if 'dG' not in header_line:
            message = 'Provided ct Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... do not begin with expected header.'
            print(message)
            raise Exception(message)
        header_items = header_line.split('\t')
        expected_seq_len = int(header_items[0])
        expected_record_lines = expected_seq_len + 1

        if len(record_lines) != expected_record_lines:
            message = 'Provided ct Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... do not match the expected %i ' % expected_seq_len
            message += 'lines from the header.'
            print(message)
            raise Exception(message)

        energy_string = header_items[1]
        energy = float(energy_string.split()[-1])
        enthalpy_string = header_items[2]
        enthalpy = float(enthalpy_string.split()[-1])
        full_name = header_items[3]
        # short_name = full_name.split()[1]

        seq = ''
        fold = ''
        seg1_highlight = ''
        seg1_fold = ''
        seg2_highlight = ''
        seg2_fold = ''
        found_seg2 = False
        last_seg_i = 1

        for i, line in enumerate(record_lines[1:], 1):
            line_split = line.strip().split('\t')
            if len(line_split) != 8:
                message = 'Provided ct Record Line:\n'
                message += line.rstrip() + '\n'
                message += '\n  ... does not have 8-column format'
                print(message)
                raise Exception(message)
            base = line_split[1]
            seq += base
            fold_i = int(line_split[4])
            if fold_i == 0:
                fold_char = '.'
            elif fold_i > i:
                fold_char = '('
            elif fold_i < i:
                fold_char = ')'
            else:
                raise Exception
            fold += fold_char

            seg_i = line_split[5]
            if int(seg_i) < i:
                found_seg2 = True
            if found_seg2 is False:
                seg1_highlight += base
                seg2_highlight += '-'
                seg1_fold += fold_char
            elif found_seg2 is True:
                seg1_highlight += '-'
                seg2_highlight += base
                seg2_fold += fold_char
            last_seg_i = seg_i

        if hybformat_file:
            seg1_fold_info, seg2_fold_info = cls._parse_hybformat_name(full_name)

        else:
            # Guess middle '-' char:
            full_name_split = full_name.split('-')
            full_name_split_len = len(full_name_split)
            slice_i = int((full_name_split_len/2))
            seg1_ref_guess = '-'.join(full_name_split[:slice_i])
            seg2_ref_guess = '-'.join(full_name_split[slice_i:])
    
            seg1_fold_info = {'ref': seg1_ref_guess,
                              'ref_start': None,
                              'ref_end': None}
            seg2_fold_info = {'ref': seg2_ref_guess,
                              'ref_start': None,
                              'ref_end': None, }
     
        seg1_fold_info.update({'highlight': seg1_highlight,
                              'seg_fold': seg1_fold})
        seg2_fold_info.update({'highlight': seg2_highlight,
                              'seg_fold': seg2_fold})
        for param in [seq, fold, seg1_highlight, seg2_highlight, (seg1_fold + seg2_fold)]:
            if len(param) != expected_seq_len:
                raise Exception('Problem in record construction, should not occur.')

        return_obj = cls(full_name, seq, fold, energy,
                         seg1_fold_info=seg1_fold_info,
                         seg2_fold_info=seg2_fold_info)
        return return_obj

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_string(cls, record_string):
        """
        Create a FoldRecord entry from a string containing an arbitrary number of lines
        corresponding to lines in the ".ct" file format.
        """
        lines = record_string.strip().split('\n')
        return cls.from_ct_lines(lines)

    # FoldRecord : Private Methods : seg_info
    def _find_segN_fold_details(self, seg_start_n, seg_end_n):
        raise NotImplementedError

    # FoldRecord : Private Methods : seg_info
    def _set_segN_fold_info(self, seg, seg_info_obj={}):
        # Create a dictionary with segment fold information.
        return_dict = {}
        segment_param_types = {'ref': str,
                               'ref_short': str,
                               'ref_start': int,
                               'ref_end': int,
                               'highlight': str,
                               'seg_fold': str,
                               }
        for param in segment_param_types.keys():
            if param in seg_info_obj:
                param_type = segment_param_types[param]
                try:
                    if seg_info_obj[param] is None:
                        return_dict[param] = None
                    else:
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

    # FoldRecord : Private Methods : seg_info
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

    # FoldRecord : Private Classmethods : Parsing : Input
    @classmethod
    def _parse_hybformat_name(cls, name, seq=None, fold=None):
        # If the provided sequence name is in the Hyb-Program output format required
        #   for parsing, find segment information using name, sequence, and fold, and return
        #   seg_fold_info dicts.

        name = name.strip()
        err_message = 'Provided Record Name:\n'
        err_message += name.rstrip() + '\n'
        err_message += '    ... is being parsed in Hyb-Program format, but is not of the '
        err_message += 'required "seq1-seq2" format, where seq1/seq2 each have the format:\n'
        err_message += '  READID_READCOUNT_GENEID_TRANSRIPTID_SEQID_SEQTYPE_REFSTART_REFEND\n'
        #err_message += '  (disallowing "-" in the sequence name)'

        #name_split = name.split('-')
        ## Check that only a single '-' exists in the name.
        #if len(name_split) != 2:
        #    print(err_message)
        #    raise Exception(err_message)

        ## Check that each segment name contains seven '_' characters (8 divisions).
        #if not all((len(seg.split('_')) == 8 for seg in name_split)):
        #    print(err_message)
        #    raise Exception(err_message)

        #seg1_string, seg2_string = name.split('-')
        #seg1_split = seg1_string.split('_')
        #seg2_split = seg2_string.split('_')

        # Check that the split name has 15 divisions (with one joint middle item)
        name_split = name.split('_')
        if not len(name_split) == 15 :
            print(err_message)
            raise Exception(err_message)

        seg1_split = name_split[:7]
        seg1_split.append(name_split[7].split('-')[0])
        seg2_split = [name_split[7].split('-')[1]]
        seg2_split += name_split[8:]

        seg1_fold_info = {}
        seg1_fold_info['ref'] = '_'.join(seg1_split[2:6])
        seg1_fold_info['ref_short'] = '_'.join(seg1_split[4:6])
        seg1_fold_info['ref_start'] = seg1_split[6]
        seg1_fold_info['ref_end'] = seg1_split[7]
        seg1_len = int(seg1_fold_info['ref_start']) - int(seg1_fold_info['ref_start'])
        if seq is not None:
            seg1_fold_info['highlight'] = seq[1:seg1_len] + ('-' * (len(seg) - seg1_len))
        if fold is not None:
            seg1_fold_info['seg_fold'] = fold[1:seg1_len]

        seg2_fold_info = {}
        seg2_fold_info['ref'] = '_'.join(seg2_split[2:6])
        seg2_fold_info['ref_short'] = '_'.join(seg2_split[4:6])
        seg2_fold_info['ref_start'] = seg2_split[6]
        seg2_fold_info['ref_end'] = seg2_split[7]
        seg2_len = int(seg2_fold_info['ref_start']) - int(seg2_fold_info['ref_start'])
        seg2_start = seg1_len
        if seq is not None:
            seg2_fold_info['highlight'] = ('-' * seg1_len) + seq[seg2_start:]
        if fold is not None:
            seg2_fold_info['seg_fold'] = fold[seg2_start:]

        if seq is not None and (len(seq) != seg1_len + seg2_len):
            message = 'Problem with Hyb-Program-Format name parsing.'
            message += 'Sum of segment lengths: %i + %i ' % (seg1_len, seg2_len)
            message += ' = %i\nDoes not equal sequence length: ' % (seg1_len + seg2_len)
            message += str(len(seq))

        if seq is not None:
            for param in [seg1_highlight, seg2_highlight, (seg1_fold + seg2_fold)]:
                if len(param) != len(seq):
                    raise Exception('Problem in record construction, should not occur.')

        return seg1_fold_info, seg2_fold_info

    # FoldRecord : Private Classmethods : Parsing : Output
    @classmethod
    def _format_seg_info(cls, seg_info, prefix='', suffix='', indent_str=''):
        raise NotImplementedError
        # Returns a formatted string of the sgement info information
        ret_string = prefix
        ret_string += indent_str + 'Map Reference:  %s\n' % seg_info['ref']
        ret_string += indent_str + 'Read Start Pos: %s\n' % seg_info['read_start']
        ret_string += indent_str + 'Read End Pos:   %s\n' % seg_info['read_end']
        ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_info['ref_start']
        ret_string += indent_str + 'Map End Pos:    %s\n' % seg_info['ref_end']
        ret_string += indent_str + 'Map Score:      %s\n' % seg_info['score']
        ret_string += suffix
        return ret_string


class ViennaFile(object):
    """
    File-object wrapper that provides abiltity to return sets of three file lines as
    FoldRecord entries.
    The Hyb Software Package contains further information in the "name" field of the
    vienna record that can be used to infer further information about the fold divisions.
    Set this value to True with hybkit.ViennaFile.settings["hybformat_file"] = True to read this
    extra information.
    """

    # ViennaFile : Class-Level Constants
    DEFAULTS = {}
    DEFAULTS['hybformat_file'] = False

    # ViennaFile : Class-Level Variables
    settings = copy.deepcopy(DEFAULTS)

    # ViennaFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Wrapper for open() function that stores resulting file."""
        self.fh = open(*args, **kwargs)

    # ViennaFile : Public Methods : Initialization / Closing
    def __enter__(self, *args, **kwargs):
        """Open "with" syntax."""
        return self

    # ViennaFile : Public Methods : Initialization / Closing
    def __exit__(self, type, value, traceback):
        """Close "with" syntax"""
        self.close()

    # ViennaFile : Public Methods : Initialization / Closing
    def __iter__(self):
        """Return an iterator."""
        return self

    # ViennaFile : Public Methods : Reading
    def __next__(self):
        """Read next three lines and return output as FoldRecord object."""
        line_1 = next(self.fh)
        line_2 = next(self.fh)
        line_3 = next(self.fh)
        return FoldRecord.from_vienna_lines((line_1, line_2, line_3), 
                                            self.settings['hybformat_file'])

    # ViennaFile : Public Methods : Reading
    def close(self):
        """Close the file."""
        self.fh.close()

    # ViennaFile : Public Methods : Reading
    def read_record(self):
        """Return next three line of vienna file as FoldRecord object."""
        return next(self)

    # ViennaFile : Public Methods : Reading
    def read_records(self):
        """Return list of vienna records in vienna file as FoldRecord objects."""
        records = []
        for record in self:
            records.append(record)
        return records

    # ViennaFile : Public Methods : Writing
    def write_record(self, write_record):
        """
        Write a FoldRecord object to file as a vienna-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        """
        self._ensure_FoldRecord(write_record)
        record_string = write_record.to_vienna_string(newline=True)
        self.fh.write(record_string)

    # ViennaFile : Public Methods : Writing
    def write_records(self, write_records):
        """
        Write a sequence of FoldRecord objects as vienna-format lines to the vienna file.
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        """
        for write_record in write_records:
            self.write_record(write_record)

    # ViennaFile : Public Classmethods : Initialization
    @classmethod
    def open(cls, *args, **kwargs):
        """Return a new ViennaFile object."""
        return cls(*args, **kwargs)

    # ViennaFile : Private Methods : Writing
    def _ensure_FoldRecord(self, record):
        if not isinstance(record, FoldRecord):
            message = 'Item: "%s" is not a FoldRecord object.' % record
            print(message)
            raise Exception(message)


class HybViennaIter(object):
    """
    This class provides an iterator to iterate through a HybFile and ViennaFile simultaneously,
    returning a tuple of hyb_record, fold_record instances on each iteration.
    """
    # HybViennaIter : Public Methods
    def __init__(self, hybfile_handle, viennafile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennafile_handle = viennafile_handle
        self.counter = 0

    # HybViennaIter : Public Methods
    def __iter__(self):
        return self

    # HybViennaIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_vienna_record = next(self.viennafile_handle)
        return (next_hyb_record, next_vienna_record)


class HybViennaCmbIter(object):
    """
    This class provides an iterator to iterate through a HybFile and ViennaFile simultaneously.
    It is presumed that each respective hyb entry corresponds to an aligned vienna entry.
    Each ViennaRecord will be added to the corresponding HybRecord.
    Only the HybRecord entry will then be returned, containing the associated ViennaRecord entry.
    """
    # HybViennaCmbIter : Public Methods
    def __init__(self, hybfile_handle, viennafile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennafile_handle = viennafile_handle
        self.counter = 0

    # HybViennaCmbIter : Public Methods
    def __iter__(self):
        return self

    # HybViennaCmbIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_vienna_record = next(self.viennafile_handle)
        try:
            next_hyb_record.set_fold_record(next_vienna_record)
        except StopIteration:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return next_hyb_record


class ViennadFile(object):
    """
    File-object wrapper that provides abiltity to return sets of six viennad file lines as
    FoldRecord entries.
    The Hyb Software Package contains further information in the "name" field of the
    viennad record that can be used to infer further information about the fold divisions.
    Set this value to True with hybkit.ViennadFile.settings['hybformat_file'] = True to read this
    extra information.
    """
    # ViennadFile : Class-Level Constants
    DEFAULTS = {}
    DEFAULTS['hybformat_file'] = False

    # ViennadFile : Class-Level Variables
    settings = copy.deepcopy(DEFAULTS)

    # ViennadFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Wrapper for open() function that stores resulting file."""
        self.fh = open(*args, **kwargs)

    # ViennadFile : Public Methods : Initialization / Closing
    def __enter__(self, *args, **kwargs):
        """Open "with" syntax."""
        return self

    # ViennadFile : Public Methods : Initialization / Closing
    def __exit__(self, type, value, traceback):
        """Close "with" syntax"""
        self.close()

    # ViennadFile : Public Methods : Initialization / Closing
    def __iter__(self):
        """Return an iterator."""
        return self

    # ViennadFile : Public Methods : Reading
    def __next__(self):
        """
        Call io.FileIO __next__ method for next three six lines and return
        output as FoldRecord object.
        """
        line_1 = next(self.fh)
        line_2 = next(self.fh)
        line_3 = next(self.fh)
        line_4 = next(self.fh)
        line_5 = next(self.fh)
        line_6 = next(self.fh)
        return FoldRecord.from_viennad_lines((line_1, line_2, line_3, line_4, line_5, line_6),
                                             hybformat_file=self.settings['hybformat_file'])

    # ViennadFile : Public Methods : Reading
    def close(self):
        """Close the file."""
        self.fh.close()

    # ViennadFile : Public Methods : Reading
    def read_record(self):
        """Return next six lines of viennad file as FoldRecord object."""
        return next(self)

    # ViennadFile : Public Methods : Reading
    def read_records(self):
        """Return list of viennad records in viennad file as FoldRecord objects."""
        records = []
        for record in self:
            records.append(record)
        return records

    # ViennadFile : Public Methods : Writing
    def write_record(self, write_record):
        """
        Write a FoldRecord object to file as a viennad-format string.
        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        """
        self._ensure_FoldRecord(write_record)
        record_string = write_record.to_viennad_string(newline=True)
        self.fh.write(record_string)

    # ViennadFile : Public Methods : Writing
    def write_records(self, write_records):
        """
        Write a sequence of FoldRecord objects as viennad-format lines to the viennad file.
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        """
        for write_record in write_records:
            self.write_record(write_record)

    # ViennadFile : Public Classmethods : Initialization
    @classmethod
    def open(cls, *args, **kwargs):
        """Return a new ViennadFile object."""
        return cls(*args, **kwargs)

    # ViennadFile : Private Methods
    def _ensure_FoldRecord(self, record):
        if not isinstance(record, FoldRecord):
            message = 'Item: "%s" is not a FoldRecord object.' % record
            print(message)
            raise Exception(message)


class HybViennadIter(object):
    """
    This class provides an iterator to iterate through a HybFile and ViennadFile simultaneously,
    returning a tuple of hyb_record, fold_record instances on each iteration.
    """

    # HybViennaIter : Public Methods
    def __init__(self, hybfile_handle, viennadfile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennadfile_handle = viennadfile_handle
        self.counter = 0

    # HybViennaIter : Public Methods
    def __iter__(self):
        return self

    # HybViennaIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_viennad_record = next(self.viennadfile_handle)
        return (next_hyb_record, next_viennad_record)


class HybViennadCmbIter(object):
    """
    This class provides an iterator to iterate through a HybFile and ViennadFile simultaneously.
    It is presumed that each respective hyb entry corresponds to an aligned viennad entry.
    Each FoldRecord will be added to the corresponding HybRecord.
    Only the HybRecord entry will then be returned, containing the associated FoldRecord entry.
    """
    # HybViennaCmbIter : Public Methods
    def __init__(self, hybfile_handle, viennadfile_handle):
        self.hybfile_handle = hybfile_handle
        self.viennadfile_handle = viennadfile_handle
        self.counter = 0

    # HybViennaCmbIter : Public Methods
    def __iter__(self):
        return self

    # HybViennaCmbIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_viennad_record = next(self.viennadfile_handle)
        try:
            next_hyb_record.set_fold_record(next_viennad_record)
        except StopIteration:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return next_hyb_record


class CtFile(object):
    """
    File-object wrapper that provides abiltity to return sets of ct file lines as
    FoldRecord entries.
    The Hyb Software Package contains further information in the "name" field of the
    ct record that can be used to infer further information about the fold divisions.
    Set this value to True with hybkit.CtFile.settings['hybformat_file'] = True to read this
    extra information.
    """
    # ViennaFile : Class-Level Constants
    DEFAULTS = {}
    DEFAULTS['hybformat_file'] = False

    # ViennaFile : Class-Level Variables
    settings = copy.deepcopy(DEFAULTS)

    # CtFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Wrapper for open() function that stores resulting file."""
        self.fh = open(*args, **kwargs)

    # CtFile : Public Methods : Initialization / Closing
    def __enter__(self, *args, **kwargs):
        """Open "with" syntax."""
        return self

    # CtFile : Public Methods : Initialization / Closing
    def __exit__(self, type, value, traceback):
        """Close "with" syntax"""
        self.close()

    # CtFile : Public Methods : Initialization / Closing
    def __iter__(self):
        """Return an iterator."""
        return self

    # CtFile : Public Methods
    def __next__(self):
        """
        Call return the first line of the next entry.
        Read the expected number of following lines in the entry, and read that number
        lines further. Return lines as FoldRecord object.
        """
        header = next(self.fh)
        record_lines = [header]
        expected_line_num = int(header.strip().split()[0])
        for i in range(expected_line_num):
            record_lines.append(next(self.fh))
        ret_record = FoldRecord.from_ct_lines(record_lines, 
                                              hybformat_file=self.settings['hybformat_file'])
        return ret_record

    # CtFile : Public Methods : Reading
    def close(self):
        """Close the file."""
        self.fh.close()

    # CtFile : Public Methods
    def read_record(self):
        """Return next lines of ct file as FoldRecord object."""
        return next(self)

    # CtFile : Public Methods
    def read_records(self):
        """Return list of records in ct file as FoldRecord objects."""
        records = []
        for record in self:
            records.append(record)
        return records

    # No write_record method is implmeneted for ct files, as the FoldRecord object does not
    #   contain the complete set of ct record information.
    # def write_record(self, write_record):

    # No write_records method is implmeneted for ct files, as the FoldRecord object does not
    #   contain the complete set of ct record information.
    # def write_records(self, write_records):

    # CtFile : Public Classmethods : Initialization
    @classmethod
    def open(cls, *args, **kwargs):
        """Return a new CtFile object."""
        return cls(*args, **kwargs)


class HybCtIter(object):
    """
    This class provides an iterator to iterate through a HybFile and CtFile simultaneously,
    returning a tuple of hyb_record, fold_record instances on each iteration.
    """
    # HybCtIter : Public Methods
    def __init__(self, hybfile_handle, ctfile_handle):
        self.hybfile_handle = hybfile_handle
        self.ctfile_handle = ctfile_handle
        self.counter = 0

    # HybCtIter : Public Methods
    def __iter__(self):
        return self

    # HybCtIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_ct_record = next(self.ctfile_handle)
        return (next_hyb_record, next_ct_record)


class HybCtCmbIter(object):
    """
    This class provides an iterator to iterate through a HybFile and CtFile simultaneously.
    It is presumed that each respective hyb entry corresponds to an aligned ct entry.
    Each FoldRecord will be added to the corresponding HybRecord.
    Only the HybRecord entry will then be returned, containing the associated FoldRecord entry.
    """
    # HybCtCmbIter : Public Methods
    def __init__(self, hybfile_handle, ctfile_handle):
        self.hybfile_handle = hybfile_handle
        self.ctfile_handle = ctfile_handle
        self.counter = 0

    # HybCtCmbIter : Public Methods
    def __iter__(self):
        return self

    # HybCtCmbIter : Public Methods
    def __next__(self):
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_ct_record = next(self.ctfile_handle)
        try:
            next_hyb_record.set_fold_record(next_viennad_record)
        except StopIteration:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return next_hyb_record
