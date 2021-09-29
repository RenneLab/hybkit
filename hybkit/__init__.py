#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
This module contains classes and methods for reading, writing, and manipulating data 
in the ".hyb" genomic sequence format.

This is primarily based on two classes for storage of 
chimeric sequence information and associated fold-information:

+---------------------+----------------------------------------------------------------------+
| :class:`HybRecord`  | Class for storage of hybrid sequence records                         |
+---------------------+----------------------------------------------------------------------+
| :class:`FoldRecord` | Minimal class for storage of predicted RNA                           |
|                     | secondary structure information for chimeric sequence reads          |
+---------------------+----------------------------------------------------------------------+
    
It also includes classes for reading, writing, and iterating over files containing that 
information:

+-------------------------+------------------------------------------------------------------+
| :class:`HybFile`        | Class for reading and writing ".hyb"-format files                |
|                         | containing chimeric RNA sequence information                     |
|                         | as :class:`HybRecord` objects                                    |
+-------------------------+------------------------------------------------------------------+
| :class:`ViennaFile`     | Class for reading and writing Vienna (.vienna)-format files      |
|                         | containing RNA secondary structure information in dot-bracket    |
|                         | format as :class:`FoldRecord` objects                            |
+-------------------------+------------------------------------------------------------------+
| :class:`CtFile`         | Class for reading Connectivity Table (.ct)-format files          |
|                         | containing predicted RNA secondary-structure information         |
|                         | as used by the RNAStructure software package as                  |
|                         | :class:`FoldRecord` objects                                      |
|                         | http://rna.urmc.rochester.edu/Text/File_Formats.html#CT          |
+-------------------------+------------------------------------------------------------------+
| :class:`HybFoldIter`    | Class for concurrent iteration over a :class:`HybFile` and a     |
|                         | :class:`ViennaFile` or :class:`CtFile`                           |
+-------------------------+------------------------------------------------------------------+

Todo:
    Add Hybrecord.to_csv_header()
    Add analysis.target_region analyses.
    Create hyb-format database download script.
    Add user-friendly individual scripts.
    Implement all sample analyses with bash workflows and individual scripts.
    Implement all sample analyses with nextflow workflows and individual scripts.
    Make decision and clean "extra" scripts.

"""

import os
import sys
import io
import types
import csv
import copy
try:
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    bio_module_error = None
except ModuleNotFoundError as bio_module_error:
    Bio, SeqIO, Seq, SeqRecord = None, None, None, None
    

# Perform *Initial* hybkit submodule imports, remainder at code end.
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import hybkit
import hybkit.settings
import hybkit.type_finder
import hybkit.region_finder
import hybkit.__about__
# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__,\
                             __email__, __license__, __maintainer__, __status__, __version__

class HybRecord(object):
    """
    Class for storing and analyzing chimeric (hybrid) RNA-seq reads in ".hyb" format.

    Hyb format entries are a GFF-related file format described by Travis, et al. 
    (see :ref:`References`)
    that contain information about a genomic sequence read identified to be a chimera by 
    anlaysis software. Each line contains 15 or 16 columns separated by tabs ("\\\\t") and provides
    annotations on each components. An example .hyb format line 
    from Gay et al. (See :ref:`References`)::
 
        2407_718\tATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC\t.\tMIMAT0000078_MirBase_miR-23a_microRNA\t1\t21\t1\t21\t0.0027\tENSG00000188229_ENST00000340384_TUBB2C_mRNA\t23\t49\t1181\t1207\t1.2e-06

    These columns are respectively described in hybkit as:

         id, seq, energy, [seg1\_]ref, [seg1\_]read_start, [seg1\_]read_end, [seg1\_]ref_start, 
         [seg1\_]ref_end, [seg1\_]score, [seg2\_]read_start, [seg2\_]read_end, [seg2\_]ref_start, 
         [seg2\_]ref_end, [seg2\_]score, [flag1=val1; flag2=val2;flag3=val3...]"

    The preferred method for reading hyb records from lines is with 
    the :func:`HybRecord.from_line` constructor::

        # line = "2407_718\tATC..."
        hyb_record = hybkit.HybRecord.from_line(line)

    This constructor parses hyb files using the 
    :class:`HybFile` wrapper class described below.
    For example, to print all hybrid identifiers in a ".hyb" file::

        with hybkit.HybFile('path/to/file.hyb', 'r') as hyb_file:
            for hyb_record in hyb_file:
                print(hyb_record.id)

    HybRecord objects can also be constructed directly. A minimum amount of data necessary 
    for a HybRecord object is the genomic sequence and its corresponding identifier. 
    
    Examples:
        ::

            hyb_record_1 = hybkit.HybRecord('1_100', 'ACTG')
            hyb_record_2 = hybkit.HybRecord('2_107', 'CTAG', '-7.3')

    Details about segments are provided via dict objects with the keys
    specific to each segment. Data can be provided either as strings or 
    as floats/integers (where relevant).
    For example, to create a HybRecord object representing the example line given above:: 

        seg1_props = {'ref_name': 'MIMAT0000078_MirBase_miR-23a_microRNA',
                     'read_start': '1',
                     'read_end': '21',
                     'ref_start': '1',
                     'ref_end': '21',
                     'score': '0.0027'}
        seg2_props = {'ref_name': 'ENSG00000188229_ENST00000340384_TUBB2C_mRNA',
                     'read_start': 23,
                     'read_end': 49,
                     'ref_start': 1181,
                     'ref_end': 1207,
                     'score': 1.2e-06}
        seq_id = '2407_718' 
        seq = 'ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC'
        energy = None
        
        hyb_record = hybkit.HybRecord(seq_id, seq, energy, seg1_props, seg2_props)
        # OR
        hyb_record = hybkit.HybRecord(seq_id, seq, seg1_props=seg1_props, seg2_props=seg2_props)

    Args:
        id (str): Identifier for the hyb record
        seq (str): Nucleotide sequence of the hyb record
        energy (str, optional): Predicted energy of sequence folding in kcal/mol
        seg1_props (dict, optional): Properties of segment 1 of the record, containing possible:
            keys: ('ref_name', 'read_start', 'read_end', 'ref_start', 'ref_end', 'score')
        seg2_props (dict, optional): Properties of segment 2 of the record, containing possible:
            keys: ('ref_name', 'read_start', 'read_end', 'ref_start', 'ref_end', 'score')
        flags (dict, optional): Dict with keys of flags for the record and their associated values.
            By default flags must be defined in :attr:`ALL_FLAGS` but custom 
            flags can be supplied in :attr:`settings['custom_flags'] <settings>`.
            This setting can also be disabled by setting 'allow_undefined_flags' 
            to :obj:`True` in :attr:`HybRecord.settings`.
        fold_record (FoldRecord, optional): Set the record's :attr:`fold_record` attribute 
            as the provided FoldRecord object using :meth:`set_fold_record` on initializtaion.

    .. _HybRecord-Attributes:

    Attributes:
        id (str): Identifier for the hyb record (Hyb format: "<read-num>_<read-count>")

        seq (str): Nucleotide sequence of the hyb record
        energy (str or None): Predicted energy of folding
        seg1_props (dict): Information on chimeric segment 1, contains keys: 
            'ref_name' (str), 'read_start' (int), 'read_end' (int), 'ref_start' (int), 
            'ref_end' (int), and 'score' (float).
        seg2_props (dict): Information on segment 2, contains keys:
            'ref_name' (str), 'read_start' (int), 'read_end' (int), 'ref_start' (int), 
            'ref_end' (int), and 'score' (float).
        flags (dict): Dict of flags with possible keys and values as defined in
            the :ref:`Flags` section of the :ref:`Hybkit Hyb File Specification`.
        mirna_props (dict or None): Link to appropriate seg1_props or seg2_props dict 
            corresponding to a record's miRNA (if present), assigned by the 
            :func:`eval_mirna` method.
        target_props (dict or None): Link to appropriate seg1_props or seg2_props dict 
            corresponding to a record's target of a miRNA (if present), assigned by the  
            :func:`eval_mirna` method.
        fold_record (FoldRecord): Information on the predicted secondary structure of the sequence
            set by :func:`set_fold_record`.
    """

    # HybRecord : Class-Level Constants
    #: Record columns 1-3 defining parameters of the overall hybrid, defined by the Hyb format
    HYBRID_COLUMNS = [
                      'id',      #: str, Hybrid read identifier, ex: "1257_12"
                      'seq',     #: str, Hybrid nucleotide sequence, ex: "ATCGGCTAATCGGTCA..."
                      'energy',  #: str(of float), Intra-hybrid folding energy, ex: "-11.33"
                     ]

    #: Record columns 4-9 and 10-15, respectively, defining annotated parameters of  
    #: seg1 and seg2 respectively, defined by the Hyb format
    SEGMENT_COLUMNS = [
                       'ref_name',    # str, Mapping Reference Identity: ex:
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
    _HYB_FLAGS = [
                 'count_total',            # str(int), total represented hybrids
                 'count_last_clustering',  # str(int), total represented hybrids at last clustring
                 'two_way_merged',         # "0" or "1", boolean representation of whether
                                           #   entries with mirrored 5' and 3' hybrids were merged
                 'seq_IDs_in_cluster',     # str, comma-separated list of all ids of hybrids
                                           #   merged into this hybrid entry.
                ]
    # Additional flag specifications utilized by hybkit
    _HYBKIT_FLAGS = [
                    'read_count',   # str(int), number of sequence reads represented by record
                                    #   if merged record, represents total for all merged entries
                    'orient',       # str, orientation of strand in relation to transcript. Options:
                                    #   "F" (Forward / "Sense"),     "IF" (Inferred Forward),
                                    #   "R" (Reverse / "Antisense"), "IR" (Inferred Reverse),
                                    #   "U" (Unknown), or "IC" (Inferred Conflicting)
                    'seg1_type',    # str, assigned type of segment 1, ex: "miRNA" or "mRNA"
                    'seg2_type',    # str, assigned type of segment 2, ex: "miRNA" or "mRNA"
                    'seg1_det',     # str, arbitrary detail about segment 1
                    'seg2_det',     # str, arbitrary detail about segment 2
                    'miRNA_seg',    # str, indicates which (if any) segment mapping is a miRNA
                                    #   options are "N" (none), "5p" (seg1), "3p" (seg2),
                                    #   "B" (both), or "U" (unknown)
                    'target_reg',   # str, assigned region of the miRNA target.
                                    #   options are "5pUTR", "coding", "3pUTR",
                                    #   "N" (none), or "U" (unknown)
                    'extended',     # int, "TRUE" or "FALSE", boolean representation of whether
                                    #   record sequences were bioinformatically extended as is
                                    #   performed by the Hyb software package.
                    'dataset',      # str, label for sequence dataset id (eg. source file), when 
                                    #   combining records from different datasets.
                   ]

    #: Flags defined by the hybkit package. Flags 1-4 are utilized by the Hyb software package.
    #: For information on flags, see the :any:`Flags` portion of the 
    #: :any:`hybkit Hyb File Specification`.
    ALL_FLAGS = _HYB_FLAGS + _HYBKIT_FLAGS

    #: Class-level settings. See :attr:`settings.HybRecord_settings` for descriptions.
    settings = hybkit.settings.HybRecord_settings

    #: Link to :class:`type_finder.TypeFinder` class for parsing sequence identifiers
    #: in assigning segment types by :func:`eval_types`.
    #:
    #: :meta hide-value:
    TypeFinder = hybkit.type_finder.TypeFinder

    #: Link to :class:`region_finder.RegionFinder` class for idenfitying target region
    #: in assigning segment types by :func:`target_region_eval`.
    #: 
    #: :meta hide-value:
    RegionFinder = hybkit.region_finder.RegionFinder

    # Placeholder for set of allowed flags filled on first use.
    _flagset = None

    # HybRecord : Public Methods : Initialization
    def __init__(self, 
                 id, 
                 seq, 
                 energy=None,
                 seg1_props={},
                 seg2_props={}, 
                 flags={},
                 read_count=None,
                 fold_record=None,
                 ):
        self.id = id
        self.seq = seq

        if energy is not None:
            self.energy = energy
        else:
            energy = self.settings['field_placeholder']

        self.seg1_props = self._make_seg_props_dict(seg1_props)
        self.seg2_props = self._make_seg_props_dict(seg2_props)
        self.flags = self._make_flags_dict(flags)

        self.mirna_props = None     # Placeholder variable for eval_mirna
        self.target_props = None    # Placeholder variable for eval_mirna

        if read_count is not None:
            if read_count in flags:
                message = '"read_count" paramater defined both in function call and flags.\n'
                message += 'Please define only once.'
                raise Exception(message)
            else:
                self.set_flag('read_count', int(read_count))

        if fold_record is not None:
            self.set_fold_record(fold_record)
        else:
            self.fold_record = None 

        self._post_init_tasks()    # Method stub for subclassing

    # HybRecord : Public Methods : flags
    def set_flag(self, flag_key, flag_val, allow_undefined_flags=None):
        """Set the value of record flag_key to flag_val.

        Args:
            flag_key (str): Key for flag to set.
            flag_val : Value for flag to set.
            allow_undefined_flags (bool or None, optional): Allow inclusion of flags not  
                defined in :attr:`ALL_FLAGS` or in :attr:`settings['custom_flags'] <settings>`.
                If None (default), uses setting in 
                :attr:`settings['allow_undefined_flags'] <settings>`.
        """

        if allow_undefined_flags is None:
            allow_undefined_flags = self.settings['allow_undefined_flags']

        if self._flagset is None:
            self._flagset = set(self.ALL_FLAGS + list(self.settings.custom_flags))

        if (not allow_undefined_flags and 
            flag_key not in self._flagset):
            message = 'Flag "%s" is not defined. Please check flag key' % flag_key
            message += ' or run with: "allow_undefined_flags=True"'
            print(message)
            raise Exception(message)

        self.flags[flag_key] = flag_val

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg1_type(self, require=False):
        """Return the "seg1_type" flag if defined, or return None.

        Args:
            require (bool, optional): If True, raise an error if seg1_type is not defined.
        """
        if require:
            return self._get_flag('seg1_type')
        else:
            return self._get_flag_or_none('seg1_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg2_type(self, require=False):
        """Return the "seg2_type" flag if defined, or return None.

        Args:
            require (bool, optional): If True, raise an error if seg2_type is not defined.
        """
        if require:
            return self._get_flag('seg2_type')
        else:
            return self._get_flag_or_none('seg2_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg_types(self, require=False):
        """Return a tuple of the ("seg1_type", "seg2_type") flags where defined, or None.

        Args:
            require (bool, optional): If True, raise an error if either flag is not defined.
        """
        if require:
            return (self._get_flag('seg1_type'), self._get_flag('seg2_type'))
        else:
            return (self._get_flag_or_none('seg1_type'), self._get_flag_or_none('seg2_type'))

    # HybRecord : Public Methods : Flag_Info : get_read_count
    def get_read_count(self, require=False):
        """Return the "read_count" flag if defined, otherwise return None.


        Args:
            require (bool, optional): If True, raise an error if the "read_count" flag 
                is not defined.
        """
        
        if require:
            ret_val = self._get_flag('read_count')
        else:
            ret_val = self._get_flag_or_none('read_count')
        
        if ret_val is None:
            return ret_val
        return int(ret_val)

    # HybRecord : Public Methods : Flag_Info : record_count
    def get_record_count(self, require=False):
        """If the "count_total" flag is defined, return it, otherwise return 1 (this record).

        Args:
            require (bool, optional): If True, raise an error if the "read_count" flag 
                is not defined.
        """
        if require:
            ret_val = self._get_flag('count_total')
        else:
            ret_val = self._get_flag_or_none('count_total')

        if ret_val is None:
            ret_val = 1
        return int(ret_val)

    # HybRecord : Public Methods : Flag_Info
    def get_count(self, count_mode, require=False):
        """Return either of :func:`get_read_count` or :func:`get_record_count`.

        Args:
            count_mode (str) : Mode for returned count: one of : {'read', 'record'}
                If 'read', require the 'read_count' flag to be defined.
                If 'record', return '1' if the 'count_total' flag is not defined.
        """

        if count_mode in {'read'}:
            ret_val = self.get_read_count(require=require)
        elif count_mode in {'record'}:
            ret_val = self.get_record_count(require=require)
        else:
            message = 'Unrecognized Count Mode: "%s"\n' % count_mode
            message += 'Allowed options are: %s' % ', '.join(['read', 'record'])
            print(message)
            raise Exception(message)
        return ret_val

    # HybRecord : Public Methods : Flag_Info : find_seg_type
    def eval_types(self, allow_unknown=None, check_complete=None):
        """Find the types of each segment using the the :class:`TypeFinder` class.

        This method provides :attr:`seg1_props` and :attr:`seg2_props` 
        to the :class:`TypeFinder` class, linked as attribute :attr:`HybRecord.TypeFinder`.
        This uses the method: :func:`TypeFinder.method`
        set by :func:`TypeFinder.set_method` or :func:`TypeFinder.set_custom_method` to set the
        :ref:`seg1_type <seg1_type>`, :ref`seg2_type <seg2_type>` flags if not already set. 

        To use a system other than the default, prepare the :class:`TypeFinder` class by 
        preparing and setting :attr:`TypeFinder.params` and using :func:`TypeFinder.set_method`. 

        Args:
            allow_unknown (bool, optional): If True, allow segment types that cannot be
                identified and set them as "unknown". Otherwise raise an error.
                If None (default), uses setting in 
                :attr:`settings['allow_unknown_seg_types'] <settings>`.
            check_complete (bool, optional): If True, check every possibility for the 
                type of a given segment (where applicable), instead of 
                stopping after finding the first type.
                If None (default), uses setting in 
                :attr:`settings['check_complete_seg_types'] <settings>`.
        """

        # If types already set, skip.
        if self.is_set('eval_types'):
            return        

        if allow_unknown is None:
            allow_unknown = self.settings['allow_unknown_seg_types']
        
        if check_complete is None:
            check_complete = self.settings['check_complete_seg_types']

        types = []
        for seg_props in [self.seg1_props, self.seg2_props]:
            seg_type = self.TypeFinder.find(seg_props, self.TypeFinder.params, check_complete)
            if seg_type is None:
                if allow_unknown:
                    types.append('unknown')
                else:
                    message = 'Cannot identify segment type for segment:\n'
                    message += self._format_seg_props(seg_props, prefix=' '*2) + '\n'
                    print(message)
                    raise Exception(message)
            else:
                types.append(seg_type)
        self.set_flag('seg1_type', types[0])
        self.set_flag('seg2_type', types[1])

    # HybRecord : Public Methods : eval_types
    def find_seg_types(self, *args, **kwargs):
        """Find_seg_types method is deprecated. Please use :func:`eval_types`"""
        message = 'find_seg_types method is deprecated. Please use eval_types.'
        print(message)
        raise Exception(message)

    # HybRecord : Public Methods : fold_record
    def set_fold_record(self, fold_record):
        """
        Check and set provided fold_record (:class:`FoldRecord`) as :obj:`fold_record`.
         
        Ensures that fold_record argument is an instance of FoldRecord and 
        has a matching sequence to this HybRecord, then set as self.fold_record.

        Args:
            fold_record (FoldRecord): :attr:`FoldRecord` instance to set as :obj:`fold_record`.
        """
        if fold_record is None or (isinstance(fold_record, tuple) and fold_record[0] is None):
            message = 'Trying to assign None object as FoldRecord.'
            print(message)
            raise Exception(message)
        if not isinstance(fold_record, FoldRecord):
            message = 'Supplied argument to fold_record: %s' % str(fold_record)
            message += '\n   is not a FoldRecord object.'
            print(message)
            raise Exception(message)
        if not self.seq == fold_record.seq:
            message = 'Disallowed mismatch between HybRecord sequence and FoldRecord sequence.'
            message += 'Hyb : %s\n' % str(self.seq)
            message += 'Fold: %s\n' % str(fold_record.seq) 
            print(message)
            raise Exception(message)
        self.fold_record = fold_record

    # HybRecord : Public Methods : eval_mirna
    def eval_mirna(self, mirna_types=None):
        """Analyze and set mstore miRNA properties from other properties in the hyb record.

        If not already done, determine whether a miRNA exists within this record and 
        set the :ref:`miRNA_seg <mirna_seg>` flag.
        This evaluation requries the :ref:`seg1_type <seg1_type>` and :ref:`seg2_type <seg2_type>` flags to 
        be populated, which can be performed by the :func:`eval_types` method.
        If the record contains a miRNA, link the :attr:`mirna_props` and :attr:`target_props` 
        dicts to the corresponding :attr:`seg1_props` / :attr:`seg2_props` dicts as appropriate.

        Args:
            mirna_types (list, tuple, or set, optional): Iterable of strings of "types" to be 
                considered as miRNA. Otherwise, the default types are used 
                from :attr:`settings['mirna_types'] <settings>`.
        """

        if mirna_types is None:
            mirna_types = self.settings['mirna_types']

        # If miRNA_seg flag is not defined, find and set flag.
        if not self.is_set('eval_mirna'):
            self._ensure_set('eval_types')
            seg_types = self.get_seg_types()

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

        if self.has_prop('mirna_not_dimer'):
            if self.has_prop('5p_mirna'):
                self.mirna_props = self.seg1_props
                self.target_props = self.seg2_props
            elif self.has_prop('3p_mirna'):
                self.mirna_props = self.seg2_props
                self.target_props = self.seg1_props


    def mirna_detail(self, detail='all', allow_mirna_dimers=False):
        """Provide a detail about the miRNA or target following :func:`eval_mirna`.

        Analyze miRNA properties within the sequence record and provide a detail as ouptut.
        Method requires record to contain a non-dimer miRNA, otherwise will produce an error.

        Args:
            detail (str): | Type of detail to return. Options include: 
                          | 'all' : (dict of all properties, default); 
                          | 'mirna_name'      : Identifier for Assigned miRNA;
                          | 'target_name'     : Identifier for Assigned Target;
                          | 'mirna_seg_type'  : Assigned seg_type of miRNA;
                          | 'target_seg_type' : Assigned seg_type of target;
                          | 'mirna_seq'       : Annotated subsequence of miRNA;
                          | 'target_seq'      : Annotated subsequence of target;
                          | 'mirna_fold'      : Annotated fold substring of miRNA (requires fold_record set);
                          | 'target_fold'     : Annotated fold substring target (requires fold_record set);
            allow_mirna_dimers (bool, optional): Allow miRNA/miRNA dimers. 
                5p-position will be assigned as miRNA.

        """

        self._ensure_set('eval_mirna')
        mirna_flag = self._get_flag('miRNA_seg')
        

        if ((not self.has_prop('has_mirna')) or
            (not allow_mirna_dimers and not self.has_prop('mirna_not_dimer'))):
            message = 'mirna_detail method requires a hybrid containing a single mirna.\n'
            message += 'hybrecord: %s does not meet this criteria ' % str(self)
            message += 'with miRNA_seg flag: %s' % mirna_flag
            print(message)
            raise Exception(message)


        ALLOWED_DETAILS = {'all', 'mirna_name', 'target_name', 'mirna_seg_type', 'target_seg_type', 
                           'mirna_seq', 'target_seq', 'mirna_fold', 'target_fold',
                          } 

        if detail not in ALLOWED_DETAILS:
            message = 'Requested miRNA detail: "%s" ' % detail
            message += 'not in allowed types: \n    %s' % ', '.join(ALLOWED_DETAILS)
            print(message)
            raise Exception(message)

        # Analyze miRNA details
        mirna_details = {}

        if mirna_flag in {'5p', 'B'}:
            mirna_details['mirna_seg_type'] = self.get_seg1_type(require=True)
            mirna_details['target_seg_type'] = self.get_seg2_type(require=True)
            mirna_props = self.seg1_props
            target_props = self.seg2_props
        elif mirna_flag == '3p':
            mirna_details['mirna_seg_type'] = self.get_seg2_type(require=True)
            mirna_details['target_seg_type'] = self.get_seg1_type(require=True)
            mirna_props = self.seg2_props
            target_props = self.seg1_props
        else:
            message = 'Problem with eval_mirna for hybrecord: %s ' % str(self)
            message += 'Undefined value: %s found for flag: miRNA_seg' % mirna_flag
            print(message)
            raise Exception(message)

        mirna_details['mirna_name'] = mirna_props['ref_name']
        mirna_details['target_name'] = target_props['ref_name']
        mirna_details['mirna_seq'] = self._get_seg_seq(mirna_props)
        mirna_details['target_seq'] = self._get_seg_seq(target_props)
        if self.fold_record is not None:
            mirna_details['mirna_fold'] = self.fold_record._get_seg_fold(miRNA_props)
            mirna_details['target_fold'] = self.fold_record._get_seg_fold(target_props)
        else:
            mirna_details['mirna_fold'] = None
            mirna_details['target_fold'] = None

        if detail == 'all':
            return mirna_details
        else:
            return mirna_details[detail]


    # HybRecord : Public Methods : eval_target
    def prep_target_region_eval(self, region_info=None, region_csv_name=None, sep=','): 
        """
        Prepare from an input csv and assign (or only assign) a dict with information 
        on coding transcript regions. This method is required to be used before
        performing a target region evaluation with :func:`target_region_eval`.
        
        Example:
            Example format of prepared/input dict object::

                region_info = {'ENST00000372098': {'cdna_coding_start':'45340255',
                                                   'cdna_coding_end':'45340388'}}

        This method requires one of the "region_info" (supply dict directly) or 
        "region_csv_name" (csv_file for preparation) arguments to be filled.
        
        If supplying the a csv file,
        the input csv must contain a header line, and must contain the columns::

            identifier,cdna_coding_start,cdna_coding_end

        Args:
            region_info (dict, optional): Dict of region information to set as 
                :attr:`region_finder.RegionFinder.region_info`.
            region_csv_name (str, optional): String of path to csv file to read information from.
            sep (str, optional): Separator for columns of input delimited file. (Default: ',')

        """
        if region_info is None and region_csv_name is None:
            message = 'One of "region_info" or "region_csv_name" arguments must be supplied '
            message += 'to prep_target_region_eval() method.\n'
            print(message)
            raise Exception(message)

        elif region_info:
            self.RegionFinder.set_region_info(region_info)
        elif region_csv_name:
            self.RegionFinder.make_set_region_info(region_csv_name, sep=sep)
        else:
            raise Exception()

    # HybRecord : Public Methods : eval_target
    def eval_target(self, coding_types=None,
                    allow_unknown_regions=None, warn_unknown_regions=None):
        """
        For miRNA/coding-target pairs, find the region of the coding transcript targeted.

        The evaluation requires a dict containing region 
        information to be set using the :func:`prep_target_region_eval` method.
        
        If the record contains an identified mirna and identified coding target, 
        find the region in which the targeted sequence resides and store the results in the 
        :ref:`target_reg <target_reg>` flag and miRNA_eval dict.
        This evaluation requries the :ref:`seg1_type <seg1_type>`, :ref`seg2_type <seg2_type>`, 
        and :ref:`miRNA_seg <mirna_seg>` flags to be populated. This can be performed 
        by sequentially using the :func:`eval_types` and :func:`eval_mirna` methods.
        If the :ref:`miRNA_seg <mirna_seg>` flag is in {"N" (None), "B" (Both)},
        the :ref:`target_reg <target_reg>` flag will be set to {"N" (None)}. 
        If the :ref:`miRNA_seg <mirna_seg>` flag == "U" (Unknown),
        the :ref:`target_reg <target_reg>` flag will be set to {"U" (Unknown)}. 
        If the :ref:`miRNA_seg <mirna_seg>` flag is in {"3p", "5p"},
        the :ref:`target_reg <target_reg>` flag will be checked if it is a coding type.
        If the target is a coding type the evaluation will be performed and the 
        :ref:`target_reg <target_reg>` flag will be set appropriately.


        Args:
            coding_types (iterable, optional): Iterable of strings representing sequence
                types to be recognized as coding.
                If None (default), uses :attr:`settings['coding_types'] <settings>`.
            allow_unknown_regions (bool, optional):
                Allow missing identifiers in evaluation by skipping sequences instead of 
                raising an error.
                If None (default), uses setting in 
                :attr:`settings['allow_unknown_regions'] <settings>`. 
            allow_unknown_regions (bool, optional):
                Warn for missing identifiers in evaluation by printing a message.
                If None (default), uses setting in 
                :attr:`settings['allow_undefined_flags'] <settings>`.
        """
        if coding_types is None:
            coding_types = self.settings['coding_types']

        if allow_unknown_regions is None:
            allow_unknown_regions = self.settings['allow_unknown_regions']

        if warn_unknown_regions is None:
            warn_unknown_regions = self.settings['warn_unknown_regions']

        if self.region_finder.region_info == {}:
            message = 'Target region information has not yet been prepared.\n'
            message += 'Please prepare with "prep_target_region_eval()" before use.'
            print(message)
            raise Exception(message)

        # Ensure eval_types and eval_mirna has been performed.
        self._ensure_set('eval_types')
        self._ensure_set('eval_mirna')
        
        # Get miRNA flag.
        mirna_flag = self._get_flag('miRNA_seg')
        # If miRNA is unknown, target is unknown
        if mirna_flag == 'U':
           target_reg = 'U'
           reason = 'miRNA status is unknown (miRNA_flag = "U").'
        # If no miRNA or both miRNA, there is no target.
        elif miRNA_flag in {'N', 'B'}:
           target_reg = 'N' 
        # If miRNA and target, identify relevant seg_props
        if mirna_flag in {'5p', '3p'}:
            target_type = self.mirna_detail('target_seg_type')
            if target_type not in coding_types:
                target_reg = 'NON'
            else:
                target_name = self.target_props['ref_name']
                target_start = self.target_props['ref_start']
                target_end = self.target_props['ref_end']
                target_reg = self.region_finder.find(target_type, target_start, target_end)
                if target_reg is None:
                    target_reg = 'U'
                    reason = 'Target ID not in supplied region info (region_info): ' + target_name
                    
        if target_reg == 'U':
            if allow_unknown_regions:
                if warn_unknown_regions:
                    message = 'WARNING: target_region_eval: ' + reason
                    print(message)
    
            else:
                message = 'Problem with target_region_eval for hybrecord: %s \n' % str(self)
                message += reason + '\n'
                print(message)
                raise Exception(message)

        self.set_flag('target_reg', target_reg)               

    # HybRecord : Public Methods : Record Properties
    def is_set(self, prop):
        """
        Return True if HybRecord property "prop" is set (if relevant) and is not None.
        Options described in :attr:`IS_SET_PROPERTIES`.

        Args:
            prop (str): Property / Analysis to check
        """

        if prop not in self.IS_SET_PROPERTIES:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available proprties are:\n' + ', '.join(self.IS_SET_PROPERTIES)
            print(message)
            raise Exception(message)

        if prop in {'energy', 'fold_record'}:
            ret_bool = (getattr(self, prop) is not None)
        elif prop == 'eval_types':
            ret_bool = all(st is not None for st in self.get_seg_types())
        elif prop == 'eval_mirna':
            ret_bool = self._get_flag_or_none('miRNA_seg') is not None
        elif prop == 'eval_target':
            ret_bool = self._get_flag_or_none('target_reg') is not None
        return ret_bool
        

    # HybRecord : Public Methods : Record Properties
    def has_prop(self, prop, prop_compare=None):
        """
        Check if HybRecord has property of "prop_type". 
 
        Check property against list of allowed properties in :attr:`PROPERTIES`.
        If query property has a comparator, provide this in prop_compare.
        Raises an error if property is not set (use :func:`is_set` to check).
        
        Args:
            prop: Property to check
            prop_compare (optional): Optional comparator to check.

        Examples:
            General Record Properties::

                # hyb_record = hybkit.HybRecord(id, seq....)
                is_id = hyb_record.has_prop('id', 'target_identifier')
                seq_is_ATCG = hyb_record.has_peroperty('seq', 'ATCG')
                seq_endswith_ATCG = hyb_record.has_prop('seq_suffix', 'ATCG')
                
            Record Type Properties::

                # hyb_record = hybkit.HybRecord(id, seq....)
                has_seg_types = hyb_record.has_prop('has_seg_types')  # -> False
                hyb_record.find_types()
                has_seg_types = hyb_record.has_prop('has_seg_types')  # -> True
                # Requires Type Analysis
                is_5p_mrna = hyb_record.has_prop('seg1_type', 'mRNA') 
                has_mRNA = hyb_record.has_prop('seg_type_contains', 'mRNA')

            miRNA Properties::

                # hyb_record = hybkit.HybRecord(id, seq....)
                # hyb_record.find_types()
                mirna_analyzed = hyb_record.has_prop('has_mirna_seg')  # -> False
                hyb_record.eval_mirna()
                mirna_analyzed = hyb_record.has_prop('has_mirna_seg')  # -> True
                # Requires mirna evaluation
                has_mirna = hyb_record.has_prop('has_mirna')  # Requires miRNA Analysis
                has_5p_mirna = hyb_record.has_prop('5p_mirna')
                
            Target Region Properties::
                # hyb_record = hybkit.HybRecord(id, seq....)
                # hyb_record.find_types()
                # hyb_record.eval_mirna()
                targets_analyzed = hyb_record.has_prop('has_target_reg')  # -> False
                hyb_record.target_region_eval()
                targets_analyzed = hyb_record.has_prop('has_target_reg')  # -> True
                has_coding_target = hyb_record.has_prop('target_coding') 
               
        """

        if prop not in self._HAS_PROPERTIES_SET:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available proprties are:\n' + ', '.join(self.HAS_PROPERTIES)
            print(message)
            raise Exception(message)

        # Check if a substring compares to a desired property string.
        if prop in self.STR_PROPERTIES:
            if not prop_compare:
                message = 'Property: %s  requires a comparison string. ' % prop
                message += 'Please provide an argument to prop_compare.'

            prop_split = prop.split('_')
            assert len(prop_split) in {2, 3, 4}
            check_attr = '_'.join(prop_split[:-1])
            check_type = prop_split[-1]

            check_info = None
            multi_check_type = None
            if check_attr in {'id', 'seq'}:
                check_info = getattr(self, check_attr)
            elif check_attr == 'seg1':
                check_info = self.seg1_props['ref_name']
            elif check_attr == 'seg2':
                check_info = self.seg2_props['ref_name']
            elif check_attr in {'any_seg', 'all_seg'}:
                check_info = {self.seg1_props['ref_name'],
                              self.seg2_props['ref_name']}
                multi_check = check_attr.split('_')[0]
            elif check_attr == 'seg1_type':
                check_info = self.get_seg1_type()
                self._ensure_set('eval_types')
            elif check_attr == 'seg2_type':
                check_info = self.get_seg2_type()
                self._ensure_set('eval_types')
            elif check_attr in {'any_seg_type', 'all_seg_type'}:
                self._ensure_set('eval_types')
                check_info = (self.get_seg1_type(), self.get_seg2_type())
                multi_check = check_attr.split('_')[0]
            else:
                raise Exception('Unknown Field')

            # Wrap single check_info value in a list, if not already.
            if not multi_check:
                check_info = {check_info}

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
            elif check_type == 'is':
                 ret_val = bool(prop_compare in check_info)
                # ret_val = any((prop_copmpare == val) for val in check_info)
            else:
                raise Exception(prop)

        # Check mirna-specific properties (requires mirna-evaluation)
        elif prop in self.MIRNA_PROPERTIES:
            self._ensure_set('eval_mirna')
            if prop == 'has_mirna':
                ret_val = self.flags['miRNA_seg'] in {'5p', '3p', 'B'}
            elif prop == 'no_mirna':
                ret_val = self.flags['miRNA_seg'] not in {'5p', '3p', 'B'}
            elif prop == 'mirna_dimer':
                ret_val = self.flags['miRNA_seg'] == 'B'
            elif prop == 'mirna_not_dimer':
                ret_val = self.flags['miRNA_seg'] in {'5p', '3p'}
            elif prop == '5p_mirna':
                ret_val = self.flags['miRNA_seg'] in {'5p', 'B'}
            elif prop == '3p_mirna':
                ret_val = self.flags['miRNA_seg'] in {'3p', 'B'}
            else:
                raise Exception(prop)

        elif prop in self.TARGET_PROPERTIES:
            self._ensure_set('eval_target')
            if prop == 'has_target':
                raise NotImplementedError()
            elif prop == 'target_none':
                ret_val = (self._get_flag('target_reg') == 'N')
            elif prop == 'target_unknown':
                ret_val = (self._get_flag('target_reg') == 'U')
            elif prop == 'target_ncrna':
                ret_val = (self._get_flag('target_reg') == 'NON')
            elif prop == 'target_5p_utr':
                ret_val = (self._get_flag('target_reg') == '5pUTR')
            elif prop == 'target_coding':
                ret_val = (self._get_flag('target_reg') == 'C')
            elif prop == 'target_3p_utr':
                ret_val = (self._get_flag('target_reg') == '3pUTR')
        return ret_val         


    # HybRecord : Public Methods : Record Parsing
    def to_line(self, newline=False, sep='\t'):
        """
        Return a hyb-format string representation of the Hyb record.

        Args:
            newline (bool, optional): If True, end the returned string with a newline.
            sep (str, optional): Default: "\\\\t", Provide a different separator 
                between fields.
        """
        line_items = []
        for item_key in self.HYBRID_COLUMNS:
            line_items.append(getattr(self, item_key, '.'))
        for seg_dict in [self.seg1_props, self.seg2_props]:
            for item_key in self.SEGMENT_COLUMNS:
                if item_key in seg_dict and seg_dict[item_key] is not None:
                    line_items.append(seg_dict[item_key])
                else:
                    line_items.append('.')

        flag_string = self._make_flag_string()

        if flag_string:
            line_items.append(flag_string)

        ret_string = sep.join((str(x) for x in line_items))
        if newline:
            ret_string += '\n'
        return ret_string

    # HybRecord : Public Methods : Record Parsing
    def to_csv(self, newline=False):
        """
        Return a comma-separated hyb-format string representation of the Hyb record.

        Args:
            newline (bool, optional): If True, end the returned string with a newline.
        """
        return self.to_line(newline, sep=',')

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_record(self, mode='hybrid', annotate=True):
        """
        Return nucleotide sequence as BioPython SeqRecord object.

        Args:
            mode (str, optional): Determines which sequence component to return. 
                                  Undefined options will return an error. Options are:
                                  'hybrid': Entire hybrid sequence (default); 
                                  'seg1': Sequence 1 (if defined); 
                                  'seg2': Sequence 2 (if defined);
                                  'miRNA': miRNA sequence of miRNA/target pair (if defined); 
                                  'target': Target sequence of miRNA/target pair (if defined);
            annotate (bool, optional): Add name of components to fasta sequence identifier
                                       if present.
        """
        if Bio is None:
            message = 'Please Install BioPython Package and ensure it can be imported.'
            print(message)
            raise bio_module_error
            
        allowed_modes = ['hybrid', 'seg1', 'seg2', 'mirna', 'target']
        allowed_modes_set = set(allowed_modes)
        if mode not in allowed_modes:
            message = 'Mode %s not allowed for parsing as fasta record.\n' % mode
            message += '    Allowed Modes: ' + ', '.join(allowed_modes)
            print(message)
            raise Exception(message)

        fasta_description = ''
        if mode == 'hybrid':
            fasta_id = self.seq_id
            fastq_seq = self.seq
            if annotate:
                if self.seg1_props['ref_name'] is not None:
                    fasta_description += ' ' + self._format_seg_props_line(self.seg1_props)
                if self.seg2_props['ref_name'] is not None:
                    fasta_description += ' ' + self._format_seg_props_line(self.seg2_props)
                fasta_description = fasta_description.lstrip()

        if mode in {'mirna', 'target'}:
            if not self._is_set('eval_mirna'):
                message = 'eval_mirna must be performed before miRNA/target fasta output'
                print(message)
                raise Exception(message)
            elif self.has_prop('no_mirna'):
                message = 'miRNA / target cannot be output as fasta because record ' + str(self)
                message += 'does not have a miRNA.'
                print(message)
                raise Exception(message)
            elif self.has_prop('mirna_dimer'):
                message = 'miRNA / target cannot be output as fasta because record ' + str(self)
                message += 'has a miRNA dimer.'
                print(message)
                raise Exception(message)

            if annotate:
                fasta_description += mode + '--'

            if self.has_prop('5p_mirna'):
                if mode == 'mirna':
                   mode = 'seg1'
                else:
                   mode = 'seg2'
            elif self.has_prop('3p_mirna'):
                if mode == 'mirna':
                    mode = 'seg2'
                else:
                    mode = 'seg1'
         
        if mode in {'seg1', 'seg2'}:
            seq_props = getattr(self, mode+'_props')
            if any (seq_props[v] == None for v in ['read_start', 'read_end']):
                message = 'Record %s, Segment: %s is missing one of ' % (str(self), mode)
                message += 'read_start/read_end and cannot be converted to FASTA.'
                print(message)
                raise Exception(message)
            read_start, read_end = seq_props['read_start'], seq_props['read_end']
            fasta_name = self.seq_id + ':%i-%i' % (read_start, read_end)
            fasta_seq = self._get_seg_seq(seg_props)
            if annotate:
                fasta_description += self._format_seg_props_line(seq_props)

        fasta_record = SeqRecord(Seq(fasta_seq), 
                                 id=fasta_name,
                                 description=fasta_description,
                                )

        return fasta_record

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_str(self, mode='hybrid', annotate=True):
        """
        Return nucleotide sequence as a fasta string.

        Args:
            mode (str, optional): Determines which sequence component to return. 
                                  Undefined options will return an error. Options are:
                                  'hybrid': Entire hybrid sequence (default); 
                                  'seg1': Sequence 1 (if defined); 
                                  'seg2': Sequence 2 (if defined);
                                  'miRNA': miRNA sequence of miRNA/target pair (if defined); 
                                  'target': Target sequence of miRNA/target pair (if defined);
            annotate (bool, optional): Add name of components to fasta sequence identifier
                                       if present.
        """
        return self.to_fasta_record(mode=mode, annotate=annotate).format('fasta')

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

    # HybRecord : Public Classmethods : Record Construction
    @classmethod
    def from_line(cls, line, hybformat_id=False, hybformat_ref=False):
        """
        Construct a HybRecord instance from a line in ".hyb" format.

        The Hyb Software Package contains further information in the "id" field of the
        line that can be used to infer read counts represented by the hyb record.
        Additionally, the Hyb Software Package also utilizes a database by default that contains 
        further information in the names of each respective reference sequence.

        Args:
            line (str): Hyb-format line containing record information.
            hybformat_id (bool, optional): Read count information from identifier in
                "<id>_<count>" format. (Default: False)
            hybformat_ref (bool, optional): Read additional record information from 
                identifier in "<gene_id>_<transcript_id>_<gene_name>_<seg_type>" format.
                (Default: False)

        Returns:
            :class:`HybRecord` instance containing record information.
        """

        line_items = line.strip().split('\t')
        # print(line_items)
        hyb_id = line_items[0]
        seq = line_items[1]
        energy = line_items[2]
        seg1_props = {}
        seg1_props['ref_name'] = line_items[3]
        seg1_props['read_start'] = line_items[4]
        seg1_props['read_end'] = line_items[5]
        seg1_props['ref_start'] = line_items[6]
        seg1_props['ref_end'] = line_items[7]
        seg1_props['score'] = line_items[8]
        seg2_props = {}
        seg2_props['ref_name'] = line_items[9]
        seg2_props['read_start'] = line_items[10]
        seg2_props['read_end'] = line_items[11]
        seg2_props['ref_start'] = line_items[12]
        seg2_props['ref_end'] = line_items[13]
        seg2_props['score'] = line_items[14]
        flags = {}
        if len(line_items) > 15:
            flags = cls._read_flags(line_items[15])

        if hybformat_id:
            read_id, read_count = cls._parse_hybformat_id(hyb_id)
            if 'read_count' not in flags:
                flags['read_count'] = read_count
                 
        if hybformat_ref:
            for i, seg_props in enumerate([seg1_props, seg2_props], start=1):
                ref = seg_props['ref']
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
              
        return_obj = cls(hyb_id, seq, energy, seg1_props, seg2_props, flags)
        return return_obj

    # HybRecord : Private Constants
    #: Properties for the ".is_set()" method.
    IS_SET_PROPERTIES = {
        'energy', 'fold_record', 'eval_types', 'eval_mirna', 'eval_target'
    }
    #: String-comparison properties for the ".has_prop()" method.
    STR_PROPERTIES = [
        'id_is', 'id_prefix', 'id_suffix', 'id_contains',
        'seq_is', 'seq_prefix', 'seq_suffix', 'seq_contains',
        'seg1_is', 'seg1_prefix', 'seg1_suffix', 'seg1_contains',
        'seg2_is', 'seg2_prefix', 'seg2_suffix', 'seg2_contains',
        'any_seg_is', 'any_seg_prefix', 'any_seg_suffix', 'any_seg_contains',
        'all_seg_is', 'all_seg_prefix', 'all_seg_suffix', 'all_seg_contains',
        'seg1_type', 'seg1_type_prefix', 'seg1_type_suffix', 'seg1_type_contains',
        'seg2_type', 'seg2_type_prefix', 'seg2_type_suffix', 'seg2_type_contains',
        'any_seg_type_is', 'any_seg_type_prefix', 'any_seg_type_suffix', 'any_seg_type_contains',
        'all_seg_type_is', 'all_seg_type_prefix', 'all_seg_type_suffix', 'all_seg_type_contains',
    ]
    #: miRNA-evaluation-related properties for the ".has_prop()" method.
    MIRNA_PROPERTIES = [
        'has_mirna', 'no_mirna', 'mirna_dimer', 'mirna_not_dimer',
        '3p_mirna', '5p_mirna',
    ]
    #: Target-evaluation-related properties for the ".has_prop()" method.
    TARGET_PROPERTIES = [
        'target_none', 'target_unknown', 'target_ncrna', 
        'target_5p_utr', 'target_3p_utr', 'target_coding', 
    ]

    #: All allowed properties for the ".has_prop()" method.
    HAS_PROPERTIES = STR_PROPERTIES + MIRNA_PROPERTIES + TARGET_PROPERTIES

    _HAS_PROPERTIES_SET = set(HAS_PROPERTIES)

    # HybRecord : Private Methods : Initialization
    def _post_init_tasks(self):
        # Stub for subclassing
        pass

    # HybRecord : Private Methods : Record Parsing
    def _format_seg_props(self, seg_props, prefix='', suffix='', indent_str=''):
        # Returns a formatted string of the segement properties
        ret_string = prefix
        ret_string += indent_str + 'Map Reference:  %s\n' % seg_props['ref_name']
        ret_string += indent_str + 'Read Start Pos: %s\n' % seg_props['read_start']
        ret_string += indent_str + 'Read End Pos:   %s\n' % seg_props['read_end']
        ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_props['ref_start']
        ret_string += indent_str + 'Map End Pos:    %s\n' % seg_props['ref_end']
        ret_string += indent_str + 'Map Score:      %s\n' % seg_props['score']
        ret_string += suffix
        return ret_string

    # HybRecord : Private Methods : Record Parsing
    def _format_seg_props_line(self, seg_props, prefix='', suffix=''):
        # Returns a single-line formatted string of the segement properties
        ret_string = prefix
        if seg_props['ref_name'] is not None:
            ret_string += seg_props['ref_name']
            read_vals = (str(seg_props[key]) if seg_props[key] is not None 
                         else '??' for key in ['read_start', 'read_end'])
            if any((read_vals[i] != '??' for i in [0, 1])):
                ret_string += ':%s-%s' % read_vals 
            ref_vals = (str(seg_props[key]) if seg_props[key] is not None 
                         else '??' for key in ['ref_start', 'ref_end'])
            if any((ref_vals[i] != '??' for i in [0, 1])):
                ret_string += '(%s-%s)' % ref_vals 
        ret_string += suffix
        return ret_string

    # HybRecord : Private Methods : Segment Parsing
    def _get_seg_seq(self, seg_props):
        if any (seg_props[v] == None for v in ['read_start', 'read_end']):
            message = 'Segement subsequence cannot be obtained for '
            message += 'Record %s, Segment %s.\n' % (str(self), seg_props['ref_name'])
            message += 'Record segment is missing one of read_start/read_end.'
            print(message)
            raise Exception(message)
        read_start, read_end = seg_props['read_start'], seg_props['read_end']
        return self.seq[read_start-1:read_end]

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
    def _make_flag_string(self):
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
        if self._flagset is None:
            self._flagset = set(self.ALL_FLAGS + list(self.settings.custom_flags))
        return_list = []
        for flag in self.ALL_FLAGS + self.settings['custom_flags']:
            if flag in self.flags:
                return_list.append(flag)
        for flag in self.flags:
            if flag not in self._flagset:
                return_list.append(flag)
        return return_list

    # HybRecord : Private Methods : flags
    def _make_flags_dict(self, flag_obj, allow_undefined_flags=None):
        #  allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        #  If argument is provided to the method, it overrides default behavior.
        #  Otherwise, the method falls back to the object-defaults.
        if allow_undefined_flags is None:
            allow_undefined_flags = self.settings['allow_undefined_flags']

        if self._flagset is None:
            self._flagset = set(self.ALL_FLAGS + list(self.settings['custom_flags']))

        if not isinstance(flag_obj, dict):
            message = '"flag_obj" argument must be a dict obj. Defined keys are:'
            message += (self.ALL_FLAGS + self._custom_flags).join(', ')
            print(message)
            raise Exception(message)

        if not allow_undefined_flags:
            for flag in flag_obj:
                if flag not in self._flagset:
                    message = 'Flag "%s" is not defined. Please check flag key' % flag
                    message += ' or run with: "allow_undefined_flags=True"\n'
                    message += 'Defined Flags are: '
                    message += ', '.join(self.ALL_FLAGS + self._custom_flags)
                    print(message)
                    raise Exception(message)
        return flag_obj

    # HybRecord : Private Methods : seg_props
    def _make_seg_props_dict(self, seg_props_obj={}):
        # Create a dictionary with mapping entries, ensuring each read data point is either a
        #   placeholder or is of the correct data type.
        return_dict = {}
        segment_column_types = {'ref_name': str,
                                'read_start': int,
                                'read_end': int,
                                'ref_start': int,
                                'ref_end': int,
                                'score': str,
                                }
        for column in self.SEGMENT_COLUMNS:
            if column in seg_props_obj:
                if seg_props_obj[column] == self.settings['hyb_placeholder']:
                    return_dict[column] = None
                else:
                    column_type = segment_column_types[column]
                    try:
                        return_dict[column] = column_type(seg_props_obj[column])
                    except TypeError:
                        message = 'Error in setting segN_props dict.\n'
                        message += 'Entry "%s" for column: %s' % (seg_props_obj[column], column)
                        message += 'could not be converted to type: %s' % str(column_type)
                        print(message)
                        raise
            else:
                return_dict[column] = None
        return return_dict

    # HybRecord : Private Methods : properties
    def _ensure_set(self, prop):
        if not self.is_set(prop):
            message = 'Problem with HybRecord instance: %s\n' % str(self)
            message += 'Method requries set attribute/evaluation: "%s" before use.' % prop
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

        if cls._flagset is None:
            cls._flagset = set(cls.ALL_FLAGS + list(cls.settings['custom_flags']))
        flag_string = flag_string.rstrip()
        flag_string = flag_string.rstrip(';')
        flag_pairs = [flag_pair.split('=') for flag_pair in flag_string.split(';')]
        flags = {}
        for flag_key, flag_value in flag_pairs:
            if not allow_undefined_flags and flag_key not in cls._flagset:
                message = 'Problem: Unidefined Flag: %s\n' % flag_key
                message += 'Defined Flags: '
                message += ', '.join(self.ALL_FLAGS + self._custom_flags)
                print(message)
                raise Exception(message)
            flags[flag_key] = flag_value
        return flags


class HybFile(object):
    """
    File-Object wrapper that provides abiltity to return file lines as HybRecord entries.
    """

    #: Class-level settings. See :obj:`settings.HybFile_settings` for descriptions.
    settings = hybkit.settings.HybFile_settings

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
    # Check if provided argument ("record") is an instance of HybRecord.
    def _ensure_HybRecord(self, record):
        if not isinstance(record, HybRecord):
            message = 'Item: "%s" is not a HybRecord object.' % record
            print(message)
            raise Exception(message)


class FoldRecord(object):
    """
    Class for storing secondary structure (folding) information for a nucleotide sequence.
    
    This class supports the following file types:
    (Data courtesy of Gay et al. [see :ref:`References`])

    .. _vienna_file_format:

    * | The .vienna file format used by the RNAStructure package (see :ref:`References`):

      Example:
          ::

              34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-...
              TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
              .....((((((.((((((......)))))).))))))   (-11.1)

    * | The .ct file format utilized by the UNAFold Software Package:

      Example:
          ::

              41	dG = -8	dH = -93.9	seq1_name-seq2_name
              1	A	0	2	0	1	0	0
              2	G	1	3	0	2	0	0
              ...
              40	G	39	41	11	17	39	41
              41	T	40	0	10	18	40	0


    A minimum amount of data necessary for a FoldRecord object is a sequence identifier,
    a genomic sequence, and its fold representaiton.

    Args:
        id (str): Identifier for record
        seq (str): Nucleotide sequence of record.
        fold (str): Fold representation of record.
        energy (str or float, optional): Energy of folding for record.

    .. _FoldRecord-Attributes:

    Attributes:
        id (str): Sequence Identifier (often seg1name-seg2name)
        seq (str): Genomic Sequence
        fold (str): Dot-bracket Fold Representation, '(', '.', and ')' characters
        energy (float or None): Predicted energy of folding
    """

    # FoldRecord : Class-Level Constants
    #: Class-level settings. See :obj:`settings.FoldRecord_settings` for descriptions.
    settings = hybkit.settings.FoldRecord_settings

    # FoldRecord : Public Methods : Initialization
    def __init__(self, id, seq, fold, energy):
        self.id = id                 # Sequence Identifier (often seg1name-seg2name)
        self.seq = seq               # Genomic Sequence
        self.fold = fold             # Fold Representation, str of '(', '.', and ')' characters
        self.energy = float(energy)  # Predicted energy of folding

    # FoldRecord : Public Methods : Parsing : Vienna
    def to_vienna_lines(self, newline=False):
        """
        Return a list of lines for the record in vienna format. 

        See (:ref:`Vienna File Format <vienna_file_format>`).

        Args:
            newline (bool, optional): If True, add newline character to the end of each
                returned line. (Default: False)
        """
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
        """
        Return a 3-line string for the record in vienna format. 

        See (:ref:`Vienna File Format <vienna_file_format>`).

        Args:
            newline (bool, optional): If True, terminate the returned string with a newline
                character. (Default: False)
        """
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_vienna_lines(newline=False)) + suffix)

    # FoldRecord : Public MagicMethods : Comparison
    def __eq__(self, other):
        """Return True if two records have matching sequences and identifiers."""
        return (self.id == other.id and self.seq == other.seq)

    # FoldRecord : Public MagicMethods : Evaluation
    def __hash__(self):
        """Return a hash of the record ".id" attribute"""       
        return hash(self.id)

    # FoldRecord : Public MagicMethods : Evaluation
    def __bool__(self):
        """Return True wherever the class is defined."""
        return True

    # FoldRecord : Public MagicMethods : Evaluation
    def __len__(self):
        """Return the length of the genomic sequence"""
        return len(self.seq)

    # FoldRecord : Public MagicMethods : Printing
    def __str__(self):
        """Print the identifier of the record."""
        return '<FoldRecord ID: %s>' % self.id

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_lines(cls, 
                          record_lines, 
                          skip_bad_fold_records=None, 
                          warn_bad_fold_records=None
                         ):
        """
        Construct instance from a list of 3 strings of Vienna-format lines.

        Args:
            record_lines (str or tuple): Iterable of 3 strings corresponding to lines of a
                vienna-format record.
            skip_bad_fold_records (bool, optional): If True, return None when parsing 
                badly-formatted entries instead of raising an error. 
                If None (default), uses setting in 
                :attr:`settings['skip_bad_fold_records'] <settings>`.
            warn_bad_fold_records (bool, optional): If True, print a warning message when 
                attempting to parse badly-formatted entries.
                If None (default), uses setting in 
                :attr:`settings['warn_bad_fold_records'] <settings>`.
        """

        if skip_bad_fold_records is None:
            skip_bad_fold_records = cls.settings['skip_bad_fold_records']
        if warn_bad_fold_records is None:
            warn_bad_fold_records = cls.settings['warn_bad_fold_records']
        fail_ret_val = (None, ''.join(record_lines))


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

        #If no fold was created, potentially due to low-complexity sequence
        if "(99" in line_3:
           if skip_bad_fold_records:
                if warn_bad_fold_records:
                    message = 'WARNING: Improper Vienna: No Fold (Energy = 99*.*)'
                return fail_ret_val
           else:
                message = 'ERROR: Improper Vienna: No Fold (Energy = 99*.*)'
                print(message)
                raise Exception(message)

        if len(line_3_split) != 2:
            message = 'Provided Vienna Record Line 3:\n'
            message += line_3.rstrip() + '\n'
            message += str(line_3_split) + '\n'
            message += '\n  ... does not have required ".(.).<tab>(-1.23)" format.'
            print(message)
            raise Exception(message)
        fold = line_3_split[0]
        energy = float(line_3_split[1].strip('()'))

        return_obj = cls(rec_id, seq, fold, energy,
                         seg1_fold_info=seg1_fold_info,
                         seg2_fold_info=seg2_fold_info)
        #print(return_obj.to_viennad_string())
        return return_obj

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_string(cls, record_string, hybformat_file=False):
        """
        Construct instance from a string representing 3 Vienna-format lines.

        The Hyb Software Package contains further information in the "name" field of the
        vienna record that can be used to infer further information about the fold divisions.

        Args:
            record_lines (str or tuple): Iterable of 3 strings corresponding to lines of a
                vienna-format record.
            hybformat_file (bool, optional): If True, extra information stored in the 
                record identifier by Hyb will be parsed.
        """
        lines = record_string.strip().split('\n')[0:3]
        return cls.from_vienna_lines(lines)

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_lines(cls, record_lines):
        """
        Create a FoldRecord entry from a list of an arbitrary number of strings
        corresponding to lines in the ".ct" file format.
        """
        header_line = record_lines[0].strip()
        if not any((x in header_line for x in ['dG', 'Energy', 'ENERGY'])):
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

        for i, line in enumerate(record_lines[1:], 1):
            line_split = line.strip().split('\t')
            if len(line_split) not in {6, 8}:
                message = 'Provided ct Record Line:\n'
                message += line.rstrip() + '\n'
                message += '\n  ... does not have 6-column or 8-column format'
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

        return_obj = cls(full_name, seq, fold, energy)
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

    # FoldRecord : Private Classmethods : Parsing : Output
    @classmethod
    def _format_seg_props(cls, seg_props, prefix='', suffix='', indent_str=''):
        raise NotImplementedError
        # Returns a formatted string of the sgement info information
        ret_string = prefix
        ret_string += indent_str + 'Map Reference:  %s\n' % seg_props['ref']
        ret_string += indent_str + 'Read Start Pos: %s\n' % seg_props['read_start']
        ret_string += indent_str + 'Read End Pos:   %s\n' % seg_props['read_end']
        ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_props['ref_start']
        ret_string += indent_str + 'Map End Pos:    %s\n' % seg_props['ref_end']
        ret_string += indent_str + 'Map Score:      %s\n' % seg_props['score']
        ret_string += suffix
        return ret_string


class FoldFile(object):
    """
    Base class for file-object wrappers that return file lines as FoldRecord objects.

    See :class:`ViennaFile` or :class:`CtFile`.

    """

    #: Class-level settings. See :attr:`settings.FoldFile_settings` for descriptions.
    settings = hybkit.settings.FoldFile_settings

    # FoldFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Stub method to be replaced by subclasses."""
        message = 'FoldFile is a base class and is not meant to be used directly.' 
        print(message)
        raise NotImplementedError(message)

    # FoldFile : Public Methods : Initialization / Closing
    def __enter__(self, *args, **kwargs):
        """Open "with" syntax."""
        return self

    # FoldFile : Public Methods : Initialization / Closing
    def __exit__(self, type, value, traceback):
        """Close "with" syntax"""
        self.close()

    # FoldFile : Public Methods : Initialization / Closing
    def __iter__(self):
        """Return an iterator."""
        return self

    # FoldFile : Public Methods : Reading
    def __next__(self):
        """Stub method to be replaced by subclasses."""
        message = 'FoldFile is a base class and is not meant to be used directly.' 
        print(message)
        raise NotImplementedError(message)

    # FoldFile : Public Methods : Reading
    def close(self):
        """Close the file handle."""
        self.fh.close()

    # FoldFile : Public Methods : Reading
    def read_record(self):
        """Return next :class:`FoldRecord` based on the appropriate file type."""
        return next(self)

    # FoldFile : Public Methods : Reading
    def read_records(self):
        """Return list of all :class:`FoldRecord` objects based on the appropriate file type."""
        records = []
        for record in self:
            records.append(record)
        return records

    # FoldFile : Public Methods : Writing
    def write_record(self, write_record):
        """
        Write a FoldRecord object as the appropriate record/file type.

        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.
        """
        self._ensure_FoldRecord(write_record)
        record_string = self._to_record_string(write_record, newline=True)
        self.fh.write(record_string)

    # FoldFile : Public Methods : Writing
    def write_records(self, write_records):
        """
        Write a sequence of FoldRecord objects as the appropraite record/file type.
        
        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.
        """
        for write_record in write_records:
            self.write_record(write_record)

    # FoldFile : Public Classmethods : Initialization
    @classmethod
    def open(cls, *args, **kwargs):
        """Return a new FoldFile object."""
        return cls(*args, **kwargs)

    # FoldFile : Private Methods
    def _post_init_tasks(self):
        """Stub for convenient subclassing"""
        pass

    # FoldFile : Private Methods
    def _ensure_FoldRecord(self, record):
        if not isinstance(record, FoldRecord):
            message = 'Item: "%s" is not a FoldRecord object.' % record
            print(message)
            raise Exception(message)

    # FoldFile : Private Methods
    def _to_record_string(self, write_record, newline):
        """Stub method to be replaced by subclasses."""
        message = 'FoldFile is a base class and is not meant to be used directly.'
        print(message)
        raise NotImplementedError(message)


class ViennaFile(FoldFile):
    """
    Vienna file wrapper that returns ".vienna" file lines as FoldRecord objects.
    """

    # ViennaFile : Public Methods : Reading
    def __next__(self):
        """Read next three lines and return output as FoldRecord object."""
        line_1 = next(self.fh)
        line_2 = next(self.fh)
        line_3 = next(self.fh)
        return FoldRecord.from_vienna_lines((line_1, line_2, line_3))

    # ViennaFile : Private Methods
    def _to_record_string(self, write_record, newline):
        """Return a :class:`Fold Record` as a Vienna-format string."""
        return write_record.to_vienna_string(newline=newline)


class CtFile(FoldFile):
    """
    Ct file wrapper that returns ".ct" file lines as FoldRecord objects.
    """

    # CtFile : Public Methods
    def __next__(self):
        """
        Return the next ct record as a :class:`FoldRecord` object.
 
        Call next(self.fh) to return the first line of the next entry.
        Determine the expected number of following lines in the entry, and read that number
        of lines further. Return lines as a FoldRecord object.
        """
        header = next(self.fh)
        record_lines = [header]
        expected_line_num = int(header.strip().split()[0])
        for i in range(expected_line_num):
            record_lines.append(next(self.fh))
        ret_record = FoldRecord.from_ct_lines(record_lines)
        return ret_record

    # CtFile : Private Methods
    def _to_record_string(self, write_record, newline):
        """Return a :class:`Fold Record` as a Ct-format string."""
        message = 'No write_record is implmeneted for ct files, as the FoldRecord '
        message += 'object does not contain the complete set of ct record information.'
        print(message)
        raise NotImplementedError(message)


class HybFoldIter(object):
    """
    Iterator for simultaneous iteration over a :class:`HybFile` and :class:`FoldFile` object.

    This class provides an iterator to iterate through a :class:`HybFile` and one of a 
    :class:`ViennaFile`, or :class:`CtFile` simultaneously.
    Each :class:`HybRecord` and :class:`FoldRecord` in the respective
    files are checked to ensure they have matching sequences.
    The obtained :class:`FoldRecord` will be set as 
    :attr:`.HybRecord.fold_record` of the returned :class:`HybRecord` object.

    Args:
        hybfile_handle (HybFile) : HybFile object for iteration
        foldfile_handle (FoldFile) : FoldFile object for iteration

    Returns:
        (:class:`HybRecord`, :class:`FoldRecord`) 
    """
    # HybFoldIter : Public Methods
    def __init__(self, hybfile_handle, foldfile_handle):
        """Please see :class:`HybFoldIter` for initialization information."""
        self.hybfile_handle = hybfile_handle
        self.foldfile_handle = foldfile_handle
        self.counter = 0

    # HybFoldIter : Public Methods
    def __iter__(self):
        """Return an iterator object."""
        return self

    # HybFoldIter : Public Methods
    def __next__(self):
        """Read and return (class:`HybRecord`, :class:`FoldRecord`)"""
        self.counter += 1
        next_hyb_record = next(self.hybfile_handle)
        next_fold_record = next(self.foldfile_handle)
        try:
            next_hyb_record.set_fold_record(next_fold_record)
        except:
            print('For %s counter iteration: %i ...' % (str(self), self.counter))
            raise
        return (next_hyb_record, next_fold_record)


# Import the remainder of hybkit code to connect.
import hybkit.analysis
import hybkit.plot
import hybkit.util

