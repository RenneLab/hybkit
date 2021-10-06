#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
This module contains classes and methods for reading, writing, and manipulating data 
in the ".hyb" genomic sequence format ([Travis2014]_).

This is primarily based on three classes for storage of 
chimeric sequence information and associated fold-information:

+----------------------------+---------------------------------------------------------------+
| :class:`HybRecord`         | Class for storage of hybrid sequence records                  |
+----------------------------+---------------------------------------------------------------+
| :class:`FoldRecord`        | Class for storage of predicted RNA                            |
|                            | secondary structure information for chimeric sequence reads   |
+----------------------------+---------------------------------------------------------------+
| :class:`DynamicFoldRecord` | Class for storage of predicted RNA                            |
|                            | secondary structure information for sequence constructed from |
|                            | aligned portions of chimeric sequence reads                   |
+----------------------------+---------------------------------------------------------------+
    
It also includes classes for reading, writing, and iterating over files containing that 
information:

+-------------------------+------------------------------------------------------------------+
| :class:`HybFile`        | Class for reading and writing ".hyb"-format files                |
|                         | [Travis2014]_ containing chimeric RNA sequence information       |
|                         | as :class:`HybRecord` objects                                    |
+-------------------------+------------------------------------------------------------------+
| :class:`ViennaFile`     | Class for reading and writing Vienna (.vienna)-format files      |
|                         | [ViennaFormat]_ containing RNA secondary structure information   |
|                         | dot-bracket format as :class:`FoldRecord` objects                |
+-------------------------+------------------------------------------------------------------+
| :class:`CtFile`         | Class for reading Connectivity Table (.ct)-format files          |
|                         | [CTFormat]_ containing predicted RNA secondary-structure         |
|                         | information as used by UNAFold_ as                               |
|                         | :class:`FoldRecord` objects                                      |
+-------------------------+------------------------------------------------------------------+
| :class:`HybFoldIter`    | Class for concurrent iteration over a :class:`HybFile` and       |
|                         | one of a :class:`ViennaFile` or :class:`CtFile`                  |
+-------------------------+------------------------------------------------------------------+

Todo:
    Add Hybrecord.to_csv_header()

"""

import os
import sys
import io
import types
import csv
import copy
from collections import Counter
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
import hybkit.__about__
# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__,\
                             __email__, __license__, __maintainer__, __status__, __version__

class HybRecord(object):
    """
    Class for storing and analyzing chimeric (hybrid) RNA-seq reads in ".hyb" format.

    Hyb file (".hyb") format entries are a GFF-related file format described by [Travis2014]_.
    that contain information about a genomic sequence read identified to be a chimera by 
    anlaysis software. Each line contains 15 or 16 columns separated by tabs ("\\\\t") and provides
    annotations on each components. An example .hyb format line 
    from [Gay2018]_::
 
        2407_718\tATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC\t.\tMIMAT0000078_MirBase_miR-23a_microRNA\t1\t21\t1\t21\t0.0027\tENSG00000188229_ENST00000340384_TUBB2C_mRNA\t23\t49\t1181\t1207\t1.2e-06

    These columns are respectively described in hybkit as:

         id, seq, energy,
         seg1_ref_name, seg1_read_start, seg1_read_end, seg1_ref_start, seg1_ref_end, seg1_score, 
         seg2_ref_name, seg2_read_start, seg2_read_end, seg2_ref_start, seg2_ref_end, seg2_score, 
         flag1=val1;flag2=val2;flag3=val3..."

    The preferred method for reading hyb records from lines is with 
    the :func:`HybRecord.from_line` constructor::

        # line = "2407_718\tATC..."
        hyb_record = hybkit.HybRecord.from_line(line)

    This constructor parses hyb files using the 
    :class:`HybFile` class.
    For example, to print all hybrid identifiers in a hyb file::

        with hybkit.HybFile('path/to/file.hyb', 'r') as hyb_file:
            for hyb_record in hyb_file:
                print(hyb_record.id)

    HybRecord objects can also be constructed directly. A minimum amount of data necessary 
    for a HybRecord object is the genomic sequence and its corresponding identifier. 
    
    Examples:
        ::

            hyb_record_1 = hybkit.HybRecord('1_100', 'ACTG')
            hyb_record_2 = hybkit.HybRecord('2_107', 'CTAG', '-7.3')

    Details about segments are provided via python dictionaries with 
    :attr:`keys <HybRecord.SEGMENT_COLUMNS>`
    specific to each segment. Data can be provided either as strings or 
    as floats/integers (where appropriate).
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
        seg1_props (dict, optional): Properties of segment 1 of the record, containing possible
            :attr:`segment column keys <HybRecord.SEGMENT_COLUMNS>`: 
            ('ref_name', 'read_start', 'read_end', 'ref_start', 'ref_end', 'score')
        seg2_props (dict, optional): Properties of segment 2 of the record, containing possible:
            :attr:`segment column keys <HybRecord.SEGMENT_COLUMNS>`: 
            ('ref_name', 'read_start', 'read_end', 'ref_start', 'ref_end', 'score')
        flags (dict, optional): Dict with keys of flags for the record and their associated values.
            By default flags must be defined in :attr:`ALL_FLAGS` but custom 
            flags can be supplied in :attr:`settings['custom_flags'] <HybRecord.settings>`.
            This setting can also be disabled by setting 'allow_undefined_flags' 
            to :obj:`True` in :attr:`HybRecord.settings`.
        fold_record (FoldRecord, optional): Set the record's :attr:`fold_record` attribute 
            as the provided FoldRecord object using :meth:`set_fold_record` on initializtaion.

    .. _HybRecord-Attributes:

    Attributes:
        id (str): Identifier for the hyb record (Hyb format: "<read-num>_<read-count>")

        seq (str): Nucleotide sequence of the hyb record
        energy (str or None): Predicted energy of folding
        seg1_props (dict): Information on chimeric segment 1, contains 
            :attr:`segment column keys <HybRecord.SEGMENT_COLUMNS>`: 
            'ref_name' (str), 'read_start' (int), 'read_end' (int), 'ref_start' (int), 
            'ref_end' (int), and 'score' (float).
        seg2_props (dict): Information on segment 2, contains
            :attr:`segment column keys <HybRecord.SEGMENT_COLUMNS>`: 
            'ref_name' (str), 'read_start' (int), 'read_end' (int), 'ref_start' (int), 
            'ref_end' (int), and 'score' (float).
        flags (dict): Dict of flags with possible 
            :attr:`flag keys <HybRecord.ALL_FLAGS>` and values as defined in
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
    #   “feature1=value1;feature2=value2;..."
    #   Flags utilized in the Hyb software package, see specification.
    _HYB_FLAGS = ['count_total',
                  'count_last_clustering',
                  'two_way_merged',
                  'seq_IDs_in_cluster',
                ]
    # Additional flag specifications utilized by hybkit, see specificiation.
    _HYBKIT_FLAGS = ['read_count',
                     'orient',
                     'seg1_type',
                     'seg2_type',
                     'seg1_det',
                     'seg2_det',
                     'miRNA_seg',
                     'target_reg',
                     'ext',
                     'dataset',
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
        """Set the value of record ``flag_key`` to ``flag_val``.

        Args:
            flag_key (str): Key for flag to set.
            flag_val : Value for flag to set.
            allow_undefined_flags (bool or None, optional): Allow inclusion of flags not  
                defined in :attr:`ALL_FLAGS` or in :attr:`settings['custom_flags'] <HybRecord.settings>`.
                If None (default), uses setting in 
                :attr:`settings['allow_undefined_flags'] <HybRecord.settings>`.
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
        """Return the :ref:`seg1_type <seg1_type>` flag if defined, or return None.

        Args:
            require (bool, optional): If ``True``, raise an error if seg1_type is not defined.
        """
        if require:
            return self._get_flag('seg1_type')
        else:
            return self._get_flag_or_none('seg1_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg2_type(self, require=False):
        """Return the :ref:`seg2_type <seg2_type>` flag if defined, or return None.

        Args:
            require (bool, optional): If ``True``, raise an error if seg2_type is not defined.
        """
        if require:
            return self._get_flag('seg2_type')
        else:
            return self._get_flag_or_none('seg2_type')

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg_types(self, require=False):
        """Return :ref:`seg1_type <seg1_type>`, :ref:`seg2_type <seg2_type>` flags, or None.

        Args:
            require (bool, optional): If ``True``, raise an error if either flag is not defined.
        """
        if require:
            return (self._get_flag('seg1_type'), self._get_flag('seg2_type'))
        else:
            return (self._get_flag_or_none('seg1_type'), self._get_flag_or_none('seg2_type'))

    # HybRecord : Public Methods : Flag_Info : get_read_count
    def get_read_count(self, require=False):
        """Return the :ref:`read_count <read_count>` flag if defined, otherwise return None.

        Args:
            require (bool, optional): If ``True``, raise an error if the "read_count" flag 
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
        """Return :ref:`count_total <count_total>` flag if defined, or return 1 (this record).

        Args:
            require (bool, optional): If ``True``, raise an error if the "count_total" flag 
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
            count_mode (str) : | Mode for returned count, one of : {'read', 'record'}
                               | ``read``   : Require the 'read_count' flag to be defined.
                               | ``record`` : Return '1' if the 'count_total' flag is not defined.
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

    # HybRecord : Public Methods : Flag_Info : eval_type
    def eval_types(self, allow_unknown=None, check_complete=None):
        """Find the types of each segment using the the :class:`TypeFinder` class.

        This method provides :attr:`seg1_props` and :attr:`seg2_props` 
        to the :class:`TypeFinder` class, linked as attribute :attr:`HybRecord.TypeFinder`.
        This uses the method: :func:`TypeFinder.method`
        set by :func:`TypeFinder.set_method` or :func:`TypeFinder.set_custom_method` to set the
        :ref:`seg1_type <seg1_type>`, :ref:`seg2_type <seg2_type>` flags if not already set. 

        To use a type-finding method other than the default, 
        prepare the :class:`TypeFinder` class by 
        preparing and setting :attr:`TypeFinder.params` and using :func:`TypeFinder.set_method`. 

        Args:
            allow_unknown (bool, optional): If ``True``, allow segment types that cannot be
                identified and set them as "unknown". Otherwise raise an error.
                If None (default), uses setting in 
                :attr:`settings['allow_unknown_seg_types'] <HybRecord.settings>`.
            check_complete (bool, optional): If ``True``, check every possibility for the 
                type of a given segment (where applicable), instead of 
                stopping after finding the first type.
                If None (default), uses setting in 
                :attr:`settings['check_complete_seg_types'] <HybRecord.settings>`.
        """

        # If types already set, skip.
        if self.is_set('eval_types'):
            return        

        if allow_unknown is None:
            allow_unknown = self.settings['allow_unknown_seg_types']
        
        if check_complete is None:
            check_complete = self.settings['check_complete_seg_types']

        if self.TypeFinder.find is None:
            self.TypeFinder.set_method('hybformat')

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
        fold_record.ensure_matches_hyb_record(self)
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
                from :attr:`settings['mirna_types'] <HybRecord.settings>`.
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
        Unless ``allow_mirna_dimers`` is ``True``, 
        this method requires record to contain a non-dimer miRNA, 
        otherwise an error will be raised.

        Args:
            detail (str): | Type of detail to return. Options include: 
                          | ``all``             : (dict of all properties, default) 
                          | ``mirna_ref``       : Identifier for Assigned miRNA
                          | ``target_ref``      : Identifier for Assigned Target
                          | ``mirna_seg_type``  : Assigned seg_type of miRNA
                          | ``target_seg_type`` : Assigned seg_type of target
                          | ``mirna_seq``       : Annotated subsequence of miRNA
                          | ``target_seq``      : Annotated subsequence of target
                          | ``mirna_fold``      : Annotated fold substring of miRNA 
                            (requires fold_record set)
                          | ``target_fold``     : Annotated fold substring target 
                            (requires fold_record set)
            allow_mirna_dimers (bool, optional): Allow miRNA/miRNA dimers. 
                The 5p-position will be assigned as the "miRNA", 
                and the 3p-position will be assigned as the "target".

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


        ALLOWED_DETAILS = {'all', 'mirna_ref', 'target_ref', 'mirna_seg_type', 'target_seg_type', 
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

        mirna_details['mirna_ref'] = mirna_props['ref_name']
        mirna_details['target_ref'] = target_props['ref_name']
        mirna_details['mirna_seq'] = self._get_seg_seq(mirna_props)
        mirna_details['target_seq'] = self._get_seg_seq(target_props)
        if self.fold_record is not None:
            mirna_details['mirna_fold'] = self.fold_record._get_seg_fold(mirna_props, self) 
            mirna_details['target_fold'] = self.fold_record._get_seg_fold(target_props, self)
        else:
            mirna_details['mirna_fold'] = None
            mirna_details['target_fold'] = None

        if detail == 'all':
            return mirna_details
        else:
            return mirna_details[detail]

    # HybRecord : Public Methods : Record Properties
    def is_set(self, prop):
        """
        Return ``True`` if HybRecord property "prop" is set (if relevant) and is not ``None``.
        Options described in :attr:`SET_PROPS`.

        Args:
            prop (str): Property / Analysis to check
        """

        if prop not in self._SET_PROPS_SET:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available proprties are:\n' + ', '.join(self.SET_PROPS)
            print(message)
            raise Exception(message)

        if prop in {'energy', 'fold_record'}:
            ret_bool = (getattr(self, prop) is not None)
        elif prop == 'full_seg_props':
            ret_bool = (all(k in self.seg1_props and k is not None 
                            for k in self.SEGMENT_COLUMNS)
                        and all(k in self.seg1_props and k is not None 
                                for k in self.SEGMENT_COLUMNS)
                       )
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
        Return ``True`` if HybRecord has property: ``prop``. 
 
        Check property against list of allowed properties in :attr:`HAS_PROPS`.
        If query property has a string comparator, provide this in prop_compare.
        Raises an error if a prerequisite field is not set 
        (use :meth:`is_set` to check whether properties are set).
     
        Specific properties available to check are described in attributes:
        
            ======================= =============================================
            :attr:`GEN_PROPS`       General Record Properties
            :attr:`STR_PROPS`       Field String Comparison Properties
            :attr:`MIRNA_PROPS`     miRNA-Associated Record Properties
            :attr:`MIRNA_STR_PROPS` miRNA-Associated String Comparison Properties
            :attr:`TARGET_PROPS`    miRNA-Target-Associated Properties
            ======================= =============================================
 
        Args:
            prop (str):                   Property to check
            prop_compare (str, optional): Comparator to check.

        """

        if prop not in self._HAS_PROPS_SET:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available proprties are:\n' + ', '.join(self.HAS_PROPS)
            print(message)
            raise Exception(message)

        # Check if a substring compares to a desired property string.
        if prop in self._GEN_PROPS_SET:
            if prop == 'has_indels':
                self._ensure_set('full_seg_props')
                has_indels = False
                for seg_props in (self.seg1_props, self.seg2_props):
                    read_len = seg_props['read_end'] - seg_props['read_start']        
                    ref_len = seg_props['ref_end'] - seg_props['ref_start']        
                    has_indels = (read_len != ref_len)
                    if has_indels:
                        break
                ret_val = has_indels

        # Check if a substring compares to a desired property string.
        elif prop in self._ALL_STR_PROPS_SET:
            if not prop_compare:
                message = 'Property: %s  requires a comparison string. ' % prop
                message += 'Please provide an argument to prop_compare.'

            if prop in self._MIRNA_STR_PROPS_SET:
                self._ensure_set('eval_mirna')

            prop_split = prop.split('_')
            assert len(prop_split) in {2, 3, 4}
            check_attr = '_'.join(prop_split[:-1])
            check_type = prop_split[-1]

            check_info = None
            multi_check = None
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
            elif check_attr == 'mirna':
                check_info = self.mirna_detail('mirna_ref', allow_mirna_dimers=True)
            elif check_attr == 'target':
                check_info = self.mirna_detail('target_ref', allow_mirna_dimers=True)
            elif check_attr == 'mirna_seg_type':
                check_info = self.mirna_detail('mirna_seg_type', allow_mirna_dimers=True)
            elif check_attr == 'target_seg_type':
                check_info = self.mirna_detail('target_seg_type', allow_mirna_dimers=True)
            else:
                raise Exception('Unknown Field: ' + check_attr)

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
        elif prop in self._MIRNA_PROPS_SET:
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

        elif prop in self._TARGET_PROPS_SET:
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
        Return a hyb format string representation of the record.

        Args:
            newline (bool, optional): If ``True``, end the returned string with a newline
            sep (str, optional): Separator between fields (Default: "\\\\t")

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
        Return a comma-separated hyb-format string representation of the record.

        Args:
            newline (bool, optional): If ``True``, end the returned string with a newline.
        """
        return self.to_line(newline, sep=',')

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_record(self, mode='hybrid', annotate=True):
        """
        Return nucleotide sequence as BioPython SeqRecord object.

        Args:
            mode (str, optional): | Determines which sequence component to return. Options:
                                  | ``hybrid``: Entire hybrid sequence (default) 
                                  | ``seg1``: Sequence 1 (if defined) 
                                  | ``seg2``: Sequence 2 (if defined)
                                  | ``miRNA``: miRNA sequence of miRNA/target pair 
                                    (if defined, else None) 
                                  | ``target``: Target sequence of miRNA/target pair 
                                    (if defined, else None)
            annotate (bool, optional): Add name of components to fasta sequence identifier
                                       if present.
        """
        if Bio is None:
            message = 'Please install BioPython package and ensure it can be imported.'
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
            fasta_id = self.id
            fasta_seq = self.seq
            if annotate:
                if self.seg1_props['ref_name'] is not None:
                    fasta_description += ' ' + self._format_seg_props_line(self.seg1_props)
                if self.seg2_props['ref_name'] is not None:
                    fasta_description += ' ' + self._format_seg_props_line(self.seg2_props)
                fasta_description = fasta_description.lstrip()

        if mode in {'mirna', 'target'}:
            if not self.is_set('eval_mirna'):
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
            seg_props = getattr(self, mode+'_props')
            if any (seg_props[v] == None for v in ['read_start', 'read_end']):
                message = 'Record %s, Segment: %s is missing one of ' % (str(self), mode)
                message += 'read_start/read_end and cannot be converted to FASTA.'
                print(message)
                raise Exception(message)
            read_start, read_end = seg_props['read_start'], seg_props['read_end']
            fasta_id = self.id + ':%i-%i' % (read_start, read_end)
            fasta_seq = self._get_seg_seq(seg_props)
            if annotate:
                fasta_description += self._format_seg_props_line(seg_props)

        fasta_record = SeqRecord(Seq(fasta_seq), 
                                 id=fasta_id,
                                 description=fasta_description,
                                )

        return fasta_record

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_str(self, mode='hybrid', annotate=True):
        """
        Return nucleotide sequence as a fasta string.

        Args:
            mode (str, optional): | Determines which sequence component to return. Options:
                                  | ``hybrid``: Entire hybrid sequence (default) 
                                  | ``seg1``: Sequence 1 (if defined) 
                                  | ``seg2``: Sequence 2 (if defined)
                                  | ``miRNA``: miRNA sequence of miRNA/target pair 
                                    (if defined, else None) 
                                  | ``target``: Target sequence of miRNA/target pair 
                                    (if defined, else None)
            annotate (bool, optional): Add name of components to fasta sequence identifier
                                       if present.
        """
        return self.to_fasta_record(mode=mode, annotate=annotate).format('fasta')

    # HybRecord : Public MagicMethods : Comparison
    def __eq__(self, other):
        """Return ``True`` if ".id" and ".seq" attributes match."""
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
        """Return ``True`` wherever the class is defined."""
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
        Construct a HybRecord instance from a single-line hyb-format string.

        The Hyb software package ([Travis2014]_) records read-count information 
        in the "id" field of the record, which can be read by setting ``hybformat_id=True``.
        Additionally, the Hyb hOH7 database contains the segment type in the
        identifier of each reference in the 4th field, which can be read by setting
        ``hybformat_ref=True``.

        Args:
            line (str): Hyb-format string containing record information.
            hybformat_id (bool, optional): Read count information from identifier in
                "<read_id>_<read_count>" format. (Default: False)
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
                ref = seg_props['ref_name']
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
    #: Properties for the :meth:`is_set` method.
    #: 
    #: * ``energy``         : record.energy is not None
    #: * ``full_seg_props`` : Each seg key is in segN_props dict and is not None 
    #: * ``fold_record``    : record.fold_record has been set
    #: * ``eval_types``     : seg1_type and seg2_type flags have been set
    #: * ``eval_mirna``     : miRNA_seg flag has been set
    #:
    SET_PROPS = [
        'energy', 'full_seg_props', 'fold_record', 'eval_types', 'eval_mirna', 'eval_target', 
    ]
    _SET_PROPS_SET = set(SET_PROPS)

    #: General record properties for the :meth:`has_prop` method.
    #:
    #: * ``has_indels`` : either seg1 or seg2 alignments has insertions/deltions,
    #:   shown by differing read/refernce length for the same alignment
    GEN_PROPS = [
        'has_indels'
    ]

    #: String-comparison properties for the :meth:`has_prop` method.
    #:
    #: * **Field Types:**
    #:  
    #:   * ``id``           : record.id
    #:   * ``seq``          : record.seq
    #:   * ``seg1``         : record.seg1_props['ref_name']
    #:   * ``seg2``         : record.seg2_props['ref_name']
    #:   * ``any_seg``      : record.seg1_props['ref_name'] OR record.seg1_props['ref_name']
    #:   * ``all_seg``      : record.seg1_props['ref_name'] AND record.seg1_props['ref_name']
    #:   * ``seg1_type``    : seg1_type flag 
    #:   * ``seg2_type``    : seg2_type flag 
    #:   * ``any_seg_type`` : seg1_type OR seg2_type flags
    #:   * ``all_seg_type`` : seg1_type AND seg2_type flags
    #:
    #: * **Comparisons:**
    #:
    #:   * ``is``       : Comparison string matches field exactly
    #:   * ``prefix``   : Comparison string matches beginning of field
    #:   * ``suffix``   : Comparison string matches end of field 
    #:   * ``contains`` : Comparison string is contained within field
    #:
    STR_PROPS = [
        'id_is', 'id_prefix', 'id_suffix', 'id_contains',
        'seq_is', 'seq_prefix', 'seq_suffix', 'seq_contains',
        'seg1_is', 'seg1_prefix', 'seg1_suffix', 'seg1_contains',
        'seg2_is', 'seg2_prefix', 'seg2_suffix', 'seg2_contains',
        'any_seg_is', 'any_seg_prefix', 'any_seg_suffix', 'any_seg_contains',
        'all_seg_is', 'all_seg_prefix', 'all_seg_suffix', 'all_seg_contains',
        'seg1_type_is', 'seg1_type_prefix', 'seg1_type_suffix', 'seg1_type_contains',
        'seg2_type_is', 'seg2_type_prefix', 'seg2_type_suffix', 'seg2_type_contains',
        'any_seg_type_is', 'any_seg_type_prefix', 'any_seg_type_suffix', 'any_seg_type_contains',
        'all_seg_type_is', 'all_seg_type_prefix', 'all_seg_type_suffix', 'all_seg_type_contains',
    ]
    #: miRNA-evaluation-related properties for the :meth:`has_prop` method.
    #: Requires :ref:`miRNA_seg <mirna_seg>` field to be set by :meth:`eval_mirna` method.
    #:
    #: * ``has_mirna``       : Seg1 or seg2 has been identified as a miRNA
    #: * ``no_mirna``        : Seg1 and seg2 have been identified as not a miRNA
    #: * ``mirna_dimer``     : Both seg1 and seg2 have been identified as a miRNA
    #: * ``mirna_not_dimer`` : Only one of seg1 or seg2 has been identifed as a miRNA
    #: * ``5p_mirna``        : Seg1 (5p) has been identifed as a miRNA
    #: * ``3p_mirna``        : Seg2 (3p) has been identifed as a miRNA
    #: 
    MIRNA_PROPS = [
        'has_mirna', 'no_mirna', 'mirna_dimer', 'mirna_not_dimer',
        '5p_mirna', '3p_mirna',
    ]

    #: miRNA-evaluation & string-comparison properties for the :meth:`has_prop` method.
    #: Requires :ref:`miRNA_seg <mirna_seg>` field to be set by :meth:`eval_mirna` method.
    #:
    #: * Field Types:
    #:  
    #:   * ``mirna``       : segN_props['ref_name'] for identified miRNA segN_props
    #:   * ``target``      : segN_props['ref_name'] for identified target segN_props
    #:   * ``mirna_type``  : segN_type for identified miRNA segN for miRNA/target hybrid
    #:   * ``target_type`` : segN_type for identified target segN for miRNA/target hybrid
    #:
    #: * Comparisons:
    #:
    #:   * ``is``       : Comparison string matches field exactly
    #:   * ``prefix``   : Comparison string matches beginning of field
    #:   * ``suffix``   : Comparison string matches end of field 
    #:   * ``contains`` : Comparison string is contained within field
    #:
    MIRNA_STR_PROPS = [
        'mirna_is', 'mirna_prefix', 'mirna_suffix', 'mirna_contains',
        'target_is', 'target_prefix', 'target_suffix', 'target_contains',
        'mirna_seg_type_is', 'mirna_seg_type_prefix', 
        'mirna_seg_type_suffix', 'mirna_seg_type_contains',
        'target_seg_type_is', 'target_seg_type_prefix', 
        'target_seg_type_suffix', 'target_seg_type_contains',
    ]
    #: Target-evaluation-related properties for the :meth:`has_prop` method.
    #: Requires :ref:`target_reg <target_reg>` field to be set.
    #:
    #: * ``target_none``    : Identified to have no miRNA target
    #: * ``target_unknown`` : Unknown whether there is a miRNA target
    #: * ``target_ncrna``   : miRNA target is identified as in a noncoding transcript
    #: * ``target_5p_utr``  : miRNA target is identified as in the 5p UnTranslated Region
    #:   of a coding transcript
    #: * ``target_3p_utr``  : miRNA target is identified as in the 5p UnTranslated Region
    #:   of a coding transcript
    #: * ``target_coding``  : miRNA target is identified as in coding region 
    #:   of a coding transcript
    #: 
    TARGET_PROPS = [
        'target_none', 'target_unknown', 'target_ncrna', 
        'target_5p_utr', 'target_3p_utr', 'target_coding', 
    ]

    #: All allowed properties for the :meth:`has_prop()` method.
    #: See :attr:`GEN_PROPS`, :attr:`STR_PROPS`, :attr:`MIRNA_PROPS`, 
    #: :attr:`MIRNA_STR_PROPS`, and :attr:`TARGET_PROPS` for details.
    HAS_PROPS = GEN_PROPS + STR_PROPS + MIRNA_PROPS + MIRNA_STR_PROPS + TARGET_PROPS

    _GEN_PROPS_SET = set(GEN_PROPS)
    _STR_PROPS_SET = set(STR_PROPS)
    _MIRNA_PROPS_SET = set(MIRNA_PROPS)
    _MIRNA_STR_PROPS_SET = set(MIRNA_STR_PROPS)
    _ALL_STR_PROPS_SET = _STR_PROPS_SET | _MIRNA_STR_PROPS_SET
    _TARGET_PROPS_SET = set(TARGET_PROPS)
    _HAS_PROPS_SET = set(HAS_PROPS)

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
            read_vals = [str(seg_props[key]) if seg_props[key] is not None 
                         else '??' for key in ['read_start', 'read_end']]
            if any([read_vals[i] != '??' for i in [0, 1]]):
                ret_string += ':%s-%s' % tuple(read_vals)
            ref_vals = [str(seg_props[key]) if seg_props[key] is not None 
                         else '??' for key in ['ref_start', 'ref_end']]
            if any([ref_vals[i] != '??' for i in [0, 1]]):
                ret_string += '(%s-%s)' % tuple(ref_vals)
        ret_string += suffix
        return ret_string

    # HybRecord : Private Methods : Segment Parsing
    def _ensure_props_read_start_end(self):
        for seg_n, seg_props in enumerate([self.seg1_props, self.seg2_props], start=1):
            for key in ('read_start', 'read_end'):
                if key not in seg_props:
                    message = '"read_start" and "read_end" keys required to be found in the '
                    message += 'segN_prop dicts for creating a "dynamic_seq" object\n'
                    message += 'Record: %s, dict: seg%i_props' % (str(self), seg_n)
                    message += 'is missing key: .' % key 
                    print(message)
                    raise Exception(message)

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

    # HybRecord : Private Methods : Segment Parsing
    def _get_dynamic_seq(self):
        self._ensure_props_read_start_end()
        seq1 = self._get_seg_seq(self.seg1_props) 
        seq2 = self._get_seg_seq(self.seg2_props) 
        return  seq1 + seq2

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
    (Data courtesy of [Gay2018]_)

    .. _vienna_file_format:

    * | The ".vienna" file format used by the ViennaRNA_ package ([ViennaFormat]_; [Lorenz2011]_):

      Example:
          ::

              34_151138_MIMAT0000076_MirBase_miR-21_microRNA_1_19-...
              TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG
              .....((((((.((((((......)))))).))))))   (-11.1)

    .. _ct_file_format:

    * | The ".ct" file format used by UNAFold_ and other packages ([CTFormat]_, [Zuker2003]_):

      Example:
          ::

              41	dG = -8	dH = -93.9	seq1_name-seq2_name
              1	A	0	2	0	1	0	0
              2	G	1	3	0	2	0	0
              ...
              ...
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
            newline (bool, optional): If ``True``, add newline character to the end of each
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
            newline (bool, optional): If ``True``, terminate the returned string with a newline
                character. (Default: False)
        """
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_vienna_lines(newline=False)) + suffix)

    # FoldRecord : Public Methods : HybRecord Comparison
    def count_hyb_record_mismatches(self, hyb_record):
        """
        Count mismatches between ``hyb_record.seq`` and ``fold_record.seq``

        Args:
            hyb_record (HybRecord): hyb_record for comparison.
        """
        if (self.seq == hyb_record.seq):
            return 0
        else:
            mismatches = 0
            for i in range(max([len(hyb_record.seq), len(self.seq)])):
                if hyb_record.seq[i:i+1] != self.seq[i:i+1]:
                    mismatches += 1
            return mismatches

    # FoldRecord : Public Methods : HybFile Comparison
    def matches_hyb_record(self, hyb_record):
        """
        Return ``True`` if self.seq == hyb_record.seq

        Args:
            hyb_record (HybRecord): hyb_record to compare.
        """
        return (self.seq == hyb_record.seq)

    # FoldRecord : Public Methods : HybFile Comparison
    def ensure_matches_hyb_record(self, hyb_record):
        """
        Ensure self.seq == hyb_record.seq

        Args:
            hyb_record (HybRecord): hyb_record to compare.
        """
        if not self.matches_hyb_record(hyb_record):
            message = 'Disallowed mismatch between HybRecord sequence and FoldRecord sequence.\n'
            message += 'Hyb : %s\n' % str(hyb_record.seq)
            message += 'Fold: %s\n' % str(self.seq) 
            print(message)
            raise Exception(message)

    # FoldRecord : Public MagicMethods : Comparison
    def __eq__(self, other):
        """Return ``True`` if two records have matching sequences and identifiers."""
        return (self.id == other.id and self.seq == other.seq)

    # FoldRecord : Public MagicMethods : Evaluation
    def __hash__(self):
        """Return a hash of the record ".id" attribute"""       
        return hash(self.id)

    # FoldRecord : Public MagicMethods : Evaluation
    def __bool__(self):
        """Return ``True`` wherever the class is defined."""
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
                          error_mode='raise',
                         ):
        """
        Construct instance from a list of 3 strings of vienna-format ([ViennaFormat]_) lines.

        Args:
            record_lines (str or tuple): Iterable of 3 strings corresponding to lines of a
                vienna-format record.
            error_mode (str, optional): | Error mode. Options: 
                                        | ``raise`` : Raise an error when 
                                          encountered and exit program
                                        | ``warn_return`` : Print a warning and return 
                                          the error_value 
                                        | ``return`` : Return the error value with no program output.
        """

        error_mode_options = {'raise', 'warn_return', 'return'}
        if error_mode not in error_mode_options:
            message = 'Provided error mode: %s is not in allowed options\n'
            message += '    ' + ', '.join(error_mode_options)
            print(message)
            raise Exception(message)

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
           if 'return' in error_mode:
                if 'warn' in error_mode:
                    message = 'WARNING: Improper Vienna: No Fold (Energy = 99*.*)'
                    print(message)
                return ('NOFOLD', ''.join(record_lines))
           else:
                message = 'ERROR: Improper Vienna: No Fold (Energy = 99*.*)'
                print(message)
                raise Exception(message)

        if len(line_3_split) != 2:
            message = 'Provided Vienna Record Line 3:\n'
            message += line_3.rstrip() + '\n'
            message += str(line_3_split) + '\n'
            message += '\n  ... does not have required "..((.))..<tab>(-1.23)" format.'
            print(message)
            raise Exception(message)
        fold = line_3_split[0]
        energy = float(line_3_split[1].strip('()'))

        return_obj = cls(rec_id, seq, fold, energy)
        return return_obj

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_string(cls, record_string, error_mode='raise'):
        """
        Construct instance from a string representing 3 vienna-format ([ViennaFormat]_) lines.

        Args:
            record_string (str or tuple): 3-line string containing 
                a vienna-format record
            error_mode (str, optional): 'string representing the error mode.
                Options: "raise": Raise an error when encountered and exit program;
                "warn_return": Print a warning and return the error_value ;
                "return": Return the error value with no program output.
                record is encountered.
        """
        lines = record_string.strip().split('\n')[0:3]
        return cls.from_vienna_lines(lines)

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_lines(cls, record_lines, error_mode='raise'):
        """
        Create a FoldRecord entry from a list of an arbitrary number of strings
        corresponding to lines in the ".ct" file format ([CTFormat]_).

        Args:
            record_lines (list or tuple): list containing lines of ct record
            error_mode (str, optional):   | Error mode. Options: 
                                          | ``raise`` : Raise an error when 
                                            encountered and exit program
                                          | ``warn_return`` : Print a warning and return 
                                            the error_value 
                                          | ``return`` : Return the error value with no program output.
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
    def from_ct_string(cls, record_string, error_mode='raise'):
        """
        Create a FoldRecord entry from a multi-line string containing information in the
        ".ct" file format ([CTFormat]_).

        Args:
            record_string (str): String containing lines of ct record
            error_mode (str, optional): | Error mode. Options: 
                                        | ``raise`` : Raise an error when 
                                          encountered and exit program
                                        | ``warn_return`` : Print a warning and return 
                                          the error_value 
                                        | ``return`` : Return the error value with no program output
        """
        lines = record_string.strip().split('\n')
        return cls.from_ct_lines(lines)

    # FoldRecord : Private Classmethods : Parsing : Output
    #@classmethod
    #def _format_seg_props(cls, seg_props, prefix='', suffix='', indent_str=''):
    #    raise NotImplementedError
    #    # Returns a formatted string of the sgement info information
    #    ret_string = prefix
    #    ret_string += indent_str + 'Map Reference:  %s\n' % seg_props['ref']
    #    ret_string += indent_str + 'Read Start Pos: %s\n' % seg_props['read_start']
    #    ret_string += indent_str + 'Read End Pos:   %s\n' % seg_props['read_end']
    #    ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_props['ref_start']
    #    ret_string += indent_str + 'Map End Pos:    %s\n' % seg_props['ref_end']
    #    ret_string += indent_str + 'Map Score:      %s\n' % seg_props['score']
    #    ret_string += suffix
    #    return ret_string

    # HybRecord : Private Methods : Segment Parsing
    def _get_seg_fold(self, seg_props, hyb_record=None):
        seg_start, seg_end = seg_props['read_start'], seg_props['read_end']
        return self.fold[seg_start-1:seg_end]


class DynamicFoldRecord(FoldRecord):
    """
    Class for storing secondary structure (folding) information for a nucleotide sequence.
    
    Instead of expecting the nucleotide sequence to match a potential :attr:`HybRecord.seq`
    attribute exactly, this type of fold record is 
    expected to be reconstructed from aligned regions of 
    a chimeric read. For chimeras with overlapping alignments, the sequence will be longer. 
    For chimeras with gapped alignments, the sequence will be shorter.

    For an example read with overlapping aligned portions::

        Orignal:
        seg1: 11111111111111111111
        seg2:                   2222222222222222222       
        seq: TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG

        Dynamic:
        seg1: 11111111111111111111
        seg2:                     2222222222222222222       
        seq:  TAGCTTATCAGACTGATGTTTTAGCTTATCAGACTGATG

    For an example read with gapped aligned portions::

        Orignal:
        seg1:  1111111111111111
        seg2:                    222222222222222222       
        seq:  TAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG

        Dynamic:
        seg1: 1111111111111111
        seg2:                 222222222222222222       
        seq:  AGCTTATCAGACTGATTAGCTTATCAGACTGATG

    This type of sequence is found in the Hyb program \\*_hybrids_ua.hyb file type. 
    This is primarily relevant in error-checking when setting a :obj:`HybRecord.fold_record`
    attribute.

    | The primary diffences in this class from the base :class:`FoldRecord` class 
        include modified versions of the methods:
    | :meth:`count_hyb_record_mismatches`
    | :meth:`matches_hyb_record`
    | :meth:`ensure_matches_hyb_record`

    Args:
        id (str): Identifier for record
        seq (str): Nucleotide sequence of record.
        fold (str): Fold representation of record.
        energy (str or float, optional): Energy of folding for record.

    Attributes:
        id (str): Sequence Identifier (often seg1name-seg2name)
        seq (str): Genomic Sequence
        fold (str): Dot-bracket Fold Representation, '(', '.', and ')' characters
        energy (float or None): Predicted energy of folding
    """

    # DynamicFoldRecord : Public Methods : HybRecord Comparison
    def count_hyb_record_mismatches(self, hyb_record):
        """
        Count mismatches between dynamic hyb_record.seq and fold_record.seq

        Args:
            hyb_record (HybRecord): hyb_record for comparison
        """
        dynamic_seq = hyb_record._get_dynamic_seq()
        if (self.seq == dynamic_seq):
            return 0
        else:
            mismatches = 0
            for i in range(max([len(dynamic_seq), len(self.seq)])):
                if dynamic_seq[i:i+1] != self.seq[i:i+1]:
                    mismatches += 1
            return mismatches

    # DynamicFoldRecord : Public Methods : HybFile Comparison
    def matches_hyb_record(self, hyb_record):
        """Calculate dynamic sequence from hyb record and compare to sequence in DynamicFoldRecord

        Args:
            hyb_record (HybRecord): hyb_record for comparison
        """
        dynamic_seq = hyb_record._get_dynamic_seq()
        if (self.seq == dynamic_seq):
            return True
        else:
            match_str = ''
            for i in range(max([len(dynamic_seq), len(self.seq)])):
                if dynamic_seq[i:i+1] == self.seq[i:i+1]:
                    match_str += '|'
                else:
                    match_str += '.'
                mismatch_count = match_str.count('.')
            return mismatch_count <= self.settings['allowed_mismatches']

    # DynamicFoldRecord : Public Methods : HybFile Comparison
    def ensure_matches_hyb_record(self, hyb_record):
        """
        Ensure the dynamic fold record sequence matches hyb_record.seq

        Args:
            hyb_record (HybRecord): hyb_record for comparison
        """
        #if True:
        #    dynamic_seq = hyb_record._get_dynamic_seq()
        #    message  = 'HybRecord Seq:         %s\n' % str(hyb_record.seq)
        #    message += 'HybRecord Dynamic Seq: %s\n' % str(dynamic_seq)
        #    message += 'DynamicFoldRecord Seq: %s\n' % str(self.seq) 
        #    print(message)
        if not self.matches_hyb_record(hyb_record):
            dynamic_seq = hyb_record._get_dynamic_seq()
            match_str = ''
            for i in range(max([len(dynamic_seq), len(self.seq)])):
                if dynamic_seq[i:i+1] == self.seq[i:i+1]:
                    match_str += '|'
                else:
                    match_str += '.'
                mismatch_count = match_str.count('.')
            message = 'Disallowed mismatch between HybRecord sequence '
            message += 'and DynamicFoldRecord sequence.\n'
            message += 'ID: %s\n' % hyb_record.id
            dataset = hyb_record._get_flag_or_none('dataset')
            if dataset is not None:
                message += 'Dataset: %s\n' % dataset
            message += 'HybRecord Seq:         %s\t(%i)\n' % (hyb_record.seq, len(hyb_record.seq))
            message += 'HybRecord Dynamic Seq: %s\t(%i)\n' % (dynamic_seq, len(dynamic_seq))
            message += '                       %s\t' % (match_str)
            message += '(%i of %i)\n' % (mismatch_count, 
                                         self.settings['allowed_mismatches'])
            message += 'DynamicFoldRecord Seq: %s\t(%i)\n' % (self.seq, len(self.seq)) 
            raise Exception(message)

    # DynamicFoldRecord : Private Methods : Segment Parsing
    def _get_seg_fold(self, seg_props, hyb_record):
        hyb_record._ensure_props_read_start_end()
        seg_start, seg_end = seg_props['read_start'], seg_props['read_end']
        seg1_start, seg1_end = hyb_record.seg1_props['read_start'], hyb_record.seg1_props['read_end']
        seg1_len = seg1_end - seg1_start + 1
        seg2_start, seg2_end = hyb_record.seg2_props['read_start'], hyb_record.seg2_props['read_end']
        seg1_fold = self.fold[0:seg1_len]
        seg2_fold = self.fold[seg1_len:]
        #print(seg1_fold, seg2_fold)
        #print(self.fold, self.seq)
        #print(len(seg1_fold + seg2_fold), len(self.seq), len(hyb_record._get_dynamic_seq()))
        assert len(seg1_fold + seg2_fold) == len(self.seq) == len(hyb_record._get_dynamic_seq())
        assert len(seg1_fold) == len(hyb_record._get_seg_seq(hyb_record.seg1_props))
        if (seg_start, seg_end) == (seg1_start, seg1_end):
            return seg1_fold 
        elif (seg_start, seg_end) == (seg2_start, seg2_end):
            return seg2_fold
        else:
            raise Exception() 


class FoldFile(object):
    """
    Base class for file-object wrappers that return file lines as FoldRecord objects.

    See :class:`ViennaFile` or :class:`CtFile`.

    """

    #: Class-level settings. See :attr:`settings.FoldFile_settings` for descriptions.
    settings = hybkit.settings.FoldFile_settings

    _foldrecord_types = {'strict': FoldRecord, 'dynamic': DynamicFoldRecord} 

    # FoldFile : Public Methods : Initialization / Closing
    def __init__(self, *args, **kwargs):
        """Wrapper for open() function that stores resulting file."""
        self.fh = open(*args, **kwargs)
        self._post_init_tasks()
        
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
        """Return :class:`FoldRecord` via :meth:`read_record` based on the appropriate file type."""
        return self.read_record()

    # FoldFile : Public Methods : Reading
    def close(self):
        """Close the file handle."""
        self.fh.close()

    # FoldFile : Public Methods : Reading
    def read_record(self):
        """Stub for implementation by subclasses."""
        message = 'FoldFile is a base class and is not meant to be used directly.' 
        print(message)
        raise NotImplementedError(message)

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

    # FoldFile : Public Methods : Writing
    def write_direct(self, write_string):
        """
        Write a string directly to the underlying file handle.
        """
        self.fh.write(write_string)

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
    def read_record(self, error_mode=None):
        """
        Read next three lines and return output as FoldRecord object.

        Args:
            error_mode (str, optional): 'string representing the error mode.
                If None, defaults to :attr:`settings['error_mode'] <ViennaFile.settings>`.
                Options: "raise": Raise an error when encountered and exit program;
                "warn_return": Print a warning and return the error_value ;
                "return": Return the error value with no program output.
                record is encountered.
        """
        line_1 = next(self.fh)
        line_2 = next(self.fh)
        line_3 = next(self.fh)
        if error_mode is None:
            error_mode = self.settings['foldfile_error_mode']
        record_type = self._foldrecord_types[self.settings['foldrecord_type']]
        record = record_type.from_vienna_lines((line_1, line_2, line_3), error_mode)
        return record


    # ViennaFile : Private Methods
    def _to_record_string(self, write_record, newline):
        """Return a :class:`Fold Record` as a Vienna-format string."""
        return write_record.to_vienna_string(newline=newline)


class CtFile(FoldFile):
    """
    Ct file wrapper that returns ".ct" file lines as FoldRecord objects.
    """

    # CtFile : Public Methods
    def read_record(self, error_mode=None):
        """
        Return the next ct record as a :class:`FoldRecord` object.

        Call next(self.fh) to return the first line of the next entry.
        Determine the expected number of following lines in the entry, and read that number
        of lines further. Return lines as a FoldRecord object.
 
        Args:
            error_mode (str, optional): 'string representing the error mode.
                If None, defaults to :attr:`settings['error_mode'] <CtFile.settings>`
                Options: "raise": Raise an error when encountered and exit program;
                "warn_return": Print a warning and return the error_value ;
                "return": Return the error value with no program output.
                record is encountered.
        """
        record_lines = [header]
        if error_mode is None:
            error_mode = self.settings['error_mode']
        header = next(self.fh)
        expected_line_num = int(header.strip().split()[0])
        for i in range(expected_line_num):
            record_lines.append(next(self.fh))
        record_type = self._foldrecord_types[self.settings['foldrecord_type']]
        record = record_type.from_ct_lines(record_lines, error_mode)
        return record

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
    :class:`ViennaFile`, or :class:`CtFile` simultaneously to return 
    a :class:`HybRecord` and :class:`FoldRecord`.
    
    Basic error checking / catching is performed based on the value of the 
    :attr:`~settings['error_mode'] <HybFoldIter.settings>` setting.

    The obtained :class:`FoldRecord` will be set as 
    :attr:`.HybRecord.fold_record` of the returned :class:`HybRecord` object.

    Args:
        hybfile_handle (HybFile) : HybFile object for iteration
        foldfile_handle (FoldFile) : FoldFile object for iteration
        combine (bool, optional) : Use HybRecord.set_fold_record(FoldRecord) 
            and return only the HybRecord.

    Returns:
        (:class:`HybRecord`, :class:`FoldRecord`) 
    """

    #: Class-level settings. See :attr:`settings.HybFoldIter_settings` for descriptions.
    settings = hybkit.settings.HybFoldIter_settings

    # HybFoldIter : Public Methods
    def __init__(self, hybfile_handle, foldfile_handle, combine=False):
        """Please see :class:`HybFoldIter` for initialization information."""
        self.hybfile_handle = hybfile_handle
        self.foldfile_handle = foldfile_handle
        self.counters = Counter()
        self.combine = combine
        self.sequential_skips = 0
        self.last_hyb_record = None
        self.last_fold_record = None

    # HybFoldIter : Public Methods
    def report(self):
        """Create a report of information from iteration"""
        ret_lines = ['HybFoldIter Iteration Report:']
        add_line = 'Combined Iteration Attempts: ' + str(self.counters['total_read_attempts'])
        ret_lines.append(add_line)
        add_line = 'Hyb Record Iteration Attempts: ' + str(self.counters['hyb_record_read_attempts'])
        ret_lines.append(add_line)
        add_line = 'Fold Record Iteration Attempts: ' + str(self.counters['fold_record_read_attempts'])
        ret_lines.append(add_line)
        #add_line = 'Total Skipped Fold-Only Records: ' + str(self.counters['fold_only_skips'])
        #ret_lines.append(add_line)
        add_line = 'Total Skipped Record Pairs: ' + str(self.counters['pair_skips'])
        ret_lines.append(add_line)
        return ret_lines

    # HybFoldIter : Public Methods
    def print_report(self):
        """Create a report of information from iteration"""
        ret_lines = self.report()
        print('\n'.join(ret_lines) + '\n')

    # HybFoldIter : Public Methods
    def __iter__(self):
        """Return an iterator object."""
        return self

    # HybFoldIter : Public Methods
    def __next__(self):
        """Read and return (class:`HybRecord`, :class:`FoldRecord`)"""
        self.counters['total_read_attempts'] += 1
        try:
            self.counters['hyb_record_read_attempts'] += 1
            next_hyb_record = self.hybfile_handle.read_record()
            self.counters['fold_record_read_attempts'] += 1
            next_fold_record = self.foldfile_handle.read_record(error_mode='return')
            error = ''
            do_skip = False
            if 'foldrecord_nofold' in self.settings['error_checks']:
                if isinstance(next_fold_record, tuple) and next_fold_record[0] =='NOFOLD':
                    error = 'Improper FoldRecord: No Fold (Energy = 99*.*)' 
                #while isinstance(next_fold_record, tuple) and next_fold_record[0] =='NOFOLD':
                #    if 'skip' in self.settings['error_mode']:
                #        if self.settings['error_mode'] == 'warn_skip':
                #            print('WARNING: SkipFold: Improper FoldRecord: No Fold (Energy = 99*.*)')  
                #        self.sequential_skips += 1
                #        self.counters['fold_only_skips'] += 1
                #        if self.sequential_skips > self.settings['max_sequential_skips']:
                #            message = 'ERROR: Skipped %i records in a row ' % self.sequential_skips
                #            message += '(max: %i)\n' % self.settings['max_sequential_skips']
                #            message += 'Check for misalignment of records, or disable setting.'
                #            print(message)
                #            raise Exception(message)
                #        self.counters['fold_record_read_attempts'] += 1
                #        next_fold_record = self.foldfile_handle.read_record(error_mode='return')
                #    else:
                #        error = 'Improper FoldRecord: No Fold (Energy = 99*.*)' 

            if not error and 'hybrecord_indel' in self.settings['error_checks']:
                if next_hyb_record.has_prop('has_indels'):
                    error = 'HybRecord: %s has InDels.' % str(next_hyb_record)
                
            if not error and 'max_mismatch' in self.settings['error_checks']:
                hyb_fold_mismatches = next_fold_record.count_hyb_record_mismatches(next_hyb_record)
                if hyb_fold_mismatches > FoldRecord.settings['allowed_mismatches']:
                    error = 'HybRecord: %s ' % str(next_hyb_record)
                    error += 'has: %i ' % hyb_fold_mismatches
                    error += 'mismatches of %i allowed ' % FoldRecord.settings['allowed_mismatches']

            if error:
                if self.settings['error_mode'] == 'raise':
                    error = 'ERROR: ' + error
                    raise Exception(error)
                elif self.settings['error_mode'] == 'warn_skip':
                    print('WARNING: SkipPair:', error)
                elif self.settings['error_mode'] == 'warn_return':
                    print('WARNING:', error)
    
                if 'skip' in self.settings['error_mode']:
                    self.sequential_skips += 1
                    self.counters['pair_skips'] += 1
                    if self.sequential_skips > self.settings['max_sequential_skips']:
                        message = 'ERROR: Skipped %i record pairs in a row ' % self.sequential_skips
                        message += '(max: %i)\n' % self.settings['max_sequential_skips']
                        message += 'Check for misalignment of records, or disable setting.'
                        print(message)
                        raise Exception(message)
                    do_skip = True

        except StopIteration:
            raise
        except:
            message = 'Error at Counter iteration: '
            message += '%s\n' % self.counters['total_read_attempts']
            message += 'Last HybRecord: %s\n' % str(self.last_hyb_record) 
            message += 'Last FoldRecord: %s\n' % str(self.last_fold_record) 
            message += 'Next HybRecord: %s\n' % str(next_hyb_record) 
            message += 'Next FoldRecord: %s\n' % str(next_fold_record) 
            print(message)
            raise

        if do_skip:
            return next(self)
             
        if self.combine:
            next_hyb_record.set_fold_record(next_fold_record)
            ret_obj = next_hyb_record
        else:
            ret_obj = (next_hyb_record, next_fold_record)

        self.last_hyb_record = next_hyb_record
        self.last_fold_record = next_fold_record
        self.sequential_skips = 0
        return ret_obj


# Import the remainder of hybkit code to connect.
import hybkit.analysis
import hybkit.plot
import hybkit.util

