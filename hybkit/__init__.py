#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

r"""
Module storing primary hybkit classes and hybkit API.

This module contains classes and methods for reading, writing, and manipulating data
in the hyb genomic sequence format ([Travis2014]_). For more information, see the
:ref:`hybkit Hyb File Specification`.

An example string of a hyb-format line from [Gay2018]_ is::

    2407_718\tATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC\t.\tMIMAT0000078_MirBase_miR-23a_microRNA\t1\t21\t1\t21\t0.0027\tENSG00000188229_ENST00000340384_TUBB2C_mRNA\t23\t49\t1181\t1207\t1.2e-06


Hybkit functionality is primarily based on classes for storage and evaluation of
chimeric genomic sequences and associated fold-information:

+----------------------------+---------------------------------------------------------------+
| :class:`HybRecord`         | Class to store a single hyb (hybrid) sequence record          |
+----------------------------+---------------------------------------------------------------+
| :class:`FoldRecord`        | Class to store predicted RNA                                  |
|                            | secondary structure information for hybrid reads              |
+----------------------------+---------------------------------------------------------------+

Also included are classes for reading, writing, and iterating over files containing hybrid
information:

+-------------------------+------------------------------------------------------------------+
| :class:`HybFile`        | Class for reading and writing hyb-format files                   |
|                         | [Travis2014]_ containing chimeric RNA sequence information       |
|                         | as :class:`HybRecord` objects                                    |
+-------------------------+------------------------------------------------------------------+
| :class:`ViennaFile`     | Class for reading and writing Vienna-format files                |
|                         | [ViennaFormat]_ containing RNA secondary structure information   |
|                         | in dot-bracket format as :class:`FoldRecord` objects             |
+-------------------------+------------------------------------------------------------------+
| :class:`CtFile`         | -BETA- Class for reading Connectivity Table (.ct)-format files   |
|                         | [CTFormat]_ containing predicted RNA secondary-structure         |
|                         | information as used by UNAFold_ as                               |
|                         | :class:`FoldRecord` objects                                      |
+-------------------------+------------------------------------------------------------------+
| :class:`HybFoldIter`    | Class for concurrent iteration over a :class:`HybFile` and       |
|                         | a :class:`ViennaFile` or :class:`CtFile`                         |
+-------------------------+------------------------------------------------------------------+

"""

import copy
import logging
import os
import sys
from collections import Counter
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Type, Union

from typing_extensions import Self

Bio, Seq, SeqRecord = None, None, None
try:
    import Bio

    # from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ModuleNotFoundError:
    pass

# ----- Linting Directives:
# ruff: noqa: D214 E402 F401 SLF001 TRY301 B028

# Perform *Initial* hybkit submodule imports, remainder at code end.
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import hybkit
import hybkit.__about__
import hybkit.settings
import hybkit.type_finder
from hybkit.__about__ import (
    __author__,
    __contact__,
    __credits__,
    __date__,
    __deprecated__,
    __email__,
    __license__,
    __maintainer__,
    __status__,
    __version__,
)
from hybkit.errors import HybkitArgError, HybkitConstructorError, HybkitIterError, HybkitMiscError

# ----- Begin Custom Types -----
SegProps = Dict[str, Union[float, int, str]]
StrOrNum = Union[float, int, str]
FlagsDict = Dict[str, Any]
Value = Union[float, int, str, bool, None]
FlexFoldRecord = Union['FoldRecord', Tuple['FoldRecord', Any]]
MirnaDetailsArg = Literal[
    'all', 'mirna_ref', 'target_ref', 'mirna_seg_type',
    'target_seg_type', 'mirna_seq', 'target_seq',
    'mirna_fold', 'target_fold'
    ]
ToFastaRecordArg = Literal['hybrid', 'seg1', 'seg2', 'mirna', 'target']
FoldSeqArg = Literal['static', 'dynamic']
ErrorModeArg = Literal['raise', 'warn_return', 'return']
IterErrorModeArg = Literal['raise', 'warn_return', 'warn_skip', 'skip', 'return']

FieldsHeaderReturn = Literal[
    'id', 'seq', 'energy', 'seg1_ref_name', 'seg1_read_start', 'seg1_read_end',
    'seg1_ref_start', 'seg1_ref_end', 'seg1_score', 'seg2_ref_name',
    'seg2_read_start', 'seg2_read_end', 'seg2_ref_start', 'seg2_ref_end',
    'seg2_score', 'flags'
    ]
CsvHeaderReturn = Literal[
    'id,seq,energy,seg1_ref_name,seg1_read_start,seg1_read_end,'
    'seg1_ref_start,seg1_ref_end,seg1_score,seg2_ref_name,'
    'seg2_read_start,seg2_read_end,seg2_ref_start,seg2_ref_end,'
    'seg2_score,flags'
    ]
FoldReturn = Union[
    Tuple[None, str],
    Tuple[Literal['NOFOLD'], str],
    Tuple[Literal['NOENERGY'], str],
    Self
    ]


# ----- Begin HybRecord Class -----
class HybRecord:
    r"""
    Class for storing and analyzing chimeric (hybrid) RNA-seq reads in hyb format.

    Hyb file (hyb) format entries are a GFF-related file format described by [Travis2014]_
    that contain information about a genomic sequence read identified to be a hybrid by
    a chimeric read caller. Each line contains 15 or 16
    columns separated by tabs ("\\t") and provides
    annotations on each component. An example hyb-format line
    from [Gay2018]_::

        2407_718\tATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC\t.\tMIMAT0000078_MirBase_miR-23a_microRNA\t1\t21\t1\t21\t0.0027\tENSG00000188229_ENST00000340384_TUBB2C_mRNA\t23\t49\t1181\t1207\t1.2e-06

    The columns are respectively described in hybkit as:

        ``id``, ``seq``, ``energy``,
        ``seg1_ref_name``, ``seg1_read_start``, ``seg1_read_end``,
        ``seg1_ref_start``, ``seg1_ref_end``, ``seg1_score``,
        ``seg2_ref_name``, ``seg2_read_start``, ``seg2_read_end``,
        ``seg2_ref_start``, ``seg2_ref_end``, ``seg2_score``,
        ``flags``


    (For more information, see the
    :ref:`hybkit Hyb File Specification`)

    The preferred method for reading hyb records from lines is with
    the :func:`HybRecord.from_line` constructor::

        # line = "2407_718\tATC..."
        hyb_record = hybkit.HybRecord.from_line(line)

    This is the constructor used by the  :class:`HybFile` class to parse hyb files.
    For example, to print all hybrid identifiers in a hyb file::

        with hybkit.HybFile('path/to/file.hyb', 'r') as hyb_file:
            # performs "hyb_record = hybkit.HybRecord.from_line(line)" for each line in file
            for hyb_record in hyb_file:
                print(hyb_record.id)

    HybRecord objects can also be constructed directly. A minimum amount of data necessary
    for a HybRecord object is the genomic sequence and its corresponding identifier.

    Examples:
        ::

            hyb_record_1 = hybkit.HybRecord('1_100', 'ACTG')
            hyb_record_2 = hybkit.HybRecord('2_107', 'CTAG', '-7.3')
            hyb_record_3 = hybkit.HybRecord('3_295', 'CTTG', energy='-10.3')

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
        energy (:obj:`str` or :obj:`float`, optional): Predicted energy of sequence
            folding in kcal/mol
        seg1_props (:obj:`dict`, optional): Properties of segment 1 of the record,
            containing possible
            :attr:`segment column <HybRecord.SEGMENT_COLUMNS>` keys:
            (``ref_name``, ``read_start``, ``read_end``, ``ref_start``, ``ref_end``, ``score``)
        seg2_props (:obj:`dict`, optional): Properties of segment 2 of the record,
            containing possible:
            :attr:`segment column <HybRecord.SEGMENT_COLUMNS>` keys:
            (``ref_name``, ``read_start``, ``read_end``, ``ref_start``, ``ref_end``, ``score``)
        flags (:obj:`dict`, optional): Dict with keys of flags for the record and their
            associated values.
            By default flags must be defined in :attr:`ALL_FLAGS` but custom
            flags can be supplied by changing
            :attr:`HybRecord.settings['custom_flags'] <HybRecord.settings>`.
            This setting can also be disabled by setting 'allow_undefined_flags'
            to :obj:`True` in :attr:`HybRecord.settings`.
        allow_undefined_flags (:obj:`bool`, optional): If :obj:`True`, allows flags
            not defined in :attr:`ALL_FLAGS` or
            :attr:`HybRecord.settings['custom_flags'] <HybRecord.settings>`
            to be added to the record. If not provided, defaults to the value in
            :attr:`HybRecord.settings['allow_undefined_flags'] <HybRecord.settings>`.


    .. _HybRecord-Attributes:

    Attributes:
        id (str): Identifier for the hyb record (Hyb format: ``<read-num>_<read-count>``)
        seq (str): Nucleotide sequence of the hyb record
        energy (str): Predicted energy of folding
        seg1_props (dict): Information on chimeric segment 1, contains
            :attr:`segment column <HybRecord.SEGMENT_COLUMNS>` keys:
            ``ref_name`` (:obj:`str`), ``read_start`` (:obj:`int`), ``read_end`` (:obj:`int`),
            ``ref_start`` (:obj:`int`), ``ref_end`` (:obj:`int`), and ``score`` (:obj:`float`).
        seg2_props (dict): Information on segment 2, contains
            :attr:`segment column <HybRecord.SEGMENT_COLUMNS>` keys:
            ``ref_name`` (:obj:`str`), ``read_start`` (:obj:`int`), ``read_end`` (:obj:`int`),
            ``ref_start`` (:obj:`int`), ``ref_end`` (:obj:`int`), and ``score`` (:obj:`float`).
        flags (dict): Dict of flags with possible
            :attr:`flag keys <HybRecord.ALL_FLAGS>` and values as defined in
            the :ref:`Flags` section of the :ref:`Hybkit Hyb File Specification`.
        fold_record (FoldRecord): Information on the predicted secondary structure of the sequence
            set by :func:`set_fold_record`.
        allow_undefined_flags (bool): Whether to allow undefined flags to be set.
    """

    # Start HybRecord Attributes
    # HybRecord : Class-Level Constants
    #: Record columns 1-3 defining parameters of the overall hybrid, defined by the Hyb format
    HYBRID_COLUMNS = (
        'id',      #: str, Hybrid read identifier, ex: "1257_12"
        'seq',     #: str (alphabetic), Hybrid seq, ex: "ATCGGCTAATCGGTCA..."
        'energy',  #: float or str(of float), hybrid folding energy, ex: "-11.33"
    )

    #: Record columns 4-9 and 10-15, respectively, defining annotated parameters of
    #: seg1 and seg2 respectively, defined by the Hyb format
    SEGMENT_COLUMNS = (
        'ref_name',    # str, Reference ID, ex: "MIMAT0000076_MirBase_miR-21_microRNA"
        'read_start',  # int, Start position of mapping within hybrid, ex: "0"
        'read_end',    # int, End position of mapping within hybrid, ex: "21"
        'ref_start',   # int, Start position of mapping within reference, ex: "1723"
        'ref_end',     # int, End position of mapping within reference, ex: "1744"
        'score',       # str(int or float) Alignment Score, can be BLAST e-score
                       #   or mapping alignment score, depending on analysis
                       #   implementation type.
    )

    # Arbitrary details included in column 16 of the hyb format in the form:
    #   “feature1=value1;feature2=value2;..."
    #   Flags utilized in the Hyb software package, see specification.
    _HYB_FLAGS = (
        'count_total',
        'count_last_clustering',
        'two_way_merged',
        'seq_IDs_in_cluster',
    )
    # Additional flag specifications utilized by hybkit, see specification.
    _HYBKIT_FLAGS = (
        'read_count',
        'orient',
        'det',
        'seg1_type',
        'seg2_type',
        'seg1_det',
        'seg2_det',
        'miRNA_seg',
        'target_reg',
        'ext',
        'dataset',
    )

    #: Flags defined by the hybkit package. Flags 1-4 are utilized by the Hyb software package.
    #: For information on flags, see the :any:`Flags` portion of the
    #: :any:`hybkit Hyb File Specification`.
    ALL_FLAGS = _HYB_FLAGS + _HYBKIT_FLAGS

    #: Class-level settings. See :attr:`settings.HybRecord_settings_info` for descriptions.
    settings = hybkit.settings.HybRecord_settings

    #: Link to :class:`type_finder.TypeFinder` class for parsing sequence identifiers
    #: in assigning segment types by :func:`eval_types`.
    #:
    #: :meta hide-value:
    TypeFinder = hybkit.type_finder.TypeFinder

    # HybRecord : Private Constants
    #: Properties for the :meth:`is_set` method.
    #:
    #: * ``energy``         : :attr:`energy` is not None
    #: * ``full_seg_props`` : Each seg key is in :ref:`segN_props <HybRecord-Attributes>`
    #:   dict and is not None
    #: * ``fold_record``    : :ref:`fold_record <HybRecord-Attributes>` has been set
    #: * ``eval_types``     : :ref:`seg1_type <HybRecord-Attributes>` and
    #:   :ref:`seg2_type <HybRecord-Attributes>` flags have been set
    #: * ``eval_mirna``     : :ref:`miRNA_seg <HybRecord-Attributes>` flag has been set
    #:
    SET_PROPS = (
        'energy', 'full_seg_props', 'fold_record',
        'eval_types', 'eval_mirna', 'eval_target',
    )

    #: General record properties for the :meth:`prop` method.
    #:
    #: * ``has_indels`` : either seg1 or seg2 alignments has insertions/deletions,
    #:   shown by differing read/reference length for the same alignment
    GEN_PROPS = (
        'has_indels',
    )

    #: String-comparison properties for the :meth:`prop` method.
    #:
    #: * **Field Types:**
    #:
    #:   * ``id``           : :ref:`record.id <HybRecord-Attributes>`
    #:   * ``seq``          : :ref:`record.seq <HybRecord-Attributes>`
    #:   * ``seg1``         : :ref:`seg1_props['ref_name'] <HybRecord-Attributes>`
    #:   * ``seg2``         : :ref:`seg2_props['ref_name'] <HybRecord-Attributes>`
    #:   * ``any_seg``      : :ref:`seg1_props['ref_name'] <HybRecord-Attributes>` OR
    #:     :ref:`seg1_props['ref_name'] <HybRecord-Attributes>`
    #:   * ``seg1_type``    : :ref:`seg1_type <seg1_type>` flag
    #:   * ``seg2_type``    : :ref:`seg2_type <seg2_type>` flag
    #:   * ``any_seg_type`` : :ref:`seg1_type <seg1_type>` OR :ref:`seg2_type <seg2_type>` flags
    #:
    #: * **Comparisons:**
    #:
    #:   * ``is``       : Comparison string matches field exactly
    #:   * ``prefix``   : Comparison string matches beginning of field
    #:   * ``suffix``   : Comparison string matches end of field
    #:   * ``contains`` : Comparison string is contained within field
    #:
    STR_PROPS = (
        'id_is', 'id_prefix', 'id_suffix', 'id_contains',
        'seq_is', 'seq_prefix', 'seq_suffix', 'seq_contains',
        'seg1_is', 'seg1_prefix', 'seg1_suffix', 'seg1_contains',
        'seg2_is', 'seg2_prefix', 'seg2_suffix', 'seg2_contains',
        'any_seg_is', 'any_seg_prefix', 'any_seg_suffix', 'any_seg_contains',
        'seg1_type_is', 'seg1_type_prefix', 'seg1_type_suffix', 'seg1_type_contains',
        'seg2_type_is', 'seg2_type_prefix', 'seg2_type_suffix', 'seg2_type_contains',
        'any_seg_type_is', 'any_seg_type_prefix', 'any_seg_type_suffix', 'any_seg_type_contains',
    )
    #: miRNA-evaluation-related properties for the :meth:`prop` method.
    #: Requires :ref:`miRNA_seg <mirna_seg>` flag to be set by :meth:`eval_mirna` method.
    #:
    #: * ``has_mirna``       : **Either or Both** Seg1 or seg2 hve been **identified as** a miRNA
    #: * ``no_mirna``        : **Both** Seg1 and seg2 have been identified as **Not** a miRNA
    #: * ``mirna_dimer``     : **Both** seg1 and seg2 have been **identified as** a miRNA
    #: * ``mirna_not_dimer`` : **One and Only One** of seg1 or seg2
    #:   has been **identified as** a miRNA
    #: * ``5p_mirna``        : Seg1 (5p) has been identified as a miRNA
    #: * ``3p_mirna``        : Seg2 (3p) has been identified as a miRNA
    #:
    MIRNA_PROPS = (
        'has_mirna', 'no_mirna', 'mirna_dimer', 'mirna_not_dimer',
        '5p_mirna', '3p_mirna',
    )

    #: miRNA-evaluation & string-comparison properties for the :meth:`prop` method.
    #: Requires :ref:`miRNA_seg <mirna_seg>` flag to be set by :meth:`eval_mirna` method.
    #:
    #: * Field Types:
    #:
    #:   * ``mirna``       : `segN_props['ref_name'] <HybRecord-Attributes>``
    #      for identified miRNA segN_props
    #:   * ``target``      : `segN_props['ref_name'] <HybRecord-Attributes>``
    #      for identified target segN_props
    #:   * ``mirna_type``  : :ref:`segN_type <seg1_type>` flag
    #      for identified miRNA segN for miRNA/target hybrid
    #:   * ``target_type`` : :ref:`segN_type <seg1_type>` flag
    #      for identified target segN for miRNA/target hybrid
    #:
    #: * Comparisons:
    #:
    #:   * ``is``       : Comparison string matches field exactly
    #:   * ``prefix``   : Comparison string matches beginning of field
    #:   * ``suffix``   : Comparison string matches end of field
    #:   * ``contains`` : Comparison string is contained within field
    #:
    MIRNA_STR_PROPS = (
        'mirna_is', 'mirna_prefix', 'mirna_suffix', 'mirna_contains',
        'target_is', 'target_prefix', 'target_suffix', 'target_contains',
        'mirna_seg_type_is', 'mirna_seg_type_prefix',
        'mirna_seg_type_suffix', 'mirna_seg_type_contains',
        'target_seg_type_is', 'target_seg_type_prefix',
        'target_seg_type_suffix', 'target_seg_type_contains',
    )

    #: All allowed properties for the :meth:`prop()` method.
    #: See :attr:`GEN_PROPS`, :attr:`STR_PROPS`, :attr:`MIRNA_PROPS`, and
    #: :attr:`MIRNA_STR_PROPS`
    HAS_PROPS = GEN_PROPS + STR_PROPS + MIRNA_PROPS + MIRNA_STR_PROPS

    # Start HybRecord Public Methods
    # HybRecord : Public Methods : Initialization
    def __init__(
            self,
            id: str,
            seq: str,
            energy: Optional[StrOrNum] = None,
            seg1_props: Optional[SegProps] = None,
            seg2_props: Optional[SegProps] = None,
            flags: Optional[FlagsDict] = None,
            read_count: Optional[int] = None,
            allow_undefined_flags: Optional[bool] = None,
            ) -> None:
        """Describe __init__ method description in class docstring."""
        if id is None or seq is None:
            message = 'HybRecord initialization requires "id" '
            message += 'and "seq" parameters to be defined.\n'
            message += f'ID: "{id}", Seq: "{seq}"'
            raise HybkitConstructorError(message)
        self.id = self._ensure_attr_types(id, 'id')
        self.seq = self._ensure_attr_types(seq, 'seq')
        self.energy = self._ensure_attr_types(energy, 'energy')

        if seg1_props is None:
            seg1_props = {}
        else:
            seg1_props = copy.deepcopy(seg1_props)
            self._ensure_attr_types(seg1_props, 'seg_props')
        self.seg1_props = self._make_seg_props_dict(seg1_props)

        if seg2_props is None:
            seg2_props = {}
        else:
            seg2_props = copy.deepcopy(seg2_props)
            self._ensure_attr_types(seg2_props, 'seg_props')
        self.seg2_props = self._make_seg_props_dict(seg2_props)

        if allow_undefined_flags is None:
            self.allow_undefined_flags = self.settings['allow_undefined_flags']
        else:
            self.allow_undefined_flags = allow_undefined_flags

        if flags is None:
            flags = {}
        else:
            flags = self._make_flags_dict(copy.deepcopy(flags))
            self._ensure_attr_types(flags, 'flags')
        self.flags = flags
        self.fold_record = None     # Placeholder variable for fold_record

        if read_count is not None:
            if 'read_count' in flags and int(flags['read_count']) != int(read_count):
                message = '"read_count" parameter defined both in function call and flags,\n'
                message += 'but values do not match. Please define only once.'
                raise HybkitConstructorError(message)
            self.set_flag('read_count', str(read_count))

        # if fold_record is not None:
        #     self.set_fold_record(fold_record)

        self._post_init_tasks()    # Method stub for subclassing

    # HybRecord : Public Methods : flags
    def set_flag(
            self,
            flag_key: str,
            flag_val: Value,
            allow_undefined_flags: Optional[bool] = None
            ) -> None:
        """
        Set the value of record ``flag_key`` to ``flag_val``.

        Args:
            flag_key (str): Key for flag to set.
            flag_val : Value for flag to set.
            allow_undefined_flags (:obj:`bool`, optional):
                Allow inclusion of flags not
                defined in :attr:`ALL_FLAGS` or in
                :attr:`settings['custom_flags'] <HybRecord.settings>`.
                If not provided, uses setting in
                :attr:`'HybRecord.allow_undefined_flags'` (Defaults to value in:
                :attr:`settings['allow_undefined_flags'] <HybRecord.settings>` ).
        """
        if allow_undefined_flags is None:
            allow_undefined_flags = self.allow_undefined_flags

        self._ensure_flagset()

        if (not allow_undefined_flags
                and flag_key not in self._flagset):
            message = 'Flag "%s" is not defined. Please check flag key, ' % flag_key
            message += 'run with: "allow_undefined_flags=True", '
            message += 'add flag to "custom_flags" setting, '
            message += 'or set "HybRecord.settings.allow_undefined_flags = True".'
            raise HybkitMiscError(message)

        self.flags[flag_key] = str(flag_val)

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg1_type(
            self,
            require: bool = False
            ) -> Optional[str]:
        """
        Return the :ref:`seg1_type <seg1_type>` flag if defined, or return None.

        Args:
            require: If ``True``, raise an error if seg1_type is not defined.
        """
        return self._get_flag('seg1_type', require)

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg2_type(
            self,
            require: bool = False
            ) -> Optional[str]:
        """
        Return the :ref:`seg2_type <seg2_type>` flag if defined, or return None.

        Args:
            require (:obj:`bool`, optional): If ``True``,
                raise an error if seg2_type is not defined.
        """
        return self._get_flag('seg2_type', require)

    # HybRecord : Public Methods : Flag_Info : seg_type
    def get_seg_types(
            self,
            require: bool = False
            ) -> Tuple[Union[str, None], Union[str, None]]:
        """
        Return "seg1_type" (or None), "seg2_type" (or None) flags.

        Return a tuple of the :ref:`seg1_type <seg1_type>` and :ref:`seg2_type <seg2_type>` flags
        for each respective flag that is defined, or None for each flag that is not.

        Args:
            require (:obj:`bool`, optional): If ``True``,
                raise an error if either flag is not defined.
        """
        ret_tuple = (
            self._get_flag('seg1_type', require),
            self._get_flag('seg2_type', require)
            )
        return ret_tuple

    # HybRecord : Public Methods : Flag_Info : get_read_count
    def get_read_count(
            self,
            require: bool = False
            ) -> Optional[int]:
        """
        Return the :ref:`read_count <read_count>` flag if defined, otherwise return None.

        Args:
            require (:obj:`bool`, optional): If ``True``, raise an error if the "read_count" flag
                is not defined.
        """
        ret_val = self._get_flag('read_count', require)

        if ret_val is None:
            return ret_val
        return int(ret_val)

    # HybRecord : Public Methods : Flag_Info : record_count
    def get_record_count(
            self,
            require: bool = False
            ) -> int:
        """
        Return :ref:`count_total <count_total>` flag if defined, or return 1 (this record).

        Args:
            require (:obj:`bool`, optional): If ``True``, raise an error if the "count_total" flag
                is not defined.
        """
        ret_val = self._get_flag('count_total', require)

        if ret_val is None:
            ret_val = 1
        return int(ret_val)

    # HybRecord : Public Methods : get_mirna_props
    def get_mirna_props(
            self,
            allow_mirna_dimers: bool = False,
            require: bool = True
            ) -> Optional[Dict]:
        """
        Return the seg_props dict corresponding to the miRNA segment, if set.

        If :meth:`eval_mirna` has been run, return the seg_props dict corresponding to the
        miRNA segment type as determined by checking the :ref:`miRNA_seg <mirna_seg>` flag,
        or :obj:`None` if the record does not contain a miRNA.

        Args:
            allow_mirna_dimers (:obj:`bool`, optional): If ``True``, consider miRNA dimers
                as a miRNA/target pair and return the 5p miRNA segment properties.
            require (:obj:`bool`, optional): If ``True``, raise an error if the read does not
                contain a miRNA-annotated segment (Default: ``True``).
        """
        self._ensure_set('eval_mirna')
        if self.prop('has_mirna'):
            if self.prop('mirna_dimer'):
                if allow_mirna_dimers:
                    return self.seg1_props
                elif require:
                    message = 'Record contains a dimer of mirna-annotated segments.\n'
                    message += 'Please use "allow_mirna_dimers=True" to allow this.'
                    raise HybkitMiscError(message)
                else:
                    return None
            elif self.prop('5p_mirna'):
                return self.seg1_props
            elif self.prop('3p_mirna'):
                return self.seg2_props
            else:
                message = 'Uncaught If Branch.'
                raise RuntimeError(message)
        elif require:
            message = 'Record does not contain a miRNA-annotated segment.'
            raise HybkitMiscError(message)
        else:
            return None

    def get_target_props(
            self,
            allow_mirna_dimers: bool = False,
            require: bool = True,
            ) -> Optional[Dict]:
        """
        Return the seg_props dict corresponding to the target segment, if set.

        If :meth:`eval_mirna` has been run, return the seg_props dict corresponding to the
        target segment type as determined by checking the :ref:`miRNA_seg <mirna_seg>` flag,
        (and returning the other segment), or :obj:`None` if the record does not contain a miRNA
        or contains two miRNAs.

        Args:
            allow_mirna_dimers (:obj:`bool`, optional): If ``True``, consider miRNA dimers
                as a miRNA/target pair and return the 3p miRNA segment properties as the
                arbitrarily-selected "target" of the dimer pair.
            require (:obj:`bool`, optional): If ``True``, raise an error if the read does not
                contain a single target-annotated segment (Default: ``True``).
        """
        self._ensure_set('eval_mirna')
        if self.prop('has_mirna'):
            if self.prop('mirna_dimer'):
                if allow_mirna_dimers:
                    return self.seg2_props
                elif require:
                    message = 'Record contains a dimer of mirna-annotated segments.\n'
                    message += 'Please use "allow_mirna_dimers=True" to allow this.'
                    raise HybkitMiscError(message)
                else:
                    return None
            elif self.prop('5p_mirna'):
                return self.seg2_props
            elif self.prop('3p_mirna'):
                return self.seg1_props
            else:
                message = 'Uncaught If Branch.'
                raise RuntimeError(message)
        elif require:
            message = 'Record does not contain a miRNA-annotated segment.'
            raise HybkitMiscError(message)
        else:
            return None

    # HybRecord : Public Methods : Flag_Info : eval_type
    def eval_types(
            self,
            allow_unknown: Optional[bool] = None
            ) -> None:
        """
        Find the types of each segment using the the :class:`TypeFinder` class.

        This method provides :attr:`HybRecord.seg1_props` and :attr:`HybRecord.seg2_props`
        to the :class:`TypeFinder` class, linked as attribute :attr:`HybRecord.TypeFinder`.
        This uses the method: :func:`TypeFinder.find <hybkit.type_finder.TypeFinder.find>`
        set by :func:`TypeFinder.set_method <hybkit.type_finder.TypeFinder.set_method>`
        or :func:`TypeFinder.set_custom_method <hybkit.type_finder.TypeFinder.set_custom_method>`
        to set the
        :ref:`seg1_type <seg1_type>`, :ref:`seg2_type <seg2_type>` flags if not already set.

        To use a type-finding method other than the default,
        prepare the :class:`TypeFinder <hybkit.type_finder.TypeFinder>` class by
        preparing and setting :attr:`TypeFinder.params <hybkit.type_finder.TypeFinder.params>`
        and using :meth:`TypeFinder.set_method <hybkit.type_finder.TypeFinder.set_method>`.

        Args:
            allow_unknown (:obj:`bool`, optional): If ``True``, allow segment types that cannot be
                identified and set them as "unknown". Otherwise raise an error.
                If not provided uses setting in
                :attr:`settings['allow_unknown_seg_types'] <HybRecord.settings>`.
        """
        # If types already set, skip.
        if self.is_set('eval_types'):
            return

        # If allow_unknown is provided, use it. Otherwise use the value in settings.
        if allow_unknown is None:
            allow_unknown = self.settings['allow_unknown_seg_types']

        # Check that TypeFinder method has been set, and if not then set the default.
        self.TypeFinder.check_set_method()

        types = []
        for seg_props in [self.seg1_props, self.seg2_props]:
            seg_type = self.TypeFinder.find(seg_props)
            if seg_type is None:
                if allow_unknown:
                    types.append('unknown')
                else:
                    message = 'Cannot identify segment type for segment:\n'
                    message += self._format_seg_props(seg_props, prefix=(' ' * 2)) + '\n'
                    raise HybkitMiscError(message)
            else:
                types.append(seg_type)
        self.set_flag('seg1_type', types[0])
        self.set_flag('seg2_type', types[1])

    # HybRecord : Public Methods : fold_record
    def set_fold_record(
            self,
            fold_record: FlexFoldRecord,
            allow_energy_mismatch: bool = False
            ) -> None:
        """
        Check and set provided fold_record (:class:`FoldRecord`) as attribute fold_record.

        Ensures that fold_record argument is an instance of FoldRecord and
        has a matching sequence to this HybRecord, then set as
        :ref:`HybRecord.fold_record <HybRecord-Attributes>`.

        Args:
            fold_record (FoldRecord): :attr:`FoldRecord` instance to set as
                :ref:`HybRecord.fold_record <HybRecord-Attributes>`.
            allow_energy_mismatch (:obj:`bool`, optional): If ``True``, allow mismatched fold_record
                and HybRecord energy. Otherwise, raise an error.

        """
        # TODO: update fold_record_reading tuple
        if fold_record is None or (isinstance(fold_record, tuple) and fold_record[0] is None):
            message = 'Trying to assign None object as FoldRecord.'
            raise HybkitMiscError(message)

        if isinstance(fold_record, FoldRecord):
            use_fold_record = fold_record
        elif isinstance(fold_record, tuple) and isinstance(fold_record[0], FoldRecord):
            use_fold_record = fold_record[0]
        else:
            message = 'Supplied argument to fold_record: %s' % str(fold_record)
            message += '\n   is not a FoldRecord object.'
            raise HybkitConstructorError(message)

        use_fold_record.ensure_matches_hyb_record(self)
        self.fold_record = use_fold_record
        if self.fold_record.energy is not None:
            if ((not allow_energy_mismatch)
                    and self.energy not in {None, '.'}
                    and str(self.fold_record.energy) != str(self.energy)):
                message = 'ERROR: HybRecord energy: %s ' % self.energy
                message += 'and FoldRecord energy: %s ' % self.fold_record.energy
                message += 'do not match!'
                raise HybkitConstructorError(message)
            self.energy = self.fold_record.energy

    # HybRecord : Public Methods : eval_mirna
    def eval_mirna(
            self,
            override: bool = False,
            mirna_types: bool = None
            ) -> None:
        """
        Analyze and set miRNA properties from type properties in the hyb record.

        If not already done, determine whether a miRNA exists within this record and
        set the :ref:`miRNA_seg <mirna_seg>` flag.
        This evaluation requires the :ref:`seg1_type <seg1_type>`
        and :ref:`seg2_type <seg2_type>` flags to
        be populated, which can be performed by the :func:`eval_types` method.

        Args:
            override (:obj:`bool`, optional): If ``True``, override existing
                :ref:`miRNA_seg <mirna_seg>` flag if present.
            mirna_types (:obj:`list`, :obj:`tuple`, or :obj:`set`, optional): Iterable
                of string representing sequence types
                considered as miRNA. Otherwise, the types are used
                from :attr:`settings['mirna_types'] <HybRecord.settings>`
                (it is suggested that this be provided as a :obj:`set` for fastest checking).
        """
        if mirna_types is None:
            mirna_types = self.settings['mirna_types']

        # If miRNA_seg flag is not defined, or override is enabled, find and set flag.
        if self.not_set('eval_mirna') or override:
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

    def mirna_details(
            self,
            detail: MirnaDetailsArg = 'all',
            allow_mirna_dimers: bool = False,
            ) -> Optional[Union[Dict, str]]:
        """
        Provide a detail about the miRNA or target following :func:`eval_mirna`.

        Analyze miRNA properties within the sequence record and provide a detail as output.
        Unless ``allow_mirna_dimers`` is ``True``,
        this method requires record to contain a non-dimer miRNA,
        otherwise an error will be raised.

        Args:
            detail (str): | Type of detail to return. Options include:
                          | ``all``             : Dict of all properties (default)
                          | ``mirna_ref``       : Identifier for Assigned miRNA
                          | ``target_ref``      : Identifier for Assigned Target
                          | ``mirna_seg_type``  : Assigned seg_type of miRNA
                          | ``target_seg_type`` : Assigned seg_type of target
                          | ``mirna_seq``       : Annotated subsequence of miRNA
                          | ``target_seq``      : Annotated subsequence of target
                          | ``mirna_fold``      : Annotated fold substring of miRNA
                            (requires fold_record set)
                          | ``target_fold``     : Annotated fold substring of target
                            (requires fold_record set)
            allow_mirna_dimers (:obj:`bool`, optional): Allow miRNA/miRNA dimers.
                The 5p-position will be assigned as the "miRNA",
                and the 3p-position will be assigned as the "target".

        """
        self._ensure_set('eval_mirna')
        mirna_flag = self._get_flag('miRNA_seg')

        if ((not self.prop('has_mirna'))
                or (not allow_mirna_dimers and not self.prop('mirna_not_dimer'))):
            message = 'mirna_detail method requires a hybrid containing a single mirna.\n'
            message += 'hyb_record: %s does not meet this criteria ' % str(self)
            message += 'with miRNA_seg flag: %s' % mirna_flag
            raise HybkitMiscError(message)

        allowed_details = {
            'all', 'mirna_ref', 'target_ref', 'mirna_seg_type', 'target_seg_type',
            'mirna_seq', 'target_seq', 'mirna_fold', 'target_fold',
        }

        if detail not in allowed_details:
            message = 'Requested miRNA detail: "%s" ' % detail
            message += 'not in allowed types: \n    %s' % ', '.join(allowed_details)
            raise HybkitArgError(message)

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
            message = 'Uncaught If Branch.'
            raise RuntimeError(message)

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
    def mirna_detail(self, *args, **kwargs):  # noqa: ANN201 ANN101 ANN002 ANN003
        """
        Deprecate, alias for :meth:`mirna_details`.

        .. deprecated:: v0.3.0
        """
        return self.mirna_details(*args, **kwargs)

    # HybRecord : Public Methods : Record Properties
    def is_set(
            self,
            prop: str,
            ) -> bool:
        """
        Return ``True`` if HybRecord property "prop" is set (if relevant) and is not ``None``.

        Options described in :attr:`SET_PROPS`.

        Args:
            prop (str): Property / Analysis to check
        """
        if prop not in self._SET_PROPS_SET:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available properties are:\n' + ', '.join(self.SET_PROPS)
            raise HybkitMiscError(message)

        if prop in {'energy', 'fold_record'}:
            ret_bool = (getattr(self, prop) is not None)
        elif prop == 'full_seg_props':
            ret_bool = (
                all(
                    k in self.seg1_props and self.seg1_props[k] is not None
                    for k in self.SEGMENT_COLUMNS
                )
                and all(
                    k in self.seg2_props and self.seg1_props[k] is not None
                    for k in self.SEGMENT_COLUMNS
                )
            )
        elif prop == 'eval_types':
            ret_bool = all(st is not None for st in self.get_seg_types())
        elif prop == 'eval_mirna':
            ret_bool = self._get_flag('miRNA_seg') is not None
        elif prop == 'eval_target':
            ret_bool = self._get_flag('target_reg') is not None
        return ret_bool

    # HybRecord : Public Methods : Record Properties
    def not_set(
            self,
            prop: str
            ) -> bool:
        """
        Return ``False`` if HybRecord property "prop" is set (if relevant) and is not ``None``.

        ( returns ``not is_set(prop)`` )

        Args:
            prop (str): Property / Analysis to check
        """
        return not self.is_set(prop)

    # HybRecord : Public Methods : Record Properties
    def prop(
            self,
            prop: str,
            prop_compare: Optional[str] = None,
            ) -> bool:
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
            ======================= =============================================

        Args:
            prop (str):                          Property to check
            prop_compare (:obj:`str`, optional): Comparator to check.

        """
        if prop not in self._HAS_PROPS_SET:
            message = 'Requested Property: %s is not defined. ' % prop
            message += 'Available properties are:\n' + ', '.join(self.HAS_PROPS)
            raise HybkitMiscError(message)

        # Check if a substring compares to a desired property string.
        if prop in self._GEN_PROPS_SET:
            if prop == 'has_indels':
                if not self.is_set('full_seg_props'):
                    message = 'Checking for indels requires full_seg_props to be set.'
                    raise HybkitMiscError(message)
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
                raise HybkitMiscError(message)

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
            elif check_attr == 'any_seg':
                check_info = {self.seg1_props['ref_name'],
                              self.seg2_props['ref_name']}
                multi_check = check_attr.split('_')[0]
            elif check_attr == 'seg1_type':
                check_info = self.get_seg1_type()
                self._ensure_set('eval_types')
            elif check_attr == 'seg2_type':
                check_info = self.get_seg2_type()
                self._ensure_set('eval_types')
            elif check_attr == 'any_seg_type':
                self._ensure_set('eval_types')
                check_info = (self.get_seg1_type(), self.get_seg2_type())
                multi_check = check_attr.split('_')[0]
            elif check_attr == 'mirna':
                check_info = self.mirna_details('mirna_ref', allow_mirna_dimers=True)
            elif check_attr == 'target':
                check_info = self.mirna_details('target_ref', allow_mirna_dimers=True)
            elif check_attr == 'mirna_seg_type':
                check_info = self.mirna_details('mirna_seg_type', allow_mirna_dimers=True)
            elif check_attr == 'target_seg_type':
                check_info = self.mirna_details('target_seg_type', allow_mirna_dimers=True)
            # else:
            #    raise RuntimeError('Unknown Field: ' + check_attr)

            # Wrap single check_info value in a list, if not already.
            if not multi_check:
                check_info = {check_info}

            for info in check_info:
                if info is None:
                    message = 'HybRecord Instance: %s does not have a ' % (str(self))
                    message += 'value for requested property: %s' % check_attr
                    raise HybkitMiscError(message)

            if check_type == 'prefix':
                ret_val = any((val.startswith(prop_compare)) for val in check_info)
            elif check_type == 'contains':
                ret_val = any((prop_compare in val) for val in check_info)
            elif check_type == 'suffix':
                ret_val = any((val.endswith(prop_compare)) for val in check_info)
            elif check_type == 'is':
                ret_val = bool(prop_compare in check_info)
                # ret_val = any((prop_compare == val) for val in check_info)
            # else:
            #    raise RuntimeError(prop)

        # Check mirna-specific properties (requires mirna-evaluation)
        elif prop in self._MIRNA_PROPS_SET:
            self._ensure_set('eval_mirna')
            if prop == 'has_mirna':
                ret_val = self.flags['miRNA_seg'].upper() in {'5P', '3P', 'B'}
            elif prop == 'no_mirna':
                ret_val = self.flags['miRNA_seg'].upper() not in {'5P', '3P', 'B'}
            elif prop == 'mirna_dimer':
                ret_val = self.flags['miRNA_seg'].upper() == 'B'
            elif prop == 'mirna_not_dimer':
                ret_val = self.flags['miRNA_seg'].upper() in {'5P', '3P'}
            elif prop == '5p_mirna':
                ret_val = self.flags['miRNA_seg'].upper() in {'5P', 'B'}
            elif prop == '3p_mirna':
                ret_val = self.flags['miRNA_seg'].upper() in {'3P', 'B'}
            # else:
            #    raise RuntimeError(prop)

        # elif prop in self._TARGET_PROPS_SET:
        #    raise NotImplementedError('Target properties not yet implemented.')
            # self._ensure_set('eval_target')
            # if prop == 'target_none':
            #     ret_val = (self._get_flag('target_reg') == 'N')
            # elif prop == 'target_unknown':
            #     ret_val = (self._get_flag('target_reg') == 'U')
            # elif prop == 'target_ncrna':
            #     ret_val = (self._get_flag('target_reg') == 'NON')
            # elif prop == 'target_5p_utr':
            #     ret_val = (self._get_flag('target_reg') == '5pUTR')
            # elif prop == 'target_coding':
            #     ret_val = (self._get_flag('target_reg') == 'C')
            # elif prop == 'target_3p_utr':
            #     ret_val = (self._get_flag('target_reg') == '3pUTR')
        return ret_val

    # HybRecord : Public Methods : Record Properties
    def has_prop(self, *args, **kwargs): # noqa: ANN201 ANN101 ANN002 ANN003
        """
        Return ``True`` if HybRecord has property: ``prop``.

        .. deprecated:: v0.3.0
           Use :meth:`prop` instead.

        """
        return self.prop(*args, **kwargs)

    # HybRecord : Public Methods : Record Parsing
    def to_line(
            self,
            newline: bool = True,
            sep: str = '\t',
            ) -> str:
        r"""
        Return a hyb-format string representation of the record.

        Args:
            newline (:obj:`bool`, optional): Terminate returned string with
                a newline (default: ``True``)
            sep (:obj:`str`, optional): Separator between fields (Default: "\\t")

        """
        line_items = []
        for item_key in self.HYBRID_COLUMNS:
            item_val = getattr(self, item_key, '.')
            if item_val is None:
                item_val = '.'
            line_items.append(item_val)
        for seg_dict in [self.seg1_props, self.seg2_props]:
            for item_key in self.SEGMENT_COLUMNS:
                if item_key in seg_dict and seg_dict[item_key] is not None:
                    line_items.append(seg_dict[item_key])
                else:
                    line_items.append('.')

        flag_string = self._make_flag_string()

        if flag_string:
            line_items.append(flag_string)

        ret_string = sep.join(str(x) for x in line_items)
        if newline:
            ret_string += '\n'
        return ret_string

    # HybRecord : Public Methods : Record Parsing
    def to_csv(
            self,
            newline: bool = False
            ) -> str:
        """
        Return a comma-separated hyb-format string representation of the record.

        Args:
            newline (:obj:`bool`, optional): If ``True``, end the returned string with a newline.
        """
        return self.to_line(newline, sep=',')

# HybRecord : Public Methods : Record Parsing
    def to_fields(
            self,
            missing_obj : Optional[Value] = None
            ) -> dict:
        """
        Return a python dictionary representation of the record.

        Returns a dictionary with keys corresponding to the fields in the hyb-format
        file, and values corresponding to the values in the record. Output can be
        used with the pandas DataFrame constructor or csv.DictWriter.

        Args:
            missing_obj (optional): Object to use for missing values. Default = :obj:`None`.

        """
        ret_dict = {}
        for item_key in self.HYBRID_COLUMNS:
            item_val = getattr(self, item_key, None)
            if item_val is None:
                item_val = copy.deepcopy(missing_obj)
            ret_dict[item_key] = item_val
        for i, seg_dict in enumerate([self.seg1_props, self.seg2_props], start=1):
            for item_key in self.SEGMENT_COLUMNS:
                out_item_key = f'seg{i}_{item_key}'
                if item_key in seg_dict and seg_dict[item_key] is not None:
                    ret_dict[out_item_key] = seg_dict[item_key]
                else:
                    ret_dict[out_item_key] = copy.deepcopy(missing_obj)
        if self.flags:
            ret_dict['flags'] = copy.deepcopy(self.flags)
        else:
            ret_dict['flags'] = copy.deepcopy(missing_obj)
        return ret_dict

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_record(
            self,
            mode: ToFastaRecordArg = 'hybrid',
            annotate: bool = True,
            allow_mirna_dimers: bool = False
        ) -> SeqRecord:
        """
        Return nucleotide sequence as BioPython SeqRecord object.

        Args:
            mode (:obj:`str`, optional): | Determines which sequence component to return. Options:
                                  | ``hybrid``: Entire hybrid sequence (default)
                                  | ``seg1``: Sequence 1 (if defined)
                                  | ``seg2``: Sequence 2 (if defined)
                                  | ``miRNA``: miRNA sequence of miRNA/target pair
                                    (if defined, else None)
                                  | ``target``: Target sequence of miRNA/target pair
                                    (if defined, else None)
            annotate (:obj:`bool`, optional): Add name of components to fasta sequence identifier
                                       if present.
            allow_mirna_dimers (:obj:`bool`, optional): | If ``True``, allow miRNA dimers to be
                                                   | returned as miRNA sequence (the 5p segment
                                                   | will be selected as the "miRNA").
        """
        global Bio  # noqa: PLW0602
        if Bio is None:
            message = 'BioPython is required for fasta output.'
            raise ModuleNotFoundError(message)
        allowed_modes = ['hybrid', 'seg1', 'seg2', 'mirna', 'target']
        allowed_modes_set = set(allowed_modes)
        mode = mode.lower()
        if mode not in allowed_modes_set:
            message = 'Mode %s not allowed for parsing as fasta record.\n' % mode
            message += '    Allowed Modes: ' + ', '.join(allowed_modes)
            raise HybkitMiscError(message)

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
            if self.not_set('eval_mirna'):
                message = 'eval_mirna must be performed before miRNA/target fasta output'
                raise HybkitMiscError(message)
            elif self.prop('no_mirna'):
                message = 'miRNA / target cannot be output as fasta because record ' + str(self)
                message += 'does not have a miRNA.'
                raise HybkitMiscError(message)
            elif self.prop('mirna_dimer') and not allow_mirna_dimers:
                message = 'miRNA / target cannot be output as fasta because record ' + str(self)
                message += 'has a miRNA dimer.'
                raise HybkitMiscError(message)

            if annotate:
                fasta_description += mode + '--'

            if self.prop('5p_mirna') or self.prop('mirna_dimer'):
                if mode == 'mirna':
                    mode = 'seg1'
                else:
                    mode = 'seg2'
            elif self.prop('3p_mirna'):
                if mode == 'mirna':
                    mode = 'seg2'
                else:
                    mode = 'seg1'

        if mode in {'seg1', 'seg2'}:
            seg_props = getattr(self, mode + '_props')
            fasta_seq = self._get_seg_seq(seg_props)  # Checks read_start and read_end
            read_start, read_end = seg_props['read_start'], seg_props['read_end']
            fasta_id = self.id + ':%i-%i' % (read_start, read_end)
            if annotate:
                fasta_id += ':' + seg_props['ref_name']
                fasta_description += self._format_seg_props_line(seg_props)

        if annotate and 'dataset' in self.flags:
            fasta_id = self.flags['dataset'] + ':' + fasta_id

        fasta_record = SeqRecord(Seq(fasta_seq),
                                 id=fasta_id,
                                 description=fasta_description,
                                 )

        return fasta_record

    # HybRecord : Public Methods : Record Parsing
    def to_fasta_str(
            self,
            mode: ToFastaRecordArg = 'hybrid',
            annotate: bool = True
        ) -> str:
        """
        Return nucleotide sequence as a fasta string.

        Args:
            mode (:obj:`str`, optional): | as with :meth:`to_fasta_record` method.
            annotate (:obj:`bool`, optional): Add name of components to fasta sequence identifier
                                       if present.
        """
        return self.to_fasta_record(mode=mode, annotate=annotate).format('fasta')

    # Start HybRecord Magic Methods
    # HybRecord : Public MagicMethods : Comparison
    def __eq__(self, other: Self) -> bool:
        """Return ``True`` if ".id" and ".seq" attributes match."""
        return (self.id == other.id and self.seq == other.seq)

    # HybRecord : Public MagicMethods : Evaluation
    def __hash__(self) -> int:
        """Return a hash of the record ".id" attribute."""
        return hash(self.id)

    # HybRecord : Public MagicMethods : Evaluation
    def __bool__(self) -> Literal[True]:
        """Return ``True`` wherever the class is defined."""
        return True

    # HybRecord : Public MagicMethods : Evaluation
    def __len__(self) -> int:
        """Return the length of the genomic sequence."""
        return len(self.seq)

    # HybRecord : Public MagicMethods : Printing
    def __str__(self) -> str:
        """Print the identifier of the record."""
        return '<HybRecord ID: %s>' % self.id

    # Start HybRecord Public Classmethods
    # HybRecord : Public Classmethods : Record Construction
    @classmethod
    def from_line(
            cls,
            line: str,
            hybformat_id: bool = False,
            hybformat_ref: bool = False,
            ) -> Self:
        """
        Construct a HybRecord instance from a single-line hyb-format string.

        The Hyb software package ([Travis2014]_) records read-count information
        in the "id" field of the record, which can be read by setting ``hybformat_id=True``.
        Additionally, the Hyb hOH7 database contains the segment type in the
        identifier of each reference in the 4th field, which can be read by setting
        ``hybformat_ref=True``.

        Args:
            line (str): hyb-format string containing record information.
            hybformat_id (:obj:`bool`, optional): If ``True``, read count
                information from identifier in
                ``<read_number>_<read_count>`` format.
            hybformat_ref (:obj:`bool`, optional): If ``True``, read
                additional record information from
                identifier in ``<gene_id>_<transcript_id>_<gene_name>_<seg_type>`` format.

        Returns:
            :class:`HybRecord` instance containing record information.
        """
        line_items = line.strip().split('\t')
        if (len(line_items) < hybkit.settings.MIN_RECORD_FIELDS
            or len(line_items) > hybkit.settings.MAX_RECORD_FIELDS
            ):
            message = 'Hyb record lines require 15 or 16 fields '
            message += 'separated by tab ("\\t") characters, '
            message += 'but only %i were found:\n"%s"' % (len(line_items), line.strip())
            message += 'Line:\n%s' % line
            raise HybkitConstructorError(message)
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
        if len(line_items) > hybkit.settings.MIN_RECORD_FIELDS:
            flags = cls._read_flags(line_items[15])

        if hybformat_id:
            # If 'seq_IDs_in_cluster' flag is set, use it to set count information
            if 'seq_IDs_in_cluster' in flags and flags['seq_IDs_in_cluster'].strip():
                if 'read_count' not in flags or 'count_total' not in flags:
                    combined_read_count = 0
                    combined_record_count = 0
                    for cluster_hyb_id in flags['seq_IDs_in_cluster'].strip(' ,').split(','):
                        cluster_id, cluster_count = cls._parse_hybformat_id(cluster_hyb_id)
                        combined_read_count += int(cluster_count)
                        combined_record_count += 1
                    if 'read_count' not in flags:
                        flags['read_count'] = combined_read_count
                    if 'count_total' not in flags:
                        flags['count_total'] = combined_record_count

            # Otherwise, read count information directly from hybrid identifier
            elif 'read_count' not in flags:
                read_id, read_count = cls._parse_hybformat_id(hyb_id)
                flags['read_count'] = read_count

        if hybformat_ref:
            for i, seg_props in enumerate([seg1_props, seg2_props], start=1):
                ref = seg_props['ref_name']
                seg_type_key = 'seg%i_type' % i
                gene_id, transcript_id, gene_name, seg_type = cls._parse_hybformat_ref(ref)
                if seg_type_key in flags and flags[seg_type_key] != seg_type:
                    message = 'Problem reading in hybformat ref for reference: %s\n' % ref
                    message += 'Inferred type: %s\n' % seg_type
                    message += 'Does not equal current type flag: %s' % flags[seg_type_key]
                    raise HybkitConstructorError(message)
                elif seg_type_key not in flags:
                    flags[seg_type_key] = seg_type

        return_obj = cls(hyb_id, seq, energy, seg1_props, seg2_props, flags)
        return return_obj

    # HybRecord : Public Classmethods : Record Construction
    @classmethod
    def from_fasta_records(
            cls,
            seg1_record: SeqRecord,
            seg2_record: SeqRecord,
            hyb_id: Optional[str] = None,
            energy: Optional[StrOrNum] = None,
            flags: Optional[FlagsDict] = None,
            ) -> Self:
        """
        Construct a HybRecord instance from two BioPython SeqRecord Objects.

        Create artificial HybRecord from two SeqRecord Objects
        For the hybrid:

            | id: ``[seg1_record.id]--[seg2_record.id]``
              (overwritten by "id" parameter if provided)
            | seq: seg1_record.seq + seg2_record

        For each segment:

            | FASTA_Sequence_ID -> segN_ref_name
            | FASTA_Description -> Flags: segN_det
              (Overwritten if segN_det flag is provided directly)

        Optional fields to add via function arguments:

            | hyb_id
            | energy
            | flags

        Args:
            seg1_record (SeqRecord): Biopython SeqRecord object containing information
                on the left/first/5p hybrid segment (seg1)
            seg2_record (SeqRecord): Biopython SeqRecord object containing information
                on the right/second/3p hybrid segment (seg2)
            hyb_id (:obj:`str`, optional): Identifier for the hyb record
                (overwrites generated id if provided)
            energy (:obj:`str` or :obj:`float`, optional): Predicted energy of
                sequence folding in kcal/mol
            flags (:obj:`dict`, optional): Dict with keys of flags for the record and their
                associated values.
                Any flags provided overwrite default-generated flags.

        Returns:
            :class:`HybRecord` instance containing record information.
        """
        if flags is None:
            flags = {}

        for segn_record in [seg1_record, seg2_record]:
            if not isinstance(segn_record, Bio.SeqRecord.SeqRecord):
                message = 'Record is not a valid SeqRecord Object:\n    '
                message += str(segn_record)
                raise HybkitConstructorError(message)

        if hyb_id is None:
            hyb_id = seg1_record.id + '--' + seg2_record.id
        seq = seg1_record.seq + seg2_record.seq
        # energy =
        seg1_len = len(seg1_record.seq)
        seg2_len = len(seg2_record.seq)
        seg1_props = {}
        seg1_props['ref_name'] = seg1_record.id
        seg1_props['read_start'] = 1
        seg1_props['read_end'] = seg1_len
        seg1_props['ref_start'] = None
        seg1_props['ref_end'] = None
        seg1_props['score'] = None
        seg2_props = {}
        seg2_props['ref_name'] = seg2_record.id
        seg2_props['read_start'] = seg1_len + 1
        seg2_props['read_end'] = seg1_len + seg2_len
        seg2_props['ref_start'] = None
        seg2_props['ref_end'] = None
        seg2_props['score'] = None
        if 'seg1_det' not in flags:
            flags['seg1_det'] = seg1_record.description
        if 'seg2_det' not in flags:
            flags['seg2_det'] = seg2_record.description

        return_obj = cls(hyb_id, seq, energy, seg1_props, seg2_props, flags)
        return return_obj

    @classmethod
    # HybRecord : Public Classmethods : Record Parsing
    def to_fields_header(cls) -> FieldsHeaderReturn:
        """
        Return a list of the fields in a :class:`HybRecord` object.

        For use with the :func:`to_fields` method.
        """
        ret_fields = copy.deepcopy(list(cls.HYBRID_COLUMNS))
        for i in range(1,3):
            for item_key in cls.SEGMENT_COLUMNS:
                out_item_key = f'seg{i}_{item_key}'
                ret_fields.append(out_item_key)
        ret_fields.append('flags')
        return tuple(ret_fields)

    @classmethod
    # HybRecord : Public Classmethods : Record Parsing
    def to_csv_header(cls, newline: bool = False) -> CsvHeaderReturn:
        """
        Return a comma-separated string representation of the fields in the record.

        For use with the :func:`to_csv` method.

        Args:
            newline (:obj:`bool`, optional): If ``True``, end the returned string with a newline.
        """
        ret_str = ','.join(cls.to_fields_header())
        if newline:
            ret_str += '\n'
        return ret_str

    # Start HybRecord Private Constants
    # HybRecord : Private Constants
    _SET_PROPS_SET = frozenset(SET_PROPS)
    _GEN_PROPS_SET = frozenset(GEN_PROPS)
    _STR_PROPS_SET = frozenset(STR_PROPS)
    _MIRNA_PROPS_SET = frozenset(MIRNA_PROPS)
    _MIRNA_STR_PROPS_SET = frozenset(MIRNA_STR_PROPS)
    _ALL_STR_PROPS_SET = frozenset(_STR_PROPS_SET | _MIRNA_STR_PROPS_SET)
    # _TARGET_PROPS_SET = set(TARGET_PROPS)
    _HAS_PROPS_SET = frozenset(HAS_PROPS)

    # Placeholder for set of allowed flags filled on first use.
    _flagset = None

    # Start HybRecord Private Methods
    # HybRecord : Private Methods : Initialization
    def _post_init_tasks(self) -> None:
        # Stub for subclassing
        pass

    # HybRecord : Private Methods : Record Parsing
    # Ensure attributes passed to constructor are valid for the respective HybRecord attributes
    def _ensure_attr_types(
            self,
            value: Any,  # noqa: ANN401
            attribute: str,
            ) -> Any:  # noqa: ANN401
        err_message = None
        # Value is checked for placeholder value: '.' and converted to None if found.
        if value == '.':
            value = None
        # Check "id" is valid str
        if attribute == 'id':
            if not isinstance(value, str):
                err_message = (
                    'id must be a string. Provided id is type '
                    f'{type(value)}: {value}'
                )
            elif not value.strip():
                err_message = ('id must be a non-empty string. '
                               'Provided id: "%s" is an empty string' % value)

        # Check "seq" is valid str
        elif attribute == 'seq':
            if not isinstance(value, str):
                err_message = (f'seq must be a string. Provided seq is type {type(value)}: {value}')
            elif not value.strip():
                err_message = (f'seq must be a non-empty string. Provided seq: "{value}" is an'
                               ' empty string')
            elif not value.isalpha():
                err_message = 'seq must be alphabetic. Provided seq is non-alphabetic: %s' % value

        # Check "energy" is valid float, int, or str,
        # if float or int, convert to str of float with one decimal.
        elif attribute == 'energy':
            if not isinstance(value, (float, int, type(None), str)):
                err_message = ('energy must be a float, int, numeric str, or None. '
                               f'Provided energy: "{value}" is type {type(value)}')
            if isinstance(value, str):
                if not value.strip():
                    err_message = ('energy must be a float, int, numeric str, or None. '
                                   'Provided energy is empty string: "%s"' % value)
                elif '_' in value:
                    err_message = ('energy must be a float, int, numeric str, or None. '
                                   'Provided energy is non-numeric string: "%s"' % value)
                else:
                    try:
                        float(value)
                    except ValueError:
                        err_message = ('energy must be a float, int, numeric str, or None. '
                                       'Provided energy is non-numeric string: "%s"' % value)
            elif isinstance(value, (int, float)):
                value = '%.1f' % float(value)

        # Check "seg_props" is dict or None
        elif attribute == 'seg_props':
            if not isinstance(value, (dict)):  # Nonetype disallowed, need empty dict.
                err_message = ('segN_props must be a dict. Provided segN_props is type '
                               f'{type(value)}: {value}')
            else:
                for prop_key, prop_val in value.items():
                    if prop_key == 'ref_name':
                        if not isinstance(prop_val, (str, type(None))):
                            err_message = ('segN_props ref_name must be a string (or None). '
                                           'Provided ref_name is type '
                                           f'{type(prop_val)}: {prop_val}'
                                           )
                            break
                        elif isinstance(prop_val, str) and not prop_val.strip():
                            err_message = ('segN_props ref_name must be a non-empty string '
                                           '(or None). Provided ref_name is empty string: '
                                           f'"{prop_val}"')
                            break
                    elif prop_key in {'read_start', 'read_end', 'ref_start', 'ref_end'}:
                        if not isinstance(prop_val, (int, str, type(None))):
                            err_message = (
                                f'segN_props {prop_key} must be an int (or None). '
                                'Provided read_start: '
                                f'"{prop_val}" is type {type(prop_val)}: {prop_val}'
                            )
                            break
                        elif isinstance(prop_val, str):
                            if not prop_val.strip():
                                err_message = (
                                    f'segN_props "{prop_key}" must be an int (or None). '
                                    f'Provided "{prop_key}" is '
                                    f'empty string: "{prop_val}"'
                                )
                                break
                            elif prop_val == '.':
                                value[prop_key] = None
                            elif '_' in prop_val:
                                err_message = (
                                    f'segN_props {prop_key} must be an int (or None). '
                                    f'Provided {prop_key} is '
                                    f'non-numeric string: "{prop_val}"'
                                )
                                break
                            else:
                                try:
                                    value[prop_key] = int(prop_val)
                                except ValueError:
                                    err_message = (
                                        f'segN_props {prop_key} must be an int (or None). '
                                        f'Provided {prop_key} '
                                        f'is non-numeric string: "{prop_val}"'
                                    )
                                    break
                    elif prop_key == 'score':
                        if not isinstance(prop_val, (float, int, str, type(None))):
                            err_message = (
                                'segN_props score must be a float, int, numeric str, or None. '
                                f'Provided score: "{prop_val}" is type {type(prop_val)}" '
                                f'{prop_val}'
                            )
                            break
                    else:
                        err_message = (
                            'segN_props must be have only keys: ref_name, read_start, read_end, '
                            'ref_start, ref_end, score. Provided segN_props has key %s' % prop_key
                        )
                        break

        # Check "flags" is dict or None
        elif attribute == 'flags':
            if not isinstance(value, dict):  # Nonetype disallowed, need empty dict.
                err_message = (
                    'flags must be a dict. Provided flags is type '
                    f'{type(value)}: {value}'
                )
            else:
                for key, val in value.items():
                    if not isinstance(val, str):
                        err_message = (
                            'flags values must be str. Provided value for flag: '
                            f'{key} is type {type(val)}: {val}'
                        )
                        break
        if err_message is not None:
            raise HybkitConstructorError(err_message)
            return None
        else:
            return value

    # HybRecord : Private Methods : Record Parsing
    def _format_seg_props(
            self,
            seg_props: SegProps,
            prefix: str = '',
            suffix: str = '',
            indent_str: str = ''
            ) -> str:
        # Returns a formatted string of the segment properties
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
    def _format_seg_props_line(
            self,
            seg_props: SegProps,
            prefix: str = '',
            suffix: str = '',
            ) -> str:
        # Returns a single-line formatted string of the segment properties
        ret_string = prefix
        if seg_props['ref_name'] is not None:
            ret_string += seg_props['ref_name']
            read_vals = [str(seg_props[key]) if seg_props[key] is not None
                         else '??' for key in ['read_start', 'read_end']]
            if any(read_vals[i] != '??' for i in [0, 1]):
                ret_string += ':{}-{}'.format(*tuple(read_vals))
            ref_vals = [str(seg_props[key]) if seg_props[key] is not None
                        else '??' for key in ['ref_start', 'ref_end']]
            if any(ref_vals[i] != '??' for i in [0, 1]):
                ret_string += '({}-{})'.format(*tuple(ref_vals))
        ret_string += suffix
        return ret_string

    # HybRecord : Private Methods : Segment Parsing
    def _ensure_props_read_start_end(self) -> None:
        for seg_n, seg_props in enumerate([self.seg1_props, self.seg2_props], start=1):
            for key in ('read_start', 'read_end'):
                if key not in seg_props or seg_props[key] is None:
                    message = '"read_start" and "read_end" keys required to be found in the '
                    message += 'segN_prop dicts for creating a "dynamic_seq" object\n'
                    message += 'Record: %s, dict: seg%i_props ' % (str(self), seg_n)
                    message += 'is missing key: %s.' % key
                    raise HybkitConstructorError(message)

    # HybRecord : Private Methods : Segment Parsing
    def _get_seg_seq(self, seg_props: SegProps) -> str:
        if any(seg_props[v] is None for v in ['read_start', 'read_end']):
            message = 'Segment subsequence cannot be obtained for '
            message += f'Record {self!s}, Segment {seg_props["ref_name"]}.\n'
            message += 'Record segment is missing one of read_start/read_end.'
            raise HybkitConstructorError(message)
        read_start, read_end = seg_props['read_start'], seg_props['read_end']
        return self.seq[(read_start - 1):read_end]

    # HybRecord : Private Methods : Segment Parsing
    def _get_dynamic_seq(self) -> str:
        """Return the record "dynamic" seq by constructing from segment starts/stops."""
        self._ensure_props_read_start_end()
        seq1 = self._get_seg_seq(self.seg1_props)
        seq2 = self._get_seg_seq(self.seg2_props)
        return seq1 + seq2

    # HybRecord : Private Methods : flags
    def _get_flag(
            self,
            flag_key: str,
            require: bool = False
            ) -> Any:  # noqa: ANN401
        if flag_key in self.flags:
            return self.flags[flag_key]
        elif require:
            message = 'Expected Flag Key: %s, but it is not present in record.' % flag_key
            raise HybkitMiscError(message)
        else:
            return None

    # HybRecord : Private Methods : flags
    def _make_flag_string(self, reorder_flags: Optional[bool] = None) -> str:
        flag_string = ''
        for flag in self._get_flag_keys(reorder_flags=reorder_flags):
            flag_string += (f'{flag}={self.flags[flag]!s};')
        return flag_string.rstrip(';')

    # HybRecord : Private Methods : flags
    def _get_flag_keys(self, reorder_flags: Optional[bool] = None) -> List[str]:
        # reorder_flags option returns flags in default ordering scheme.
        #  If reorder-flags argument provided, it overrides default behavior.
        #  Otherwise, the method falls back to the object-default.
        return_list = []
        if reorder_flags is None:
            reorder_flags = self.settings['reorder_flags']
        if reorder_flags:
            return_list = self._get_ordered_flag_keys()
        else:
            return_list = [*self.flags.keys()]
        return return_list

    # HybRecord : Private Methods : flags
    def _get_ordered_flag_keys(self) -> List[str]:
        self._ensure_flagset()
        return_list = []

        for flag in self.ALL_FLAGS + tuple(self.settings['custom_flags']):
            if flag in self.flags:
                return_list.append(flag)
        for flag in self.flags:
            if flag not in self._flagset:
                return_list.append(flag)
        return return_list

    # HybRecord : Private Methods : flags
    def _make_flags_dict(self, flag_obj: FlagsDict) -> FlagsDict:
        #  allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        #  If argument is provided to the method, it overrides default behavior.
        #  Otherwise, the method falls back to the object-defaults.
        allow_undefined_flags = self.allow_undefined_flags

        self._ensure_flagset()
        if not allow_undefined_flags:
            for flag in flag_obj:
                if flag not in self._flagset:
                    message = 'Flag "%s" is not defined. Please check flag key' % flag
                    message += ' or run with: "allow_undefined_flags=True"\n'
                    message += 'Defined Flags are: '
                    message += ', '.join(self.ALL_FLAGS + tuple(self.settings['custom_flags']))
                    raise HybkitMiscError(message)

        return {k: str(v) for k, v in flag_obj.items()}

    # HybRecord : Private Methods : seg_props
    def _make_seg_props_dict(
            self,
            seg_props_obj: Optional[SegProps] = None
            ) -> SegProps:
        # Create a dictionary with mapping entries for all segment properties.
        return_dict = {}
        if seg_props_obj is None:
            seg_props_obj = {}
        for prop_key in self.SEGMENT_COLUMNS:
            if prop_key in seg_props_obj:
                return_dict[prop_key] = seg_props_obj[prop_key]
            else:
                return_dict[prop_key] = None
        return return_dict

    # HybRecord : Private Methods : properties
    def _ensure_set(
            self,
            prop: str,
            ) -> Literal[True]:
        if self.not_set(prop):
            message = 'Problem with HybRecord instance: %s\n' % str(self)
            message += 'Method requires set attribute/evaluation: "%s" before use.' % prop
            raise HybkitMiscError(message)
        return True

    # Start HybRecord Private Classmethods
    # HybRecord : Private Classmethods : hybformat record parsing
    @classmethod
    def _parse_hybformat_id(cls, record_id: str) -> Tuple[str, str]:
        # Parse id in format: "48_50002" into read_id, read_count
        split_id = tuple(record_id.split('_'))
        if len(split_id) != 2:  # noqa: PLR2004
            message = 'Failed attempt to parse record id: %s in hyb format.\n' % record_id
            message += 'Hyb-Program format record ids have form: <read_id>_<read_count>'
            raise HybkitConstructorError(message)

        return split_id

    # HybRecord : Private Classmethods : hybformat record parsing
    @classmethod
    def _parse_hybformat_ref(cls, seg_ref: str) -> Tuple[str, str, str, str]:
        # Parse reference sequence identifier in format:
        # "ENSG00000146425_ENST00000367089_DYNLT1_mRNA"
        # into <gene_id>_<transcript_id>_<gene_name>_<seg_type> information.
        split_ref = seg_ref.split('_')
        if len(split_ref) != 4:  # noqa: PLR2004
            message = 'Failed attempt to parse segment reference id: "%s"' % seg_ref
            message += ' in hyb format.\n'
            message += 'Hyb-Program format record ids have form:\n'
            message += '    <gene_id>_<transcript_id>_<gene_name>_<seg_type>'
            raise HybkitConstructorError(message)

        return (split_ref[0], split_ref[1], split_ref[2], split_ref[3])

    # HybRecord : Private Classmethods : flags
    @classmethod
    def _read_flags(
            cls,
            flag_string: str,
            allow_undefined_flags: Optional[bool] = None,
            ) -> FlagsDict:
        # allow_undefined_flags allows the inclusion of flags not defined in hybkit.
        # undefined flags allowed in this method by default, to allow the object-level setting to
        # take precedence
        if allow_undefined_flags is None:
            allow_undefined_flags = cls.settings['allow_undefined_flags']

        cls._ensure_flagset()
        flag_string = flag_string.rstrip()
        flag_string = flag_string.rstrip(';')
        flag_pairs = [flag_pair.split('=') for flag_pair in flag_string.split(';')]
        flags = {}
        for flag_key, flag_value in flag_pairs:
            if not allow_undefined_flags and flag_key not in cls._flagset:
                message = 'Problem: Undefined Flag: %s\n' % flag_key
                message += 'Defined Flags: '
                message += ', '.join(cls.ALL_FLAGS + tuple(cls.settings['custom_flags']))
                raise HybkitConstructorError(message)

            flags[flag_key] = str(flag_value)
        return flags

    # HybRecord : Private classmethods : flags
    @classmethod
    def _ensure_flagset(cls: Self) -> None:
        """Ensure the flagset has been created."""
        if cls._flagset is None:
            cls._flagset = frozenset(cls.ALL_FLAGS + tuple(cls.settings['custom_flags']))


# ----- Begin HybFile Class -----
class HybFile:
    r"""
    Wrapper for a hyb-format text file which returns entries (lines) as HybRecord objects.

    Args:
        path (str): Path to text file to open as hyb-format file.
        *args: Arguments passed to :func:`open` function to open a text file for reading/writing.
        hybformat_id (:obj:`bool`, optional): If ``True``, during parsing of lines read count
            information from identifier in ``<read_number>_<read_count>`` format.
            Defaults to value in :attr:`settings['hybformat_id'] <HybFile.settings>`.
        hybformat_ref (:obj:`bool`, optional): If ``True``, during parsing of lines read
            additional record information from
            identifier in ``<gene_id>_<transcript_id>_<gene_name>_<seg_type>`` format.
            Defaults to value in :attr:`settings['hybformat_ref'] <HybFile.settings>`.
        from_file_like (:obj:`bool`, optional): If ``True``, the first argument is treated as a
            file-like object (such as io.StringIO or gzip.GzipFile) and the remaining positional
            arguments are ignored. (Default False``)
        **kwargs: Keyword arguments passed to :func:`open` function to open a text file for
            reading/writing.

    Attributes:
        hybformat_id (bool): Read count information from identifier during line parsing
        hybformat_ref (bool): Read type information from reference name
            during line parsing
        fh (file): Underlying file handle for the HybFile object.

    """

    #: Class-level settings. See :obj:`hybkit.settings.HybFile_settings_info` for descriptions.
    settings = hybkit.settings.HybFile_settings

    # Start HybFile Public Methods
    # HybFile : Public Methods : Initialization / Closing
    def __init__(
            self,
            path: str,
            *args: Any, # noqa: ANN401
            hybformat_id: Optional[bool] = None,
            hybformat_ref: Optional[bool] = None,
            from_file_like: bool = False,
            **kwargs: Any, # noqa: ANN401
            ) -> None:
        """Describe __init__ method description in class docstring."""
        if from_file_like:
            self.fh = path
        else:
            self.fh = open(path, *args, **kwargs)  # noqa: SIM115
        if hybformat_id is None:
            self.hybformat_id = self.settings['hybformat_id']
        else:
            self.hybformat_id = hybformat_id
        if hybformat_ref is None:
            self.hybformat_ref = self.settings['hybformat_ref']
        else:
            self.hybformat_ref = hybformat_ref

    # HybFile : Public Methods : Initialization / Closing
    def __enter__(self, *args: Any, **kwargs: Any) -> Self: # noqa: ANN401
        """Open "with" syntax."""
        return self

    # HybFile : Public Methods : Initialization / Closing
    def __exit__(self, etype, value, traceback) -> None:  # noqa: ANN001
        """Close "with" syntax."""
        self.close()

    # HybFile : Public Methods : Initialization / Closing
    def __iter__(self) -> Self:
        """Return an iterator."""
        return self

    # HybFile : Public Methods : Reading
    def __next__(self) -> Self:
        """Return next line as HybRecord object."""
        next_line = self.fh.__next__()
        return HybRecord.from_line(
            next_line,
            hybformat_id=self.hybformat_id,
            hybformat_ref=self.hybformat_ref,
        )

    # HybFile : Public Methods : Reading
    def close(self) -> None:
        """Close the file."""
        self.fh.close()

    # HybFile : Public Methods : Reading
    def read_record(self) -> str:
        """Return next line of hyb file as HybRecord object."""
        return next(self)

    # HybFile : Public Methods : Reading
    def read_records(self) -> List[str]:
        """Return list of all (remaining) records in hyb file as HybRecord objects."""
        records = []
        for record in self:
            records.append(record)  # noqa: PERF402
        return records

    # HybFile : Public Methods : Writing
    def write_record(self, write_record: HybRecord) -> None:
        """
        Write a HybRecord object to file as a Hyb-format string.

        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.

        Args:
            write_record (HybRecord): Record to write.
        """
        self._ensure_hybrecord(write_record)
        record_string = write_record.to_line(newline=True)
        self.fh.write(record_string)

    # HybFile : Public Methods : Writing
    def write_records(self, write_records: Iterable[HybRecord]) -> None:
        """
        Write a sequence of HybRecord objects as hyb-format lines to the Hyb file.

        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.

        Args:
            write_records (list): List of :class:`HybRecord` objects to write.
        """
        for write_record in write_records:
            self.write_record(write_record)

    # HybFile : Public Methods : Writing
    def write_fh(self:Self, *args, **kwargs) -> None:  # noqa: ANN002 ANN003
        """Write directly to the underlying file handle."""
        self.fh.write(*args, **kwargs)

    # HybFile : Public Methods : Writing
    def write(self, *_args, **_kwargs) -> None:  # noqa: ANN002 ANN003
        """
        Implement no-op / error for "write" method to catch errors.

        Use :meth:`write_record` or :meth:`write_fh()` instead.
        """
        message = 'HybFile.write() is not implemented.\n'
        message += 'Please Use HybFile.write_record() or HybFile.write_fh() instead.'
        raise NotImplementedError(message)

    # Start HybFile Public Classmethods
    # HybFile : Public Classmethods : Initialization
    @classmethod
    def open(
            cls,
            path: str,
            *args: Any,  # noqa: ANN401
            hybformat_id: Optional[bool] = None,
            hybformat_ref: Optional[bool] = None,
            **kwargs: Any,  # noqa: ANN401
            ) -> Self:
        """
        Open a path to a text file using :func:`open` and return a HybFile object.

        Arguments match those of the Python3 built-in :func:`open` function and are
        passed directly to it.

        This method is provided as a convenience function for drop-in replacement of the
        built-in :func:`open` function.

        Specific keyword arguments are provided for HybFile-specific settings:

        Args:
            path (str): Path to file to open.
            hybformat_id (:obj:`bool`, optional): If ``True``, during parsing of lines read count
                information from identifier in ``<read_number>_<read_count>`` format.
                Defaults to value in :attr:`settings['hybformat_id'] <HybFile.settings>`.
            hybformat_ref (:obj:`bool`, optional): If ``True``, during parsing of lines read
                additional record information from
                identifier in ``<gene_id>_<transcript_id>_<gene_name>_<seg_type>`` format.
                Defaults to value in :attr:`settings['hybformat_ref'] <HybFile.settings>`.

        Example usage:
            ::

                with HybFile.open('path/to/file.hyb', 'r') as hyb_file:
                    for record in hyb_file:
                        print(record)

        Args:
            *args: Passed directly to :func:`open`.
            **kwargs: Passed directly to :func:`open`.

        Returns:
            :class:`HybFile` object.
        """
        return cls(
            path,
            *args,
            hybformat_id=hybformat_id,
            hybformat_ref=hybformat_ref,
            from_file_like=False,
            **kwargs,
        )

    # HybFile : Private Methods
    # Check if provided argument ("record") is an instance of HybRecord.
    def _ensure_hybrecord(self, record: HybRecord) -> None:
        if not isinstance(record, HybRecord):
            raise HybkitMiscError('Item: "%s" is not a HybRecord object.' % record)


# ----- Begin FoldRecord Class -----
class FoldRecord:
    r"""
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
    a genomic sequence, and its fold representation.

    .. _fold_record_types:

    Two types of FoldRecord objects are supported, 'static' and 'dynamic'. Static
    FoldRecord objects are those where the 'seq' attribute
    matches exactly to the corresponding
    :attr:`HybRecord.seq` attribute (where applicable). Dynamic FoldRecord
    objects are those where :attr:`FoldRecord.seq` is reconstructed from aligned
    regions of a :attr:`HybRecord.seq` chimeric read: Longer for chimeras with
    overlapping alignments, shorter for chimeras with gapped alignments.

    Overlapping Alignment Example::

        Static:
        seg1: 1111111111111111111111
        seg2:                   222222222222222222222
        seq:  TAGCTTATCAGACTGATGTTTTAGCTTATCAGACTGATG

        Dynamic:
        seg1: 1111111111111111111111
        seg2:                       222222222222222222222
        seq:  TAGCTTATCAGACTGATGTTTTTTTTAGCTTATCAGACTGATG

    Gapped Alignment Example::

        Static:
        seg1:   1111111111111111
        seg2:                     222222222222222222
        seq:  TTAGCTTATCAGACTGATGTTAGCTTATCAGACTGATG

        Dynamic:
        seg1: 1111111111111111
        seg2:                 222222222222222222
        seq:  AGCTTATCAGACTGATTAGCTTATCAGACTGATG

    Dynamic sequences are found in the Hyb program \*_hybrids_ua.hyb file type.
    This is primarily relevant in error-checking when setting the :meth:`HybRecord.set_fold_record`
    method.

    When the 'static' FoldRecord type is used, the following methods are used for
    :obj:`HybRecord.fold_record` error-checking:

      * :meth:`static_count_hyb_record_mismatches`

    When the 'dynamic' FoldRecord type is used, the following methods are used for
    :obj:`HybRecord.fold_record` error-checking:

      * :meth:`dynamic_count_hyb_record_mismatches`

    Args:
        id (str): Identifier for record
        seq (str): Nucleotide sequence of record.
        fold (str): Fold representation of record.
        energy (:obj:`str` or :obj:`float`, optional): Energy of folding for record.
        seq_type (:obj:`str`, optional): Expect sequence to be
            'static' (match exactly to corresponding
            :ref:`HybRecord.seq <HybRecord-Attributes>`), or
            'dynamic' (construct from pieces of
            :ref:`HybRecord.seq <HybRecord-Attributes>`).
            if not provided, defaults to
            :attr:`~settings['seq_type'] <FoldRecord.settings>` setting.
            See :obj:`hybkit.settings.FoldRecord_settings_info` for descriptions.

    .. _FoldRecord-Attributes:

    Attributes:
        id (str): Sequence Identifier (often seg1name-seg2name)
        seq (str): Genomic Sequence
        fold (str): Dot-bracket Fold Representation, '(', '.', and ')' characters
        energy (str): Predicted energy of folding
        seq_type (str): Whether sequence is 'static' or 'dynamic'
            (Default: 'static'; see Args for details)
    """

    # FoldRecord : Class-Level Constants
    #: Class-level settings. See :obj:`hybkit.settings.FoldRecord_settings_info` for descriptions.
    settings = hybkit.settings.FoldRecord_settings

    _seq_type_choices = frozenset(
        hybkit.settings.FoldRecord_settings_info['seq_type'][4]['choices'])
    _error_mode_choices = frozenset(
        hybkit.settings.FoldRecord_settings_info['error_mode'][4]['choices'])

    # Start FoldRecord Public Methods
    # FoldRecord : Public Methods : Initialization
    def __init__(
            self,
            id: str,
            seq: str,
            fold: str,
            energy: StrOrNum = None,
            seq_type: Optional[FoldSeqArg] = None
            ) -> None:
        """Describe __init__ method description in class docstring."""
        # Sequence Identifier (often seg1name-seg2name)
        self.id = self._ensure_attr_types(id, 'id')
        # Genomic Sequence
        self.seq = self._ensure_attr_types(seq, 'seq')
        # Fold Representation, str of '(', '.', and ')' characters
        self.fold = self._ensure_attr_types(fold, 'fold')
        # Predicted energy of folding
        self.energy = self._ensure_attr_types(energy, 'energy')

        if seq_type is not None:
            if seq_type.lower() in self._seq_type_choices:
                self.seq_type = seq_type.lower()
            else:
                message = (f'seq_type must be one of: {self._seq_type_choices}\n'
                           f'Provided: {seq_type}')
                raise HybkitConstructorError(message)
        else:
            self.seq_type = self.settings['seq_type']

    # FoldRecord : Public Methods : Parsing : Vienna
    def to_vienna_lines(self, newline: bool = True) -> List[str]:
        """
        Return a list of lines for the record in vienna format.

        See (:ref:`Vienna File Format <vienna_file_format>`).

        Args:
            newline (:obj:`bool`, optional): Add newline character to the end of each
                returned line. (Default: True)
        """
        ret_lines = []
        suffix = ''
        if newline:
            suffix = '\n'
        ret_id = self.id
        if not ret_id.startswith('>'):
            ret_id = '>' + ret_id
        ret_lines.append(ret_id + suffix)   # Add line 1, id
        ret_lines.append(self.seq + suffix)  # Add line 2, sequence

        # If necessary, create formatted energy string without decimal places for integers
        if isinstance(self.energy, type(None)):
            energy_str = '.'
        elif isinstance(self.energy, str):
            energy_str = self.energy
        elif abs(self.energy - round(self.energy)) > 0.00001:  # noqa: PLR2004
            energy_str = ('%.5f' % self.energy).rstrip('0')
        else:
            energy_str = '%i' % int(round(self.energy))

        line_3 = f'{self.fold}\t({energy_str})'
        ret_lines.append(line_3 + suffix)    # Add line 3, fold representation and energy
        return ret_lines

    # FoldRecord : Public Methods : Parsing : Vienna
    def to_vienna_string(self, newline: bool = True) -> str:
        """
        Return a 3-line string for the record in vienna format.

        See (:ref:`Vienna File Format <vienna_file_format>`).

        Args:
            newline (:obj:`bool`, optional): Terminate
                the returned string with a newline
                character. (Default: True)
        """
        if newline:
            suffix = '\n'
        else:
            suffix = ''
        return ('\n'.join(self.to_vienna_lines(newline=False)) + suffix)

    # FoldRecord : Public Methods : HybRecord Comparison
    def count_hyb_record_mismatches(self, hyb_record: HybRecord) -> int:
        """
        Count mismatches between ``hyb_record.seq`` and ``fold_record.seq``.

        Uses :meth:`static_count_hyb_record_mismatches` if :attr:`seq_type` is ``static``, or
        :meth:`dynamic_count_hyb_record_mismatches` if :attr:`seq_type` is ``dynamic``.

        Args:
            hyb_record (HybRecord): hyb_record for comparison.
        """
        if self.seq_type == 'static':
            return self.static_count_hyb_record_mismatches(hyb_record)
        elif self.seq_type == 'dynamic':
            return self.dynamic_count_hyb_record_mismatches(hyb_record)
        else:
            message = 'seq_type must be one of: %s\n' % str(self._seq_type_choices)
            message += 'Provided: %s' % str(self.seq_type)
            raise HybkitArgError(message)

    # FoldRecord : Public Methods : HybRecord Comparison
    def static_count_hyb_record_mismatches(self, hyb_record: HybRecord) -> int:
        """
        Count mismatches between ``hyb_record.seq`` and ``fold_record.seq``.

        Args:
            hyb_record (HybRecord): hyb_record for comparison.
        """
        if (self.seq == hyb_record.seq):
            return 0
        else:
            mismatches = 0
            for i in range(max([len(hyb_record.seq), len(self.seq)])):
                if hyb_record.seq[i:(i + 1)] != self.seq[i:(i + 1)]:
                    mismatches += 1
            return mismatches

    # FoldRecord : Public Methods : HybRecord Comparison
    def dynamic_count_hyb_record_mismatches(self, hyb_record: HybRecord) -> int:
        """
        Count mismatches between hyb_record.seq and dynamic fold_record.seq.

        Args:
            hyb_record (HybRecord): hyb_record for comparison
        """
        dynamic_seq = hyb_record._get_dynamic_seq()
        if (self.seq == dynamic_seq):
            return 0
        else:
            mismatches = 0
            for i in range(max([len(dynamic_seq), len(self.seq)])):
                if dynamic_seq[i:(i + 1)] != self.seq[i:(i + 1)]:
                    mismatches += 1
            return mismatches

    # FoldRecord : Public Methods : HybFile Comparison
    def matches_hyb_record(
            self,
            hyb_record: HybRecord,
            allowed_mismatches: Optional[int] = None
            ) -> bool:
        """
        Return ``True`` if self.seq and hyb_record.seq mismatches are <= allowed_mismatches.

        Args:
            hyb_record (HybRecord): hyb_record to compare.
            allowed_mismatches (:obj:`int`, optional): Number of mismatches allowed
                for a match. If not provided, defaults to the option in
                :attr:`settings['allowed_mismatches'] <FoldRecord.settings>`.
        """
        if allowed_mismatches is None:
            allowed_mismatches = self.settings['allowed_mismatches']
        mismatches = self.count_hyb_record_mismatches(hyb_record)
        return (mismatches <= allowed_mismatches)

    # FoldRecord : Public Methods : HybFile Comparison
    def ensure_matches_hyb_record(
            self,
            hyb_record: HybRecord,
            allowed_mismatches: Optional[int] = None
            ) -> None:
        """
        Ensure self.seq matches hyb_record.seq, else raise an error.

        Args:
            hyb_record (HybRecord): hyb_record to compare.
            allowed_mismatches (:obj:`int`, optional): Number of mismatches allowed
                for a match. If not provided, defaults to the option in
                :attr:`settings['allowed_mismatches'] <FoldRecord.settings>`.
        """
        if not self.matches_hyb_record(hyb_record, allowed_mismatches):
            if self.seq_type == 'static':
                message = 'Disallowed mismatch between HybRecord sequence '
                message += 'and FoldRecord sequence.\n'
                message += 'Hyb : %s\n' % str(hyb_record.seq)
                message += 'Fold: %s\n' % str(self.seq)
            elif self.seq_type == 'dynamic':
                dynamic_seq = hyb_record._get_dynamic_seq()
                match_str, mismatch_count = self._get_seq_mismatch_string(dynamic_seq, self.seq)
                message = 'Disallowed mismatch between HybRecord sequence '
                message += 'and dynamic FoldRecord sequence.\n'
                message += 'ID: %s\n' % hyb_record.id
                dataset = hyb_record._get_flag('dataset')
                if dataset is not None:
                    message += 'Dataset: %s\n' % dataset
                message += 'FoldRecord Seq:        %s\t(%i)\n' % (self.seq, len(self.seq))
                message += 'HybRecord Seq:         %s\t(%i)\n' % (hyb_record.seq,
                                                                  len(hyb_record.seq))
                message += 'HybRecord Dynamic Seq: %s\t(%i)\n' % (dynamic_seq, len(dynamic_seq))
                message += '                       %s\t' % (match_str)
                message += '(%i of %i)\n' % (mismatch_count,
                                             self.settings['allowed_mismatches'])
                message += 'dynamic FoldRecord Seq: %s\t(%i)\n' % (self.seq, len(self.seq))
            raise HybkitMiscError(message)

    # Start FoldRecord Magic Methods
    # FoldRecord : Public MagicMethods : Comparison
    def __eq__(self, other: Self) -> bool:
        """Return ``True`` if two records have matching sequences and identifiers."""
        return (self.id == other.id and self.seq == other.seq)

    # FoldRecord : Public MagicMethods : Evaluation
    def __hash__(self) -> int:
        """Return a hash of the record ".id" attribute."""
        return hash(self.id)

    # FoldRecord : Public MagicMethods : Evaluation
    def __bool__(self) -> Literal[True]:
        """Return ``True`` wherever the class is defined."""
        return True

    # FoldRecord : Public MagicMethods : Evaluation
    def __len__(self) -> int:
        """Return the length of the genomic sequence."""
        return len(self.seq)

    # FoldRecord : Public MagicMethods : Printing
    def __str__(self) -> str:
        """Print the identifier of the record."""
        return '<FoldRecord ID: %s>' % self.id

    # Start FoldRecord Public Classmethods
    # FoldRecord : Public Classmethods : Construction
    _ERROR_MODE_ARGS_SUFFIX = (
        """
            error_mode (:obj:`str`, optional):
                | Error mode. Options:
                | ``raise`` : Raise an error when
                  encountered and exit program
                | ``warn_return`` : Print a warning and return
                  the error_value
                | ``return`` : Return the error value
                  with no program output.
                | If not provided, uses the value in
                | :attr:`~settings['error_mode'] <HybRecord.settings>`.
            seq_type (:obj:`str`, optional):
                | Sequence type. Options:
                | ``static`` : Static sequence
                | ``dynamic`` : Dynamic sequence
                | If not provided, uses the value in
                  :attr:`~settings['seq_type'] <HybRecord.settings>`.
        """
    )

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_lines(
            cls: Self,
            record_lines: List[str],
            error_mode: Optional[ErrorModeArg] = None,
            seq_type: Optional[FoldSeqArg] = None,
            ) -> FoldReturn:
        """
        Construct instance from a list of 3 strings of vienna-format ([ViennaFormat]_) lines.

        See :ref:`Vienna File Format <vienna_file_format>` for more details.

        Args:
            record_lines (:obj:`list` or :obj:`tuple`): Iterable of 3
                strings corresponding to lines of a vienna-format record.
        """
        # See ERROR_MODE_ARGS_SUFFIX for error_mode and seq_type args

        if error_mode is None:
            error_mode = cls.settings['error_mode']
        elif error_mode not in cls._error_mode_choices:
            message = 'Provided error mode: %s is not in allowed options\n' % error_mode
            message += '    ' + ', '.join(cls._error_mode_choices)
            raise HybkitArgError(message)

        fail_ret_val = (None, ''.join(record_lines))

        if len(record_lines) != 3:  # noqa: PLR2004
            error = 'Provided Vienna Record Lines:\n'
            error += '\n'.join([line.rstrip() for line in record_lines])
            error += '\n  ... are not in required 3-line format.'
            if 'return' in error_mode:
                if 'warn' in error_mode:
                    logging.warning(error)
                return fail_ret_val
            else:
                raise HybkitConstructorError('ERROR: ' + error)

        rec_id = record_lines[0].strip().lstrip('>')
        seq = record_lines[1].strip()
        line_3 = record_lines[2].strip()
        line_3_split = line_3.split('\t')

        # If no fold was created, potentially due to low-complexity sequence
        if '(99' in line_3:
            error = 'Improper Vienna: No Fold (Energy = 99*.*)'
            if 'return' in error_mode:
                if 'warn' in error_mode:
                    logging.warning(error)
                return ('NOFOLD', ''.join(record_lines))
            else:
                raise HybkitConstructorError('ERROR: ' + error)

        if len(line_3_split) != 2 or not line_3_split[1].strip('()'):  # noqa: PLR2004
            error = 'Provided Vienna Record Line 3:\n'
            error += line_3.rstrip() + '\n'
            error += str(line_3_split) + '\n'
            error += '\n  ... does not have required "..((.))..<tab>(-1.23)" format.'
            if 'return' in error_mode:
                if 'warn' in error_mode:
                    logging.warning(error)
                return ('NOENERGY', ''.join(record_lines))
            else:
                raise HybkitConstructorError('ERROR: ' + error)

        fold = line_3_split[0]
        energy = line_3_split[1].strip('()')

        return_obj = cls(rec_id, seq, fold, energy, seq_type=seq_type)
        return return_obj

    # Add error_mode and seq_type to docstring
    from_vienna_lines.__doc__ += _ERROR_MODE_ARGS_SUFFIX

    # FoldRecord : Public Classmethods : Construction : Vienna
    @classmethod
    def from_vienna_string(
            cls,
            record_string: str,
            error_mode: Optional[ErrorModeArg] = None,
            seq_type: Optional[FoldSeqArg] = None,
            ) -> FoldReturn:
        """
        Construct instance from a string representing 3 vienna-format ([ViennaFormat]_) lines.

        See :ref:`Vienna File Format <vienna_file_format>` for more details.

        Args:
            record_string (str or tuple): 3-line string containing
                a vienna-format record
        """
        # See ERROR_MODE_ARGS_SUFFIX for error_mode and seq_type args

        lines = record_string.strip().split('\n')[0:3]
        return cls.from_vienna_lines(lines, error_mode=error_mode, seq_type=seq_type)

    # Add error_mode and seq_type to docstring
    from_vienna_string.__doc__ += _ERROR_MODE_ARGS_SUFFIX

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_lines(
            cls,
            record_lines: List[str],
            error_mode: Optional[ErrorModeArg] = None,
            seq_type: Optional[FoldSeqArg] = None,
            ) -> FoldReturn:
        """
        Create a FoldRecord from a list of record lines in ".ct" format ([CTFormat]_).

        See :ref:`CT File Format <ct_file_format>` for more details.

        Warning:
            This method is in beta stage, and is not well-tested.

        Args:
            record_lines (list or tuple): list containing lines of ct record
        """
        # See ERROR_MODE_ARGS_SUFFIX for error_mode and seq_type args

        if error_mode is None:
            error_mode = cls.settings['error_mode']
        elif error_mode not in cls._error_mode_choices:
            message = 'Provided error mode: %s is not in allowed options\n' % error_mode
            message += '    ' + ', '.join(cls._error_mode_choices)
            raise HybkitArgError(message)

        header_line = record_lines[0].strip()
        if not any(x in header_line for x in ['dG', 'Energy', 'ENERGY']):
            message = 'Provided ct Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... do not begin with expected header.'
            raise HybkitConstructorError(message)

        header_items = header_line.split('\t')
        expected_seq_len = int(header_items[0])
        expected_record_lines = expected_seq_len + 1

        if len(record_lines) != expected_record_lines:
            message = 'Provided ct Record Lines:\n'
            message += '\n'.join([line.rstrip() for line in record_lines])
            message += '\n  ... do not match the expected %i ' % expected_seq_len
            message += 'lines from the header.'
            raise HybkitConstructorError(message)

        energy_string = header_items[1]
        check_energies = ['=99%i.' % i for i in range(10)]
        if any(e in energy_string for e in check_energies):
            message = 'Improper CT: No Fold (Energy = 99*.*)'
            if 'return' in error_mode:
                if 'warn' in error_mode:
                    logging.warning(message)
                return ('NOFOLD', ''.join(record_lines))
            else:
                raise HybkitConstructorError('ERROR: ' + message)

        energy = float(energy_string.split()[-1])
        # enthalpy_string = header_items[2]
        # enthalpy = float(enthalpy_string.split()[-1])
        full_name = header_items[3]

        seq = ''
        fold = ''
        for i, line in enumerate(record_lines[1:], 1):
            line_split = line.strip().split('\t')
            if len(line_split) not in {6, 8}:
                message = 'Provided ct Record Line:\n'
                message += line.rstrip() + '\n'
                message += '\n  ... does not have 6-column or 8-column format'
                raise HybkitConstructorError(message)

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
                raise RuntimeError
            fold += fold_char

        if not len(fold):
            message = 'Improper CT: No Fold (Len = 0)'
            if 'return' in error_mode:
                if 'warn' in error_mode:
                    logging.warning(message)
                return ('NOFOLD', ''.join(record_lines))
            else:
                raise HybkitConstructorError('ERROR: ' + message)

        return_obj = cls(full_name, seq, fold, energy, seq_type=seq_type)
        return return_obj

    # Add error_mode and seq_type to docstring
    from_ct_lines.__doc__ += _ERROR_MODE_ARGS_SUFFIX

    # FoldRecord : Public Classmethods : Construction : Ct
    @classmethod
    def from_ct_string(
            cls,
            record_string: str,
            error_mode: Optional[ErrorModeArg] = None,
            seq_type: Optional[FoldSeqArg] = None,
            ) -> FoldReturn:
        """
        Create a FoldRecord entry from a multi-line string from ".ct" format ([CTFormat]_).

        See :ref:`CT File Format <ct_file_format>` for more details.

        Warning:
            This method is in beta stage, and is not well-tested.

        Args:
            record_string (str): String containing lines of ct record
        """
        # See ERROR_MODE_ARGS_SUFFIX for error_mode and seq_type args

        lines = record_string.strip().split('\n')
        return cls.from_ct_lines(lines, error_mode=error_mode, seq_type=seq_type)

    # Add error_mode and seq_type to docstring
    from_ct_string.__doc__ += _ERROR_MODE_ARGS_SUFFIX

    # FoldRecord : Private Methods : Record Parsing
    # Ensure attributes passed to constructor are valid for the respective FoldRecord attributes
    def _ensure_attr_types(
            self,
            value: Any,  # noqa: ANN401
            attribute: Any,  # noqa: ANN401
            ) -> Any:  # noqa: ANN401
        err_message = None
        # Value is checked for placeholder value: '.' and converted to None if found.
        if value == '.':
            value = None
        if attribute == 'id':
            if not isinstance(value, str):
                err_message = (
                    'id must be a string. Provided id is type '
                    f'{type(value)}: {value}'
                )
            elif not value.strip():
                err_message = ('id must be a non-empty string. '
                               f'Provided id: "{value}" is an empty string')
        elif attribute == 'seq':
            if not isinstance(value, str):
                err_message = ('seq must be a string. Provided seq is type '
                               f'{type(value)}: {value}')
            elif not value.strip():
                err_message = ('seq must be a non-empty string. Provided seq: "%s" is an empty '
                               'string' % value)
            elif not value.isalpha():
                err_message = 'seq must be alphabetic. Provided seq is non-alphabetic: %s' % value
        elif attribute == 'energy':
            if not isinstance(value, (float, int, type(None), str)):
                err_message = ('energy must be a float, int, numeric str, or None. '
                               f'Provided energy: "{value}" is type {type(value)}')
            if isinstance(value, str):
                if not value.strip():
                    err_message = ('energy must be a float, int, numeric str, or None. '
                                   f'Provided energy is empty string: "{value}"')
                elif '_' in value:
                    err_message = ('energy must be a float, int, numeric str, or None. '
                                   f'Provided energy is non-numeric string: "{value}"')
                else:
                    try:
                        float(value)
                    except ValueError:
                        err_message = ('energy must be a float, int, numeric str, or None. '
                                       f'Provided energy is non-numeric string: "{value}"')
            # Store energy value as a string
            elif isinstance(value, (float, int)):
                value = '%.1f' % value
        elif attribute == 'fold':
            if not isinstance(value, str):
                err_message = ('fold must be a string. Provided fold is type '
                               f'{type(value)}: {value}')
            elif not value.strip():
                err_message = ('fold must be a non-empty string. Provided fold: "%s" is an empty '
                               'string' % value)
            elif not all(c in '().-' for c in value):
                err_message = ('fold must be a string of characters in "().-". '
                               'Provided fold: "%s" contains invalid characters' % value)
        if err_message is not None:
            raise HybkitConstructorError(err_message)
            return None
        else:
            return value

    # FoldRecord : Private Classmethods : Parsing : Output
    # @classmethod
    # def _format_seg_props(cls, seg_props, prefix='', suffix='', indent_str=''):
    #    raise NotImplementedError
    #    # Returns a formatted string of the segment info information
    #    ret_string = prefix
    #    ret_string += indent_str + 'Map Reference:  %s\n' % seg_props['ref']
    #    ret_string += indent_str + 'Read Start Pos: %s\n' % seg_props['read_start']
    #    ret_string += indent_str + 'Read End Pos:   %s\n' % seg_props['read_end']
    #    ret_string += indent_str + 'Map Start Pos:  %s\n' % seg_props['ref_start']
    #    ret_string += indent_str + 'Map End Pos:    %s\n' % seg_props['ref_end']
    #    ret_string += indent_str + 'Map Score:      %s\n' % seg_props['score']
    #    ret_string += suffix
    #    return ret_string

    # Start FoldRecord Private Methods
    # FoldRecord : Private Methods : Segment Parsing
    def _get_seg_fold(
            self,
            seg_props: SegProps,
            hyb_record: Optional[HybRecord] = None,
            ) -> str:
        if self.seq_type == 'static':
            return self._static_get_seg_fold(seg_props, hyb_record)
        elif self.seq_type == 'dynamic':
            return self._dynamic_get_seg_fold(seg_props, hyb_record)
        else:
            message = 'seq_type must be one of: %s\n' % str(self._seq_type_choices)
            message += 'Provided: %s' % str(self.seq_type)
            raise HybkitArgError(message)

    # FoldRecord : Private Methods : Segment Parsing
    def _static_get_seg_fold(
            self,
            seg_props: SegProps,
            hyb_record: Optional[HybRecord] = None
            ) -> str:
        seg_start, seg_end = seg_props['read_start'], seg_props['read_end']
        return self.fold[(seg_start - 1):seg_end]

    # FoldRecord : Private Methods : Segment Parsing
    def _dynamic_get_seg_fold(
            self,
            seg_props: SegProps,
            hyb_record: Optional[HybRecord] = None,
            ) -> str:
        hyb_record._ensure_props_read_start_end()
        seg_start, seg_end = seg_props['read_start'], seg_props['read_end']
        seg1_start = hyb_record.seg1_props['read_start']
        seg1_end = hyb_record.seg1_props['read_end']
        seg1_len = seg1_end - seg1_start + 1
        seg2_start = hyb_record.seg2_props['read_start']
        seg2_end = hyb_record.seg2_props['read_end']
        seg1_fold = self.fold[0:seg1_len]
        seg2_fold = self.fold[seg1_len:]
        assert len(seg1_fold + seg2_fold) == len(self.seq) == len(hyb_record._get_dynamic_seq())
        assert len(seg1_fold) == len(hyb_record._get_seg_seq(hyb_record.seg1_props))
        if (seg_start, seg_end) == (seg1_start, seg1_end):
            return seg1_fold
        elif (seg_start, seg_end) == (seg2_start, seg2_end):
            return seg2_fold
        else:
            raise RuntimeError

    # FoldRecord : Private Methods : Seq Comparison
    def _get_seq_mismatch_string(self, seq1: str, seq2: str) -> str:
        match_str = ''
        for i in range(max([len(seq1), len(seq2)])):
            if seq1[i:(i + 1)] == seq2[i:(i + 1)]:
                match_str += '|'
            else:
                match_str += '.'
            mismatch_count = match_str.count('.')
        return match_str, mismatch_count


# ----- Begin FoldFile Class -----
FOLD_FILE_COMMON_ARGS_ATTRS = (
    """

    Args:
        seq_type (:obj:`str`, optional): Type of FoldRecord to return:
            ``static``, or ``dynamic``
            (if not provided, uses
            :attr:`FoldRecord.settings['seq_type'] <FoldRecord.settings>`).
        error_mode (:obj:`str`, optional): String representing the error mode.
            If None, defaults to the value set in
            :attr:`settings['error_mode'] <FoldRecord.settings>`.
            Options:
            "raise": Raise an error when encountered and exit program;
            "warn_return": Print a warning and return the error_value;
            "return": Return the error value with no warnings.
        from_file_like (:obj:`bool`, optional): If True, treat the first argument
            as a file-like object (such as io.StringIO or gzip.GzipFile) and the
            remaining positional arguments are ignored (Default ``False``).
        *args: Passed to :func:`open()`.
        **kwargs: Passed to :func:`open()`.

    Attributes:
        fh (:obj:`file`): File handle for the file being wrapped.
        foldrecord_seq_type (str): Type of FoldRecord to return (see Args)
        error_mode (str): Mode for error catching (see Args)

    Warning:
        Occasionally fold files can be poorly-formatted. In that case, this iterator
        attempts error-catching but this is not always successful so verbose error modes
        are encouraged.
    """
)


class FoldFile:
    """
    #Base class for file-object wrappers that return file lines as FoldRecord objects.

    See :class:`ViennaFile` or :class:`CtFile`.
    """

    # Args/Attrs Description in FOLD_FILE_COMMON_ARGS_ATTRS

    #: Class-level settings. See :attr:`hybkit.settings.FoldFile_settings_info` for descriptions.
    settings = hybkit.settings.FoldFile_settings

    _foldrecord_settings_info = hybkit.settings.FoldRecord_settings_info
    _foldrecord_seq_type_choices = frozenset(_foldrecord_settings_info['seq_type'][4]['choices'])
    _error_mode_choices = frozenset(_foldrecord_settings_info['error_mode'][4]['choices'])

    # Start FoldFile Public Methods
    # FoldFile : Public Methods : Initialization / Closing
    def __init__(
            self,
            *args: Any,  # noqa: ANN401
            seq_type: Optional[FoldSeqArg] = None,
            error_mode: Optional[ErrorModeArg] = None,
            from_file_like: bool = False,
            **kwargs: Any,  # noqa: ANN401
            ) -> None:
        """Wrap for open() function that stores resulting file."""
        if from_file_like:
            self.fh = args[0]
        else:
            self.fh = open(*args, **kwargs)  # noqa: SIM115

        # Set foldrecord_type
        if seq_type is None:
            self.foldrecord_seq_type = None  # Use default value in FoldRecord class
        elif seq_type in self._foldrecord_seq_type_choices:
            self.foldrecord_seq_type = seq_type
        else:
            message = 'Invalid foldrecord_type: {}. Allowed values: {}'.format(
                seq_type, self._foldrecord_seq_type_choices)
            raise HybkitArgError(message)

        # Set error_mode
        if error_mode is None:
            self.error_mode = None  # Use default value in FoldRecord class
        elif error_mode in self._error_mode_choices:
            self.error_mode = error_mode
        else:
            message = (f'Invalid error_mode: {error_mode}. '
                       f'Allowed values: {self._error_mode_choices!s}')
            raise HybkitArgError(message)
        self._post_init_tasks()  # Throws error on base class

    # FoldFile : Public Methods : Initialization / Closing
    def __enter__(self, *args: Any, **kwargs: Any) -> Self:  # noqa: ANN401
        """Open "with" syntax."""
        return self

    # FoldFile : Public Methods : Initialization / Closing
    def __exit__(self, etype: Any, value: Any, traceback: Any) -> None:  # noqa: ANN401
        """Close "with" syntax."""
        self.close()

    # FoldFile : Public Methods : Initialization / Closing
    def __iter__(self) -> Self:
        """Return an iterator."""
        return self

    # FoldFile : Public Methods : Reading
    def __next__(self) -> FoldRecord:
        """Return :class:`FoldRecord` via :meth:`read_record` for this file type."""
        return self.read_record()

    # FoldFile : Public Methods : Reading
    def close(self) -> None:
        """Close the file handle."""
        self.fh.close()

    # FoldFile : Public Methods : Reading
    #: Read the next record from the FoldFile. (Stub for implementation by subclasses)
    read_record = None

    # FoldFile : Public Methods : Reading
    def read_records(self) -> List[FoldRecord]:
        """Return list of all :class:`FoldRecord` objects for this file type."""
        records = []
        for record in self:
            records.append(record)  # noqa: PERF402
        return records

    # FoldFile : Public Methods : Writing
    def write_record(self, write_record: FoldRecord) -> None:
        """
        Write a FoldRecord object for this file type.

        Unlike the file.write() method, this method will add a newline to the
        end of each written record line.

        Args:
            write_record (:class:`FoldRecord`): :class:`FoldRecord` objects to write.
        """
        self._ensure_foldrecord(write_record)
        record_string = self._to_record_string(write_record, newline=True)
        self.fh.write(record_string)

    # FoldFile : Public Methods : Writing
    def write_records(self, write_records: Iterable[FoldRecord]) -> None:
        """
        Write a sequence of FoldRecord objects for this file type.

        Unlike the file.writelines() method, this method will add a newline to the
        end of each written record line.

        Args:
            write_records (list): List of :class:`FoldRecord` objects to write.
        """
        for write_record in write_records:
            self.write_record(write_record)

    # FoldFile : Public Methods : Writing
    def write_fh(self, *args: Any, **kwargs: Any) -> None:  # noqa: ANN401
        """Write directly to the underlying file handle."""
        self.fh.write(*args, **kwargs)

    # FoldFile : Public Classmethods : Initialization
    @classmethod
    def open(
        cls,
        path: str,
        *args: Any,  # noqa: ANN401
        **kwargs: Any,  # noqa: ANN401
        ) -> Self:
        """
        Open a path to a text file using :func:`open` and return relevant file object.

        Arguments match those of the Python3 built-in :func:`open` function and are
        passed directly to it.

        This method is provided as a convenience function for drop-in replacement of the
        built-in :func:`open` function.

        Specific keyword arguments are provided for fold-file-specific settings:

        Args:
            path (str): Path to file to open.
            seq_type (:obj:`str`, optional): Type of FoldRecord to return:
                ``static``, or ``dynamic``
                (if not provided, uses
                :attr:`FoldRecord.settings['seq_type'] <FoldRecord.settings>`).
            error_mode (:obj:`str`, optional): String representing the error mode.
                If None, defaults to the value set in
                :attr:`settings['error_mode'] <FoldRecord.settings>`.
                Options:
                "raise": Raise an error when encountered and exit program;
                "warn_return": Print a warning and return the error_value;
                "return": Return the error value with no warnings.
            *args: Passed directly to :func:`open`.
            **kwargs: Passed directly to :func:`open`.

        Returns:
            :class:`HybFile` object.
        """
        return cls(path, *args, **kwargs, from_file_like=False)

    # Start FoldFile Private Methods
    # FoldFile : Private Methods
    def _post_init_tasks(self) -> None:
        """Stub for subclassing. Raise error on base class."""
        message = 'FoldFile is a base class and is not meant to be used directly.'
        raise NotImplementedError(message)

    # FoldFile : Private Methods
    def _ensure_foldrecord(self, record: Self) -> None:
        if not isinstance(record, FoldRecord):
            raise HybkitMiscError('Item: "%s" is not a FoldRecord object.' % record)

    # Stub method to be replaced by subclasses.
    _to_record_string = None


# Append common FoldFile Args and Attributes description to docstring:
FoldFile.__doc__ += FOLD_FILE_COMMON_ARGS_ATTRS


# ----- Begin ViennaFile Class -----
class ViennaFile(FoldFile):
    """
    Vienna file wrapper that returns vienna-format file lines as FoldRecord objects.

    See :ref:`Vienna File Format <vienna_file_format>` for more information.

    .. _ViennaFile-Attributes:
    """

    # Args/Attrs Description in FOLD_FILE_COMMON_ARGS_ATTRS

    # Start ViennaFile Methods
    # ViennaFile : Public Methods : Reading
    def read_record(self, override_error_mode: Optional[ErrorModeArg] = None) -> FoldReturn:
        """
        Read next three lines and return output as FoldRecord object.

        Args:
            override_error_mode (str): Override the error_mode set in the
                :class:`ViennaFile` object. See the
                :ref:`ViennaFile Constructor <ViennaFile-Attributes>` for more
                information on allowed error modes.
        """
        line_1 = next(self.fh)
        line_2 = next(self.fh)
        line_3 = next(self.fh)
        if override_error_mode is None:
            use_error_mode = self.error_mode
        else:
            use_error_mode = override_error_mode
        record = FoldRecord.from_vienna_lines(
            (line_1, line_2, line_3),
            error_mode=use_error_mode,
            seq_type=self.foldrecord_seq_type,
        )

        return record

    # ViennaFile : Private Methods
    def _post_init_tasks(self) -> None:
        """Hold place for any post-initialization tasks."""
        pass

    # ViennaFile : Private Methods
    def _to_record_string(self, write_record: FoldRecord, newline: bool) -> str:
        """Return a :class:`Fold Record` as a Vienna-format string."""
        return write_record.to_vienna_string(newline=newline)


# Append common FoldFile Args and Attributes description to docstring:
ViennaFile.__doc__ += FOLD_FILE_COMMON_ARGS_ATTRS


# ----- Begin CtFile Class -----
class CtFile(FoldFile):
    """
    Ct file wrapper that returns ".ct" file lines as FoldRecord objects.

        See :ref:`CT File Format <ct_file_format>` for more information.

    Warning:
            This class is in beta stage, and is not well-tested.
    """

    # Args/Attrs Description in FOLD_FILE_COMMON_ARGS_ATTRS

    # Start CtFile Methods
    # CtFile : Public Methods
    def read_record(self) -> FoldReturn:
        """
        Return the next CT record as a :class:`FoldRecord` object.

        Call next(self.fh) to return the first line of the next entry.
        Determine the expected number of following lines in the entry, and read that number
        of lines further. Return lines as a FoldRecord object.
        """
        header = next(self.fh)
        record_lines = [header]
        expected_line_num = int(header.strip().split()[0])
        for _ in range(expected_line_num):
            record_lines.append(next(self.fh))
        record_type = self._foldrecord_types[self.settings['foldrecord_type']]
        record = record_type.from_ct_lines(
            record_lines,
            error_mode=self.error_mode)
        return record

    # CtFile : Disable Record Writing Methods
    #: CtFile Record Writing Not Implemented
    write_record = None

    # CtFile : Disable Record Writing Methods
    #: CtFile Record Writing Not Implemented
    write_records = None

    # CtFile : Private Methods
    def _post_init_tasks(self) -> None:
        """Hold place for any post-initialization tasks."""
        pass

    # CtFile : Private Methods
    def _to_record_string(self, write_record: FoldRecord, newline: bool) -> None:
        """(NOOP / Not-Implemented) Return a :class:`Fold Record` as a CT-format string."""
        message = 'No write_record is implemented for ct files, as the FoldRecord '
        message += 'object does not contain the complete set of ct record information.'
        raise NotImplementedError(message)


# Append common FoldFile Args and Attributes description to docstring:
CtFile.__doc__ += FOLD_FILE_COMMON_ARGS_ATTRS


# ----- Begin HybFoldIter Class -----
class HybFoldIter:
    """
    Iterator for simultaneous iteration over a :class:`HybFile` and ``FoldFile`` object.

    This class provides an iterator to iterate through a :class:`HybFile` and one of a
    :class:`ViennaFile`, or :class:`CtFile` simultaneously to return
    a :class:`HybRecord` and :class:`FoldRecord`.

    Basic error checking / catching is performed based on the value of the
    :attr:`~settings['error_mode'] <HybFoldIter.settings>` setting.

    Args:
        hybfile_handle (HybFile) : HybFile object for iteration
        foldfile_handle (:class:`ViennaFile` or :class:`CtFile`) : :class:`ViennaFile` or
            :class:`CtFile` object for iteration
        combine (:obj:`bool`, optional) : Use HybRecord.set_fold_record(FoldRecord)
            and return only the HybRecord.
        iter_error_mode (str, optional) : Error mode to use for reading :class:`FoldRecord`
            objects. If not set, defaults to the value in
            :attr:`settings['iter_error_mode'] <HybFoldIter.settings>`.

    Returns:
        (:class:`HybRecord`, :class:`FoldRecord`)
    """

    #: Class-level settings. See :attr:`settings.HybFoldIter_settings_info` for descriptions.
    settings = hybkit.settings.HybFoldIter_settings

    _iter_error_modes = frozenset(
        hybkit.settings.HybFoldIter_settings_info['iter_error_mode'][4]['choices'])

    # Start HybFoldIter Methods
    # HybFoldIter : Public Methods
    def __init__(
            self,
            hybfile_handle: HybFile,
            foldfile_handle: FoldFile,
            combine: bool = False,
            iter_error_mode: Optional[IterErrorModeArg] = None
            ) -> None:
        """Please see :class:`HybFoldIter` for initialization information."""
        if not isinstance(hybfile_handle, HybFile):
            message = 'hybfile_handle must be an instance of a HybFile object.'
            raise HybkitIterError(message)
        if not isinstance(foldfile_handle, FoldFile):
            message = 'foldfile_handle must be an instance of a FoldFile object.'
            raise HybkitIterError(message)
        if iter_error_mode is None:
            self.iter_error_mode = self.settings['iter_error_mode']
        elif iter_error_mode not in self._iter_error_modes:
            message = 'iter_error_mode must be one of the following: '
            message += ', '.join(self._iter_error_modes)
            raise HybkitIterError(message)
        else:
            self.iter_error_mode = iter_error_mode
        self.hybfile_handle = hybfile_handle
        self.foldfile_handle = foldfile_handle
        self.counters = Counter()
        self.combine = combine
        self.sequential_skips = 0
        self.last_hyb_record = None
        self.last_fold_record = None

    # HybFoldIter : Public Methods
    def report(self) -> List[str]:
        """Return a report of information from iteration."""
        ret_lines = ['HybFoldIter Iteration Report:']
        add_line = 'Combined Iteration Attempts: '
        add_line += str(self.counters['total_read_attempts'])
        ret_lines.append(add_line)
        add_line = 'Hyb Record Iteration Attempts: '
        add_line += str(self.counters['hyb_record_read_attempts'])
        ret_lines.append(add_line)
        add_line = 'Fold Record Iteration Attempts: '
        add_line += str(self.counters['fold_record_read_attempts'])
        ret_lines.append(add_line)
        # add_line = 'Total Skipped Fold-Only Records: ' + str(self.counters['fold_only_skips'])
        # ret_lines.append(add_line)
        add_line = 'Total Skipped Record Pairs: ' + str(self.counters['pair_skips'])
        ret_lines.append(add_line)
        return ret_lines

    # HybFoldIter : Public Methods
    def print_report(self) -> None:
        """Print a report of information from iteration."""
        ret_lines = self.report()
        print('\n'.join(ret_lines) + '\n')

    # HybFoldIter : Public Methods
    def __iter__(self) -> Self:
        """Return an iterator object."""
        return self

    # HybFoldIter : Public Methods
    def __next__(
        self
        ) -> Union[
            HybRecord,
            Tuple[HybRecord, FoldRecord],
            Tuple[HybRecord, FoldRecord, str],
            ]:
        """Read and return (:class:`HybRecord`, :class:`FoldRecord`)."""
        self.counters['total_read_attempts'] += 1
        next_hyb_record = None
        next_fold_record = None
        iter_error_mode = self.iter_error_mode
        try:
            # Attempt to read next HybRecord
            self.counters['hyb_record_read_attempts'] += 1
            next_hyb_record = self.hybfile_handle.read_record()
            # Attempt to read next FoldRecord
            self.counters['fold_record_read_attempts'] += 1
            next_fold_record = self.foldfile_handle.read_record(override_error_mode='return')
            # Set initial loop variables
            error = ''
            do_skip = False
            # Check for "NoFold" error
            if ('foldrecord_nofold' in self.settings['error_checks']
                    and isinstance(next_fold_record, tuple)
                    and next_fold_record[0] == 'NOFOLD'
                ):
                error = 'Improper FoldRecord: No Fold (Energy = 99*.*)'

            # Check for "NoEnergy" error
            if (not error and isinstance(next_fold_record, tuple)
                    and next_fold_record[0] == 'NOENERGY'):
                error = 'Improper FoldRecord: No Energy (no <Tab> in 3rd line)'

            # Check for "InDel" errors
            if (not error
                    and 'hybrecord_indel' in self.settings['error_checks']
                    and next_hyb_record.prop('has_indels')
                ):
                error = 'HybRecord: %s has InDels.' % str(next_hyb_record)

            # Check for "Mismatch" errors
            if not error and 'max_mismatch' in self.settings['error_checks']:
                hyb_fold_mismatches = next_fold_record.count_hyb_record_mismatches(next_hyb_record)
                if hyb_fold_mismatches > FoldRecord.settings['allowed_mismatches']:
                    error = 'HybRecord: %s ' % str(next_hyb_record)
                    error += 'has: %i ' % hyb_fold_mismatches
                    error += 'mismatches of '
                    error += '%i allowed ' % FoldRecord.settings['allowed_mismatches']

            # Check for "EnergyMismatch" errors
            if (not error
                    and 'energy_mismatch' in self.settings['error_checks']
                    and next_fold_record.energy is not None
                    and next_hyb_record.energy not in {None, '.'}
                    and str(next_fold_record.energy) != str(next_hyb_record.energy)
                ):
                    if not error:
                        error = 'HybRecord: %s ' % str(next_hyb_record)
                    error += 'has hyb-record / fold-record energy mismatch: '
                    error += '%s / ' % str(next_hyb_record.energy)
                    error += '%s\n' % str(next_fold_record.energy)
                    error += next_fold_record.to_vienna_string()

            # If an error exists, deal with it depending on the value of iter_error_mode
            if error:
                if iter_error_mode == 'raise':
                    raise HybkitIterError('ERROR: ' + error)
                elif iter_error_mode == 'warn_skip':
                    logging.warning(error)
                elif iter_error_mode == 'warn_return':
                    logging.warning(error)

                if 'skip' in iter_error_mode:
                    self.sequential_skips += 1
                    self.counters['pair_skips'] += 1
                    if self.sequential_skips > self.settings['max_sequential_skips']:
                        message = 'ERROR: Skipped %i ' % self.sequential_skips
                        message += 'record pairs in a row '
                        message += '(max: %i)\n' % self.settings['max_sequential_skips']
                        message += 'Check for misalignment of records, or disable setting.'
                        raise HybkitIterError(message)
                    do_skip = True

        except StopIteration:
            raise
        except BaseException:
            message = 'Error at Counter iteration: '
            message += '%s\n' % self.counters['total_read_attempts']
            message += 'Last HybRecord: %s\n' % str(self.last_hyb_record)
            message += 'Last FoldRecord: %s\n' % str(self.last_fold_record)
            if next_hyb_record is not None:
                message += 'Next HybRecord: %s\n' % str(next_hyb_record)
            if next_fold_record is not None:
                message += 'Next FoldRecord: %s\n' % str(next_fold_record)
            message += '\n' + '\n'.join(self.report()) + '\n'
            logging.warning(message)
            raise

        if do_skip:
            return next(self)

        if self.combine:
            next_hyb_record.set_fold_record(next_fold_record, allow_energy_mismatch=True)
            ret_obj = next_hyb_record
        else:
            ret_obj = (next_hyb_record, next_fold_record)

        if iter_error_mode == 'warn_return':
            if isinstance(ret_obj, tuple):
                ret_obj = (*ret_obj, error)
            else:
                ret_obj = (ret_obj, error)

        self.last_hyb_record = next_hyb_record
        self.last_fold_record = next_fold_record
        self.sequential_skips = 0
        return ret_obj

# Import the remainder of hybkit code to connect.
import hybkit.analysis
import hybkit.plot
import hybkit.util
