
hybkit Hyb File Specification
==================================

    Version: |spec_version|

    The ".hyb" (hyb file) format is described by [Travis2014]_ along with the Hyb software package
    as a "gff-related format that contains sequence identifiers, read sequences, 1-based
    mapping coordinates, and annotation information for each chimera".

    Each line in a hyb file (a hyb "record") contains information about an
    RNA sequence read identified as a chimera by an RNA hybridization analysis.
    Each line contains 15 or 16 columns separated by tab characters ("\\t") and provides
    information on each of the two aligned segments identified within the sequence read.
    The columns are described as follows by [Travis2014]_:

        | Column 1, unique sequence identifier.
        | Column 2, read sequence [...].
        | Column 3, predicted binding energy in kcal/mol.
        | Columns 4–9, mapping information for first fragment of read:
          name of matched transcript, coordinates in
          read, coordinates in transcript, mapping score.
        | Columns 10–15, mapping information for second fragment of read.
        | Column 16 (optional, [...]), list of annotations in the format:
          ‘‘feature1=value1; feature2=value2;..."

    The hybkit project uses an extended version of this description, including assigning
    columns reference names, and defining allowed flags.


Columns
-------------

    .. table::
        :widths: auto

        == =============== ===================================================================
        #  Name            Description
        == =============== ===================================================================
        1  id              Hybrid Read Identifier
        2  seq             Read Nucleotide Sequence
        3  energy          Predicted Gibbs Free-Energy of Intra-Hybrid Folding
        4  seg1_ref_name   Segment 1 Mapping Reference Identity
        5  seg1_read_start Segment 1 Mapping Start on Read
        6  seg1_read_end   Segment 1 Mapping End on Read
        7  seg1_ref_start  Segment 1 Mapping Start on Reference
        8  seg1_ref_end    Segment 1 Mapping End on Reference
        9  seg1_score      Segment 1 Mapping Score
        10 seg2_ref_name   Segment 2 Mapping Reference Identity
        11 seg2_read_start Segment 2 Mapping Start on Read
        12 seg2_read_end   Segment 2 Mapping End on Read
        13 seg2_ref_start  Segment 2 Mapping Start on Reference
        14 seg2_ref_end    Segment 2 Mapping End on Reference
        15 seg2_score      Segment 2 Mapping Score
        16 flags           Hybrid Read Analysis Details
        == =============== ===================================================================


.. toctree::
   :maxdepth: 2
   :caption: Contents:


.. _Flags:

Flags
-----

    Hyb Flags:
        The following four flags are used by the Hyb software package ([Travis2014]_).
        The definitions provided describe how these flags are used in the hybkit package.

        .. _count_total:

        :obj:`count_total` - Integer: Total represented hybrid records, if records have been
        combined.

        .. _count_last_clustering:

        :obj:`count_last_clustering` - Integer: Total represented hybrid records at last
        clustering.

        .. _two_way_merged:

        :obj:`two_way_merged` - {"0" or "1"} Boolean representation of whether
        entries with mirrored 5' and 3' hybrids were merged if the record is a combined record.

        .. _seq_IDs_in_cluster:

        :obj:`seq_IDs_in_cluster` -  String: Comma-separated list of all record IDs of hybrids
        merged into this hybrid entry.

    hybkit Flags:
        The following flags are used by hybkit.

        .. _read_count:

        :obj:`read_count` -  Integer: Number of sequence reads represented by this record.
        If the record is combined, this represents the total read count for all merged entries.

        .. _orient:

        :obj:`orient` -  String: Orientation of strand. Options:
        "F" (Forward), "IF" (Inferred Forward),
        "R" (Reverse), "IR" (Inferred Reverse),
        "U" (Unknown), or "IC" (Inferred Conflicting).

        .. _seg1_type:

        :obj:`seg1_type` - String: Assigned segment type of segment 1, ex: "miRNA" or "mRNA".

        .. _seg2_type:

        :obj:`seg2_type` - String: Assigned segment type of segment 2, ex: "miRNA" or "mRNA".

        .. _seg1_det:

        :obj:`seg1_det` -  String: Arbitrary detail about segment 1.

        .. _seg2_det:

        :obj:`seg2_det` -  String: Arbitrary detail about segment 2.

        .. _mirna_seg:

        :obj:`miRNA_seg` -  String: Indicates which (if any) segment mapping is a miRNA.
        Options are "N" (none), "5p" (seg1), "3p" (seg2),
        "B" (both), or "U" (unknown).

        .. _target_reg:

        :obj:`target_reg` -  String: Assigned region of the miRNA target.
        Options are "5pUTR", "C" ([C]oding), "3pUTR",
        "NON" ([NON]coding), "N" ([N]one), or "U" ([U]nknown).

        .. _ext:

        :obj:`ext` -  Integer: "0" or "1", Boolean representation of whether
        record sequences were bioinformatically extended as is
        performed by the Hyb software package.

        .. _dataset:

        :obj:`dataset` -  String: Label for sequence dataset id (eg. source file), when
        combining records from different datasets.



Other Details
-------------

    .. table::
        :widths: auto

        ============ ==================================================================================
        Item         Role
        ============ ==================================================================================
        "\\t" (tab)  Column Delimiter
        "."          Missing Data Placeholder (equivalent to None)
        ".hyb"       File Suffix
        ".hyb.gz"    gzipped File Suffix
        Disallowed   Header Lines
        Disallowed   In-file Comments
        ============ ==================================================================================


Example
-------

An example .hyb format line (courtesy of [Gay2018]_)::

    2407_718	ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC	.	MIMAT0000078_MirBase_miR-23a_microRNA	1	21	1	21	0.0027	ENSG00000188229_ENST00000340384_TUBB2C_mRNA	23	49	1181	1207	1.2e-06
