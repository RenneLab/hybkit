
hybkit Specification
==================================

hybkit Specification, Version: |version|

Hybkit uses an extended version of the specification provided for ".hyb" format-files by 

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Flags
-----

    .. _count_total:
    
    'count_total',            # str(int), total represented hybrids

    .. _count_last_clustering:

    'count_last_clustering',  # str(int), total represented hybrids at last clustring
                 
    .. _two_way_merged:
    
    'two_way_merged',         # "0" or "1", boolean representation of whether
    #   entries with mirrored 5' and 3' hybrids were merged
    
    .. _seq_IDs_in_cluster:

    'seq_IDs_in_cluster',     # str, comma-separated list of all ids of hybrids
    #   merged into this hybrid entry.

    .. _read_count:

    read_count,   # str(int), number of sequence reads represented by record
    #   if merged record, represents total for all merged entries

    .. _orient:

    'orient',       # str, orientation of strand. Options:
    #   "F" (Forward), "IF" (Inferred Forward),
    #   "R" (Reverse), "IR" (Inferred Reverse),
    #   "U" (Unknown), or "IC" (Inferred Conflicting)

    .. _seg1_type:

    :obj:`seg1_type`    # str, assigned type of segment 1, ex: "miRNA" or "mRNA"

    .. _seg2_type:

    'seg2_type',    # str, assigned type of segment 2, ex: "miRNA" or "mRNA"

    .. _seg1_det:
 
    'seg1_det',     # str, arbitrary detail about segment 1

    .. _seg2_det:

    'seg2_det',     # str, arbitrary detail about segment 2

    .. _miRNA_seg:

    'miRNA_seg',    # str, indicates which (if any) segment mapping is a miRNA
    #   options are "N" (none), "3p" (seg1), "5p" (seg2),
    #   "B" (both), or "U" (unknown)

    .. _target_reg:

    'target_reg',   # str, assigned region of the miRNA target.
    #   options are "5pUTR", "coding", "3pUTR",
    #   "N" (none), or "U" (unknown)

    .. _ext:
  
    'ext',          # int, "0" or "1", boolean representation of whether
    #   record sequences were bioinformatically extended as is
    #   performed by the Hyb software package.

    .. _source:

    'source',       # str, label for sequence source id (eg. file), when 
    #   combining records from different sources.


