hybkit.analysis Module
======================

.. automodule:: hybkit.analysis


Type Analysis
-------------

The type analysis provides an analysis of what segment types 
are included in the analyzed hyb files. 

Before using the analysis, the :ref:`seg1_type <seg1_type>` and 
:ref:`seg2_type <seg2_type>` flags must be set for the record, 
as is done by :func:`HybRecord.find_seg_types`. 
A count is added to 
the analysis dict for each hybrid type (Ex: "miRNA-mRNA") with 
segments placed in sorted order for non-redundant typing.
The analysis additionally reports the number of individual segment
types.

.. autofunction:: type_dict
.. autofunction:: combine_type_dicts
.. autofunction:: addto_type
.. autofunction:: format_type
.. autofunction:: write_type


miRNA Count Analysis
--------------------

The mirna_count analysis determines what type each record is 
with regard to mirna and counts them accordingly.
This includes:

    5p_mirna_hybrids: Hybrids with a |5p| miRNA.
    3p_mirna_hybrids: Hybrids with a |3p| miRNA.
    mirna_dimer_hybrids: Hybrids with both a |5p| and |3p| miRNA.
    no_mirna_hybrids: Hybrids with no miRNA.
    (And additionally includes:)
    all_mirna_hybrids: Hybrids that fall into the first three categories.

Before using the analysis, the :ref:`mirna_seg <mirna_seg>` flag 
must be set for each record as can be done by sequential use of the 
:func:`HybRecord.find_seg_types` and :func:`HybRecord.mirna_analysis` 
methods.

.. autofunction:: mirna_count_dict
.. autofunction:: combine_mirna_count_dicts
.. autofunction:: addto_mirna_count
.. autofunction:: format_mirna_count
.. autofunction:: write_mirna_count


Summary Analysis
----------------

.. autofunction:: summary_dict
.. autofunction:: combine_summary_dicts
.. autofunction:: addto_summary
.. autofunction:: write_summary


miRNA Target Analysis
---------------------

.. autofunction:: mirna_target_dict
.. autofunction:: combine_mirna_target_dicts
.. autofunction:: addto_mirna_target
.. autofunction:: process_mirna_target
.. autofunction:: format_mirna_target
.. autofunction:: write_mirna_target


miRNA Fold Analysis
-------------------

.. autofunction:: mirna_fold_dict
.. autofunction:: combine_mirna_fold_dicts
.. autofunction:: addto_mirna_fold
.. autofunction:: process_mirna_fold
.. autofunction:: format_mirna_fold
.. autofunction:: write_mirna_fold

