#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Example fold analysis performed as a python workflow.

Provided as an example of direct usage of hybkit functions.

See: 'README.rst' for this analysis for more information.
"""

import os

import hybkit

# Set mirna types as custom to include KSHV-miRNAs
hybkit.util.set_setting('mirna_types', set_value=['miRNA', 'KSHV-miRNA'])
# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.util.set_setting('hybformat_id', set_value=True)
# Set FoldRecord seq_type as "dynamic" to work with Hyb-format *_hybrid_ua.vienna files
#   which have a modified sequence.
hybkit.util.set_setting('seq_type', 'dynamic')

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
input_hyb_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb')
input_vienna_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.vienna')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')
out_dir = os.path.join(analysis_dir, 'output_python')
out_hyb_name = os.path.join(out_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua_qc.hyb')
data_label = 'WT_BR1'
analysis_basename = out_hyb_name.replace('.hyb', '')

# Begin Analysis
print('\nPerforming Pattern Analysis...')
if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Using Input Files:')
print(f'    {input_hyb_name}\n    {input_vienna_name}\n')

# Set hybrid segment types to remove as part of quality control (QC)
remove_types = ['rRNA', 'mitoch-rRNA']

# Set the method of finding segment type
match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(match_legend_file)
hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)

# Initialize Fold Analysis:
hyb_analysis = hybkit.analysis.Analysis(analysis_types=['fold', 'energy'], name='WT_BR1')

# Use the combined iterator to iterate over the hyb and vienna files simultaneously,
#   returning hyb records containing their associated fold record.
in_file_label = os.path.basename(input_hyb_name).replace('.hyb', '')
with hybkit.HybFile.open(input_hyb_name, 'r') as input_hyb,\
     hybkit.ViennaFile.open(input_vienna_name, 'r') as input_vienna:
    # hybkit.HybFile.open(out_hyb_name, 'w') as out_hyb:

    hyb_fold_iter = hybkit.HybFoldIter(input_hyb, input_vienna, combine=True)
    for _, hyb_record in enumerate(hyb_fold_iter):
        # Find Segment types and miRNA information
        hyb_record.eval_types()
        hyb_record.eval_mirna()

        # Set dataset flag of record
        hyb_record.set_flag('dataset', in_file_label)

        # Determine if record has type that is excluded
        use_record = True
        for remove_type in remove_types:
            if hyb_record.has_prop('any_seg_type_is', remove_type):
                use_record = False
                break

        # If record has an excluded type, continue to next record without analyzing.
        if not use_record:
            continue

        # If the record contains a non-duplex miRNA, then analyze folding.
        # Equivalent to ('has_mirna' and not 'has_mirna_dimer')
        if hyb_record.has_prop('mirna_not_dimer'):
            # Perform miRNA Fold Analysis
            hyb_analysis.add_hyb_record(hyb_record)
            # Write the record to the output hyb file.
            # out_hyb.write_record(hyb_record)

# Print iterator report after Iteration
print()
hyb_fold_iter.print_report()
print()

# Write pattern analysis for input file to outputs.
print('Outputting Analyses to:\n    %s\n' % analysis_basename)
hyb_analysis.write_analysis_results_special(out_basename=analysis_basename)
hyb_analysis.plot_analysis_results(out_basename=analysis_basename)

print('\nDone!\n')
