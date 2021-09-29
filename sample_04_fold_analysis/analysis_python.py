#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Analysis for sample_fold_analysis performed as a python workflow.

Provided as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
For further detail, see "fold_analysis_notes.rst".
"""

import os
import sys
import hybkit

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
input_hyb_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb')
input_viennad_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.viennad')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')
out_dir = os.path.join(analysis_dir, 'output')
out_hyb_name = os.path.join(out_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua_coding.hyb')
data_label = 'WT_BR1'
out_analysis_basename = out_hyb_name.replace('.hyb', '')

# Begin Analysis
print('\nPerforming Fold Analysis...')
if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Using Input Files:')
print('    ' + '\n    '.join([input_hyb_name, input_viennad_name]) + '\n')

# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.HybFile.settings['hybformat_id'] = True

# Tell the FoldRecord to allow (by skipping) poorly-formatted viennad entries, instead of 
#   raising an error.
#hybkit.FoldRecord.settings['skip_bad'] = True

# Create a variable mirna-types for use in the miRNA analysis, that includes kshv mirna.
mirna_types = list(hybkit.HybRecord.MIRNA_TYPES) + ['kshv_microRNA']

# Set hybrid segment types to remove as part of quality control (QC)
remove_types = ['rRNA', 'mitoch_rRNA']

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)

# Prepare fold_analysis dict:
analysis_dict = hybkit.analysis.mirna_fold_dict()

# Use the combined iterator to iterate over the hyb and viennad files simultaneously, 
#   returning hyb records containing their associated fold record.
#   NOTE: This dataset has been checked to make sure each viennad-record has 6 lines,
#     And matches the corresponding hyb record line. If this is not the case, this 
#     Analysis may fail or return poor results.
in_file_label = os.path.basename(input_hyb_name).replace('.hyb', '')
with hybkit.HybFile.open(input_hyb_name, 'r') as input_hyb,\
     hybkit.ViennadFile.open(input_viennad_name, 'r') as input_viennad,\
     hybkit.HybFile.open(out_hyb_name, 'w') as out_hyb:

    for hyb_record in hybkit.HybFoldIter(input_hyb, input_viennad, combine=True):
        # Find Segment types
        hyb_record.find_seg_types()

        # Determine if record has type that is excluded
        use_record = True
        for remove_type in remove_types:
            if hyb_record.has_property('seg_type', remove_type):
                use_record = False
                break

        # If record has an excluded type, continue to next record without analyzing.
        if not use_record:
            continue

        # Perform record analysis
        hyb_record.mirna_analysis(mirna_types=mirna_types)

        # If the record contains a non-duplex miRNA, then analyze folding.
        # Equivalent to 'has_mirna' and not 'has_mirna_dimer'
        if hyb_record.has_property('has_mirna_not_dimer'):
            # Perform miRNA-Fold Analysis
            hybkit.analysis.addto_mirna_fold(hyb_record, 
                                             analysis_dict,
                                             skip_no_fold_record=True)

            # Write the record to the output hyb file.
            out_hyb.write_record(hyb_record)


# Write mirna_fold analysis for input file to outputs.
print('Outputting Analyses to:\n    %s\n' % out_analysis_basename)
analysis_dict = hybkit.analysis.process_mirna_fold(analysis_dict)
hybkit.analysis.write_mirna_fold(out_analysis_basename,
                                 analysis_dict,
                                 multi_files=True,
                                 name=data_label,
                                 )
                            
print('Done!') 
