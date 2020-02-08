#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Analysis for sample_fold_analysis performed as a python workflow.

Provided as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
"""

import os
import sys
import datetime

# Ensure hybkit is accessible
analysis_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(analysis_dir, '..'))
import hybkit

SHORT_CHECK = True  # DEBUG
SHORT_CHECK = False  # DEBUG

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
start_time = datetime.datetime.now()  # DEBUG

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Using Input Files:')
print('    ' + '\n    '.join([input_hyb_name, input_viennad_name]) + '\n')

# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.HybFile.settings['hybformat_id'] = True

# Tell the FoldRecord to allow (by skipping) poorly-formatted viennad entries, instead of 
#   raising an error.
hybkit.FoldRecord.settings['skip_bad'] = True

# Create a variable mirna-types for use in the miRNA analysis, that includes kshv mirna.
mirna_types = list(hybkit.HybRecord.MIRNA_TYPES) + ['kshv_microRNA']

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)
#hybkit.HybRecord.set_find_type_params(params)

count = 0  # DEBUG

# Prepare fold_analysis dict:
analysis_dict = hybkit.analysis.mirna_fold_dict()

# Use the combined iterator to iterate over the hyb and viennad files simultaneously, 
#   returning hyb records containing their associated fold record.
in_file_label = os.path.basename(input_hyb_name).replace('.hyb', '')
with hybkit.HybFile.open(input_hyb_name, 'r') as input_hyb,\
     hybkit.ViennadFile.open(input_viennad_name, 'r') as input_viennad,\
     hybkit.HybFile.open(out_hyb_name, 'w') as out_hyb:
    combined_iter = hybkit.HybViennadCmbIter(input_hyb, input_viennad)
    for hyb_record in combined_iter:
        #print(hyb_record)

        hyb_record.find_seg_types()

        if SHORT_CHECK: # DEBUG
            if count > 10000:
                #break
                break
            count += 1

        # Perform record analysis
        hyb_record.mirna_analysis(mirna_types=mirna_types)

        # Equivalent to 'has_mirna' and not 'has_mirna_dimer'
        if hyb_record.has_property('has_mirna_not_dimer'):
            hybkit.analysis.running_mirna_folds(hyb_record, 
                                                analysis_dict,
                                                skip_no_fold_record=True)
            out_hyb.write_record(hyb_record)


# Write mirna_fold analysis for input file to outputs.
print('Outputting Analyses to:\n    %s\n' % out_analysis_basename)
analysis_dict = hybkit.analysis.process_mirna_folds(analysis_dict)
hybkit.analysis.write_mirna_folds(out_analysis_basename,
                                  analysis_dict,
                                  multi_files=True,
                                  name=data_label,
                                 )
                             
print('Time taken: %s\n' % str(datetime.datetime.now() - start_time)) # DEBUG
sys.stdout.flush()  # DEBUG
