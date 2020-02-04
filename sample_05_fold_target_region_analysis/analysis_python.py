#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Analysis for sample_fold_analysis performed as a python workflow, as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
'''

import os
import sys
import datetime

# Ensure hybkit is accessible
analysis_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(analysis_dir, '..'))
import hybkit

SHORT_CHECK = True  # DEBUG
#SHORT_CHECK = False  # DEBUG

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
input_hyb_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.hyb')
input_viennad_name = os.path.join(analysis_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua.viennad')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')
region_info_csv = os.path.join(hybkit.package_dir, 'databases', 'Human_mRNAs_mod.csv')
out_dir = os.path.join(analysis_dir, 'output')
out_hyb_name = os.path.join(out_dir, 'WT_BR1_comp_hOH7_KSHV_hybrids_ua_coding.hyb')

# Begin Analysis
print('\nPerforming Analysis')
start_time = datetime.datetime.now()  # DEBUG

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Analyzing Files:')
print('    ' + '\n    '.join([input_hyb_name, input_viennad_name]) + '\n')

# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.HybFile.settings['hybformat_id'] = True

# Tell the FoldRecord to allow (by skipping) poorly-formatted viennad entries, instead of 
#   raising an error.
hybkit.FoldRecord.settings['skip_bad'] = True

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)
#hybkit.HybRecord.set_find_type_params(params)

# Read the csv file containing coding sequence region info:
hybkit.HybRecord.make_set_region_info(region_info_csv)

count = 0  # DEBUG

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
        hyb_record.mirna_analysis()
        hyb_record.target_region_analysis()
        if hyb_record.has_property('has_target'):
            print(hyb_record.flags['target_reg'])
            out_hyb.write_record(hyb_record)            

        

sys.exit()
# Write mirna_analysis for input file to outputs. 
analysis_file_basename = out_file_path.replace('.hyb', '')
print('Outputting Analyses to:\n    %s\n' % analysis_file_basename)
hybkit.analysis.write_full(analysis_file_basename, 
                           analysis_dict_by_record, 
                           multi_files=True, 
                           name=in_file_label)

analysis_dicts_by_record.append(analysis_dict_by_record)

sys.stdout.flush()  # DEBUG

combined_analysis_dict_by_record = hybkit.analysis.combine_full_dicts(analysis_dicts_by_record) 
combined_analysis_file_basename = os.path.join(out_dir, 'combined_analysis')
print('Outputting Combined Analysis to:\n    %s\n' % combined_analysis_file_basename)
hybkit.analysis.write_full(combined_analysis_file_basename, 
                       combined_analysis_dict_by_record, 
                       multi_files=True,
                       name='Combined')

print('Time taken: %s' % str(datetime.datetime.now() - start_time)) # DEBUG
sys.stdout.flush()  # DEBUG
