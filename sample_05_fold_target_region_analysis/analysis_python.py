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
region_info_csv = os.path.join(hybkit.__about__.data_dir, 'hybkit_coding_ref_combined.csv')
out_dir = os.path.join(analysis_dir, 'output')
data_label = 'WT_BR1'
out_base_name = 'WT_BR1_comp_hOH7_KSHV'
out_base = os.path.join(out_dir, out_base_name)

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

# Create a variable mirna-types for use in the miRNA analysis, that includes kshv mirna.
mirna_types = list(hybkit.HybRecord.MIRNA_TYPES) + ['kshv_microRNA']

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)
#hybkit.HybRecord.set_find_type_params(params)

# Read the csv file containing coding sequence region info:
hybkit.HybRecord.make_set_region_info(region_info_csv)

count = 0  # DEBUG

# Prepare 5 catgories for 2 classes in output_categories dict
output_classes = {'cellular': 'Cellular', 'kshv': 'KSHV'}
record_types = {'5pUTR': '5pUTR',
                'coding': 'Coding', 
                '3pUTR': '3pUTR',
                'unknown': 'Unknown',
                'noncoding': 'Noncoding'}
output_categories = {}
for out_class in output_classes.keys():
    output_categories.update({(out_class + '_' + rec_type):[] for rec_type in record_types})

# Create an analysis dict and open an output file for each category.
for cat in output_categories.keys():
    output_categories[cat].append(hybkit.analysis.mirna_fold_dict())
    file_name = out_base + '_' + cat + '.hyb'
    output_categories[cat].append(hybkit.HybFile.open(file_name, 'w'))
    cat_split = cat.split('_')
    pretty_name = output_classes[cat_split[0]] + ' ' + record_types[cat_split[1]]
    output_categories[cat].append(pretty_name)

# Use the combined iterator to iterate over the hyb and viennad files simultaneously, 
#   returning hyb records containing their associated fold record.
in_file_label = os.path.basename(input_hyb_name).replace('.hyb', '')
with hybkit.HybFile.open(input_hyb_name, 'r') as input_hyb,\
     hybkit.ViennadFile.open(input_viennad_name, 'r') as input_viennad:

    for hyb_record in hybkit.HybFoldIter(input_hyb, input_viennad, combine=True):
        if SHORT_CHECK: # DEBUG
            if count > 10000:
                break
            count += 1

        # Perform record analysis
        hyb_record.find_seg_types()
        hyb_record.mirna_analysis(mirna_types=mirna_types)

        # Find miRNA-containing records.
        # Equivalent to 'has_mirna' and not 'has_mirna_dimer'
        if hyb_record.has_property('has_mirna_not_dimer'):
            # Assign as KSHV miRNA or Cellular miRNA
            if hyb_record.has_property('seg_contains', 'kshv'):
                label_prefix = 'kshv_'
            else:
                label_prefix = 'cellular_'

            hyb_record.target_region_analysis(allow_unknown_regions=True)
            if hyb_record.has_property('has_target'):
                # Set coding target region labels
                region = hyb_record.mirna_details['target_reg']
                if region not in {'5pUTR', 'coding', '3pUTR', 'unknown'}:
                    raise Exception('Unknown Region: %s' % region) 
                record_params = output_categories[label_prefix + region]
                [record_analysis_dict, record_out_file, pretty_name] = record_params
            else:
                category = label_prefix + 'noncoding'
                [record_analysis_dict, record_out_file, pretty_name] = output_categories[category]

            hybkit.analysis.running_mirna_folds(hyb_record, 
                                                record_analysis_dict,
                                                skip_no_fold_record=True)
            record_out_file.write_record(hyb_record)

# Write mirna_fold analysis for each catetory.
print('\nOutputting Analyses with prefix:\n    %s' % out_base)
for category in output_categories.keys():
    print('    Writing analyses for %s.' % category)
    analysis_name = out_base + '_' + category
    analysis_dict, out_hyb_file, pretty_name = output_categories[category]
    analysis_dict = hybkit.analysis.process_mirna_folds(analysis_dict)
    hybkit.analysis.write_mirna_folds(analysis_name,
                                      analysis_dict,
                                      multi_files=True,
                                      name=data_label + ', ' + pretty_name
                                     )
                             
print('Time taken: %s\n' % str(datetime.datetime.now() - start_time)) # DEBUG
sys.stdout.flush()  # DEBUG
