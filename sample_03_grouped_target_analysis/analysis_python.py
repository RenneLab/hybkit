#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Analysis for sample_grouped_target_analysis performed as a python workflow.

Provided as an example of direct usage of hybkit functions. 
File names are hardcoded, and functions are accessed directly.
"""

import os
import sys
import datetime

# Ensure hybkit is accessible
analysis_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(analysis_dir, '..'))
import hybkit

SHORT_CHECK = True   # DEBUG
SHORT_CHECK = False  # DEBUG
ONE_CHECK = True     # DEBUG
ONE_CHECK = False    # DEBUG

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
analysis_label = 'KSHV_Hyb_Combined'
out_dir = os.path.join(analysis_dir, 'output')
input_files = [
    os.path.join(analysis_dir, 'GSM2720020_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720021_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720022_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720023_D11_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720024_D11_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720025_D11_BR1.hyb')
]
out_file_path = os.path.join(analysis_dir, 'output', 'KSHV_Hyb_Combined.hyb')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')

# Begin Analysis
print('\nPerforming Grouped Target Analysis...')
start_time = datetime.datetime.now()  # DEBUG

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Analyzing Files:')
print('    ' + '\n    '.join(input_files) + '\n')

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)

# Create custom list of miRNA types for analysis
mirna_types = list(hybkit.HybRecord.MIRNA_TYPES) + ['kshv_microRNA']

print('Outputting KSHV-Specific Hybrids to:\n    %s\n' % out_file_path)

# Open out-file for writing, one output file will be used for all inputs.
with hybkit.HybFile(out_file_path, 'w') as out_kshv_file:
    # Iterate through input files:
    for in_file_path in input_files:
        in_file_name = os.path.basename(in_file_path)
        in_file_label = in_file_name.split('_')[0]
        print('Analyzing:\n    %s' % in_file_path)
        with hybkit.HybFile(in_file_path, 'r') as in_file:
        
            count = 0 # DEBUG

            # Iterate over each record of the input file
            for hyb_record in in_file:

                if SHORT_CHECK: # DEBUG
                    count += 1
                    if count > 10000:
                        break

                # Analyze and output only sequences where a segment identifier contains the string "kshv"
                if hyb_record.has_property('seg_contains', 'kshv'):
                    # Find the segment types of each record
                    hyb_record.set_flag('source', in_file_label)
                    hyb_record.find_seg_types()
                    hyb_record.mirna_analysis(mirna_types=mirna_types)
        
                    # Write the records to the output file
                    out_kshv_file.write_record(hyb_record)

        if ONE_CHECK:  # DEBUG
            break

sys.stdout.flush()  # DEBUG
print('\nPerforming Combined Target Analysis...\n')
# Reiterate over output HybFile and perform target analysis.
with hybkit.HybFile(out_file_path, 'r') as out_kshv_file:
    # Prepare target analysis dict:
    target_dict = {}

    for hyb_record in out_kshv_file:
        # Repeat .mirna_analysis, so that hyb_record object has associated metadata
        hyb_record.mirna_analysis(mirna_types=mirna_types)

        # Perform target-analysis of mirna within kshv-associated data.
        hybkit.analysis.running_mirna_targets(hyb_record, target_dict, 
                                              count_mode=count_mode,
                                              double_count_duplexes=True, # Includes mirna duplexes
                                              # Limits output to KSHV miRNA:
                                              mirna_contains='kshv')
        
    # Process and sort dictionary of miRNA and targets
    results = hybkit.analysis.process_mirna_targets(target_dict)
    (sorted_target_dict,       # Contains same data as target_dict, but with keys sorted by count
     counts_dict,              # dict with keys: mirna_id, values: total mirna-specific hybrids
     target_type_counts_dict,  # dict with keys: mirna_id, values: dict of targeted type counts
     total_count               # integer: total number of mirna found.
     ) = results

    # Write target information to output file
    # Set analysis basename without ".hyb" extension
    analysis_basename = out_file_path.replace('.hyb','')
    print('Writing Analyses to files named after:\n    %s\n' % analysis_basename)

    hybkit.analysis.write_mirna_targets(analysis_basename, 
                                        sorted_target_dict,
                                        counts_dict,
                                        target_type_counts_dict,
                                        name=analysis_label,
                                        multi_files=True,    # Default
                                        sep=',',             # Default
                                        file_suffix='.csv',  # Default
                                        spacer_line=True,    # Default
                                        make_plots=True,     # Default
                                        max_mirna=25)        # Default is 10

print('Total time: %s' % str(datetime.datetime.now() - start_time))  # DEBUG

print('Done\n')

sys.stdout.flush()  # DEBUG
