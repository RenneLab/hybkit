#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Analysis for summary_analysis pipeline performed as a python workflow.

Provided as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
See: "summary_analysis_notes.rst" for more information.
"""

import os
import sys
import hybkit

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)

hybkit.settings.Analysis_settings['count_mode'] = count_mode 

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.join(analysis_dir, 'output')
input_files = [
    os.path.join(analysis_dir, 'GSM2720017_UI_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720018_UI_BR2.hyb'),
    os.path.join(analysis_dir, 'GSM2720019_UI_BR3.hyb'),
    os.path.join(analysis_dir, 'GSM2720020_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720021_WT_BR2.hyb'),
    os.path.join(analysis_dir, 'GSM2720022_WT_BR3.hyb'),
    os.path.join(analysis_dir, 'GSM2720023_D11_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720024_D11_BR2.hyb'),
    os.path.join(analysis_dir, 'GSM2720025_D11_BR3.hyb')
]
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')

# Begin Analysis

print('\nPerforming QC & Summary Analysis...')

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Analyzing Files:')
print('    ' + '\n    '.join(input_files) + '\n')

# Set the method of finding segment type
match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(match_legend_file)
hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)

# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.HybFile.settings['hybformat_id'] = True
hybkit.HybFile.settings['hybformat_record'] = True

# Initialize listto store independent analyses.
summary_analyses = [] 

# Set hybrid segment types to remove as part of quality control (QC)
remove_types = ['rRNA', 'mitoch-rRNA']

# Iterate over each input file, find the segment types, and save the output 
#   in the output directory.
for in_file_path in input_files:
    in_file_name = os.path.basename(in_file_path)
    in_file_label = in_file_name.replace('.hyb', '')
    out_file_name = in_file_name.replace('.hyb', '_qc.hyb')
    out_file_path = os.path.join(out_dir, out_file_name)

    print('Analyzing:\n    %s' % in_file_path)
    print('Outputting to:\n    %s\n' % out_file_path)
    sys.stdout.flush()

    # Initialize file-specific analysis 
    file_summary_analysis = hybkit.analysis.SummaryAnalysis(name=in_file_label)

    # Open one HybFile entry for reading, and one for writing
    with hybkit.HybFile(in_file_path, 'r') as in_file, \
         hybkit.HybFile(out_file_path, 'w') as out_file:

        # Iterate over each record of the input file
        for hyb_record in in_file:
            # Find the segments type of each record
            hyb_record.eval_types()

            # Determine if record has type that is excluded
            use_record = True
            for remove_type in remove_types:
                if hyb_record.has_prop('any_seg_type_is', remove_type):
                    use_record = False
                    break

            # If record has an excluded type, continue to next record without analyzing.
            if not use_record:
                continue 

            # Perform record analysis
            hyb_record.eval_mirna()

            # Add record details to SummaryAnalysis
            file_summary_analysis.add(hyb_record)

            # Write the modified record to the output file  
            out_file.write_record(hyb_record)

    # Write analysis for input file to outputs. 
    analysis_file_basename = out_file_path.replace('.hyb', '')
    print('Outputting Analyses to:\n    %s\n' % analysis_file_basename)
    sys.stdout.flush()
    file_summary_analysis.write(analysis_file_basename) 
    file_summary_analysis.plot(analysis_file_basename) 
    summary_analyses.append(file_summary_analysis)

# Create a combined summary analysis
combined_summary_analysis = hybkit.analysis.SummaryAnalysis(name='combined_analysis')
for file_analysis in summary_analyses:
    combined_summary_analysis.update(file_analysis)

# Write combined summary analysis to files.
combined_analysis_file_basename = os.path.join(out_dir, 'combined_analysis')
print('Outputting Combined Analysis to:\n    %s\n' % combined_analysis_file_basename)
sys.stdout.flush()
combined_summary_analysis.write(combined_analysis_file_basename) 
combined_summary_analysis.plot(combined_analysis_file_basename) 

print('Done!')  
