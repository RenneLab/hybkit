#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Analysis for sample_hyb_analysis performed as a python workflow, as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
'''

import os
import sys
import datetime

# Ensure hybkit is accessible
analysis_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(analysis_dir, '..'))
import hybkit

SHORT_CHECK = True # DEBUG
#SHORT_CHECK = False # DEBUG

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.join(analysis_dir, 'output')
input_files = [
    os.path.join(analysis_dir, 'GSM2720017_UI_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720018_UI_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720019_UI_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720020_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720021_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720022_WT_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720023_D11_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720024_D11_BR1.hyb'),
    os.path.join(analysis_dir, 'GSM2720025_D11_BR1.hyb')
]
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')

# Begin Analysis

print('\nPerforming Analysis')
start_time = datetime.datetime.now()

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Analyzing Files:')
print('    ' + '\n    '.join(input_files) + '\n')

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)
#hybkit.HybRecord.set_find_type_params(params)

# Initialize Analysis Dict Object
analysis_dict = hybkit.analysis.full_analysis_dict()

# Iterate over each input file, find the segment types, and save the output 
#   in the output directory.
for in_file_path in input_files:
    in_file_name = os.path.basename(in_file_path)
    out_file_name = in_file_name.replace('.hyb', '_typed.hyb')
    out_file_path = os.path.join(out_dir, out_file_name)

    print('Analyzing:\n    %s' % in_file_path)
    print('Outputting to:\n    %s\n' % out_file_path)

    # Open one HybFile entry for reading, and one for writing
    with hybkit.HybFile(in_file_path, 'r') as in_file, \
         hybkit.HybFile(out_file_path, 'w') as out_file:

        count = 0 # DEBUG

        # Iterate over each record of the input file
        for hyb_record in in_file:
            # Find the segments type of each record
            hyb_record.find_seg_types()

            if SHORT_CHECK: # DEBUG
                if count > 10000:
                    break
                count += 1

            # Perform record analysis
            hyb_record.mirna_analysis()

            # Add mirna_analysis details to mirna_analysis_dict
            hybkit.analysis.running_full(hyb_record, analysis_dict)
 
            # Write the modified record to the output file  
            out_file.write_record(hyb_record)

    # Write mirna_analysis for input file to outputs. 
    analysis_file_basename = out_file_path.replace('.hyb', '')
    print('Outputting Analyses to:\n    %s\n' % analysis_file_basename)
    hybkit.analysis.write_full(analysis_file_basename, analysis_dict, multi_files=True)

    sys.stdout.flush()  # DEBUG
    if SHORT_CHECK:
        break # DEBUG
 
print('Time taken: %s' % str(datetime.datetime.now() - start_time)) # DEBUG
sys.stdout.flush()  # DEBUG
