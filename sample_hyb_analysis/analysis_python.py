#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Analysis for Sample_analysis performed as a python workflow, as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
'''

import os
import sys
import hybkit
import datetime

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

print('Performing Analysis\n')
print('Starting at: %s' % str(datetime.datetime.now()))
print('Creating Output Directory:\n    %s\n' % out_dir)
os.mkdir(out_dir)

print('Analyzing Files:')
print('    ' + '\n    '.join(input_files) + '\n')

# Set the method of finding segment type
match_parameters = hybkit.HybRecord.make_string_match_parameters(match_legend_file)
hybkit.HybRecord.select_find_type_method('string_match', match_parameters)
#hybkit.HybRecord.set_find_type_params(params)

# Create an empty list to store output file names:
out_file_paths = []

# Initialize Analysis Information Dictionaries
hyb_analysis_dict = hybkit.HybRecord.running_hyb_analysis_dict()
mirna_analysis_dict = hybkit.HybRecord.running_mirna_analysis_dict()

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

        # Iterate over each record of the input file
        for hyb_record in in_file:
            # Find the segments type of each record
            hyb_record.find_seg_types()

            # Add seg_type analysis details to hyb_analysis_dict
            hyb_record.running_hyb_analysis(hyb_analysis_dict)

            # Perform microRNA analysis
            hyb_record.mirna_analysis()

            # Add mirna_analysis details to mirna_analysis_dict
            hyb_record.running_mirna_analysis(mirna_analysis_dict)
 
            # Write the modified record to the output file  
            out_file.write_record(hyb_record)

            #print(hyb_analysis_dict)
            #print(mirna_analysis_dict)
            #input('\nContinue?...\n')

    # Add the completed output file to the list of output files:
    out_file_paths.append(out_file_path)
 
    print('HybRecord Analysis:')
    for key in hyb_analysis_dict.keys():
        print(key, hyb_analysis_dict[key])
    print('\nmiRNA Analysis:')
    for key in mirna_analysis_dict.keys():
        print(key, mirna_analysis_dict[key])
    break
         
print('Ending At: %s' % str(datetime.datetime.now()))
sys.stdout.flush()
