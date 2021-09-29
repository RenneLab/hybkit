#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Analysis for sample_hyb_analysis performed as a python workflow.

Provided as an example of direct 
usage of hybkit functions. File names are hardcoded, and functions are accessed directly.
See "target_analysis_notes.rst" for more information.
"""

import os
import sys
import hybkit

# Set count_mode:
# count_mode = 'read'    # Count reads represented by each record, instead of number of records.
count_mode = 'record'  # Count each record/line as one, unless record is combined.
                       #   (Default count mode, but specified here for readability)
hybkit.settings.Analysis_settings['count_mode'] = count_mode 

# Set mirna types as custom to include KSHV-miRNAs
hybkit.settings.HybRecord_settings['mirna_types'] = ['miRNA', 'microRNA', 'kshv-miRNA']

# Allow mirna/mirna dimers in analysis.
hybkit.settings.Analysis_settings['allow_mirna_dimers'] = True

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.join(analysis_dir, 'output')

in_file_label = 'GSM2720020_WT_BR1'
in_file_path = os.path.join(analysis_dir, 'GSM2720020_WT_BR1.hyb')
out_file_path = os.path.join(analysis_dir, 'output', 'GSM2720020_WT_BR1_KSHV_only.hyb')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')

# Begin Analysis
print('\nPerforming Target Analysis...')

if not os.path.isdir(out_dir):
    print('Creating Output Directory:\n    %s\n' % out_dir)
    os.mkdir(out_dir)

print('Analyzing File:')
print('    ' + in_file_path + '\n')

# Set the method of finding segment type
match_params = hybkit.HybRecord.TypeFinder.make_string_match_params(match_legend_file)
hybkit.HybRecord.TypeFinder.set_method('string_match', match_params)

# Set hybrid segment types to remove as part of quality control (QC)
remove_types = ['rRNA', 'mitoch-rRNA']

# Begin Analysis
print('Outputting KSHV-Specific Hybrids to:\n    %s\n' % out_file_path)

# Open one HybFile entry for reading, and one for writing
with hybkit.HybFile(in_file_path, 'r') as in_file, \
     hybkit.HybFile(out_file_path, 'w') as out_kshv_file:

    # Iterate over each record of the input file
    for hyb_record in in_file:

        # Analyze and output only sequences where a segment identifier contains the string "kshv"
        if hyb_record.has_prop('any_seg_contains', 'kshv'):
            # Find the segment types of each record
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

            hyb_record.eval_mirna()

            # Write the records to the output file
            out_kshv_file.write_record(hyb_record)

print('Performing Target Analysis...\n')
# Reiterate over output HybFile and perform target analysis.
with hybkit.HybFile(out_file_path, 'r') as out_kshv_file:
    # Prepare target analysis dict:
    target_analysis = hybkit.analysis.TargetAnalysis(name=in_file_label)

    for hyb_record in out_kshv_file:
        mirna_name = hyb_record.mirna_detail('mirna_name',
            allow_mirna_dimers=hybkit.settings.Analysis_settings['allow_mirna_dimers']
        )
        if 'kshv' in mirna_name:
            target_analysis.add(hyb_record)
    

    # results = target_analysis.results() 

    # Write target information to output file
    # Set analysis basename without ".hyb" extension
    analysis_basename = out_file_path.replace('.hyb', '')
    print('Writing Individual Analysis Files to Name Base:\n    %s' % analysis_basename)
    target_analysis.write_individual(analysis_basename)
    target_analysis.plot_individual(analysis_basename)
    print('Writing Combined Analysis Files to Name Base:\n    %s' % analysis_basename)
    target_analysis.write(analysis_basename)
    target_analysis.plot(analysis_basename)

print('Done\n')
