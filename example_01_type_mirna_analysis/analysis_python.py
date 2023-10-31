#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Analysis for type/mirna analysis performed as a python workflow.

Provided as an example of direct use of the Hybkit API.
See: "README.rst" for this analysis for more information.
"""

import os
import sys

import hybkit

# Set mirna types as custom to include KSHV-miRNAs
hybkit.settings.HybRecord_settings['mirna_types'] = ['miRNA', 'KSHV-miRNA']

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.join(analysis_dir, 'output_python')
input_files = sorted(f for f in os.listdir(analysis_dir) if f.endswith('.hyb'))
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

# Set hybrid segment types to remove as part of quality control (QC)
remove_types = ['rRNA', 'mitoch-rRNA']

# Initialize Combined Analysis
combined_analysis = hybkit.analysis.Analysis(
    analysis_types=['type', 'mirna'], name='Combined Analysis'
)

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
    file_analysis = hybkit.analysis.Analysis(
        analysis_types=['type', 'mirna'], name=in_file_label
    )

    # Open one HybFile entry for reading, and one for writing
    with hybkit.HybFile(in_file_path, 'r', hybformat_id=True) as in_file, \
         hybkit.HybFile(out_file_path, 'w') as out_file:

        # Track last record identifier, to use only one record per read
        last_record_id = None

        # Iterate over each record of the input file
        for hyb_record in in_file:
            hyb_record.eval_types()  # Find segment types
            hyb_record.eval_mirna()  # Find miRNA details

            # Determine if record has type that is excluded
            use_record = True
            for remove_type in remove_types:
                if hyb_record.has_prop('any_seg_type_is', remove_type):
                    use_record = False
                    break

            # If record has an excluded type, continue to next record without analyzing.
            if not use_record:  # noqa: SIM114
                continue

            # If record is a duplicate, skip it.
            elif hyb_record.id == last_record_id:
                continue
            last_record_id = hyb_record.id

            # Set the dataset label for the record
            hyb_record.set_flag('dataset', in_file_label)

            # Add record details to analyses
            file_analysis.add_hyb_record(hyb_record)
            combined_analysis.add_hyb_record(hyb_record)

            # Write the modified record to the output file
            out_file.write_record(hyb_record)

    # Write analysis for input file to outputs.
    analysis_file_basename = out_file_path.replace('.hyb', '')
    print('Outputting Analyses to:\n    %s\n' % analysis_file_basename)
    sys.stdout.flush()
    file_analysis.write_analysis_results_special(analysis_file_basename)
    file_analysis.plot_analysis_results(analysis_file_basename)

# Write combined summary analysis to files.
combined_analysis_file_basename = os.path.join(out_dir, 'combined_analysis')
print('Outputting Combined Analysis to:\n    %s\n' % combined_analysis_file_basename)
sys.stdout.flush()
combined_analysis.write_analysis_results_special(combined_analysis_file_basename)
combined_analysis.plot_analysis_results(combined_analysis_file_basename)

print('Done!')
