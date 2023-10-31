#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Example grouped target analysis performed as a python workflow.

Provided as an example of direct usage of hybkit functions.

See: 'README.rst' for this analysis for more information.
"""

import os

import hybkit

# Set mirna types as custom to include KSHV-miRNAs
hybkit.util.set_setting('mirna_types', set_value=['miRNA', 'KSHV-miRNA'])
# Tell hybkit that identifiers are in Hyb-Program standard format.
hybkit.util.set_setting('hybformat_id', set_value=True)

# Set script directories and input file names.
analysis_dir = os.path.abspath(os.path.dirname(__file__))
analysis_label = 'kshv_combined'
analysis_name = 'KSHV Combined'
out_dir = os.path.join(analysis_dir, 'output_python')
input_files = [f for f in os.listdir(analysis_dir) if f.endswith('.hyb')]
input_files.sort()
out_file_path = os.path.join(analysis_dir, 'output_python', analysis_label + '.hyb')
match_legend_file = os.path.join(analysis_dir, 'string_match_legend.csv')

# Begin Analysis
print('\nPerforming Grouped Target Analysis...')
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

# Begin Analysis
print('Outputting KSHV-Specific Hybrids to:\n    %s\n' % out_file_path)

# Open out-file for writing, one output file will be used for all inputs.
with hybkit.HybFile(out_file_path, 'w') as out_kshv_file:
    # Iterate through input files:
    for in_file_path in input_files:
        in_file_name = os.path.basename(in_file_path)
        in_file_label = in_file_name.replace('.hyb', '')
        print('Analyzing:\n    %s' % in_file_path)

        with hybkit.HybFile(in_file_path, 'r') as in_file:
            # Track last record identifier, to use only one record per read
            last_record_id = None

            # Iterate over each record of the input file
            for hyb_record in in_file:

                # Analyze only sequences where any segment identifier contains the string "kshv"
                if not hyb_record.has_prop('any_seg_contains', 'kshv'):
                    continue

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

                # Evaluate whether hybrid contains a miRNA.
                hyb_record.eval_mirna()

                # If assigned 'miRNA' does not contain string 'kshv', skip.
                if (not hyb_record.has_prop('has_mirna')  # noqa: SIM114
                        or not hyb_record.has_prop('mirna_contains', 'kshv')):

                    continue
                # If record is a duplicate, skip it.
                elif hyb_record.id == last_record_id:
                    continue
                last_record_id = hyb_record.id

                # Set dataset flag of record
                hyb_record.set_flag('dataset', in_file_label)

                # Write the records to the output file
                out_kshv_file.write_record(hyb_record)

print('\nPerforming Combined Target Analysis...\n')
# Repeat iteration over output HybFile and perform target analysis.
with hybkit.HybFile(out_file_path, 'r') as out_kshv_file:
    # Initialize target analysis:
    hyb_analysis = hybkit.analysis.Analysis(analysis_types='target', name=analysis_name)

    # Perform target-analysis of mirna within kshv-associated data.
    hyb_analysis.add_hyb_records(out_kshv_file)

    # Write target information to output file
    # Set analysis basename without ".hyb" extension
    analysis_basename = out_file_path.replace('.hyb', '')
    print('Writing Analysis Files to Name Base:\n    %s' % analysis_basename)
    hyb_analysis.write_analysis_results_special(out_basename=analysis_basename)
    hyb_analysis.plot_analysis_results(out_basename=analysis_basename)

print('Done\n')
