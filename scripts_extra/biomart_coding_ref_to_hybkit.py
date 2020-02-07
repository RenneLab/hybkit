#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Convert coding sequence information from biomart TSV to hybkit-format.

Convert file(s) output by the get_biomart_coding.py or get_biomart_coding_archives.py
scripts to a format readable by hybkit.
This entails reading in each successive file (that should be provided in order from
most recent to least recent). An "identifier" column is added, in the fashion provided
by hybkit with format:
    <gene_id>_<transcript_id>_<gene_name>_<seg_type>. 
In this case, this uses columns:
    For releases >= 89:
    'Gene stable ID', 'Transcript stable ID', 'Gene name', and 'Gene type'

    For release 88:
    'Gene ID', 'Transcript ID', 'Associated Gene Name', 'Gene type' 

    For releases 79 to 87:
    'Ensembl Gene ID', 'Ensembl Transcript ID', 'Associated Gene Name', 'Gene type;

    For releases 54 to 78:
    'Ensembl Gene ID', 'Ensembl Transcript ID', 'Associated Gene Name', 'Gene Biotype'

"""

import datetime
import textwrap
import csv
import sys
import os

try:
    import hybkit
except ModuleNotFoundError:
    message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
    message += 'on your $PYTHONPATH'
    print(message)
    raise


READ_DIALECT = 'excel-tab'
WRITE_DIALECT = 'excel'
TOTAL_INDENT = 2
BODY_INDENT = 2
PROPERTY_INDENT = 4
FINAL_WIDTH = 80

wd = os.path.join(hybkit.__about__.package_dir, 'working')

if not os.path.isdir(wd):
    print('Creating Output Directory:\n    %s\n' % wd)
    os.mkdir(wd)

data_dir = os.path.join(wd)

out_file = os.path.join(wd, 'hybkit_coding_ref_combined.csv')
out_file_detail_name = out_file.replace('.csv', '.notes.txt')
input_files = [os.path.join(data_dir, 'Ensembl_coding_details_99.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_98.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_97.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_96.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_95.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_94.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_93.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_92.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_91.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_90.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_89.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_88.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_87.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_86.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_85.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_84.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_83.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_82.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_81.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_80.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_79.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_78.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_77.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_75.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_67.tsv'),
               os.path.join(data_dir, 'Ensembl_coding_details_54.tsv'),
               ]

write_columns = ['identifier', 'cdna_coding_start', 'cdna_coding_end']
id_fields = [['Gene stable ID', 'Gene ID', 'Ensembl Gene ID'],
             ['Transcript stable ID', 'Transcript ID', 'Ensembl Transcript ID'],
             ['Gene name', 'Associated Gene Name'],
             ['Gene type', 'Gene Biotype']]

written_identifiers = {}
all_in_file_details = []
print('\nWriting to combined output file:\n    %s\n' % out_file)

with open(out_file, 'w', newline='') as out_file_obj:
    writer = csv.DictWriter(out_file_obj, dialect=WRITE_DIALECT,
                            fieldnames=write_columns, extrasaction='ignore')
    writer.writeheader()
    for in_file in input_files:

        in_file_detail = in_file.replace('.tsv', '.notes.txt')
        print('Reading input file details:\n    %s' % in_file_detail)
        with open(in_file_detail, 'r') as in_file_detail_obj:
            all_in_file_details.append(in_file_detail_obj.read().rstrip())

        print('Reading input file:\n    %s' % in_file)
        with open(in_file, 'r', newline='') as in_file_obj:
            reader = csv.DictReader(in_file_obj, dialect=READ_DIALECT)
            for line_info in reader:
                #print(line_info)
                if 'Gene type' in line_info and line_info['Gene type'] == 'protein_coding':
                    line_info['Gene type'] = 'mRNA'  # Change to match Hyb Program nomenclature
                elif 'Gene Biotype' in line_info and line_info['Gene Biotype'] == 'protein_coding':
                    line_info['Gene Biotype'] = 'mRNA'  # Change to match Hyb Program nomenclature
                id_items = []
                for field_options in id_fields:
                    for option in field_options:
                        if option in line_info:
                            id_items.append(line_info[option])
                            break
                if len(id_items) != 4:
                    message = 'Problem Construction ID for file:\n    %s\n' % in_file
                    message += 'And line:\n    %s\n' % str(line_info)
                    print(message)
                    raise Exception(message)
                identifier = '_'.join(id_items)

                if (identifier not in written_identifiers
                    and 'cDNA coding start' in line_info and 'cDNA coding end' in line_info
                    and line_info['cDNA coding start'] and line_info['cDNA coding end']):
                    cds_length = abs(int(line_info['cDNA coding end'])
                                     - int(line_info['cDNA coding start']))
                    if cds_length < 50:
                        #print('WARNING: Bad CDS length: %i' % cds_length)
                        #print(line_info)
                        #input()
                        continue
                    line_info['identifier'] = identifier
                    line_info['cdna_coding_start'] = line_info['cDNA coding start']
                    line_info['cdna_coding_end'] = line_info['cDNA coding end']
                    writer.writerow(line_info)
                    written_identifiers[identifier] = None
              

    print('Writing details to file: %s\n' % out_file_detail_name)
    with open(out_file_detail_name, 'w') as out_detail_file:
        total_prefix = (' ' * TOTAL_INDENT)
        body_prefix = (' ' * BODY_INDENT)
        property_prefix = (' ' * PROPERTY_INDENT)
        total_width = FINAL_WIDTH - (TOTAL_INDENT)
        body_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT)
        property_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT + PROPERTY_INDENT)

        header = '-- %s -- ' % os.path.basename(out_file)
        description = 'File: "%s" ' % os.path.basename(out_file)
        description += 'was created by %s ' % (os.path.basename(__file__))
        description += '(part of the hybkit project, release %s)' % hybkit.__about__.__version__
        description += ' on %s.\n' % str(datetime.datetime.now().isoformat())
        description += 'Sequences downloaded from the biomart interface '
        description += 'were processed into columns: '
        description += ','.join(write_columns)
        description = textwrap.fill(description, width=body_width) + '\n'

        body = textwrap.indent('\n'.join([description]), body_prefix)
        full_text = textwrap.indent('\n'.join([header, body]), total_prefix)

        out_detail_file.write('\n\n\n'.join([full_text.rstrip()] + all_in_file_details))

print('\nDone!\n') 
