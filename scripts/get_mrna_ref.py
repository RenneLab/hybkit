#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Acquire list of ensembl protein-coding human genome transcripts
by using the "biomart" python module, available at:
https://pypi.org/project/biomart/
'''

import datetime
import textwrap
import os
try:
    import biomart
except ImportError:
    message = 'This script requires the "biomart" python module, available at:\n    '
    message += 'https://pypi.org/project/biomart/'
    print(message)
    raise

try:
    import hybkit
except ImportError:
    message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
    message += 'on your $PYTHONPATH'
    print(message)
    raise

# server_url = 'http://www.biomart.org/biomart'
server_url = 'http://www.ensembl.org/biomart'

out_file_base = 'Ensembl_coding_details'
out_file_name = out_file_base + '.tsv'
out_file_detail_name = out_file_base + '.notes.txt'

TOTAL_INDENT = 2
BODY_INDENT = 2
PROPERTY_INDENT = 4
FINAL_WIDTH = 80

print('\nBeginning Acquision of Coding Sequence Information from: %s\n' % server_url)
server = biomart.BiomartServer(server_url)
# Debugging settings and options:
server.verbose = True

# print('\nEnsembl Databases:\n')
# server.show_databases()

# print('\nEnsembl Datasets:\n')
# server.show_datasets()

#h_genes = server.datasets['hsapiens_gene_ensembl']
#h_genes.show_filters()
#h_genes.show_attributes()

h_genes = server.datasets['hsapiens_gene_ensembl'] 
print(dir(h_genes))

features_attributes = [
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'external_gene_name',
    'external_transcript_name',
    'gene_biotype',
    'transcript_biotype',
    'mirbase_id',
    'refseq_mrna',
    'refseq_ncrna',
    ]

sequence_attributes = [
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'external_gene_name',
    '5_utr_start',
    '5_utr_end',
    'cdna_coding_start',
    'cdna_coding_end',
    '3_utr_start',
    '3_utr_end',
    'ensembl_peptide_id',
    'strand',
    'transcript_start',
    'transcript_end',
    'transcription_start_site',
    'transcript_length',
    # 'transcript_exon_intron'
    ] 

with_sequence = False
if with_sequence:
    sequence_attributes.append('transcript_exon_intron')

use_attributes = sequence_attributes

filters = {
           # 'chromosome_name': #'1', #,'2'],#'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y',
           # 'transcript_biotype': 'protein_coding',
           # 'biotype': 'miRNA',
           'ensembl_gene_id': 'ENSG00000139618'
          } 

use_filters = filters

search_settings = {
    'filters': use_filters,
    'attributes': use_attributes,
    }

# for raw_line in response.iter_lines():
#     line = raw_line.decode('utf-8')
#     print(line)

#response = h_genes.search(search_settings, header=1)
#with open(out_file_name, 'wb') as out_file:
#    for raw_line in response.iter_lines():
#        out_file.write(raw_line)

print('Writing details to file: %s\n' % out_file_detail_name)

with open(out_file_detail_name, 'w') as out_detail_file:
    total_prefix = (' ' * TOTAL_INDENT)
    body_prefix = (' ' * BODY_INDENT)
    property_prefix = (' ' * PROPERTY_INDENT)
    total_width = FINAL_WIDTH - (TOTAL_INDENT)
    body_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT)
    property_width = FINAL_WIDTH - (TOTAL_INDENT + BODY_INDENT + PROPERTY_INDENT)

    header = '-- %s -- ' % out_file_name
    description = 'File: "%s" ' % out_file_name
    description += 'was created by %s ' % (os.path.basename(__file__))
    description += '(part of the hybkit project, release %s)' % hybkit.__about__.__version__
    description += ' on %s.\n' % str(datetime.datetime.now().isoformat())
    description += 'Details for each sequence were downloaded from the biomart interface '
    description += 'to the Ensembl database, accessed via: "%s"\n' % server_url
    description = textwrap.fill(description, width=body_width) + '\n'

    properties = 'The query was made with the following properties:\nFilters:\n'
    filter_text = ''
    for f_key in use_filters:
        filter_text += '%s: %s\n' % (f_key, str(use_filters[f_key]))
    properties += textwrap.indent(filter_text, property_prefix)
    properties += '\nAttributes:\n'
    properties += textwrap.indent('\n'.join(use_attributes), property_prefix)
    properties += '\n\n'

    body = textwrap.indent('\n'.join([description, properties]), body_prefix)
    full_text = textwrap.indent('\n'.join([header, body]), total_prefix)

    out_detail_file.write(full_text)

print('Done.\n')
