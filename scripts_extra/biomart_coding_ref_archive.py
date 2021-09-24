#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Acquire lists of deatils of archive Ensembl protein-coding human genome transcripts.

WARNING: This script is not a part of the official hybkit package, but is a part of the 
preparation of the databases included with hybkit. It is included here as a record and
reference of how the data was prepared, and can be adapted and modified as needed for
users' specific purposes.

This module utilizes a modified version of the pybiomart package to allow utilization of 
archive Ensembl databases. The original biomart package is available at:
    https://pypi.org/project/biomart/
The modified version is available at 
    http://www.github.com/dstrib/pybiomart

The biomart module also depends on the pandas package for data manipulation.
"""

import datetime
import textwrap
import sys
import os

sys.path.append('/home/ds/programs/pybiomart/src')
try:
    import pybiomart
except ImportError:
    message = 'This script requires the "pybiomart" python module, available at:\n    '
    message += 'https://pypi.org/project/pybiomart/'
    print(message)
    raise

try:
    import hybkit
except ModuleNotFoundError:
    message = 'The "hybkit" module cannot be found. Please ensure this is accessible '
    message += 'on your $PYTHONPATH'
    print(message)
    raise

TOTAL_INDENT = 2
BODY_INDENT = 2
PROPERTY_INDENT = 4
FINAL_WIDTH = 80


# Release list was last updated on 2020-02-06.
release_list = [
                ('99', 'Current', 'http://www.ensembl.org'),
                ('98', 'Sep 2019', 'http://Jul2019.archive.ensembl.org'),
                ('97', 'Jul 2019', 'http://Apr2019.archive.ensembl.org'),
                ('96', 'Apr 2019', 'http://Jan2019.archive.ensembl.org'),
                ('95', 'Jan 2019', 'http://Oct2018.archive.ensembl.org'),
                ('94', 'Oct 2018', 'http://Jul2018.archive.ensembl.org'),
                ('93', 'Jul 2018', 'http://Apr2018.archive.ensembl.org'),
                ('92', 'Apr 2018', 'http://Dec2017.archive.ensembl.org'),
                ('91', 'Dec 2017', 'http://Aug2017.archive.ensembl.org'),
                ('90', 'Aug 2017', 'http://May2017.archive.ensembl.org'),
                ('89', 'May 2017', 'http://Mar2017.archive.ensembl.org'),
                ('88', 'Mar 2017', 'http://Dec2016.archive.ensembl.org'),
                ('87', 'Dec 2016', 'http://Oct2016.archive.ensembl.org'),
                ('86', 'Oct 2016', 'http://Jul2016.archive.ensembl.org'),
                ('85', 'Jul 2016', 'http://Mar2016.archive.ensembl.org'),
                ('84', 'Mar 2016', 'http://Dec2015.archive.ensembl.org'),
                ('83', 'Dec 2015', 'http://Sep2015.archive.ensembl.org'),
                ('82', 'Sep 2015', 'http://Jul2015.archive.ensembl.org'),
                ('81', 'Jul 2015', 'http://May2015.archive.ensembl.org'),
                ('80', 'May 2015', 'http://Mar2015.archive.ensembl.org'),
                ('79', 'Mar 2015', 'http://Dec2014.archive.ensembl.org'),
                ('78', 'Dec 2014', 'http://Oct2014.archive.ensembl.org'),
                ('77', 'Oct 2014', 'http://Feb2014.archive.ensembl.org'),
                ('75', 'Feb 2014', 'http://May2012.archive.ensembl.org'),
                ('67', 'May 2012', 'http://May2009.archive.ensembl.org'),
                ('54', 'May 2009', 'http://may2009.archive.ensembl.org'),
               ]

releases = {info[0]: info for info in release_list}

print('\n --- Beginning Archive Database Acquision ---')
print(__doc__)

wd = os.path.join(hybkit.__about__.package_dir, 'working')
if not os.path.isdir(wd):
    print('Creating Output Directory:\n    %s\n' % wd)
    os.mkdir(wd)

print('Changing working directory to:\n    %s\n' % wd)
os.chdir(wd)

for release, release_date, server_url in releases.values(): 
    print('Working on release: %s' % ', '.join([release, release_date, server_url])) 
    
    out_file_base = 'Ensembl_coding_details_' + release
    out_file_name = out_file_base + '.tsv'
    out_file_detail_name = out_file_base + '.notes.txt'
    
    print('\nBeginning Acquision of Coding Sequence Information from: %s\n' % server_url)
    sys.stdout.flush()

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
        'gene_biotype',
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
        'cds_start',
        'cds_end',
        # 'transcript_exon_intron'
        ] 
    
    sequence_attributes_78 = []
    with_sequence = False
    if with_sequence:
        sequence_attributes.append('transcript_exon_intron')
    
    use_attributes = sequence_attributes
    
    filters = {
               # 'chromosome_name': ['1','2'], #'1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y',
               'transcript_biotype': 'protein_coding',
               # 'link_ensembl_gene_id': 'ENSG00000139618'
              } 
    
    use_filters = filters
   

    # server = pybiomart.Server(host=server_url)
    # marts = server.list_marts()
    # print(marts)
 
    #server = pybiomart.Server(host='http://www.ensembl.org')
    #print(server.list_marts())
    #sys.exit()

    dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',
                                host=server_url)

    if int(release) <= 78:
        use_attributes = []
        for attribute in sequence_attributes:
            if attribute in dataset.attributes:
                use_attributes.append(attribute)
            elif attribute == 'external_gene_name':
                use_attributes.append('external_gene_id')
                print('Replacing "external_gene_name" with "external_gene_id"')
            else:
                print('Removing Attribute: %s' % attribute)

        use_filters = {'biotype':filters['transcript_biotype']}

    #print(dataset.list_filters())
    #print(dataset.list_attributes())
    #sys.exit()

    query = dataset.query(attributes=use_attributes,
                          filters=use_filters)


    #import npandas as pd
    #pd.options.display.float_format = '{:,.0f}'.format

    #print('Data Acquired.')
    #print(query.iloc[1])
    #print(query.astype(int)) 

    query.to_csv(out_file_name, sep='\t', index=False, float_format='%.0f')

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
        description += 'to the Ensembl database, release %s ' % release
        description += '(%s), accessed via: "%s"\n' % (release_date, server_url)
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
