#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Helper data objects for automatic testing of hybkit code.
"""

import copy
import os
import hybkit
import pytest
from contextlib import nullcontext as does_not_raise

# Real data examples for testing.
# Example hyb record string for testing.
HYB_STR_1 = (
    '695_804	ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC	.	'
    'MIMAT0000078_MirBase_miR-23a_microRNA	1	21	1	21	0.0027	'
    'ENSG00000188229_ENST00000340384_TUBB2C_mRNA	23	49	1181	1207	1.2e-06	dataset=test;'
)

# Example Vienna record string for testing.
VIENNA_STR_1 = (
    '>695_804_MIMAT0000078_MirBase_miR-23a_microRNA_1_21-'
    '695_804_ENSG00000188229_ENST00000340384_TUBB2C_mRNA_1181_1207\n'
    'ATCACATTGCCAGGGATTTCCATCCCCAACAATGTGAAAACGGCTGTC\n'
    '.((((((((...(((((....)))))...))))))))...........	(-15)\n'
)
# Example Vienna record string with sequence mismatch for testing.
VIENNA_STR_1_MIS = (
    '>695_804_MIMAT0000078_MirBase_miR-23a_microRNA_1_21-'
    '695_804_ENSG00000188229_ENST00000340384_TUBB2C_mRNA_1181_1207\n'
    'ATCACATAGCCAGGGATTTCCATCCCCAACAATGTGAAAACGGCTGTC\n'
    '.((((((((...(((((....)))))...))))))))...........	(-15)\n'
)

# Artificial Hyb record strings and parameters for eval_type and eval_mirna testing.
ART_HYB_PROPS_1 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_microRNA	21	40	21	40	0.001	'
        'dataset=artificial'
    ),
    'seg1_type': 'microRNA',
    'seg2_type': 'microRNA',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_microRNA',
    'miRNA_seg': 'B',
    'true_prop_argsets': [
        ('has_mirna',),
        ('mirna_dimer',),
        ('5p_mirna',),
        ('3p_mirna',),
    ],
    'false_prop_argsets': [
        ('no_mirna',),
        ('mirna_not_dimer',),
    ],
    'true_is_set_argsets' : [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets' : [
        ('eval_target',),
        ('fold_record',),
    ],
    'one_mirna_error': True,
    'two_mirna_error': False,
    'mirna_seg_type': 'microRNA',
    'target_seg_type': 'microRNA',
    'mirna_ref': 'ARTSEG1_SOURCE_NAME_microRNA',
    'target_ref': 'ARTSEG2_SOURCE_NAME_microRNA',
    'mirna_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'target_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'mirna_fasta_id': '1_1000:1-20',
    'mirna_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'target_fasta_id': '1_1000:21-40',
    'target_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_microRNA',
}
    
ART_HYB_PROPS_2 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_coding	21	40	21	40	0.001	'
        'dataset=artificial'
    ),
    'seg1_type': 'microRNA',
    'seg2_type': 'coding',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_coding',
    'miRNA_seg': '5p',
    'true_prop_argsets': [
        ('has_mirna',),
        ('5p_mirna',),
        ('mirna_not_dimer',),
    ],
    'false_prop_argsets': [
        ('no_mirna',),
        ('3p_mirna',),
        ('mirna_dimer',),
    ],
    'true_is_set_argsets' : [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets' : [
        ('eval_target',),
        ('fold_record',),
    ],
    'one_mirna_error': False,
    'two_mirna_error': False,
    'mirna_seg_type': 'microRNA',
    'target_seg_type': 'coding',
    'mirna_ref': 'ARTSEG1_SOURCE_NAME_microRNA',
    'target_ref': 'ARTSEG2_SOURCE_NAME_coding',
    'mirna_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'target_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'mirna_fasta_id': '1_1000:1-20',
    'mirna_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'target_fasta_id': '1_1000:21-40',
    'target_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_coding',
}
ART_HYB_PROPS_3 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_coding	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_microRNA	21	40	21	40	0.001	'
        'dataset=artificial'
    ),
    'seg1_type': 'coding',
    'seg2_type': 'microRNA',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_coding',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_microRNA',
    'miRNA_seg': '3p',
    'true_prop_argsets': [
        ('has_mirna',),
        ('3p_mirna',),
        ('mirna_not_dimer',),
    ],
    'false_prop_argsets': [
        ('no_mirna',),
        ('5p_mirna',),
        ('mirna_dimer',),
    ],
    'true_is_set_argsets' : [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets' : [
        ('eval_target',),
        ('fold_record',),
    ],
    'one_mirna_error': False,
    'two_mirna_error': False,
    'mirna_seg_type': 'microRNA',
    'target_seg_type': 'coding',
    'mirna_ref': 'ARTSEG2_SOURCE_NAME_microRNA',
    'target_ref': 'ARTSEG1_SOURCE_NAME_coding',
    'mirna_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'target_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'mirna_fasta_id': '1_1000:21-40',
    'mirna_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_microRNA',
    'target_fasta_id': '1_1000:1-20',
    'target_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_coding',
}
ART_HYB_PROPS_4 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_coding	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_coding	21	40	21	40	0.001	'
        'dataset=artificial'
    ),
    'seg1_type': 'coding',
    'seg2_type': 'coding',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_coding',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_coding',
    'miRNA_seg': 'N',
    'true_prop_argsets': [
        ('no_mirna',),
    ],
    'false_prop_argsets': [
        ('has_mirna',),
        ('5p_mirna',),
        ('3p_mirna',),
        ('mirna_not_dimer',),
        ('mirna_dimer',),
    ],
    'true_is_set_argsets' : [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets' : [
        ('eval_target',),
        ('fold_record',),
    ],
    'one_mirna_error': True,
    'two_mirna_error': True,
    'mirna_seg_type': None,
    'target_seg_type': None,
    'mirna_ref': None,
    'target_ref': None,
    'mirna_seq': None,
    'target_seq': None,
    'mirna_fasta_id': None,
    'mirna_fasta_id_annotate': None,
    'target_fasta_id': None,
    'target_fasta_id_annotate': None,
}

ART_HYB_PROPS_ALL = [ART_HYB_PROPS_1, ART_HYB_PROPS_2, ART_HYB_PROPS_3, ART_HYB_PROPS_4]

# Add additional properties to check string propertier
ART_HYB_STR_PROPS = copy.deepcopy(ART_HYB_PROPS_2)
ART_HYB_STR_PROPS['true_prop_argsets'] = [
    ('id_is', '1_1000'),
    ('id_prefix', '1_'),
    ('id_suffix', '000'),
    ('id_contains', '1_100'),
    ('seq_is', 'AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG'),
    ('seq_prefix', 'AAAAAAAAAAAAAAAAAAAA'),
    ('seq_suffix', 'GGGGGGGGGGGGGGGGGGGG'),
    ('seq_contains', 'AAAAAAAGGGGGGGGGGG'),
    ('seg1_is', 'ARTSEG1_SOURCE_NAME_microRNA'),
    ('seg1_prefix', 'ARTSEG1_SOURCE_NAME'),
    ('seg1_suffix', 'microRNA'),
    ('seg1_contains', '_SOURCE_NAME_m'),
    ('seg2_is', 'ARTSEG2_SOURCE_NAME_coding'),
    ('seg2_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('seg2_suffix', 'coding'),
    ('seg2_contains', '_SOURCE_NAME_c'),
    ('any_seg_is', 'ARTSEG2_SOURCE_NAME_coding'),
    ('any_seg_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('any_seg_suffix', 'coding'),
    ('any_seg_contains', '_SOURCE_NAME_c'),
    ('seg1_type_is', 'microRNA'),
    ('seg1_type_prefix', 'micro'),
    ('seg1_type_suffix', 'RNA'),
    ('seg1_type_contains', 'icro'),
    ('seg2_type_is', 'coding'),
    ('seg2_type_prefix', 'cod'),
    ('seg2_type_suffix', 'ing'),
    ('seg2_type_contains', 'od'),
    ('any_seg_type_is', 'coding'),
    ('any_seg_type_prefix', 'cod'),
    ('any_seg_type_suffix', 'ing'),
    ('any_seg_type_contains', 'od'),
    ('mirna_is', 'ARTSEG1_SOURCE_NAME_microRNA'),
    ('mirna_prefix', 'ARTSEG1_SOURCE_NAME'),
    ('mirna_suffix', 'microRNA'),
    ('mirna_contains', '_SOURCE_NAME_m'),
    ('target_is', 'ARTSEG2_SOURCE_NAME_coding'),
    ('target_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('target_suffix', 'coding'),
    ('target_contains', '_SOURCE_NAME_c'),
    ('mirna_seg_type_is', 'microRNA'),
    ('mirna_seg_type_prefix', 'micro'),
    ('mirna_seg_type_suffix', 'RNA'),
    ('mirna_seg_type_contains', 'icro'),
    ('target_seg_type_is', 'coding'),
    ('target_seg_type_prefix', 'cod'),
    ('target_seg_type_suffix', 'ing'),
    ('target_seg_type_contains', 'od'),
]

ART_HYB_STR_PROPS['false_prop_argsets'] = [ \
    (vals[0], vals[1]+'XXX') for vals in ART_HYB_STR_PROPS['true_prop_argsets'] if len(vals) == 2
]

# Add bad hyb stirngs with missing name or seq
ART_BAD_HYB_STRS = []
for i in range(2):
    source_str = ART_HYB_PROPS_1['hyb_str']
    items = source_str.split('\t')
    items[i] = '.'
    bad_hyb_str = '\t'.join(items)
    ART_BAD_HYB_STRS.append(bad_hyb_str)

# Add bad hyb strings with missing columns
source_str = ART_HYB_PROPS_1['hyb_str']
items = source_str.split('\t')
del(items[4])
bad_hyb_str = '\t'.join(items)
ART_BAD_HYB_STRS.append(bad_hyb_str)

# Example Artificial Hyb / Viennarecord strings and parameters for eval_type and eval_mirna testing.
ART_HYB_VIENNA_PROPS_1 = {
    'hyb_str': (
        '1_1000	GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_coding	21	40	21	40	0.001	'
        'dataset=artificial'
    ),
    'vienna_str': (
        '>1_1000_ARTSEG1_SOURCE_NAME_microRNA-ARTSEG2_SOURCE_NAME_coding\n'
        'GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA\n'
        '...((((((((((((((......))))))))))))))...	(-15)'
    ),
    'overlapping': False,
    'seg1_seq':  'GGGCCCCCCCCCCCCCCGGG',
    'seg2_seq':  'AAAGGGGGGGGGGGGGGAAA',
    'seg1_fold': '...((((((((((((((...',
    'seg2_fold': '...))))))))))))))...',
    'true_prop_argsets': [
    ],
    'false_prop_argsets': [
    ],
    'true_is_set_argsets' : [
    ],
    'false_is_set_argsets' : [
    ],
}

ART_HYB_VIENNA_PROPS_2 = {
    'hyb_str': (
        '1_1000	GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	24	1	24	0.001	'
        'ARTSEG2_SOURCE_NAME_coding	17	40	17	40	0.001	'
        'dataset=artificial'
    ),
    'vienna_str': (
        '>1_1000_ARTSEG1_SOURCE_NAME_microRNA-ARTSEG2_SOURCE_NAME_coding\n'
        'GGGCCCCCCCCCCCCCCGGGAAAGCGGGAAAGGGGGGGGGGGGGGAAA\n'
        '...((((((((((((((......()......))))))))))))))...	(-15)'
    ),
    'overlapping': True,
    'seg1_seq': 'GGGCCCCCCCCCCCCCCGGGAAAG',
    'seg2_seq': 'CGGGAAAGGGGGGGGGGGGGGAAA',
    'seg1_fold': '...((((((((((((((......(',
    'seg2_fold': ')......))))))))))))))...',
    'true_prop_argsets': [
    ],
    'false_prop_argsets': [
    ],
    'true_is_set_argsets' : [
    ],
    'false_is_set_argsets' : [
    ],
}

ART_BAD_HYB_VIENNA_PROPS_1 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_1['hyb_str'] = ART_BAD_HYB_VIENNA_PROPS_1['hyb_str'].replace('1_1000','1_1000XXX')

ART_BAD_HYB_VIENNA_PROPS_2 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_2['hyb_str'] = ART_BAD_HYB_VIENNA_PROPS_2['hyb_str'].replace('CCCGGGA','CCCGGGAXXX')


# Empty seg_props dictionary for testing.
EMPTY_SEG_PROPS = {
    'ref_name': None,
    'read_start': None,
    'read_end': None,
    'ref_start': None,
    'ref_end': None,
    'score': None
}

# List of data objects for testing
TEST_HYBID_STR = '123_456'
TEST_SEQ_STR = 'ATGCATGCATGC'
TEST_FOLD_STR = '..((..))...'
TEST_ENERGY_STR = '-7.5'
TEST_SEG_PROPS = {
    'ref_name': 'test_seg',
    'read_start': 1,
    'read_end': 21,
    'ref_start': 1,
    'ref_end': 21,
    'score': 'e-10'
}
TEST_SEG_PROPS_STR = {p: str(TEST_SEG_PROPS[p]) for p in TEST_SEG_PROPS.keys()}
TEST_FLAGS_OBJ = {
    'count_total': 10,
    'count_last_clustering': 11,
    'two_way_merged': 'TRUE',
    'seq_IDs_in_cluster': 'test_id_1,test_id_2',
    'read_count': 4,
    'orient': 'F',
    'det': 'Arbitrary_Test_Detail',
    'seg1_type': 'microRNA',
    'seg2_type': 'mRNA',
    'seg1_det': 'Arbitrary_Test_Seg1_Detail',
    'seg2_det': 'Arbitrary_Test_Seg2_Detail',
    'miRNA_seg': '5p',
    'target_reg': '3p',
    'ext': 'FALSE',
    'dataset': 'test_dataset',
}
TEST_READ_COUNT = TEST_FLAGS_OBJ['read_count']
TEST_READ_COUNT_STR = str(TEST_READ_COUNT)
TEST_RECORD_COUNT = TEST_FLAGS_OBJ['count_total']
TEST_RECORD_COUNT_STR = str(TEST_RECORD_COUNT)
TEST_FLAGS_STR = {flag: str(TEST_FLAGS_OBJ[flag]) for flag in TEST_FLAGS_OBJ.keys()}
TEST_FLAGS_STR_LINE = (
    'count_total=10;count_last_clustering=11;two_way_merged=TRUE;seq_IDs_in_cluster=test_id_1,test_id_2;'
    + 'read_count=4;orient=F;det=Arbitrary_Test_Detail;seg1_type=microRNA;seg2_type=mRNA;'
    + 'seg1_det=Arbitrary_Test_Seg1_Detail;seg2_det=Arbitrary_Test_Seg2_Detail;miRNA_seg=5p;'
    + 'target_reg=3p;ext=FALSE;dataset=test_dataset'
)
TEST_FLAGS_STR_REVERSED = {flag: TEST_FLAGS_STR[flag] for flag in reversed([*TEST_FLAGS_STR.keys()])}
TEST_OBJECTS = {
    'None': None,
    'id_str': TEST_HYBID_STR,
    'test_seq_str': TEST_SEQ_STR,
    'dot_str': '.',
    'fold_str': TEST_FOLD_STR,
    'int_str': '473',
    'float_str': '473.0',
    'empty_str': '',
    'space_str': ' ',
    'tab_str': '\t',
    'int': 473,
    'float': 473.0,
    'list': [TEST_HYBID_STR],
    'tuple': (TEST_HYBID_STR,),
    'set': {TEST_HYBID_STR},
    'empty_dict': {},
    'dict_rand_key': {TEST_HYBID_STR: True},
    'seg_props_dict': copy.deepcopy(EMPTY_SEG_PROPS),
}
NONSPACE_STRS = ['id_str', 'int_str', 'float_str', 'test_seq_str', 'fold_str']
ALL_STRS = ['empty_str', 'space_str', 'tab_str', *NONSPACE_STRS]
NONE_OBJS = {'None', 'dot_str'}
ID_ALLOWED_TYPES = {*NONSPACE_STRS}
SEQ_ALLOWED_TYPES = {'test_seq_str'}
ENERGY_ALLOWED_TYPES = {*NONE_OBJS, 'int_str', 'int', 'float_str', 'float'}
FOLD_ALLOWED_TYPES = {'fold_str'}
SEG_PROPS_ALLOWED_TYPES = {'seg_props_dict', 'empty_dict', 'None'}
REF_NAME_ALLOWED_TYPES = {*NONE_OBJS, *NONSPACE_STRS}
READ_START_ALLOWED_TYPES = {*NONE_OBJS, 'int_str', 'int'}
READ_END_ALLOWED_TYPES = {*NONE_OBJS, 'int_str', 'int'}
REF_START_ALLOWED_TYPES = {*NONE_OBJS, 'int_str', 'int'}
REF_END_ALLOWED_TYPES = {*NONE_OBJS, 'int_str', 'int'}
SCORE_ALLOWED_TYPES = {*NONE_OBJS, *ALL_STRS, 'int', 'float'}

# ----- Set Testing File Paths -----
auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
test_out_dir = os.path.join(auto_tests_dir, 'output_autotest')
hyb_file_name = os.path.join(test_out_dir, 'test_hybrid_py.hyb')
vienna_file_name = os.path.join(test_out_dir, 'test_hybrid_py.vienna')
hyb_autotest_file_name = os.path.join(test_out_dir, 'test_hybrid_py_autotest.hyb')
out_basename = hyb_file_name.replace('.hyb', '')
out_flag_string_name = os.path.join(test_out_dir, 'flag_string.txt')
match_legend_file_name = os.path.join(auto_tests_dir, 'test_string_match.csv')
bad1_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_1.csv')
bad2_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_2.csv')
bad3_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_3.csv')
id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_id_map.csv')
bad1_id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_id_map_1.csv')