#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project: https://www.github.com/RenneLab/hybkit

"""
Helper data objects for automatic testing of hybkit code.
"""

import copy
import os

from hybkit.errors import (
    HybkitArgError,
    HybkitConstructorError,
    HybkitError,
    HybkitIterError,
    HybkitMiscError,
)

# from contextlib import nullcontext as does_not_raise

# ----- Start Data Constants -----
NEG_10 = -10.0
ZERO = 0.0
HALF = 0.5
ONE = 1.0
TWO = 2.0
TWO_I = 2
FIVE = 5.0
FIVE_I = 5
TEN = 10.0
TEN_I = 10
FORTY_I = 40

# ----- Start Error Constants -----
ERROR_TYPE_STRINGS = {
    'HybkitError': HybkitError,
    'HybkitMiscError': HybkitMiscError,
    'HybkitArgError': HybkitArgError,
    'HybkitConstructorError': HybkitConstructorError,
    'HybkitIterError': HybkitIterError,
}
for init_str in list(ERROR_TYPE_STRINGS):
    ERROR_TYPE_STRINGS[init_str.lower()] = ERROR_TYPE_STRINGS[init_str]

# ----- Start Real Hyb/Vienna Data Examples -----
# Example hyb record string for testing.
HYB_STR_1 = (
    '695_804	ATCACATTGCCAGGGATTTCCAATCCCCAACAATGTGAAAACGGCTGTC	.	'
    'MIMAT0000078_MirBase_miR-23a_microRNA	1	21	1	21	0.0027	'
    'ENSG00000188229_ENST00000340384_TUBB2C_mRNA	23	49	1181	1207	1.2e-06	dataset=test;\n'
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

# ----- Start Artificial Hyb Data Examples -----
ART_HYB_PROPS_1 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_microRNA	21	40	21	40	0.001	'
        'dataset=artificial\n'
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
        ('has_indels',),
    ],
    'true_is_set_argsets': [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets': [
        ('eval_target',),
        ('fold_record',),
    ],
    'has_mirna': True,
    'has_one_mirna': False,
    'mirna_seg_props': 'seg1_props',
    'target_seg_props': 'seg2_props',
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
        'ARTSEG2_SOURCE_NAME_mRNA	21	40	21	40	0.001	'
        'dataset=artificial\n'
    ),
    'seg1_type': 'microRNA',
    'seg2_type': 'mRNA',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_mRNA',
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
        ('has_indels',),
    ],
    'true_is_set_argsets': [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets': [
        ('eval_target',),
        ('fold_record',),
    ],
    'has_mirna': True,
    'has_one_mirna': True,
    'mirna_seg_props': 'seg1_props',
    'target_seg_props': 'seg2_props',
    'mirna_seg_type': 'microRNA',
    'target_seg_type': 'mRNA',
    'mirna_ref': 'ARTSEG1_SOURCE_NAME_microRNA',
    'target_ref': 'ARTSEG2_SOURCE_NAME_mRNA',
    'mirna_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'target_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'mirna_fasta_id': '1_1000:1-20',
    'mirna_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_microRNA',
    'target_fasta_id': '1_1000:21-40',
    'target_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_mRNA',
}
ART_HYB_PROPS_3 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_mRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_microRNA	21	40	21	40	0.001	'
        'dataset=artificial\n'
    ),
    'seg1_type': 'mRNA',
    'seg2_type': 'microRNA',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_mRNA',
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
        ('has_indels',),
    ],
    'true_is_set_argsets': [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets': [
        ('eval_target',),
        ('fold_record',),
    ],
    'has_mirna': True,
    'has_one_mirna': True,
    'mirna_seg_props': 'seg2_props',
    'target_seg_props': 'seg1_props',
    'mirna_seg_type': 'microRNA',
    'target_seg_type': 'mRNA',
    'mirna_ref': 'ARTSEG2_SOURCE_NAME_microRNA',
    'target_ref': 'ARTSEG1_SOURCE_NAME_mRNA',
    'mirna_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'target_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'mirna_fasta_id': '1_1000:21-40',
    'mirna_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_microRNA',
    'target_fasta_id': '1_1000:1-20',
    'target_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_mRNA',
}
ART_HYB_PROPS_4 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_mRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_mRNA	21	40	21	40	0.001	'
        'dataset=artificial\n'
    ),
    'seg1_type': 'mRNA',
    'seg2_type': 'mRNA',
    'seg1_seq': 'AAAAAAAAAAAAAAAAAAAA',
    'seg2_seq': 'GGGGGGGGGGGGGGGGGGGG',
    'hybrid_fasta_id': '1_1000',
    'hybrid_fasta_id_annotate': 'artificial:1_1000',
    'seg1_fasta_id': '1_1000:1-20',
    'seg1_fasta_id_annotate': 'artificial:1_1000:1-20:ARTSEG1_SOURCE_NAME_mRNA',
    'seg2_fasta_id': '1_1000:21-40',
    'seg2_fasta_id_annotate': 'artificial:1_1000:21-40:ARTSEG2_SOURCE_NAME_mRNA',
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
        ('has_indels',),
    ],
    'true_is_set_argsets': [
        ('energy',),
        ('full_seg_props',),
        ('eval_types',),
        ('eval_mirna',),
    ],
    'false_is_set_argsets': [
        ('eval_target',),
        ('fold_record',),
    ],
    'has_mirna': False,
    'has_one_mirna': False,
    'mirna_seg_props': None,
    'target_seg_props': None,
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
    ('seg2_is', 'ARTSEG2_SOURCE_NAME_mRNA'),
    ('seg2_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('seg2_suffix', 'mRNA'),
    ('seg2_contains', '_SOURCE_NAME_mR'),
    ('any_seg_is', 'ARTSEG2_SOURCE_NAME_mRNA'),
    ('any_seg_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('any_seg_suffix', 'mRNA'),
    ('any_seg_contains', '_SOURCE_NAME_mR'),
    ('seg1_type_is', 'microRNA'),
    ('seg1_type_prefix', 'micro'),
    ('seg1_type_suffix', 'RNA'),
    ('seg1_type_contains', 'icro'),
    ('seg2_type_is', 'mRNA'),
    ('seg2_type_prefix', 'mR'),
    ('seg2_type_suffix', 'NA'),
    ('seg2_type_contains', 'RN'),
    ('any_seg_type_is', 'mRNA'),
    ('any_seg_type_prefix', 'mR'),
    ('any_seg_type_suffix', 'NA'),
    ('any_seg_type_contains', 'RNA'),
    ('mirna_is', 'ARTSEG1_SOURCE_NAME_microRNA'),
    ('mirna_prefix', 'ARTSEG1_SOURCE_NAME'),
    ('mirna_suffix', 'microRNA'),
    ('mirna_contains', '_SOURCE_NAME_m'),
    ('target_is', 'ARTSEG2_SOURCE_NAME_mRNA'),
    ('target_prefix', 'ARTSEG2_SOURCE_NAME'),
    ('target_suffix', 'mRNA'),
    ('target_contains', '_SOURCE_NAME_mR'),
    ('mirna_seg_type_is', 'microRNA'),
    ('mirna_seg_type_prefix', 'micro'),
    ('mirna_seg_type_suffix', 'RNA'),
    ('mirna_seg_type_contains', 'icro'),
    ('target_seg_type_is', 'mRNA'),
    ('target_seg_type_prefix', 'mR'),
    ('target_seg_type_suffix', 'NA'),
    ('target_seg_type_contains', 'RN'),
]

ART_HYB_STR_PROPS['false_prop_argsets'] = [
    (vals[0], (vals[1] + 'XXX')) for vals
    in ART_HYB_STR_PROPS['true_prop_argsets'] if len(vals) == TWO
]

# Add bad hyb strings with missing name or seq
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
del items[4:7]
bad_hyb_str = '\t'.join(items)
ART_BAD_HYB_STRS.append(bad_hyb_str)

# ----- Start Artificial Vienna Data -----
ART_HYB_VIENNA_PROPS_1 = {
    'hyb_str': (
        '1_1000	GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_mRNA	21	40	21	40	0.001	'
        'dataset=artificial\n'
    ),
    'vienna_str': (
        '>1_1000_ARTSEG1_SOURCE_NAME_microRNA-ARTSEG2_SOURCE_NAME_mRNA\n'
        'GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA\n'
        '...((((((((((((((......))))))))))))))...	(-10.0)\n'
    ),
    'overlapping': False,
    'mismatches': 0,
    'seg1_seq': 'GGGCCCCCCCCCCCCCCGGG',
    'seg2_seq': 'AAAGGGGGGGGGGGGGGAAA',
    'seg1_fold': '...((((((((((((((...',
    'seg2_fold': '...))))))))))))))...',
    'true_prop_argsets': [
    ],
    'false_prop_argsets': [
    ],
    'true_is_set_argsets': [
    ],
    'false_is_set_argsets': [
    ],
}

ART_HYB_VIENNA_PROPS_2 = {
    'hyb_str': (
        '1_1000	GGGCCCCCCCCCCCCCCGGGAAAGGGGGGGGGGGGGGAAA	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	24	1	24	0.001	'
        'ARTSEG2_SOURCE_NAME_mRNA	17	40	17	40	0.001	'
        'dataset=artificial\n'
    ),
    'vienna_str': (
        '>1_1000_ARTSEG1_SOURCE_NAME_microRNA-ARTSEG2_SOURCE_NAME_mRNA\n'
        'GGGCCCCCCCCCCCCCCGGGAAAGCGGGAAAGGGGGGGGGGGGGGAAA\n'
        '...((((((((((((((......()......))))))))))))))...	(-10.0)\n'
    ),
    'overlapping': True,
    'mismatches': 0,
    'seg1_seq': 'GGGCCCCCCCCCCCCCCGGGAAAG',
    'seg2_seq': 'CGGGAAAGGGGGGGGGGGGGGAAA',
    'seg1_fold': '...((((((((((((((......(',
    'seg2_fold': ')......))))))))))))))...',
    'true_prop_argsets': [
    ],
    'false_prop_argsets': [
    ],
    'true_is_set_argsets': [
    ],
    'false_is_set_argsets': [
    ],
}

# Static Seq mismatched
ART_BAD_HYB_VIENNA_PROPS_1 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_1['hyb_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_1['hyb_str'].replace('CCCGGGA', 'CCCGGGAXXX')
)
ART_BAD_HYB_VIENNA_PROPS_1['mismatches'] = 3

# Correct dynamic but flagged as static
ART_BAD_HYB_VIENNA_PROPS_2 = copy.deepcopy(ART_HYB_VIENNA_PROPS_2)
ART_BAD_HYB_VIENNA_PROPS_2['overlapping'] = False
ART_BAD_HYB_VIENNA_PROPS_2['mismatches'] = 3

# Dynamic with Mismatched Seq
ART_BAD_HYB_VIENNA_PROPS_3 = copy.deepcopy(ART_HYB_VIENNA_PROPS_2)
ART_BAD_HYB_VIENNA_PROPS_3['hyb_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_3['hyb_str'].replace('CCCGGGA', 'CCCGGGAXXX')
)
ART_BAD_HYB_VIENNA_PROPS_3['mismatches'] = 3

# Mismatched Energy Values
ART_BAD_HYB_VIENNA_PROPS_4 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_4['vienna_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_4['vienna_str'].replace('(-10.0)', '(-15.0)')
)

# Missing vienna string
ART_BAD_HYB_VIENNA_PROPS_5 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_5['vienna_str'] = '\n\n\n'

# Bad vienna string
ART_BAD_HYB_VIENNA_PROPS_6 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_6['vienna_str'] = 'bad\nvienna\nstring\n'

# Energy 99 Error
ART_BAD_HYB_VIENNA_PROPS_7 = copy.deepcopy(ART_HYB_VIENNA_PROPS_1)
ART_BAD_HYB_VIENNA_PROPS_7['vienna_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_7['vienna_str'].replace('(-10.0)', '(99.0)')
)

# Insertion in read
ART_BAD_HYB_VIENNA_PROPS_8 = copy.deepcopy(ART_HYB_VIENNA_PROPS_2)
ART_BAD_HYB_VIENNA_PROPS_8['hyb_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_8['hyb_str'].replace(
        '1	24	1	24',
        '1	24	1	22'
    )
)
# Deletion in read
ART_BAD_HYB_VIENNA_PROPS_9 = copy.deepcopy(ART_HYB_VIENNA_PROPS_2)
ART_BAD_HYB_VIENNA_PROPS_9['hyb_str'] = (
    ART_BAD_HYB_VIENNA_PROPS_9['hyb_str'].replace(
        '1	24	1	24',
        '1	24	1	26'
    )
)


# ----- Start Hyb Properties -----
EMPTY_SEG_PROPS = {
    'ref_name': None,
    'read_start': None,
    'read_end': None,
    'ref_start': None,
    'ref_end': None,
    'score': None
}

# List of data objects for testing
TEST_HYB_ID_STR = '123_456'
TEST_SEQ_STR = 'ATGCATGCATGC'
TEST_HYB_MINIMAL_STRING = '\t'.join([TEST_HYB_ID_STR, TEST_SEQ_STR] + (['.'] * 13))
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
TEST_SEG_PROPS_STR = {p: str(TEST_SEG_PROPS[p]) for p in TEST_SEG_PROPS}
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
TEST_FLAGS_STR = {flag: str(TEST_FLAGS_OBJ[flag])
                  for flag in TEST_FLAGS_OBJ}
TEST_FLAGS_STR_LINE = (
    'count_total=10;count_last_clustering=11;two_way_merged=TRUE;'
    + 'seq_IDs_in_cluster=test_id_1,test_id_2;'
    + 'read_count=4;orient=F;det=Arbitrary_Test_Detail;seg1_type=microRNA;seg2_type=mRNA;'
    + 'seg1_det=Arbitrary_Test_Seg1_Detail;seg2_det=Arbitrary_Test_Seg2_Detail;miRNA_seg=5p;'
    + 'target_reg=3p;ext=FALSE;dataset=test_dataset'
)
TEST_FLAGS_STR_REVERSED = {flag: TEST_FLAGS_STR[flag] for flag
                           in reversed([*TEST_FLAGS_STR.keys()])}
TEST_OBJECTS = {
    'None': None,
    'id_str': TEST_HYB_ID_STR,
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
    'list': [TEST_HYB_ID_STR],
    'tuple': (TEST_HYB_ID_STR,),
    'set': {TEST_HYB_ID_STR},
    'empty_dict': {},
    'dict_rand_key': {TEST_HYB_ID_STR: True},
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

# ----- Start Testing File Paths -----
auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
test_out_dir = os.path.join(auto_tests_dir, 'output_autotest')
test_data_dir = os.path.join(auto_tests_dir, 'test_data_files')
test_hyb_file_name = os.path.join(test_data_dir, 'test_hybrid.hyb')
test_vienna_file_name = os.path.join(test_data_dir, 'test_hybrid.vienna')
test_ct_file_name = os.path.join(test_data_dir, 'test_hybrid_placeholder.ct')
# hyb_autotest_file_name = os.path.join(test_out_dir, 'test_hybrid_py_autotest.hyb')
# vienna_autotest_file_name = os.path.join(test_out_dir, 'test_vienna_py_autotest.vienna')
# make_params_autotest_file_name = os.path.join(test_out_dir, 'make_params_autotest.csv')
# out_basename = hyb_file_name.replace('.hyb', '')
# out_flag_string_name = os.path.join(test_out_dir, 'flag_string.txt')
# id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_id_map.csv')
# bad1_id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_id_map_1.csv')

# ----- Start TypeFinder Testing Data -----
STRING_MATCH_PARAMS_1 = {
    'params_str': (
        """
        # Hyb Software Reference Search Terms:
        endswith,_microRNA,microRNA
        contains,_microRNA,microRNA
        startswith,MIMAT,microRNA
        startswith,Matchtype2,MatchType
        matches,MIMAT0000078_MirBase_miR-23a_microRNA,microRNA
        matches,Matchtype1,MatchType
        endswith,_mRNA,mRNA
        """
    ),
    'params_dict': {
        'endswith': [('_microRNA', 'microRNA'), ('_mRNA', 'mRNA')],
        'contains': [('_microRNA', 'microRNA')],
        'startswith': [
            ('MIMAT', 'microRNA'),
            ('Matchtype2', 'MatchType')
        ],
        'matches': [
            ('MIMAT0000078_MirBase_miR-23a_microRNA', 'microRNA'),
            ('Matchtype1', 'MatchType')
        ],
    }
}

BAD_STRING_MATCH_PARAMS_1 = {
    'params_str': (
        """
        # Hyb Software Reference Search Terms:
        equalto,_microRNA,microRNA
        endswith,_microRNA,microRNA
        contains,_microRNA,microRNA
        startswith,MIMAT,microRNA
        startswith,Matchtype2,MatchType
        matches,MIMAT0000078_MirBase_miR-23a_microRNA,microRNA
        matches,Matchtype1,MatchType
        endswith,_mRNA,mRNA
        """
    ),
    'params_dict': {
        'equalto': [('_microRNA', 'microRNA')],
        'endswith': [('_microRNA', 'microRNA'), ('_mRNA', 'mRNA')],
        'contains': [('_microRNA', 'microRNA')],
        'startswith': [
            ('MIMAT', 'microRNA'),
            ('Matchtype2', 'MatchType')
        ],
        'matches': [
            ('MIMAT0000078_MirBase_miR-23a_microRNA', 'microRNA'),
            ('Matchtype1', 'MatchType')
        ],
    }
}

BAD_STRING_MATCH_PARAMS_2 = {
    'params_str': (
        """
        # Hyb Software Reference Search Terms:
        endswith,_microRNA
        contains,_microRNA,microRNA
        startswith,MIMAT,microRNA
        startswith,Matchtype2,MatchType
        matches,MIMAT0000078_MirBase_miR-23a_microRNA,microRNA
        matches,Matchtype1,MatchType
        endswith,_mRNA,mRNA
        """
    ),
    'params_dict': {
        'endswith': [('_microRNA', None), ('_mRNA', 'mRNA')],
        'contains': [('_microRNA', 'microRNA')],
        'startswith': [
            ('MIMAT', 'microRNA'),
            ('Matchtype2', 'MatchType')
        ],
        'matches': [
            ('MIMAT0000078_MirBase_miR-23a_microRNA', 'microRNA'),
            ('Matchtype1', 'MatchType')
        ],
    }
}

ID_MAP_PARAMS_1 = {
    'params_str': (
        """
        #seg_id,seg_type
        MIMAT0000078_MirBase_miR-23a_microRNA,microRNA
        ENSG00000188229_ENST00000340384_TUBB2C_mRNA,mRNA
        ARTSEG1_SOURCE_NAME_microRNA,microRNA
        ARTSEG2_SOURCE_NAME_microRNA,microRNA
        ARTSEG2_SOURCE_NAME_mRNA,mRNA
        ARTSEG1_SOURCE_NAME_mRNA,mRNA
        Matchtype1,MatchType
        Matchtype22,MatchType
        """
    ),
    'params_dict': {
        'MIMAT0000078_MirBase_miR-23a_microRNA': 'microRNA',
        'ENSG00000188229_ENST00000340384_TUBB2C_mRNA': 'mRNA',
        'ARTSEG1_SOURCE_NAME_microRNA': 'microRNA',
        'ARTSEG2_SOURCE_NAME_microRNA': 'microRNA',
        'ARTSEG2_SOURCE_NAME_mRNA': 'mRNA',
        'ARTSEG1_SOURCE_NAME_mRNA': 'mRNA',
        'Matchtype1': 'MatchType',
        'Matchtype22': 'MatchType',
    }
}

BAD_ID_MAP_PARAMS_1 = copy.deepcopy(ID_MAP_PARAMS_1)
BAD_ID_MAP_PARAMS_1['params_str'] = (
    BAD_ID_MAP_PARAMS_1['params_str'].replace(
        'MIMAT0000078_MirBase_miR-23a_microRNA,microRNA',
        'MIMAT0000078_MirBase_miR-23a_microRNA'
    )
)

BAD_ID_MAP_PARAMS_2 = copy.deepcopy(ID_MAP_PARAMS_1)
BAD_ID_MAP_PARAMS_2['params_dict'] = {}
BAD_ID_MAP_PARAMS_2['params_str'] = (
    BAD_ID_MAP_PARAMS_2['params_str'].replace(
        'Matchtype1,MatchType',
        'Matchtype1,MatchType\nMatchtype1,MatchTypeBAD',
    )
)

ART_HYB_MATCHTYPE_PROPS = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	.	'
        'Matchtype1	.	.	.	.	.	'
        'Matchtype22	.	.	.	.	.'
    ),
    'seg1_type': 'MatchType',
    'seg2_type': 'MatchType'
}

ART_HYB_NOTYPE_PROPS = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	.	'
        'Notype1	.	.	.	.	.	'
        'Notype2	.	.	.	.	.'
    ),
    'seg1_type': None,
    'seg2_type': None
}
