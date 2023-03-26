#!/usr/bin/env python3
# Daniel Stribling  |  ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# Hybkit Project : https://www.github.com/RenneLab/hybkit

"""
Automatic testing of hybkit code.
"""

import os
import sys
import copy
from contextlib import nullcontext as does_not_raise
import argparse
import pytest
import hybkit

# ----- Set Testing Variables -----
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
# Example Artificial Hyb record strings and parameters for eval_type and eval_mirna testing.
ART_HYB_PROPS_1 = {
    'hyb_str': (
        '1_1000	AAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG	-10.0	'
        'ARTSEG1_SOURCE_NAME_microRNA	1	20	1	20	0.001	'
        'ARTSEG2_SOURCE_NAME_microRNA	21	40	21	40	0.001	'
        'dataset=artificial;'
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
        'dataset=artificial;'
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
        'dataset=artificial;'
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
        'dataset=artificial;'
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
NONSPACE_STRS = ['id_str', 'dot_str', 'int_str', 'float_str', 'test_seq_str']
ALL_STRS = ['empty_str', 'space_str', 'tab_str', *NONSPACE_STRS]
id_allowed_types = {*NONSPACE_STRS}
seq_allowed_types = {'test_seq_str'}
energy_allowed_types = {'float_str', 'float', 'dot_str'}
seg_props_allowed_types = {'seg_props_dict', 'empty_dict'}
ref_name_allowed_types = {'None', *NONSPACE_STRS}
read_start_allowed_types = {'None', 'int_str', 'int'}
read_end_allowed_types = {'None', 'int_str', 'int'}
ref_start_allowed_types = {'None', 'int_str', 'int'}
ref_end_allowed_types = {'None', 'int_str', 'int'}
score_allowed_types = {*ALL_STRS}

# ----- Set Testing File Paths -----
auto_tests_dir = os.path.abspath(os.path.dirname(__file__))
test_out_dir = os.path.join(auto_tests_dir, 'output_autotest')
hyb_file_name = os.path.join(test_out_dir, 'test_hybrid_py.hyb')
vienna_file_name = os.path.join(test_out_dir, 'test_hybrid_py.vienna')
out_basename = hyb_file_name.replace('.hyb', '')
out_flag_string_name = os.path.join(test_out_dir, 'flag_string.txt')
match_legend_file_name = os.path.join(auto_tests_dir, 'test_string_match.csv')
bad1_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_1.csv')
bad2_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_2.csv')
bad3_match_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_string_match_3.csv')
id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_id_map.csv')
bad1_id_map_legend_file_name = os.path.join(auto_tests_dir, 'test_bad_id_map_1.csv')

# hybkit.settings.HybFile_settings['hybformat_id'] = True
# hybkit.settings.HybFile_settings['hybformat_ref'] = True
# hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
# hybkit.settings.Analysis_settings['allow_mirna_dimers'] = True

# ----- Test Assistance Functions -----
# Return a copy of empty seg_props dictionary for testing.
def get_empty_seg_props():
    return copy.deepcopy(EMPTY_SEG_PROPS)

# Generate HybRecord and FoldRecord objects for testing.
def default_hyb_records():
    """Generate two HybRecord objects for tests"""
    hyb_record_1 = hybkit.HybRecord.from_line(
        HYB_STR_1,
        hybformat_id=True,
        hybformat_ref=True,
    )
    fold_record = None
    #fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)

    hyb_record_2 = hybkit.HybRecord(
        id=hyb_record_1.id,
        seq=hyb_record_1.seq,
        seg1_props=copy.deepcopy(hyb_record_1.seg1_props),
        seg2_props=copy.deepcopy(hyb_record_1.seg2_props),
        fold_record=fold_record,
    )
    return hyb_record_1, hyb_record_2, fold_record

# Get expected result string for exception testing.
def get_expected_result_string(is_allowed=False):
    """Return string identifying expected pass/error result"""
    if is_allowed:
        return 'Pass'
    else:
        return 'Raise'
    
# Get expected result context for exception testing.
def get_expected_result_context(expect_str, error_types = (TypeError, RuntimeError)):
    """Return context for testing allowed types."""
    if expect_str == 'Pass':
        return does_not_raise()
    elif expect_str == 'Raise':
        return pytest.raises((*error_types))

# ----- HybRecord Constructor Tests - Minimal -----
def test_hybrecord_constructor_minimal():
    """Test construction of HybRecord class with minimal information."""
    # Test HybRecord Minimal Constructor:
    test_record = hybkit.HybRecord(id=TEST_HYBID_STR, seq=TEST_SEQ_STR)
    # Test "seq_id" attribute 
    assert test_record.id == TEST_HYBID_STR
    # Test "seq" attribute
    assert test_record.seq == TEST_SEQ_STR  
    # Test "seg1_props" attribute
    assert test_record.seg1_props == EMPTY_SEG_PROPS  
    # Test "seg2_props" attribute
    assert test_record.seg2_props == EMPTY_SEG_PROPS 
    # Test "flags" attribute
    assert test_record.flags == {}  
    # Test "fold_record" attribute
    assert test_record.fold_record is None  
    # Test "get_seg1_type" method
    assert test_record.get_seg1_type() == None  
    # Test "get_seg1_type" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg1_type(require=True)  
    # Test "get_seg2_type" method
    assert test_record.get_seg2_type() == None  
    # Test "get_seg2_type" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg2_type(require=True)
    # Test "get_seg_types" method
    assert test_record.get_seg_types() == (None, None)  
    # Test "get_seg_types" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_seg_types(require=True)
    # Test "get_read_count" method
    assert test_record.get_read_count() is None  
    # Test "get_read_count" method empty error 
    with pytest.raises(RuntimeError):
        test_record.get_read_count(require=True)
    # Test "get_record_count" method
    assert test_record.get_record_count() == 1
    # Test "get_record_count" method empty error
    with pytest.raises(RuntimeError):
        test_record.get_record_count(require=True)
    # Test HybRecord Minimal Constructor with missing id
    with pytest.raises(RuntimeError):
        hybkit.HybRecord(id=None, seq=TEST_SEQ_STR)
    # Test HybRecord Minimal Constructor with missing seq
    with pytest.raises(RuntimeError):
        hybkit.HybRecord(id=TEST_HYBID_STR, seq=None)


# ----- HybRecord Constructor Tests - Full -----
test_parameters = [
    ('full_as_type', TEST_HYBID_STR, TEST_SEQ_STR, TEST_ENERGY_STR, TEST_SEG_PROPS, 
     TEST_SEG_PROPS, TEST_FLAGS_OBJ, TEST_READ_COUNT),
    ('full_as_str', TEST_HYBID_STR, TEST_SEQ_STR, TEST_ENERGY_STR, TEST_SEG_PROPS_STR, 
     TEST_SEG_PROPS_STR, TEST_FLAGS_STR, TEST_READ_COUNT_STR),
    ('full_as_str_rev_flags', TEST_HYBID_STR, TEST_SEQ_STR, TEST_ENERGY_STR, TEST_SEG_PROPS_STR, 
     TEST_SEG_PROPS_STR, TEST_FLAGS_STR_REVERSED, TEST_READ_COUNT_STR),
]
@pytest.mark.parametrize("test_name,test_id,test_seq,test_energy,test_seg1_props,test_seg2_props,test_flags,test_read_count",[*test_parameters])
def test_hybrecord_constructor_full(test_name, test_id, test_seq, test_energy, test_seg1_props, 
                                    test_seg2_props, test_flags, test_read_count):
    """Test construction of HybRecord class with full complement of information."""
    # Test HybRecord Minimal Constructor:
    test_record = hybkit.HybRecord(
        id=test_id,
        seq=test_seq,
        energy=test_energy,
        seg1_props=test_seg1_props,
        seg2_props=test_seg2_props,
        flags=test_flags,
        read_count=test_read_count,
    )
    # Test "seq_id" attribute 
    assert test_record.id == TEST_HYBID_STR
    # Test "seq" attribute
    assert test_record.seq == TEST_SEQ_STR  
    # Test "energy" attribute
    assert test_record.energy == TEST_ENERGY_STR
    # Test "seg1_props" attribute
    assert test_record.seg1_props == TEST_SEG_PROPS 
    # Test "seg2_props" attribute
    assert test_record.seg2_props == TEST_SEG_PROPS
    # Test "flags" attribute
    assert test_record.flags == TEST_FLAGS_STR  
    # Test "fold_record" attribute
    assert test_record.fold_record is None  
    # Test "get_seg1_type" method
    assert test_record.get_seg1_type() == TEST_FLAGS_STR['seg1_type'] 
    assert test_record.get_seg1_type(require=True) == TEST_FLAGS_STR['seg1_type'] 
    # Test "get_seg2_type" method
    assert test_record.get_seg2_type() == TEST_FLAGS_STR['seg2_type']
    assert test_record.get_seg2_type(require=True) == TEST_FLAGS_STR['seg2_type']
    # Test "get_seg_types" method
    assert test_record.get_seg_types() == (TEST_FLAGS_STR['seg1_type'], TEST_FLAGS_STR['seg2_type'])
    assert test_record.get_seg_types(require=True) == (TEST_FLAGS_STR['seg1_type'], TEST_FLAGS_STR['seg2_type'])
    # Test "get_read_count" method
    assert test_record.get_read_count() == TEST_READ_COUNT
    assert test_record.get_read_count(require=True) == TEST_READ_COUNT
    # Test "get_record_count" method
    assert test_record.get_record_count() == TEST_RECORD_COUNT
    assert test_record.get_record_count(require=True) == TEST_RECORD_COUNT
    for flag in TEST_FLAGS_STR.keys():
        assert test_record.flags[flag] == TEST_FLAGS_STR[flag]
        assert test_record._get_flag(flag) == TEST_FLAGS_STR[flag]
        assert test_record._get_flag(flag, require=True) == TEST_FLAGS_STR[flag]
    assert test_record._make_flag_string(reorder_flags=True) == TEST_FLAGS_STR_LINE
    test_record._make_flag_string(reorder_flags=False) == TEST_FLAGS_STR_LINE

    # Test "_ensure_set" private function
    assert test_record._ensure_set('energy')

    #with open(out_flag_string_name, 'w') as test_out:
    #    test_out.write(test_record._make_flag_string(reorder_flags=True) + '\n')
    #    test_out.write(test_record._make_flag_string(reorder_flags=False) + '\n')

# ----- HybRecord Constructor Failure Tests -----
default_params = [
    TEST_HYBID_STR, TEST_SEQ_STR, TEST_ENERGY_STR, TEST_SEG_PROPS_STR, 
    TEST_SEG_PROPS_STR, TEST_FLAGS_STR, TEST_READ_COUNT_STR
]
test_parameters = []
#TEST_FLAGS_STR is 7th parameter
# Test undefined flag in constructor.
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'badflag': True}
test_parameters.append(tuple(error_params))
# Test mismatched read_count and read_count_flag
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'read_count': TEST_READ_COUNT + 1}
test_parameters.append(tuple(error_params))
# Test disallowed seg1_type flag in constructor.
error_params = [''] + copy.deepcopy(default_params)
error_params[7] = {'seg1_type': 'badtype'}
test_parameters.append(tuple(error_params))

@pytest.mark.parametrize("test_name,test_id,test_seq,test_energy,test_seg1_props,test_seg2_props,test_flags,test_read_count",[*test_parameters])
def test_hybrecord_constructor_errors(test_name, test_id, test_seq, test_energy, test_seg1_props, 
                                      test_seg2_props, test_flags, test_read_count):
    """Test construction of HybRecord class with full complement of information."""
    # Test HybRecord Minimal Constructor:
    with pytest.raises((RuntimeError, TypeError)):
        test_record = hybkit.HybRecord(
            id=test_id,
            seq=test_seq,
            energy=test_energy,
            seg1_props=test_seg1_props,
            seg2_props=test_seg2_props,
            flags=test_flags,
            read_count=test_read_count,
        )
        

# ----- HybRecord Type Tests - Main Attributes -----
test_seg_props = get_empty_seg_props()
default_constructor_args = {
    'id': TEST_HYBID_STR,
    'seq': TEST_SEQ_STR,
    'energy': TEST_ENERGY_STR,
    'seg1_props': test_seg_props,
    'seg2_props': test_seg_props,
}
field_allowed_types = {
    'id': id_allowed_types,
    'seq': seq_allowed_types,
    'energy': energy_allowed_types,
    'seg1_props': seg_props_allowed_types,
    'seg2_props': seg_props_allowed_types,
}
test_parameters = []
# Setup test constructor types for each attribute in default_constructor_args:
for constructor_field in default_constructor_args.keys():
    # Setup testing of each possible data type for each field
    for test_name, test_object in TEST_OBJECTS.items():
        # Get types allowed for this field
        allowed_types = field_allowed_types[constructor_field]
        # Setup constructor arguments for this test
        constructor_args = copy.deepcopy(default_constructor_args)
        # Set test data for this field
        constructor_args[constructor_field] = test_object
        # Determine Error vs. Null Context for this test
        expect_result = get_expected_result_string(test_name in allowed_types)
        test_param_set = (
            constructor_field, 
            test_name, 
            expect_result,
            constructor_args, 
        )
        test_parameters.append(test_param_set)

@pytest.mark.parametrize("test_field,test_name,expect_str,test_input",[*test_parameters])
def test_hybrecord_obj_types(test_field, test_name, expect_str, test_input):
    expect_context = get_expected_result_context(expect_str)
    with expect_context:
        print(test_input)
        assert hybkit.HybRecord(**test_input) is not None

# ----- HybRecord Type Tests - Segment Property Attributes -----
# Setup test constructor types for seg_props attributes:
test_seg_props = get_empty_seg_props()
default_seg_props = {
    'ref_name': 'test_ref',
    'read_start': 1,
    'read_end': 21,
    'ref_start': 1,
    'ref_end': 21,
    'score': 1.0
}
default_constructor_args = {
    'id': TEST_HYBID_STR,
    'seq': TEST_SEQ_STR,
    'energy': TEST_ENERGY_STR,
    'seg1_props': default_seg_props,
    'seg2_props': default_seg_props,
}
props_allowed_types = {
    'ref_name': ref_name_allowed_types,
    'read_start': read_start_allowed_types,
    'read_end': read_end_allowed_types,
    'ref_start': ref_start_allowed_types,
    'ref_end': ref_end_allowed_types,
    'score': score_allowed_types,
}
test_parameters = []
# Setup test constructor types for each segN_props dict
for prop_set in ['seg1_props', 'seg2_props']:
    for prop_field in props_allowed_types.keys():
        # Setup testing of each possible data type for each field
        for test_name, test_object in TEST_OBJECTS.items():
            # Get types allowed for this field
            allowed_types = field_allowed_types[constructor_field]
            # Setup constructor arguments for this test
            constructor_args = copy.deepcopy(default_constructor_args)
            seg_args = copy.deepcopy(default_seg_props)
            seg_args[prop_field] = test_object
            constructor_args[prop_set] = seg_args
            # Determine Error vs. Null Context for this test
            expect_result = get_expected_result_string(test_name in allowed_types)
            test_param_set = (
                constructor_field, 
                test_name, 
                expect_result,
                constructor_args, 
            )
            test_parameters.append(test_param_set)

@pytest.mark.parametrize("test_field,test_name,expect_str,test_input",[*test_parameters])
def test_hybrecord_obj_types_seg_props(test_field, test_name, expect_str, test_input):
    expect_context = get_expected_result_context(expect_str)
    with expect_context:
        print(test_input)
        assert hybkit.HybRecord(**test_input) is not None


# ----- HybRecord eval_type(), eval_mirna(), prop(), and mirna_detail() tests -----
test_parameters = [
    ('miRNA-miRNA', ART_HYB_PROPS_1),
    ('miRNA-coding', ART_HYB_PROPS_2),
    ('coding-miRNA', ART_HYB_PROPS_3),
    ('coding-coding', ART_HYB_PROPS_4),
]

@pytest.mark.parametrize("test_name,test_params",[*test_parameters])
def test_hybrecord_type_mirna(test_name, test_params):
    """Test Hybrecord type_eval(), eval_mirna(), mirna-associated prop(), and mirna_detail() values"""    
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
        hybformat_ref=True,
    )
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
    )
    with pytest.raises((ValueError, RuntimeError)):
        test_record.mirna_detail(detail='all', allow_mirna_dimers=True)

    test_record.eval_types()
    test_record.eval_mirna()
    # Check type setting
    assert test_record.get_seg1_type() == test_params['seg1_type']
    assert test_record.get_seg2_type() == test_params['seg2_type']

    # Check fasta-generation properties - hybrid
    all_fasta_record = test_record.to_fasta_record(mode='hybrid', annotate=False)
    assert str(all_fasta_record.seq) == test_params['seg1_seq'] + test_params['seg2_seq']
    assert str(all_fasta_record.id) == test_params['hybrid_fasta_id']
    all_fasta_record = test_record.to_fasta_record(mode='hybrid', annotate=True)
    assert str(all_fasta_record.id) == test_params['hybrid_fasta_id_annotate']
    all_fasta_record_string = test_record.to_fasta_str(mode='hybrid', annotate=True)
    assert all_fasta_record_string == all_fasta_record.format('fasta')

    # Check fasta-generation properties - seg1
    seg1_fasta_record = test_record.to_fasta_record(mode='seg1', annotate=False)
    assert str(seg1_fasta_record.seq) == test_params['seg1_seq']
    assert str(seg1_fasta_record.id) == test_params['seg1_fasta_id']
    seg1_fasta_record = test_record.to_fasta_record(mode='seg1', annotate=True)
    assert str(seg1_fasta_record.id) == test_params['seg1_fasta_id_annotate']

    # Check fasta-generation properties - seg2
    seg2_fasta_record = test_record.to_fasta_record(mode='seg2', annotate=False)
    assert str(seg2_fasta_record.seq) == test_params['seg2_seq']
    assert str(seg2_fasta_record.id) == test_params['seg2_fasta_id']
    seg2_fasta_record = test_record.to_fasta_record(mode='seg2', annotate=True)
    assert str(seg2_fasta_record.id) == test_params['seg2_fasta_id_annotate']

    # Check properties of hybrids
    for prop in test_params['true_prop_argsets']:
        assert test_record.prop(*prop)
    for prop in test_params['false_prop_argsets']:
        assert not test_record.prop(*prop)
    for prop in test_params['true_is_set_argsets']:
        assert test_record.is_set(*prop)
    for prop in test_params['false_is_set_argsets']:
        assert test_record.not_set(*prop)

    # Check mirna-based properties
    assert test_record.flags['miRNA_seg'] == test_params['miRNA_seg']
    if test_params['one_mirna_error']:
        with pytest.raises((ValueError, RuntimeError)):
            test_record.mirna_detail()
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record(mode='miRNA')
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record(mode='target')
    else:
        # Check fasta-generation properties - mirna
        mirna_fasta = test_record.to_fasta_record('miRNA', annotate=False)
        assert str(mirna_fasta.seq) == test_params['mirna_seq']
        assert str(mirna_fasta.id) == test_params['mirna_fasta_id']
        mirna_fasta = test_record.to_fasta_record('miRNA', annotate=True)
        assert str(mirna_fasta.id) == test_params['mirna_fasta_id_annotate']
        # Check fasta-generation properties - target
        target_fasta = test_record.to_fasta_record('target', annotate=False)
        assert str(target_fasta.seq) == test_params['target_seq']
        assert str(target_fasta.id) == test_params['target_fasta_id']
        target_fasta = test_record.to_fasta_record('target', annotate=True)
        assert str(target_fasta.id) == test_params['target_fasta_id_annotate']

    if test_params['two_mirna_error']:
        with pytest.raises((ValueError, RuntimeError)):
            test_record.mirna_detail(detail='all', allow_mirna_dimers=True)
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record('miRNA')
        with pytest.raises((ValueError, RuntimeError)):
            test_record.to_fasta_record('target')
    else:
        mirna_detail_dict = test_record.mirna_detail(detail='all', allow_mirna_dimers=True)
        assert mirna_detail_dict['mirna_ref'] == test_params['mirna_ref']
        assert mirna_detail_dict['target_ref'] == test_params['target_ref']
        assert mirna_detail_dict['mirna_seg_type'] == test_params['mirna_seg_type']
        assert mirna_detail_dict['target_seg_type'] == test_params['target_seg_type']
        assert mirna_detail_dict['mirna_seq'] == test_params['mirna_seq']
        assert mirna_detail_dict['target_seq'] == test_params['target_seq']

# ----- HybRecord properties tests -----
test_parameters = [
    ('miRNA-coding-props', ART_HYB_STR_PROPS),
]
@pytest.mark.parametrize("test_name,test_params",[*test_parameters])
def test_hybrecord_props(test_name, test_params):
    """Test Hybrecord type_eval(), eval_mirna(), mirna-associated prop(), and mirna_detail() values"""    
    test_record = hybkit.HybRecord.from_line(
        line=test_params['hyb_str'],
        hybformat_id=True,
    )
    test_record.eval_types()
    test_record.eval_mirna()
    # Check properties of hybrids
    for prop in test_params['true_prop_argsets']:
        assert test_record.prop(*prop)
    for prop in test_params['false_prop_argsets']:
        assert not test_record.prop(*prop)
    for prop in test_params['true_is_set_argsets']:
        assert test_record.is_set(*prop)
    for prop in test_params['false_is_set_argsets']:
        assert test_record.not_set(*prop)

# ----- HybRecord Magic Methods tests -----
def test_hybrecord_magic_methods():
    test_record_1 = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'], 
        hybformat_id=True, 
        hybformat_ref=True
    )
    test_record_2 = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_2['hyb_str'], 
        hybformat_id=True, 
        hybformat_ref=True
    )
    test_record_2.id = 'NewID'
    print(str(test_record_1))
    assert test_record_1 == test_record_1
    assert not (test_record_1 != test_record_1)
    assert test_record_1 != test_record_2
    assert bool(test_record_1)
    assert str(test_record_1)
    assert hash(test_record_1)
    assert len(test_record_1) == 40

# ----- HybRecord misc disallowed option tests -----
def test_hybrecord_misc_disallowed():
    test_record = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'], 
        hybformat_id=True, 
        hybformat_ref=True
    )

    with pytest.raises(RuntimeError):
        test_record.to_fasta_record('notallowed')
    with pytest.raises(RuntimeError):
        test_record.is_set('badprop')
    with pytest.raises(RuntimeError):
        test_record.prop('badprop')
    with pytest.raises(RuntimeError):
        test_record.prop('any_seg_type_contains', None)
    with pytest.raises(RuntimeError):
        test_record.set_fold_record(None)
    with pytest.raises(RuntimeError):
        test_record.set_fold_record('not_fold_record')
    with pytest.raises(RuntimeError):
        test_record.mirna_detail('disallowed_detail')

    #TODO
    #assert not test_record.prop('has_indels')

   # Test Private Methods
    test_record.seg1_props['read_start'] = None
    with pytest.raises(RuntimeError):
        test_record._ensure_props_read_start_end()
    with pytest.raises(RuntimeError):
        test_record._get_seg_seq(test_record.seg1_props)
    with pytest.raises(RuntimeError):
        test_record._get_flag('fake_flag', require=True)
    with pytest.raises(RuntimeError):
        test_record._make_flags_dict('not_dict')
    with pytest.raises(RuntimeError):
        test_record._make_flags_dict({'bad_flag': True})
    with pytest.raises(RuntimeError):
        test_record._parse_hybformat_id('bad_id_name_continues_on')
    with pytest.raises(RuntimeError):
        test_record._parse_hybformat_ref('bad_ref_name_continues_on')
    with pytest.raises(RuntimeError):
        test_record._read_flags('bad_flag=B;')

# ----- HybRecord misc private_function_tests -----
def test_hybrecord_misc_private():
    test_record = hybkit.HybRecord.from_line(
        ART_HYB_PROPS_1['hyb_str'], 
        hybformat_id=True, 
        hybformat_ref=True
    )
    read_flags = test_record._read_flags('bad_flag=B;', allow_undefined_flags=True)
    assert read_flags['bad_flag'] == 'B'

# ----- Begin Old Section -----

def old_test_hybfile():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)

    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        with pytest.raises(RuntimeError):
            hyb_file.write_record('NotRecord')
        hyb_file.write_record(hyb_record)
        hyb_file.write_records([hyb_record, hyb_record])
    assert hybkit.util.hyb_exists(hyb_file_name)
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file:
        hyb_record_read = hyb_file.read_records()[0]
    assert hyb_record == hyb_record_read
    for test_prop in ['id', 'seq', 'energy', 'seg1_props', 'seg2_props', 'flags']:
        print('Comparing:', test_prop,
              getattr(hyb_record, test_prop), getattr(hyb_record_read, test_prop))
        assert getattr(hyb_record, test_prop) == getattr(hyb_record_read, test_prop)


def old_test_foldrecord():
    for TestClass in [hybkit.FoldRecord, hybkit.DynamicFoldRecord]:
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_string(VIENNA_STR_1, 'bad_mode')
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['', '', ''], 'bad_mode')
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['', ''])
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', 'singleval'])
        with pytest.raises(RuntimeError):
            fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', '(99'])
        fold_record = hybkit.FoldRecord.from_vienna_lines(['abc', '123', '(99'],
                                                          error_mode='warn_return')
        fold_record = TestClass.from_vienna_string(VIENNA_STR_1)
        print(str(fold_record))
        assert fold_record == fold_record
        assert not (fold_record != fold_record)
        assert bool(fold_record)
        hash(fold_record)
        print(len(fold_record))
        if not isinstance(fold_record, hybkit.DynamicFoldRecord):
            print(fold_record._get_seg_fold({'read_start': 1, 'read_end': 10}))
        print(fold_record.to_vienna_lines())
        print(fold_record.to_vienna_lines(True))
        print(fold_record.to_vienna_string())
        print(fold_record.to_vienna_string(True))


def old_test_viennafile():
    fold_record = hybkit.FoldRecord.from_vienna_string(VIENNA_STR_1)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_direct('testval')
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        with pytest.raises(RuntimeError):
            vienna_file.write_record('NotRecord')
        vienna_file.write_record(fold_record)
        vienna_file.write_records([fold_record, fold_record])
    assert hybkit.util.vienna_exists(vienna_file_name)
    assert hybkit.util.fold_exists(vienna_file_name)
    with hybkit.FoldFile.open(vienna_file_name, 'r') as fold_file:
        with pytest.raises(RuntimeError):
            fold_file.read_record()
        with pytest.raises(RuntimeError):
            fold_file._to_record_string(fold_record, False)
    with hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        fold_record_read = vienna_file.read_records()[0]

    assert fold_record == fold_record_read
    for test_prop in ['id', 'seq', 'fold', 'energy']:
        print('Comparing:', test_prop,
              getattr(fold_record, test_prop), getattr(fold_record_read, test_prop))
        assert getattr(fold_record, test_prop) == getattr(fold_record_read, test_prop)


def old_test_hybfolditer():
    if not os.path.isdir(test_out_dir):
        os.mkdir(test_out_dir)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    with hybkit.HybFile.open(hyb_file_name, 'w') as hyb_file:
        hyb_file.write_records([hyb_record] * 30)
    fold_record = hybkit.FoldRecord.from_vienna_string(VIENNA_STR_1)
    dyn_fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)
    dyn_fold_record_mis = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1_MIS)
    with hybkit.ViennaFile.open(vienna_file_name, 'w') as vienna_file:
        vienna_file.write_records([fold_record] * 30)

    # Compare FoldRecord and HybRecord
    assert not fold_record.matches_hyb_record(hyb_record)
    fold_record.count_hyb_record_mismatches(hyb_record)
    with pytest.raises(RuntimeError):
        fold_record.ensure_matches_hyb_record(hyb_record)

    # Compare DynamicFoldRecord and HybRecord
    assert dyn_fold_record.matches_hyb_record(hyb_record)
    dyn_fold_record.count_hyb_record_mismatches(hyb_record)
    dyn_fold_record.ensure_matches_hyb_record(hyb_record)

    # Compare DynamicFoldRecord With Mismatch and HybRecord
    assert not dyn_fold_record_mis.matches_hyb_record(hyb_record)
    dyn_fold_record_mis.count_hyb_record_mismatches(hyb_record)
    with pytest.raises(RuntimeError):
        dyn_fold_record_mis.ensure_matches_hyb_record(hyb_record)

    hybkit.settings.FoldFile_settings['foldrecord_type'] = 'strict'
    hybkit.settings.HybFoldIter_settings['error_mode'] = 'raise'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        with pytest.raises(RuntimeError):
            for ret_item in hf_iter:
                print(ret_item)

    hybkit.settings.HybFoldIter_settings['error_mode'] = 'warn_skip'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        with pytest.raises(RuntimeError):
            for ret_item in hf_iter:
                print(ret_item)

    hybkit.settings.HybFoldIter_settings['error_mode'] = 'warn_return'
    with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
         hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
        hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file)
        for ret_item in hf_iter:
            print(ret_item)

    hybkit.settings.FoldFile_settings['foldrecord_type'] = 'dynamic'
    for combine in [True, False]:
        with hybkit.HybFile.open(hyb_file_name, 'r') as hyb_file, \
             hybkit.ViennaFile.open(vienna_file_name, 'r') as vienna_file:
            hf_iter = hybkit.HybFoldIter(hyb_file, vienna_file, combine=combine)
            for ret_item in hf_iter:
                print(ret_item)
            hf_iter.report()
            hf_iter.print_report()

    hyb_record, fold_record = ret_item


def old_test_analyses():
    analysis_classes = {
        'type': hybkit.analysis.TypeAnalysis,
        'mirna': hybkit.analysis.MirnaAnalysis,
        'summary': hybkit.analysis.SummaryAnalysis,
        'target': hybkit.analysis.TargetAnalysis,
    }
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    with pytest.raises(RuntimeError):
        hybkit.analysis.BaseAnalysis()
    last_analysis = None
    for analysis_type in analysis_classes:
        Analysis = analysis_classes[analysis_type]
        analysis_basename = out_basename + '_' + analysis_type
        analysis = Analysis(name='Test')
        analysis2 = Analysis()
        analysis2.plot(analysis_basename)
        analysis2.add(hyb_record)
        analysis.update(analysis2)
        if last_analysis is not None:
            with pytest.raises(RuntimeError):
                analysis.update(last_analysis)
        print(analysis.results())
        print(analysis.results(None, True))
        analysis.write(analysis_basename)
        analysis.plot(analysis_basename)
        if hasattr(analysis, 'write_individual'):
            analysis.write_individual(analysis_basename)
            analysis.plot_individual(analysis_basename)
        last_analysis = analysis

    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    type_analysis = hybkit.analysis.TypeAnalysis()
    mirna_analysis = hybkit.analysis.MirnaAnalysis()
    target_analysis = hybkit.analysis.TargetAnalysis()
    type_analysis.add(hyb_record)
    hyb_record.set_flag('miRNA_seg', '3p')
    type_analysis.add(hyb_record)
    mirna_analysis.add(hyb_record)
    hyb_record.set_flag('miRNA_seg', 'B')
    type_analysis.add(hyb_record)
    mirna_analysis.add(hyb_record)
    init_val = target_analysis.settings['allow_mirna_dimers']
    target_analysis.settings['allow_mirna_dimers'] = True
    target_analysis.add(hyb_record)
    target_analysis.settings['allow_mirna_dimers'] = init_val
    hyb_record.set_flag('miRNA_seg', 'N')
    mirna_analysis.add(hyb_record)


def old_test_fold_analysis():
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record.eval_mirna()
    fold_record = hybkit.DynamicFoldRecord.from_vienna_string(VIENNA_STR_1)
    analysis = hybkit.analysis.FoldAnalysis(name='Test')
    analysis2 = hybkit.analysis.FoldAnalysis()
    hyb_record.set_fold_record(fold_record)
    analysis_basename = out_basename + '_fold'
    analysis2.plot(analysis_basename)
    analysis2.add(hyb_record)
    analysis.update(analysis2)
    print(analysis.results())
    print(analysis.results(None, True))
    analysis.write(analysis_basename)
    analysis.plot(analysis_basename)

    hyb_record.fold_record.fold = ')))' + hyb_record.fold_record.fold[3:]
    hyb_record.set_flag('miRNA_seg', 'B')
    init_val = analysis.settings['allow_mirna_dimers']
    analysis.settings['allow_mirna_dimers'] = True
    analysis.add(hyb_record)
    analysis.settings['allow_mirna_dimers'] = init_val


def old_test_util():
    original_abspath = hybkit.settings._USE_ABSPATH
    hybkit.settings._USE_ABSPATH = True
    assert hybkit.util._bool_from_string(True)
    assert hybkit.util._bool_from_string('yes')
    assert not hybkit.util._bool_from_string('no')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util._bool_from_string('invalid')
    assert hybkit.util.dir_exists('~')
    assert hybkit.util.dir_exists('${PWD}')
    assert hybkit.util.file_exists(__file__)
    assert hybkit.util.out_path_exists(hyb_file_name)

    hybkit.util.make_out_file_name(hyb_file_name, name_suffix='out',
                                   in_suffix='.hyb', out_suffix='.new',
                                   out_dir='', seg_sep='_'
                                   )
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.dir_exists('nonexistent_dir')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.file_exists('nonexistent_file')
    with pytest.raises(argparse.ArgumentTypeError):
        hybkit.util.hyb_exists(match_legend_file_name)
    use_namespace = argparse.Namespace()
    setattr(use_namespace, 'hybformat_id', True)
    hybkit.util.set_settings(use_namespace, verbose=True)
    hybkit.settings._USE_ABSPATH = original_abspath

    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.out_hybs_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', hyb_file_name, hyb_file_name, '-o', hyb_file_name]
    )
    with pytest.raises(SystemExit):
        hybkit.util.validate_args(args, script_parser)
    parser_components = [
        hybkit.util.in_hybs_parser,
        hybkit.util.in_folds_parser,
    ]
    script_parser = argparse.ArgumentParser(parents=parser_components)
    args = script_parser.parse_args(
        ['-i', hyb_file_name, hyb_file_name, '-f', vienna_file_name]
    )
    with pytest.raises(SystemExit):
        hybkit.util.validate_args(args, script_parser)


def old_test_type_finder():
    # Generic Tests
    def do_nothing(*args, **kwargs):
        pass
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    hybkit.type_finder.TypeFinder.set_custom_method(do_nothing)

    # Defualt Hybformat
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hybkit.type_finder.TypeFinder.set_method('hybformat')
    hyb_record.eval_types()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.set_method('bad_method')
    # Non-Defualt String-Match
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params('badfile')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad1_match_legend_file_name)
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_string_match_params(bad2_match_legend_file_name)
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(
        bad3_match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.settings['check_complete_seg_types'] = True
    with pytest.raises(RuntimeError):
        hyb_record.eval_types()
    # match_params = {'startswith':[('MIMAT', 'MIMAT')], 'endswith':[('microRNA', 'miRNA')]}
    # hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    # with pytest.raises(RuntimeError):
    #     hyb_record.eval_types()
    match_params = hybkit.type_finder.TypeFinder.make_string_match_params(match_legend_file_name)
    hybkit.type_finder.TypeFinder.set_method('string_match', match_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    # Non-Defualt ID-Map
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params()
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params('wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params(type_file_pairs='wrong_id_map_type')
    with pytest.raises(RuntimeError):
        hybkit.type_finder.TypeFinder.make_id_map_params([bad1_id_map_legend_file_name])
    with pytest.raises(RuntimeError):
        id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
            type_file_pairs=[('seq1type', id_map_legend_file_name),
                             ('seq2type', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params(
        type_file_pairs=[('seqtype', id_map_legend_file_name)])
    id_map_params = hybkit.type_finder.TypeFinder.make_id_map_params([id_map_legend_file_name])
    hybkit.type_finder.TypeFinder.set_method('id_map', id_map_params)
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.eval_types()
    hyb_record = hybkit.HybRecord.from_line(HYB_STR_1)
    hyb_record.seg1_props['ref_name'] = 'not_real_name'
    with pytest.raises(RuntimeError):
        hyb_record.eval_types(allow_unknown=False)

# if __name__ == '__main__':
#     test_type_finder()
