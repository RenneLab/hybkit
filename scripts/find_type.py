#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Executable script for assignment of types for segments of ".hyb" genomic sequence.

Records are stored in the HybRecord class.
"""

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import os
import hybkit

find_seg_type_hyb = hybkit.HybRecord.find_seg_type_hyb
find_seg_type_string_match = hybkit.HybRecord.find_seg_type_string_match
make_string_match_parameters = hybkit.HybRecord.make_string_match_parameters


def assign_segment_types(in_hyb, out_hyb, **args_info):
    """Stub funciton."""
    print('stub function')


if __name__ == '__main__':

    legend = 'find_type_string_match.csv'

    mode = 'string_match'

    find_seg_type_method = find_seg_type_methods[mode]
    find_param_method = find_seg_type_param_methods[mode]

    find_type_params = find_param_method(legend)

    hybkit.HybRecord.set_find_type_method(find_seg_type_method, find_type_params)

    with hybkit.HybFile(test_file, 'r') as hf, hybkit.HybFile(test_out, 'w') as of:
        for record in hf:
            record.find_seg_types()
            of.write_record(record)
    # Implement Argparse here.
    #    args_info = parse_args()
    #    if not args_info:
    #        sys.exit()
    #    assign_segment_types(**args_info)
