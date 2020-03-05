#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Read a '.hyb' format file and check for errors.
"""

import sys
import os 
import argparse
import hybkit

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

# Create Command-line Argument Parser
parser_components = [
                     hybkit.util.in_hyb_parser,
                     hybkit.util.gen_opts_parser,
                     hybkit.util.hybrecord_parser,
                    ]

hyb_check_parser = argparse.ArgumentParser(
                                           parents=parser_components,
                                           description=__doc__,
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                          )

# Define main script function.
def hyb_check(in_hyb_file, verbose=False, silent=False):

    if not silent:
        print('\nPerforming Error Checking of Hyb File...')
   
    if verbose:
        print('Input File:\n    ' + in_hyb_file) 

    with hybkit.HybFile.open(in_hyb_file, 'r') as in_hyb:
        for record in in_hyb:
            pass

    if verbose:
        print('\nDone.\n')


# Execute the script function
if __name__ == '__main__':
    args = hyb_check_parser.parse_args() 
    # set settings
    hyb_check(args.in_hyb, verbose=True)


