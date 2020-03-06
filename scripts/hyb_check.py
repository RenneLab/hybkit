#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Read one or more '.hyb' format files and check for errors.
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
                     hybkit.util.positional_in_hybs_parser,
                     hybkit.util.gen_opts_parser,
                     hybkit.util.hybrecord_parser,
                    ]

hyb_check_parser = argparse.ArgumentParser(
                                           parents=parser_components,
                                           description=__doc__,
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                          )

# Define main script function.
def hyb_check(in_hyb_files, verbose=False, silent=False):

    if not silent:
        print('\nPerforming Error Checking of Hyb Files...')
   
    for in_hyb_file in in_hyb_files:
        if verbose:
            print('Checking File:\n    ' + in_hyb_file) 
    
        with hybkit.HybFile.open(in_hyb_file, 'r') as in_hyb:
            for record in in_hyb:
                pass

    if verbose:
        print('\nError checking complete. No errors found.\n')


# Execute the script function
if __name__ == '__main__':
    args = hyb_check_parser.parse_args() 
    verbose = True
    hybkit.HybRecord.set_namespace_settings(args, verbose=verbose)
    hyb_check(args.in_hyb, verbose=verbose)


