#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

'''
Functions and executable script for assignment of types for segments of ".hyb" genomic sequence 
records stored in the HybRecord class.
'''

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __depreciated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

import os
import hybkit

find_seg_type_hyb = hybkit.HybRecord.find_seg_type_hyb

def find_seg_type_string_match(self, seg_info, find_type_params={}):
    '''
    Return the type of the provided segment, or return None if the segment cannot be
    identified. 
    This method attempts to find a string matching a specific pattern within the identifier
    of the aligned segment. Search options include "prefix", "contains", "suffix", and "matches"
    The required find_types_param dict should contain a key for each desired search type, 
    with a list of 2-tuples for each search-string with assigned-type. For example:
    find_type_params = {'suffix': [('_miR', 'microRNA'),
                                   ('_trans', 'mRNA')   ]}
    This dict can be generated with the associated make_string_match_parameters()
    method and an associated csv legend file with format:
        #commentline
        #search_type,search_string,seg_type
        suffix,_miR,microRNA
        suffix,_trans,mRNA    
    '''
  
    seg_name = seg_info['ref']
    found_types = []
    if 'prefix' in find_type_params:
        for search_string, search_type in find_type_params['prefix']:
            if seg_name.startswith(search_string):
                found_types.append(search_type)
    if 'contains' in find_type_params:
        for search_string, search_type in find_type_params['contains']:
            if search_string in seg_name:
                found_types.append(search_type)
    if 'suffix' in find_type_params:
        for search_string, search_type in find_type_params['suffix']:
            if seg_name.endswith(search_string):
                found_types.append(search_type)
    if 'matches' in find_type_params:
        for search_string, search_type in find_type_params['prefix']:
            if search_string == seg_name:
                found_types.append(search_type)

    if not found_types:
        return None
    elif len(found_types) == 1:
        return found_types[0]
    elif len(found_types) > 1:
        #If multiple types found, check if they are the same.
        if all(((found_type == found_types[0]) for found_type in found_types)):
            return found_types[0]
        else:
            message = 'Multiple sequence types found for item: %s' % seg_name
            message += '  ' + ', '.join(found_types)
            print(message)
            raise Exception(message)
       

def return_blank_dict():
    return {}

def make_string_match_parameters(legend_file):
    '''
    Read csv file provided in legend_file, and return a dict of search parameters 
    for use with the find_seg_type_string_match function. 
    The my_legend.csv file should have the format:
        #commentline
        #search_type,search_string,seg_type
        suffix,_miR,microRNA
        suffix,_trans,mRNA   
    search_type options include "prefix", "contains", "suffix", and "matches"
    The produced dict object contains a key for each search type, with a list of 
    2-tuples for each search-string and associated segment-type. For example:
      {'suffix': [('_miR', 'microRNA'),
                  ('_trans', 'mRNA')   ]}
    '''
    ALLOWED_SEARCH_TYPES = {'prefix':'', 'contains':'', 'suffix':'', 'matches':''}
    return_dict = {}
    with open(legend_file, 'r') as legend_file_obj:
        for line in legend_file_obj:
            # Skip Blank Lines
            if not line.split():
                continue
            # Skip Commented Lines
            if line.lstrip().startswith('#'):
                continue
            line = line.rstrip()
            split_line = line.split(',')
            if len(split_line) != 3:
                message = 'Error reading legend line: \n%s\n%s\n' % (str(line), str(split_line))
                message += 'Three comma-separated entries expected.'
                print(message)
                raise Exception(message)
            search_type = split_line[0]
            search_string = split_line[1]
            seg_type = split_line[2]
            if split_line[0] not in ALLOWED_SEARCH_TYPES:
                message = 'Read Search type: "%s"\n' % split_line[0]
                message += 'Not in allowed types: %s\n' % ', '.join(ALLOWED_SEARCH_TYPES.keys())
                message += 'For legend line: \n%s\n' % (str(line))
                print(message)
                raise Exception(mesage)

            if search_type not in return_dict:
                return_dict[search_type] = []

            return_dict[search_type].append((search_string, seg_type))
    return return_dict

find_seg_type_methods = {
                         'default': find_seg_type_hyb,
                         'hyb': find_seg_type_hyb,
                         'string_match': find_seg_type_string_match,
                         # 'user_custom': user_custom_function,
                        }

find_seg_type_param_methods = {
                               'default': return_blank_dict,
                               'hyb': return_blank_dict,
                               'string_match': make_string_match_parameters,
                               # 'user_custom': make_user_custom_parameters,
                              }


def assign_segment_types(in_hyb, out_hyb, **args_info):
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
    
    
