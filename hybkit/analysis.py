#!/usr/bin/env python3
# Daniel B. Stribling
# Renne Lab, University of Florida
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""Methods for analyzing HybRecord and FoldRecord Objects.

Todo:
    Expand Documentation
    Implement "warn" paramater for target_region_database, fold parameter.
    Change database implementation to full hyb-style identifiers by default.
"""

import copy
from collections import Counter
import hybkit

# Import module-level dunder-names:
from hybkit.__about__ import __author__, __contact__, __credits__, __date__, __deprecated__, \
                             __email__, __license__, __maintainer__, __status__, __version__

# Public Constants
TYPE_ANALYSIS_KEYS = ['hybrid_type_counts', 'seg1_types', 'seg2_types', 'all_seg_types']
MIRNA_COUNT_ANALYSIS_KEYS = ['5p_mirna_hybrids', '3p_mirna_hybrids', 'mirna_dimer_hybrids', 
                             'all_mirna_hybrids', 'no_mirna_hybrids']
DEFAULT_COUNT_MODE = 'record'
DEFAULT_HYBRID_TYPE_SEP = '-'
DEFAULT_ENTRY_SEP = ','
DEFAULT_FILE_SUFFIX = '.csv'
DEFAULT_WRITE_MULTI_FILES = False
DEFAULT_TARGET_SPACER_LINE = True
DEFAULT_MAKE_PLOTS = True
DEFAULT_MAX_MIRNA = 10

# Public Methods : HybRecord Analysis Preparation : Type Analysis
def type_dict():
    """Create a dictionary with keys of counter objects for running type analyses."""
    ret_dict = {key: Counter() for key in TYPE_ANALYSIS_KEYS}
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : miRNA Count Analysis
def mirna_count_dict():
    """Create a dictionary with keys for running miRNA analyses."""
    ret_dict = {key: 0 for key in MIRNA_COUNT_ANALYSIS_KEYS}
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : miRNA Count Analysis
def full_analysis_dict():
    """Create a dictionary with keys for running both type and miRNA analyses."""
    ret_dict = type_dict()
    ret_dict.update(mirna_count_dict())
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : Type Analysis
def combine_type_dicts(analysis_dicts):
    """Combine a list/tuple of dictionaries created from running type analyses."""
    # Check that method input is formatted correctly:
    if (not (isinstance(analysis_dicts, list) or isinstance(analysis_dicts, tuple))
        or  (len(analysis_dicts) < 2)
        or  (not all(isinstance(item, dict) for item in analysis_dicts))):
        message = 'Input to "combine_type_analysis_dicts" method must be a list/tuple of two '
        message += 'or more dicts.\n'
        message += 'Current input:\n    %s' % str(analysis_dicts)
        print(message)
        raise Exception(message)

    ret_dict = copy.deepcopy(analysis_dicts[0])
    for add_dict in analysis_dicts[1:]:
        for category in TYPE_ANALYSIS_KEYS:
            for entry in add_dict[category].keys():
                # Entry checking not needed for counter objects.
                # if entry not in ret_dict[category]:
                #     ret_dict[category][entry] = 0
                ret_dict[category][entry] += add_dict[category][entry]
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : miRNA Count Analysis
def combine_mirna_count_dicts(analysis_dicts):
    """Combine a list of two or more dicts created from running miRNA count analyses."""
    # Check that method input is formatted correctly:
    if (not (isinstance(analysis_dicts, list) or isinstance(analysis_dicts, tuple))
        or  (len(analysis_dicts) < 2)
        or  (not all(isinstance(item, dict) for item in analysis_dicts))):   
        message = 'Input to "combine_type_analysis_dicts" method must be a list/tuple of two '
        message += 'or more dicts.\n'
        message += 'Current input:\n    %s' % str(analysis_dicts)
        print(message)
        raise Exception(message)

    ret_dict = copy.deepcopy(analysis_dicts[0])
    for add_dict in analysis_dicts[1:]:
        for key in MIRNA_COUNT_ANALYSIS_KEYS:
            ret_dict[key] += add_dict[key]
    return ret_dict


# Public Methods : HybRecord Analysis
def running_types(record, analysis_dict, 
                  count_mode=DEFAULT_COUNT_MODE,
                  type_sep=DEFAULT_HYBRID_TYPE_SEP, 
                  mirna_centric_sorting=True):
    """Add information regarding various properties from the HybRecord object provided in
    "record" to the dictionary provided in the analysis_dict argument.
    """
    if not record.has_property('has_seg_types'):
        message = 'seg_type flag is required for record analysis.'
        print(message)
        raise Exception(message)

    count = record.count(count_mode)
    seg1_type = record.seg1_type()
    seg2_type = record.seg2_type()

    if mirna_centric_sorting:
        mirna_types = hybkit.HybRecord.MIRNA_TYPES
        if seg1_type in mirna_types:
            join1, join2 = seg1_type, seg2_type
        elif seg2_type in mirna_types:
            join1, join2 = seg2_type, seg1_type
        else:
            join1, join2 = sorted((seg1_type, seg2_type))
        hybrid_type = type_sep.join([join1, join2])
    else:
        hybrid_type = type_sep.join(record.seg_types_sorted())

    # _add_count(analysis_dict['hybrid_type_counts'], hybrid_type)
    # _add_count(analysis_dict['seg1_types'], seg1_type)
    # _add_count(analysis_dict['seg2_types'], seg2_type)
    # Entry Checking not necessary with counter objects.
    analysis_dict['hybrid_type_counts'][hybrid_type] += count
    analysis_dict['seg1_types'][seg1_type] += count
    analysis_dict['seg2_types'][seg2_type] += count

    for seg_type in seg1_type, seg2_type:
    #    _add_count(analysis_dict['all_seg_types'], seg_type)
        analysis_dict['all_seg_types'][seg_type] += 1

# Public Methods : HybRecord Analysis
def running_mirna_counts(record, analysis_dict,
                         count_mode=DEFAULT_COUNT_MODE):
    """Add information regarding various properties to the dictionary provided in
    the analysis_dict argument.
    """
    record._ensure_mirna_analysis()
    count = record.count(count_mode)

    # Add mirna-analysis specific details
    if record.has_property('has_mirna_dimer'):
        analysis_dict['mirna_dimer_hybrids'] += count
        analysis_dict['all_mirna_hybrids'] += count
    elif record.has_property('3p_mirna'):
        analysis_dict['3p_mirna_hybrids'] += count
        analysis_dict['all_mirna_hybrids'] += count
    elif record.has_property('5p_mirna'):
        analysis_dict['5p_mirna_hybrids'] += count
        analysis_dict['all_mirna_hybrids'] += count
    else:
        analysis_dict['no_mirna_hybrids'] += count


# Public Methods : HybRecord Analysis
def running_full(record, analysis_dict,
                 count_mode=DEFAULT_COUNT_MODE,
                 type_sep=DEFAULT_HYBRID_TYPE_SEP,
                 mirna_centric_sorting=True):
    """Add information regarding various properties to the dictionary provided in
    the analysis_dict argument.
    """
    running_types(record, analysis_dict, count_mode, type_sep, mirna_centric_sorting)
    running_mirna_counts(record, analysis_dict, count_mode)


# Public Methods : HybRecord Type Analysis Parsing
def format_types(analysis_dict, sep=DEFAULT_ENTRY_SEP):
    """Return the results of a type_analysis in a list of sep-delimited lines."""
    ret_lines = []
    ret_lines += _format_hybrid_type_counts(analysis_dict, sep)
    ret_lines.append('')
    ret_lines += _format_all_seg_types(analysis_dict, sep)
    return ret_lines


# Public Methods : HybRecord Type Analysis Writing
def write_types(file_name_base, analysis_dict,
                name=None,
                multi_files=DEFAULT_WRITE_MULTI_FILES,
                sep=DEFAULT_ENTRY_SEP, 
                file_suffix=DEFAULT_FILE_SUFFIX,
                make_plots=DEFAULT_MAKE_PLOTS):
    """Write the results of the type-analysis to a file or series of files with names based
    on file_name_base.
    """
    analyses = [
                ('types_hybrids', _format_hybrid_type_counts),
                ('types_segs', _format_all_seg_types),
                ]

    if multi_files:
        for analysis_name, analysis_method in analyses:
            analysis_file_name = file_name_base + '_' + analysis_name + file_suffix
            with open(analysis_file_name, 'w') as out_file:
                out_file.write('\n'.join(analysis_method(analysis_dict, sep)))
    else:
        write_lines = []
        analysis_file_name =  file_name_base + '_type_combined' + file_suffix
        for analysis_name, analysis_method in analyses:
            write_lines += analysis_method(analysis_dict, sep)
        with open(analysis_file_name, 'w') as out_file:
                out_file.write('\n'.join(write_lines))

    if make_plots:
        hybkit.plot.hybrid_type_counts(analysis_dict, file_name_base + '_types_hybrids', name=name)
        hybkit.plot.all_seg_types(analysis_dict, file_name_base + '_types_seg', name=name)

# Public Methods : HybRecord miRNA Count Analysis Parsing
def format_mirna_counts(analysis_dict, sep=DEFAULT_ENTRY_SEP):
    """Return the results of mirna analysis in a list of sep-delimited lines."""
    ret_lines = ['miRNA_type' + sep + 'count']
    for key in MIRNA_COUNT_ANALYSIS_KEYS:
        ret_lines.append('%s%s%i' % (key, sep, analysis_dict[key]))
    return ret_lines

# Public Methods : HybRecord miRNA Count Analysis Writing
def write_mirna_counts(file_name, analysis_dict,
                       name=None,
                       multi_files=DEFAULT_WRITE_MULTI_FILES, # Currently a dummy option
                       sep=DEFAULT_ENTRY_SEP,
                       file_suffix=DEFAULT_FILE_SUFFIX):
    """Write the results of the mirna-analysis to the file provided in file_name."""
    use_file_name = file_name
    if not use_file_name.endswith(file_suffix):
        use_file_name += file_suffix
    with open(use_file_name, 'w') as out_file:
        out_file.write('\n'.join(format_mirna_counts(analysis_dict, sep)))

    if make_plots:
        hybkit.plot.mirna_counts(analysis_dict, file_name + '_mirna_counts', name=name)
    

# Public Methods : HybRecord Analysis Preparation : Type Analysis
def combine_full_dicts(analysis_dicts):
    """Combine a list/tuple of dictionaries created from running full analyses."""
    # Check that method input is formatted correctly:
    if (not (isinstance(analysis_dicts, list) or isinstance(analysis_dicts, tuple))
        or  (len(analysis_dicts) < 2)
        or  (not all(isinstance(item, dict) for item in analysis_dicts))):
        message = 'Input to "combine_full_dicts" method must be a list/tuple of two '
        message += 'or more dicts.\n'
        message += 'Current input:\n    %s' % str(analysis_dicts)
        print(message)
        raise Exception(message)

    ret_dict = copy.deepcopy(analysis_dicts[0])
    for add_dict in analysis_dicts[1:]:
        for category in TYPE_ANALYSIS_KEYS:
            for entry in add_dict[category].keys():
                ret_dict[category][entry] += add_dict[category][entry]

    for add_dict in analysis_dicts[1:]:
        for key in MIRNA_COUNT_ANALYSIS_KEYS:
            ret_dict[key] += add_dict[key]

    return ret_dict


# Public Methods : HybRecord Multiple Analysis Writing
def write_full(file_name_base, analysis_dict,
               name=None,
               multi_files=DEFAULT_WRITE_MULTI_FILES,             
               sep=DEFAULT_ENTRY_SEP,
               file_suffix=DEFAULT_FILE_SUFFIX,
               make_plots=DEFAULT_MAKE_PLOTS):
    """Write the results of the full analysis to a file or series of files with names based
    on file_name_base.
    """
    analyses = [
                ('types_hybrids', _format_hybrid_type_counts),
                ('types_segs', _format_all_seg_types),
                ('mirna_counts', format_mirna_counts),
               ]
    
    if multi_files:
        for analysis_name, analysis_method in analyses:
            analysis_file_name = file_name_base + '_' + analysis_name + file_suffix
            with open(analysis_file_name, 'w') as out_file:
                out_file.write('\n'.join(analysis_method(analysis_dict, sep)))
    else:
        analysis_file_name = file_name_base + '_analysis' + file_suffix
        write_lines = []
        for analysis_name, analysis_method in analyses:
            write_lines += analysis_method(analysis_dict, sep)
        with open(analysis_file_name, 'w') as out_file:
            out_file.write('\n'.join(write_lines))

    if make_plots:
        hybkit.plot.hybrid_type_counts(analysis_dict, file_name_base + '_types_hybrids', name=name)
        hybkit.plot.all_seg_types(analysis_dict, file_name_base + '_types_seg', name=name)
        hybkit.plot.mirna_counts(analysis_dict, file_name_base + '_mirna_counts', name=name)

# Public Methods : HybRecord Analysis Preparation : miRNA Target Analysis
def mirna_target_dict():
    """Create a dictionary with keys for running target analyses."""
    # Created for parallel methods, potential future expansion.
    ret_dict = {}
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : miRNA Target Analysis
def combine_mirna_target_dicts(analysis_dicts):
    """Combine a list/tuple of two or more dictionaries created from running miRNA target analyses."""
    # Check that method input is formatted correctly:
    if (not (isinstance(analysis_dicts, list) or isinstance(analysis_dicts, tuple))
        or  (len(analysis_dicts) < 2)
        or  (not all(isinstance(item, dict) for item in analysis_dicts))):
        message = 'Input to "combine_mirna_target_dicts" method must be '
        message += 'a list/tuple of two or more dicts.\n'
        message += 'Current input:\n    %s' % str(analysis_dicts)
        print(message)
        raise Exception(message)

    ret_dict = copy.deepcopy(analysis_dicts[0])
    for add_dict in analysis_dicts[1:]:
        for key in add_dict:
            # Existence checking not required with counter objects.
            # if key not in ret_dict:
            #     ret_dict[key] = 0
            ret_dict[key] += add_dict[key]
    return ret_dict


# Public Methods : HybRecord Analysis
def running_mirna_targets(record, analysis_dict, 
                          count_mode=DEFAULT_COUNT_MODE,
                          double_count_duplexes=False,
                          mirna_contains=None, mirna_matches=None,
                          target_contains=None, target_matches=None):
    """Add information regarding mirna/target properties to the dictionary provided in
    the analysis_dict argument.
    If double_count_duplexes is provided as True, miRNA-miRNA duplexes will be counted 
    in both orientations. If False, the 3p miRNA will be considered as the target.
    """
    record._ensure_mirna_analysis()
    count = record.count(count_mode)

    check_pairs = []
    if record.has_property('has_mirna'):
        mirna_target_pair = ((record.mirna_info['ref'], record.mirna_details['mirna_seg_type']),
                            (record.target_info['ref'], record.mirna_details['target_seg_type']))
        check_pairs.append(mirna_target_pair)

        if double_count_duplexes and record.has_property('has_mirna_dimer'):

            check_pairs.append((mirna_target_pair[1], mirna_target_pair[0]))
    
    for ((mirna, mirna_seg_type), (target, target_seg_type)) in check_pairs:
        if any(check is not None for check in (mirna_contains, mirna_matches, 
                                               target_contains, target_matches)):
            use_pair = False
            if ((mirna_contains is not None and mirna_contains in mirna)
                or (mirna_matches is not None and mirna_matches == mirna)
                or (target_contains is not None and target_contains in target)
                or (target_matches is not None and target_matches == target)):
                use_pair = True      
        else:
            use_pair = True

        if use_pair:
            mirna_id = (mirna, mirna_seg_type)
            target_id = (target, target_seg_type) 
            if mirna_id not in analysis_dict:
            #    analysis_dict[mirna] = {}
                analysis_dict[mirna_id] = Counter()
            # Existence checking not required for counter objects.
            # _add_count(analysis_dict[mirna], target)
            analysis_dict[mirna_id][target_id] += count


# Public Methods : HybRecord miRNA Target Analysis Parsing
def process_mirna_targets(analysis_dict):
    """process and sort the results of target analysis"""
    counts = {} 
    target_type_counts = {}
    ret_dict = {}
    for mirna_id in sorted(analysis_dict.keys()):
        # Get count of all targets of a particular mirna
        counts[mirna_id] = sum(analysis_dict[mirna_id].values())
        #for key, count in analysis_dict[mirna_id].most_common():
        #    print(key, count)
        target_type_counts[mirna_id] = Counter()
        for target_id in analysis_dict[mirna_id].keys():
            target, target_seg_type = target_id
            target_type_counts[mirna_id][target_seg_type] += analysis_dict[mirna_id][target_id]
        # For each miRNA, sort targets by target count.
        # targets_by_count = sorted(analysis_dict[mirna_id].items(), key=lambda item: item[1])
        # Use Counter Implementation:
        targets_by_count = analysis_dict[mirna_id].most_common() 
        # Add sorted target dict to ret_dict
        mirna_ret_dict = {target_id: target_count for target_id, target_count in targets_by_count}
        ret_dict[mirna_id] = Counter(mirna_ret_dict)

    total_count = sum(counts.values())
    return (ret_dict, counts, target_type_counts, total_count)


# Public Methods : HybRecord miRNA Target Analysis Parsing
def format_mirna_targets(analysis_dict, counts=None, 
                         sep=DEFAULT_ENTRY_SEP, 
                         spacer_line=DEFAULT_TARGET_SPACER_LINE):
    """Return the results of target analysis in a list of sep-delimited lines."""
    header_items = ['mirna', 'mirna_type', 'target', 'target_type', 'count']
    ret_lines = []
    ret_lines.append(sep.join(header_items))
    for mirna_id in analysis_dict:
        mirna, mirna_seg_type = mirna_id
        if counts is not None:
            line_items = ['mirna', '-', 'total', '-', str(counts[mirna_id])]
            ret_lines.append(sep.join(line_items))
        for target_id in analysis_dict[mirna_id]:
            target, target_seg_type = target_id
            count = str(analysis_dict[mirna_id][target_id])
            line_items = [mirna, mirna_seg_type, target, target_seg_type, count]
            ret_lines.append(sep.join(line_items))
        if spacer_line:
            ret_lines.append('')
    return ret_lines


# Public Methods : HybRecord miRNA Target Analysis Writing
def write_mirna_targets(file_name_base, analysis_dict,
                        counts_dict=None, 
                        target_types_count_dict=None,
                        name=None,
                        multi_files=DEFAULT_WRITE_MULTI_FILES,
                        sep=DEFAULT_ENTRY_SEP,
                        file_suffix=DEFAULT_FILE_SUFFIX,
                        spacer_line=DEFAULT_TARGET_SPACER_LINE,
                        make_plots=DEFAULT_MAKE_PLOTS,
                        max_mirna=DEFAULT_MAX_MIRNA):
    """Write the results of the mirna_target to a file or series of files with names based
    on file_name_base.
    """
    if multi_files and len(analysis_dict) > max_mirna:
        message = 'miRNA-specific individual output files not supported for > 10 miRNA, with '
        message += '%i miRNA to be written in this case.\n' % len(analysis_dict)
        message += 'Please write fewer output files, or write a combined output with '
        message += 'multi_files=False.'
        print(message)
        raise Exception(message)

    if multi_files:
        for mirna_id in analysis_dict:
            mirna, mirna_seg_type = mirna_id
            analysis_file_name = file_name_base + '_' + mirna + file_suffix
            one_mirna_dict = Counter({mirna_id:analysis_dict[mirna_id]})
            if counts_dict is not None:
                one_counts_dict = {mirna_id:counts_dict[mirna_id]}
            else:
                one_counts_dict = None
            write_lines = format_mirna_targets(one_mirna_dict, one_counts_dict,
                                               sep=sep, spacer_line=False)
            with open(_sanitize_name(analysis_file_name), 'w') as out_file:
                out_file.write('\n'.join(write_lines))

            if make_plots:
                target_plot_name = _sanitize_name(file_name_base + '_' + mirna)
                hybkit.plot.mirna_targets(mirna, 
                                          analysis_dict[mirna_id], 
                                          target_plot_name,
                                          name=name,
                                          )

                target_type_plot_name = _sanitize_name(file_name_base + '_' + mirna + '_types')
                hybkit.plot.mirna_target_types(mirna,
                                               target_types_count_dict[mirna_id],
                                               target_type_plot_name,
                                               name=name
                                              )

    else:
        analysis_file_name = _sanitize_name(file_name_base + '_' + 'mirna' + file_suffix)
        analysis = format_mirna_targets(analysis_dict, counts_dict, sep, spacer_line)
        with open(analysis_file_name, 'w') as out_file:
            out_file.write('\n'.join(analysis))
        
        if make_plots:
            print('Plotting Not Supported for combined miRNA output')


# Public Methods : HybRecord Analysis Preparation : miRNA Fold Analysis
def mirna_fold_dict():
    """Create a dictionary with keys for running fold analyses."""
    ret_dict = {}
    ret_dict['base_counts'] = Counter()
    for i in range(1,28):
        ret_dict['base_counts'][i] = 0
    ret_dict['base_percent'] = None
    ret_dict['all_evaluated'] = 0
    ret_dict['all_mirna'] = 0
    ret_dict['no_mirna'] = 0
    ret_dict['all_folds'] = 0
    ret_dict['no_folds'] = 0
    return ret_dict


# Public Methods : HybRecord Analysis Preparation : miRNA Fold Analysis
def combine_mirna_fold_dicts(analysis_dicts):
    """Combine a list/tuple of two or more counter objects created from running miRNA fold analyses."""
    # Check that method input is formatted correctly:
    if (not (isinstance(analysis_dicts, list) or isinstance(analysis_dicts, tuple))
        or  (len(analysis_dicts) < 2)
        or  (not all(isinstance(item, dict) for item in analysis_dicts))):
        message = 'Input to "combine_mirna_fold_dicts" method must be '
        message += 'a list/tuple of two or more dicts.\n'
        message += 'Current input:\n    %s' % str(analysis_dicts)
        print(message)
        raise Exception(message)

    ret_dict = copy.deepcopy(analysis_dicts[0])
    ret_dict['base_fractions'] = None
    for add_dict in analysis_dicts[1:]:
        for key in analysis_dict['base_counts']:
            ret_dict['base_counts'][key] = add_dict['base_counts'][key]
        for key in ['all_evaluated', 'all_mirna', 'no_mirna', 'all_folds', 'no_folds']:
            ret_dict[key] += add_dict[key]
    return ret_dict


# Public Methods : HybRecord Analysis : miRNA Fold Analysis
def running_mirna_folds(record, analysis_dict, 
                        count_mode=DEFAULT_COUNT_MODE,
                        allow_duplexes=False,
                        skip_no_fold_record=False):
    """Add miRNA fold analysis to provided dict.

    If allow_duplexes is provided as True, miRNA-miRNA duplexes will be counted 
    considering the 5p miRNA as the "miRNA" in the hybrid.
    """
    record._ensure_mirna_analysis()
    count = record.count(count_mode, as_int=True)

    analysis_dict['all_evaluated'] += count

    if not record.has_property('has_mirna'):
        analysis_dict['no_mirna'] += count
    elif record.has_property('has_mirna'):  # For Readability
        analysis_dict['all_mirna'] += count

        if not record.has_property('has_fold_record'):
            if not skip_no_fold_record:
                message = 'Problem performing running miRNA fold analysis for record:\n'
                message += '    %s\n' % str(record)
                message += 'Record has no fold-record attribute.'
                print(message)
                raise Exception(message)
            else:
                analysis_dict['no_folds'] += count            
                return

        # Presume that record.mirna_details['mirna_fold'] is populated, as it should be:
        mirna_fold = record.mirna_details['mirna_fold']  

        # Segments should not have both fold directions, or neither fold direction.
        if (('(' in mirna_fold and ')' in mirna_fold)               
            or ('(' not in mirna_fold and ')' not in mirna_fold)
            or len(mirna_fold) < 5):
            analysis_dict['no_folds'] += count
            message = 'WARNING: running_mirna_fold: Record: %s, ' % str(record)
            message += 'Bad Fold: %s' % mirna_fold
            print(message)
            return 
        
        analysis_dict['all_folds'] += count
        for str_i in range(len(mirna_fold)):
            seq_i = str_i + 1
            if mirna_fold[str_i] in {'(', ')'}:
                analysis_dict['base_counts'][seq_i] += count


# Public Methods : HybRecord miRNA Fold Analysis Parsing
def process_mirna_folds(analysis_dict):
    """Process results of mirna fold analysis."""
    total_count = analysis_dict['all_folds']
    if total_count == 0:
        message = 'Problem with mirna fold analysis, total mirna counted is 0.'
        print(message)
        raise Exception(message)
    test_sums = [(('all_mirna', 'no_mirna'), 'all_evaluated'),
                 (('all_folds', 'no_folds'), 'all_mirna')]

    for (test_1, test_2), ex_sum in test_sums:
        if sum([analysis_dict[test_1], analysis_dict[test_2]]) != analysis_dict[ex_sum]:
            message = 'Problem with processing of mirna_folds.'
            message += 'Check of items: %s (%i) + ' % (test_1, analysis_dict[test_1])
            message += '%s (%i) != ' % (test_2, analysis_dict[test_2])
            message += '%s (%i)' % (ex_sum, analysis_dict[ex_sum])
            print(message)
            raise Exception(message)

    analysis_dict['base_fractions'] = {}
    for i, count in analysis_dict['base_counts'].items():
        analysis_dict['base_fractions'][i] = count/total_count
    return analysis_dict


# Public Methods : HybRecord miRNA Fold Analysis Parsing
def format_mirna_folds(analysis_dict,
                         sep=DEFAULT_ENTRY_SEP):
    """Return the results of mirna_fold analysis in a list of sep-delimited lines."""
    header_items = ['data_type', 'count']
    ret_lines = []
    ret_lines.append(sep.join(header_items))
    for count_type in ['all_evaluated', 'all_mirna', 'no_mirna', 'all_folds', 'no_folds']:
        ret_lines.append(sep.join([count_type, str(analysis_dict[count_type])]))

    ret_lines.append('')
    line_values = ['index']
    line_values +=  [str(key) for key in analysis_dict['base_counts'].keys()]
    ret_lines.append(sep.join(line_values))
    line_values = ['i_count']
    line_values += [str(val) for val in analysis_dict['base_counts'].values()]
    ret_lines.append(sep.join(line_values))
    
    if analysis_dict['base_fractions'] is not None:
        line_values = ['i_fraction']
        line_values += [str(val) for val in analysis_dict['base_fractions'].values()]
        ret_lines.append(sep.join(line_values))

    return ret_lines


# Public Methods : HybRecord miRNA Fold Analysis Writing
def write_mirna_folds(file_name_base, analysis_dict,
                        name=None,
                        multi_files=DEFAULT_WRITE_MULTI_FILES,
                        sep=DEFAULT_ENTRY_SEP,
                        file_suffix=DEFAULT_FILE_SUFFIX,
                        spacer_line=DEFAULT_TARGET_SPACER_LINE,
                        make_plots=DEFAULT_MAKE_PLOTS):
    """Write the results of the mirna_fold analysis to a file or series of files with names based
    on file_name_base.
    """

    write_lines = format_mirna_folds(analysis_dict, sep=sep)
    if not multi_files:
        fold_details_name = file_name_base + '_fold_info' + file_suffix
        with open(fold_details_name, 'w') as fold_details_file:
            for line in write_lines:
                fold_details_file.write(line + '\n')

    elif multi_files:
        fold_counts_name = file_name_base + '_fold_counts' + file_suffix
        with open(fold_counts_name, 'w') as fold_counts_file:
            for line in write_lines[:6]:
                fold_counts_file.write(line + '\n')

        fold_base_counts_name = file_name_base + '_fold_base_counts' + file_suffix
        with open(fold_base_counts_name, 'w') as fold_base_counts_file:
            for line in write_lines[7:9]:
                fold_base_counts_file.write(line + '\n')

        if analysis_dict['base_fractions'] is not None:
            fold_base_fractions_name = file_name_base + '_fold_base_fractions' + file_suffix
            with open(fold_base_fractions_name, 'w') as fold_base_fractions_file:
                fold_base_fractions_file.write(write_lines[7] + '\n')
                fold_base_fractions_file.write(write_lines[9] + '\n')

    if make_plots:
        hybkit.plot.mirna_folds(analysis_dict, 
                                file_name_base,
                                name=name,
                                )


# Private Methods : Utility
def _add_count(count_dict, key):
    if key not in count_dict:
        count_dict[key] = 0
    count_dict[key] += 1


def _sanitize_name(file_name):
    for char, replace in [('*', 'star'), (',','com')]:
        file_name = file_name.replace(char, replace)
    return file_name

# Private Methods : HybRecord Type Analysis Parsing
def _format_hybrid_type_counts(analysis_dict, sep=DEFAULT_ENTRY_SEP):
    """Return the results of hybrid_type_counts in a list of sep-delimited lines."""
    # Sort by count in descending order
    ret_lines = ['hybrid_type' + sep + 'count']
    # sorted_pairs = sorted(analysis_dict['hybrid_type_counts'].items(),
    #                       key=lambda item: item[1], reverse=True)
    # Use Counter sorting:
    sorted_pairs =  analysis_dict['hybrid_type_counts'].most_common()
    ret_lines += ['%s%s%i' % (key, sep, count) for (key, count) in sorted_pairs]
    return ret_lines


# Private Methods : HybRecord Type Analysis Parsing
def _format_all_seg_types(analysis_dict, sep=DEFAULT_ENTRY_SEP):
    """Return the results of all_seg_types in a list of sep-delimited lines."""
    # Sort by count in descending order
    #sorted_pairs = sorted(analysis_dict['all_seg_types'].items(),
    #                      key=lambda item: item[1], reverse=True)
    # Use Counter sorting:
    sorted_pairs =  analysis_dict['all_seg_types'].most_common()

    ret_lines = ['seg_type' + sep + 'count']
    ret_lines += ['%s%s%i' % (key, sep, count) for (key, count) in sorted_pairs]
    return ret_lines
