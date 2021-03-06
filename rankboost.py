#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

from collections import namedtuple

import argparse
import logging
import math
import os
import pandas
import re
import time
import yaml

_log = logging.getLogger(__name__)

# Set this such that printing of alleles is not truncated
pandas.options.display.max_colwidth = 300

# A named tuple used during the boosting process. It contains all the necessary values for boosting.
Boost = namedtuple('Boost', (
    # The number of peptides defining the IAR.
    'npa',
    # The number of high-ranking peptides defining the IAR.
    'nph',
    # The number of MHCs reacting to peptides within the IAR.
    'nMHC',
    # The median and maximum number of MHCs stimulated by expressed IARs in the patient
    'med_max_MHC',
    # The expression of the isoform producing the IAR.
    'TPM',
    # The median number of MHCs stimulated by expressed IARs in the patient
    'medianTPM',
    # The number of 9-mers overlapping with 10-mers in the IAR.
    'overlap',
    # The number of tumor score >= wt score + delta events
    'tndelta'
    ))
# Update this every time a new argument is added to boost
base_boost = {'npa': 0, 'nph': 0, 'nMHC': 0, 'med_max_MHC': 0, 'TPM': 0, 'medianTPM': 0,
              'overlap': 0, 'tndelta': 0}


def read_fasta(pepfile):
    """
    Reads a fasta and returns it as a dictionary.
    :param pepfile: Input peptides file
    :return: Dictionary with peptide names as keys and sequence and value.
    :rtype: dict
    """
    _log.info('Reading peptide fasta from "%s".', pepfile)
    peptides = {}
    name = None
    with open(pepfile, 'r') as pfh:
        for line in pfh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
            else:
                peptides[name] = line
    _log.info('Processed %s peptides.', len(peptides))
    return peptides


class EmptyInputException(Exception):
    __module__ = Exception.__module__

    def __init__(self, filetype, filename):
        self.filetype = filetype
        message = 'Input %s file (%s) was empty.' % (filetype, filename)
        super(EmptyInputException, self).__init__(message)


class NoExpressedActionableGenesException(Exception):
    __module__ = Exception.__module__

    def __init__(self, filtering, filter_level=None, fusions=None):
        if filtering:
            message = ('After filtering for only transcripts expressed at %s%% or more than the '
                       'expressed transcripts in the patient, we could find actionable genes.' %
                       (filter_level * 100.0))
        else:
            message = ('There were no expressed genes found in the input peptides file before any '
                       'filters were applied. ')
        if fusions:
            message += (' Also, a fusion gene should have been expressed since we caught it in '
                        'the RNA.')
        message += ('You can lower the threshold on the command line using '
                    '--expression_filter_threshold (or -e) and  try again. ''This lack of genes '
                    'could either mean there truly are no expressed actionable genes in the '
                    'patient\'s tumor, or that the RNA-seq quality was poor.')
        super(NoExpressedActionableGenesException, self).__init__(message)


def get_ref_pos_alt_aa(aa_change):
    """
    Determine the reference AA(s), the position in the protein, and the alternate AA(s) for a given
    mutation in the form A123B

    :param aa_change: A String describing a mutation in protein space
    :return: reference AA(s), position, Alternamte AA(s)
    :rtype: tuple(str, int, str)
    """
    # We expect the mutation to be 3-parts, non-numeric, numeric, non-numeric or NA
    expected_isdigit = False
    parts = ['']
    for char in aa_change:
        if char.isdigit() is expected_isdigit:
            parts[-1] += char
        else:
            parts.append(char)
            expected_isdigit = not expected_isdigit
    parts[1] = int(parts[1])
    if len(parts) == 2:
        parts.append('')
    return parts


def get_overlap(df):
    """
    Get the overlap
    :param pandas.core.frame.DataFrame df:
    :return: number of overlap events
    :rtype: int

    >>> df = pandas.DataFrame({'allele': ['A', 'A', 'B', 'A', 'B', 'A', 'B'], \
                                          'peptide_seq': ['ABCD', 'ABCE', 'ZYXW', \
                                                          'ABCDE',         'YXWVU', \
                                                                'ABABCE', 'ZYXWVU' ], \
                                          'peplen': [4, 4, 4, 5, 5, 6, 6]})
    >>> get_overlap(df)
    4
    >>> df = pandas.DataFrame({'allele': ['A', 'A', 'B'], \
                                          'peptide_seq': ['ABCD', 'ABCE', 'ZYXW'], \
                                          'peplen': [4, 4, 4]})
    >>> get_overlap(df)
    0
    >>> df = pandas.DataFrame({'allele': ['A', 'A', 'B'], \
                                          'peptide_seq': ['ABCD', 'ABCE', 'ZYXWV'], \
                                          'peplen': [4, 4, 5]})
    >>> get_overlap(df)
    0
    """
    # These two line swill ensure that the first group is the largest
    df = df.sort_values(by='peplen', ascending=False)
    in_stats = df.groupby('peplen', sort=False)
    if len(in_stats.groups) == 1:
        return 0
    overlaps = 0
    for outer_group_name, outer_group_vals in in_stats:
        for inner_group_name, inner_group_vals in in_stats:
            if inner_group_name >= outer_group_name:
                # Only process groups smaller than the current largest (the outer one)
                continue
            overlaps += sum(inner_group_vals.apply(lambda x:
                                                   len([y for y in outer_group_vals.itertuples() if
                                                        x['peptide_seq'] in y.peptide_seq and
                                                        x['allele'] == y.allele]),
                                                   axis=1))
    return overlaps


def get_vaf(mutations):
    """
    Given the list of mutations in the form
        <TRANSCRIPT_1>_X123Y#0.56,<TRANSCRIPT_2>_X456Y#0.91,etc  -> for SNVS

        <5'TRANSCRIPT_1>-<3'TRANSCRIPT_1>_FUSION_Junction:X-Spanning:Y,\
        <5'TRANSCRIPT_2>-<3'TRANSCRIPT_2>_FUSION_Junction:A-Spanning:B,etc  -> for FUSIONS

    return it
    in the form 0.XX for an SNV or 0.0 for a FUSION

    :param str mutations: The mutations covered by the input IAR
    :return: The VAF
    :rtype: float

    >>> get_vaf('ENST1231.1_S123K#0.56')
    0.56
    >>> get_vaf('ENST1231.1_S123K#0.56,ENST1211.1_S143K#0.61')
    0.56
    >>> get_vaf('ENST1231.1_S123K#0.43_K124S#0.61,ENST1211.1_S143K#0.43_K144S#0.61')
    0.43
    >>> get_vaf('ENST1231.1-ENST4564.2_FUSION_Junction:5-Spanning:10')
    0.0
    """
    muts = mutations.split(',')
    vaf = 1.0
    for mutation in muts:
        if 'FUSION' in mutation:
            return 0.0
        mutation = mutation.split('_')[1:]

        for mut in mutation:
            vaf = min(vaf, float(mut.split('#')[1]))
    return vaf


def get_mutations(mutations):
    """
    Given the list of mutations in the form
        <TRANSCRIPT_1>_X123Y#0.56,<TRANSCRIPT_2>_X456Y#0.91,etc  -> for SNVS

        <5'TRANSCRIPT_1>-<3'TRANSCRIPT_1>_FUSION_Junction:X-Spanning:Y,\
        <5'TRANSCRIPT_2>-<3'TRANSCRIPT_2>_FUSION_Junction:A-Spanning:B,etc  -> for FUSIONS

    return it
    in the form X>Y or FUSION

    :param str mutations: The mutations covered by the input IAR
    :return: The mutations in the reduced form
    :rtype: str

    >>> get_mutations('ENST1231.1_S123K#0.56')
    'S>K'
    >>> get_mutations('ENST1231.1_S123K#0.56,ENST1211.1_S143K#0.56')
    'S>K'
    >>> get_mutations('ENST1231.1_S123K#0.56_K124S#0.56,ENST1211.1_S143K#0.56_K144S#0.56')
    'S>K+K>S'
    >>> get_mutations('ENST1231.1_S123K#0.56_K125S#0.56,ENST1211.1_S143K#0.56_K145S#0.56')
    'S>K+X+K>S'
    >>> get_mutations('ENST1231.1_S123K#0.56_K126S#0.56,ENST1211.1_S143K#0.56_K146S#0.56')
    'S>K+2X+K>S'
    >>> get_mutations('ENST1231.1_S123K#0.56,ENST1211.1_S143K#0.56_K146S#0.56')
    'S>K+2X+K>S'
    >>> get_mutations('ENST1231.1-ENST4564.2_FUSION_Junction:5-Spanning:10')
    'FUSION'
    >>> get_mutations('ENST1231.1_S123K#0.56_K125S#0.56,ENST1211.1_S143K#0.56_K147S#0.56')
    Traceback (most recent call last):
    ...
    AssertionError

    """
    muts = mutations.split(',')
    out_muts = []
    for mutation in muts:
        if 'FUSION' in mutation:
            return 'FUSION'
        mutation = mutation.split('_')[1:]
        temp_mutation = []
        prev_pos = 0
        for mut in mutation:
            ref, curr_pos, alt = get_ref_pos_alt_aa(mut.split('#')[0])
            if temp_mutation:
                distance = (curr_pos - prev_pos - 1)
                append_string = ('' if distance == 0 else
                                 'X+' if distance == 1 else
                                 str(distance) + 'X+')
                temp_mutation.append('+%s' % append_string)
            temp_mutation.append(ref + '>' + alt)
            prev_pos = curr_pos
        out_muts.append(''.join(temp_mutation))
    # In case of neoepitope collapse, there might be multiple versions of mutations due to
    # mutations near splicing locations
    # In such case, we have to check if all mutations are subsets of the longest
    if len(set(out_muts)) > 1:
        longest_fusion = max(out_muts, key=len)
        assert (all(s in longest_fusion for s in out_muts))
        # If all mutations are subsets, return the longest mutation
        collapsed_mut = longest_fusion
    else:
        collapsed_mut = out_muts[0]
    return collapsed_mut


def get_expression(exp_type, rsem_table, transcripts):
    """
    Given the list of mutations in the form
       <TRANSCRIPT_1>_X123Y,<TRANSCRIPT_2>_X456Y,etc.    -> for SNVs
        <5'TRANSCRIPT_1>-<3'TRANSCRIPT_1>_FUSION_Junction:X-Spanning:Y,\
        <5'TRANSCRIPT_2>-<3'TRANSCRIPT_2>_FUSION_Junction:A-Spanning:B,etc  -> for FUSIONS

    return the combined expression in TPM or FPKM

    :param str exp_type: return TPM or FPKM?
    :param pandas.core.frame.DataFrame rsem_table: The isoform level expression table from rsem
    :param str transcripts: The mutations covered by the input IAR
    :return: The expression in TPM or FPKM
    :rtype: float

    >>> rsem_table=pandas.DataFrame({'transcript_id': ['A', 'B'], 'TPM': [1, 2], 'FPKM': [4, 5]})
    >>> get_expression('TPM', rsem_table, 'A_S123K')
    1.0
    >>> get_expression('FPKM', rsem_table, 'A_S123K,B_S123K')
    9.0
    >>> get_expression('TPM', rsem_table, 'A-B_FUSION_Junction:1:spanning=2')
    1.0
    >>> get_expression('FPKM', rsem_table, 'D_S123K')
    0.0

    """
    assert exp_type in ('TPM', 'FPKM')
    transcripts = transcripts.split(',')
    if 'FUSION' in transcripts[0]:
        # TODO Implement a better estimate for expression of the fusion product
        splitter = '-'
    else:
        splitter = '_'
    transcripts = [x.split(splitter)[0] for x in transcripts]
    expression = 0.0
    for transcript in transcripts:
        expn_table = rsem_table[rsem_table['transcript_id'] == transcript][exp_type]
        if expn_table.empty:
            continue
        else:
            expression += expn_table[expn_table.first_valid_index()]
    return expression


def get_minimal_iar_seq(iar_name, peptides, iars):
    """
    Get the minimal IAR sequence for the peptides from the input IAR from the input IARs file

    :param str iar_name: The name of the iar
    :param pandas.core.series.Series peptides: The mutations covered by the input IAR
    :param dict iars: The input IARs
    :return: The sequence of the IAR
    :rtype: str

    >>> df = pandas.Series(['ABCD', 'ABCE', 'CDEA'])
    >>> get_minimal_iar_seq('seq1', df, {'seq1': 'DABCDEABCEEEE'})
    'ABCDEABCE'
    >>> get_minimal_iar_seq('seq1', df, {'seq1': 'DABCDEABCE'})
    'ABCDEABCE'
    >>> get_minimal_iar_seq('seq1', df, {'seq1': 'ABCDEABCE'})
    'ABCDEABCE'
    >>> df = pandas.Series(['ABCD', 'ABCDE'])
    >>> get_minimal_iar_seq('seq1', df, {'seq1': 'ABCDEABCE'})
    'ABCDE'
    """
    iar = iars[iar_name]
    first = len(iar)
    last = 0
    end = 0
    for pept in peptides:
        pos = iar.find(pept)
        if pos < first:
            first = pos
        if pos >= last:
            last = pos
            peplen = len(pept)
            end = pos + peplen if pos + peplen > end else end
    return iar[first:end]


def get_normal_comparisons(group_vals, delta=1.5):
    """
    Given combinations of tumor and wt peptide binding scores, how many are significantly different
    (i.e. are in at least delta percent apart)

    :param pandas.core.frame.DataFrame group_vals: The mutations covered by the input IAR
    :param float delta: The difference between tumor and wt that is considered as significant
    :return: The number of tum > normal events satisfying the given delta

    >>> df = pandas.DataFrame({'normal_binding_score': [2,3.5,5,10], \
                               'binding_score': [1,2,3,pandas.np.nan]})
    >>> get_normal_comparisons(df, 1)
    3
    >>> get_normal_comparisons(df)
    2
    >>> get_normal_comparisons(df, 2)
    1
    >>> get_normal_comparisons(df, 10)
    0
    >>> get_normal_comparisons(df, 0)
    3
    """
    result = 0
    for row in group_vals.itertuples():
        if row.normal_binding_score - row.binding_score >= delta:
            result += 1
    return result


def prepare_dataframe(mhc, predfile, expfile, pepfile, efilt_threshold, boost):
    """
    Prepares an input dataframe for rank_boost using the data in in `predfile` , `expfile` and
    `pepfile`.

    :param str mhc: MHCI or MHCII?
    :param str predfile: The input predictions file.
    :param str expfile: The input isoform-level expression file.
    :param str pepfile: The input peptides file.
    :param float efilt_threshold: The expression filter threshold.
    :param Boost boost: The relative importance of different boosting criteria.
    :return: None
    """
    assert mhc in ('mhci', 'mhcii')
    _log.info('Preparing the dataframe for boosting ranks.')
    peptides = read_fasta(pepfile)
    if len(peptides) == 0:
        raise EmptyInputException('peptides', pepfile)
    _log.info('Reading in the merged MHC:peptide binding predictions from "%s".', predfile)
    try:
        preds = pandas.read_table(predfile, nrows=1)
    except pandas.io.common.EmptyDataError:
        raise EmptyInputException('predictions', predfile)
    if len(preds.columns) == 9:
        names = ['allele', 'peptide_seq', 'peptide_name', 'peptide_core', 'IC_50', 'binding_score',
                 'ENSEMBL_gene', 'HUGO_gene', 'mutations']
        normals = False
    else:
        if len(preds.columns) == 11:
            names = ['allele', 'peptide_seq', 'normal_seq', 'peptide_name', 'peptide_core', 'IC_50',
                     'binding_score', 'normal_binding_score', 'ENSEMBL_gene', 'HUGO_gene',
                     'mutations']
        elif len(preds.columns) == 10:
            names = ['allele', 'peptide_seq', 'peptide_name', 'peptide_core', 'IC_50',
                     'binding_score', 'normal_binding_score', 'ENSEMBL_gene', 'HUGO_gene',
                     'mutations']
        else:
            raise RuntimeError('The input predictions file is formatted badly.')
        normals = True
    preds = pandas.read_table(predfile, names=names)
    _log.info('Reading in the expression data from "%s".', expfile)
    try:
        expn = pandas.read_table(expfile)
    except pandas.io.common.EmptyDataError:
        raise EmptyInputException('expression', expfile)
    # Add a column to preds for the length of the peptide
    preds['peplen'] = [len(seq) for seq in preds['peptide_seq']]
    # Split the input predictions by the IAR
    all_data = preds.groupby('peptide_name')
    # Populate the stats data frame
    # Add minimum binding score
    stats = all_data.aggregate({'binding_score': min})
    # Add HUGO gene name
    stats['HUGO_gene'] = all_data.aggregate(
        {'HUGO_gene': lambda h: h[h.first_valid_index()]})
    # Add ENSEMBL gene name
    stats['ENSEMBL_gene'] = all_data.aggregate(
        {'ENSEMBL_gene': lambda e: e[e.first_valid_index()]})
    # Add total number of peptides in the IAR
    stats['num_pept'] = all_data.aggregate({'peptide_seq': len})
    # Add total number of MHCs binding to peptides from the IAR
    stats['num_MHC'] = all_data.aggregate({'allele': lambda a: len(set(a))})
    # Add the names of the binding MHCs
    stats['binding_MHCs'] = all_data.aggregate({'allele': lambda bm: ','.join(sorted(set(bm)))})
    if mhc == 'mhci':
        # Add total number of overlap events
        # can't use df.apply because of a bug
        # See https://github.com/pandas-dev/pandas/issues/13568
        overlap_array = []
        for group_name, group_vals in all_data['peptide_seq', 'allele', 'peplen']:
            overlap_array.append(get_overlap(group_vals))
        stats['overlap'] = overlap_array
    else:
        stats['overlap'] = [0] * len(stats)
    if normals:
        tn_array = []
        for group_name, group_vals in all_data['binding_score', 'normal_binding_score']:
            tn_array.append(get_normal_comparisons(group_vals))
        stats['T/N'] = tn_array
    else:
        stats['T/N'] = [0] * len(stats)
    # Add reduced mutation information for the IAR
    #stats['mutations'] = all_data.aggregate({'mutations': get_mutations})
    # Add the sequence of the IAR
    iars = []
    for pep_name, pep_seq in all_data['peptide_seq']:
        iars.append(get_minimal_iar_seq(pep_name, pep_seq, peptides))
    stats['IAR_sequence'] = iars
    # Add mutation information for the IAR
    stats['mutations'] = all_data.aggregate({'mutations': lambda m: m[m.first_valid_index()]})
    # Collapse duplicates
    stats = collapse_iars(stats)
    # Add expression info
    stats['TPM'] = stats['mutations'].apply(lambda m: get_expression('TPM', expn, m))
    # stats['FPKM'] = all_data.aggregate({'mutations':
    #                                        lambda x: get_expression('FPKM', expn, x)})
    # Add reduced mutation information for the IAR
    stats['mutations'] = stats['mutations'].apply(get_mutations)
    if stats[stats['TPM'] > 0.0].empty:
        raise NoExpressedActionableGenesException(filtering=False,
                                                  fusions=not stats[stats.mutations == 'FUSION'
                                                                    ].empty)
    median_expressed_transcript_tpm = expn[expn['TPM'] > 0.0]['TPM'].median()
    # Keep only the transcripts expressed at higher than 10% of the median expressed-transcript
    # expression.
    # TODO: This needs to be tweaked
    temp_stats = stats[(stats.mutations != 'FUSION') &
                      (stats.TPM > efilt_threshold * median_expressed_transcript_tpm)]
    # Fusions will never be filtered
    stats = temp_stats.append(stats[stats.mutations == 'FUSION'])
    if stats[stats['TPM'] > 0.0].empty:
        raise NoExpressedActionableGenesException(filtering=True, filter_level=efilt_threshold,
                                                  fusions=not stats[stats.mutations == 'FUSION'
                                                                    ].empty)
    # Order the final dataframe before returning it
    stats['pept/mhc'] = [round(x, 2) for x in stats['num_pept']/stats['num_MHC']]
    stats = stats.sort_values(by=['binding_score', 'pept/mhc', 'TPM'], ascending=[1, 0, 0])
    _log.info('Prepared the dataframe for boosting ranks.')
    median_mhc = int(stats['num_MHC'].median())
    max_mhc = stats['num_MHC'].max()
    tmp_boost = Boost(npa=0, nph=0, nMHC=0, med_max_MHC=(median_mhc, max_mhc), TPM=0,
                      medianTPM=median_expressed_transcript_tpm, overlap=0, tndelta=0)
    boost = Boost(*[sum(x) if type(x[1]) != tuple else x[1] for x in zip(boost, tmp_boost)])
    _log.info('Finished creating the data frame for rank boosting..')
    return all_data, stats, boost


def collapse_iars(stats):
    """
    Collapse duplicated IARs.

    :param pandas.core.frame.DataFrame stats: The stats dataframe generated from parsing the input
           calls
    :return: A collapsed dataframe where overlapping IARS are merged into a single call
    :rtype: pandas.core.frame.DataFrame

    >>> df = pandas.DataFrame({'IAR_sequence': ['ABCDE', 'FGHI',\
                                                'ABCD',         '12345', \
                                                'ABCDE', 'FGHIJ',\
                                                'ABC',                  '7891011'], \
                                'mutations': ['A,B', 'G', 'A', '1', 'A,C', 'G,H', 'A,D', '8']}, \
                                 index=['ne1', 'ne2', 'ne3', 'ne4', 'ne5', 'ne6', 'ne7', 'ne8'])
    >>> tst = collapse_iars(df)
    >>> tst.index
    Index([u'ne1', u'ne2', u'ne4', u'ne8'], dtype='object')
    >>> tst.loc['ne1', 'mutations']
    'A,B,C,D'
    >>> tst.loc['ne2', 'mutations']
    'G,H'
    >>> tst.loc['ne4', 'mutations']
    '1'
    >>> tst.loc['ne8', 'mutations']
    '8'
    >>> tst.loc['ne1', 'epitope_name(s)']
    'ne1,ne3,ne5,ne7'
    >>> tst.loc['ne2', 'epitope_name(s)']
    'ne2,ne6'
    >>> tst.loc['ne4', 'epitope_name(s)']
    'ne4'
    >>> tst.loc['ne8', 'epitope_name(s)']
    'ne8'
    """
    # Add new column for identical neoepitopes and for tagging rows to drop
    stats['epitope_name(s)'] = stats.index
    stats['collapsed'] = False
    for epitope1 in stats.index:
        # Get the IAR to compare against other IARs
        iar_to_check = stats.loc[epitope1]['IAR_sequence']
        # Create sets for extending names and mutations
        joined_names = list()
        joined_names.append(stats.loc[epitope1, 'epitope_name(s)'])
        joined_mutations = set(stats.loc[epitope1, 'mutations'].split(','))
        first = True
        for epitope2 in stats[epitope1:].index:
            # First epitope1 == epitope2, skip it
            if first:
                first = False
                continue
            # If the IAR is the same or is a part of a different IAR while not being the same
            # neoepitope, it has to be collapsed
            if (not stats.loc[epitope2, 'collapsed'] and (iar_to_check in
                    stats.loc[epitope2]['IAR_sequence'] or stats.loc[epitope2][
                    'IAR_sequence'] in iar_to_check)):
                # Add the neoepitope name to the extended group
                joined_names.append(stats.loc[epitope2, 'epitope_name(s)'])
                # Add the mutation to the extended group
                joined_mutations.update(stats.loc[epitope2, 'mutations'].split(','))
                # If one of the IARs was shorter than the other, add it to the set of neoepitopes to
                # drop
                if len(iar_to_check) < len(stats.loc[epitope2]['IAR_sequence']):
                    # Update IAR_sequence to the longer one
                    stats.loc[epitope1, 'IAR_sequence'] = stats.get_value(epitope2, 'IAR_sequence')
                stats.loc[epitope2, 'collapsed'] = True
        # Update the extended group and mutations
        joined_names = ','.join(sorted(joined_names))
        stats.loc[epitope1, 'epitope_name(s)'] = joined_names
        joined_mutations = ','.join(sorted(joined_mutations))
        stats.loc[epitope1, 'mutations'] = joined_mutations
    # Drop the tagged neoepitopes
    stats = stats[stats.collapsed != True]
    return stats


def boost_nph(group, nph):
    """
    Return a fraction of boost based on total number of high-binding peptides per IAR.

    :param pandas.core.frame.DataFrame group: The IAR under review.
    :param float nph: The total boost amount allotted to this criterion.
    :returns: The percent of boost provided by this criterion.
    :rtype: float

    >>> df = pandas.DataFrame({'binding_score': [0.1, 0.4, 1.5, 0.04, 0.44, 1.75, 1.0]})
    >>> boost_nph(df, 1)
    1.0
    >>> df = pandas.DataFrame({'binding_score': [0.1, 0.4, 1.5, 1.75, 1.0]})
    >>> boost_nph(df, 1)
    0.9
    >>> df = pandas.DataFrame({'binding_score': [0.1, 0.4, 1.5, 1.75]})
    >>> boost_nph(df, 1)
    0.7
    >>> df = pandas.DataFrame({'binding_score': [0.1, 1.5]})
    >>> boost_nph(df, 1)
    0.4
    >>> df = pandas.DataFrame({'binding_score': [1.5]})
    >>> boost_nph(df, 1)
    0.0
    """
    n = len(group[group['binding_score'] <= 1.0])
    return round(nph * ((n >= 1) * 0.4 + (n >= 2) * 0.3 + (n >= 3) * 0.2 + (n >= 4) * 0.1), 2)


def boost_npa(group, npa):
    """
    Return a fraction of boost based on total number of peptides per IAR.

    :param pandas.core.frame.DataFrame group: The IAR under review.
    :param float npa: The total boost amount allotted to this criterion.
    :returns: The percent of boost provided by this criterion.
    :rtype: float

    >>> df = pandas.DataFrame({'binding_score': [0.1, 0.4, 1.5, 0.04, 0.44, \
                                                 0.1, 0.4, 1.5, 0.04, 0.44]})
    >>> boost_npa(df, 1)
    0.4
    >>> df = pandas.DataFrame({'binding_score': [0.1, 0.4, 1.5, 0.04, 0.44, \
                                                 0.1, 0.4, 1.5, 0.04, 0.44,  \
                                                 0.1, 0.4, 1.5, 0.04, 0.44, \
                                                 0.1, 0.4, 1.5, 0.04, 0.44]})
    >>> boost_npa(df, 1)
    0.9
    """
    n = len(group)
    return round(npa * ((n >= 10) * 0.4 + (n >= 15) * 0.3 + (n >= 20) * 0.2 + (n >= 30) * 0.1), 2)


def boost_nmhc(group, nmhc, med_max_mhc):
    """
    Return a fraction of boost based on total number of MHCs stimulated by peptides from the IAR.

    :param group::param pandas.core.frame.DataFrame group: The IAR under review.
    :param float nmhc: The total boost amount allotted to this criterion.
    :param tuple med_max_mhc: A tuple of median and maximum number of stimulated MHCs over candidate
          IARs.
    :returns: The percent of boost provided by this criterion.
    :rtype: float

    >>> df = pandas.DataFrame({'allele': ['A', 'A', 'B', 'A', 'B', 'A', 'B']})
    >>> boost_nmhc(df, 0.5, (3, 5))
    0.0
    >>> df = pandas.DataFrame({'allele': ['A', 'A', 'B', 'A', 'B', 'A', 'C']})
    >>> boost_nmhc(df, 1, (3, 5))
    0.17
    >>> boost_nmhc(df, 1, (2, 4))
    0.5
    >>> boost_nmhc(df, 1, (3, 3))
    1.0
    >>> boost_nmhc(df, 1, (2, 3))
    1.0
    """
    boost_thresholds = [x - med_max_mhc[0] + 1 for x in range(med_max_mhc[0], med_max_mhc[1] + 1)]
    boost_thresholds = [round(1.0 * x/sum(boost_thresholds), 2) for x in boost_thresholds]
    mhcs_stimulated = len(set(group.allele))
    return round(nmhc * sum([val for ind, val in enumerate(boost_thresholds) if
                             mhcs_stimulated - med_max_mhc[0] >= ind]), 2)


def boost_tpm(group_tpm, tpm, med_tpm):
    """
    Return a fraction of boost based on the expression of the IAR.

    :param float group_tpm: The TPM of the IAR.
    :param float tpm: The total boost amount allotted to this criterion.
    :param float med_tpm: The median TPM for all expressed transcripts in the patient.
    :returns: The percent of boost provided by this criterion.
    :rtype: float

    >>> boost_tpm(10, 1, 20)
    0.0
    >>> boost_tpm(10, 1, 10)
    0.35
    >>> boost_tpm(10, 1, 5)
    0.6
    >>> boost_tpm(10, 1, 2)
    0.75
    >>> boost_tpm(10, 1, 1)
    1.0
    """
    return round(tpm * ((group_tpm >= med_tpm) * 0.35 +
                        (group_tpm >= 2 * med_tpm) * 0.25 +
                        (group_tpm >= 5 * med_tpm) * 0.15 +
                        (group_tpm >= 10 * med_tpm) * 0.25), 2)


def boost_overlap(group_overlap, overlap):
    """
    Return a fraction of boost based on the number of overlap events in the IAR.

    :param int group_overlap: The number of overla events in the IAR.
    :param overlap: The total boost amount allotted to this criterion
    :returns: The percent of boost provided by this criterion
    :rtype: float

    >>> boost_overlap(0, 1)
    0.0
    >>> boost_overlap(1, 1)
    0.3
    >>> boost_overlap(2, 1)
    0.55
    >>> boost_overlap(3, 1)
    0.75
    >>> boost_overlap(4, 1)
    0.9
    >>> boost_overlap(5, 1)
    1.0
    """
    return round(overlap * ((group_overlap >= 1) * 0.3 +
                            (group_overlap >= 2) * 0.25 +
                            (group_overlap >= 3) * 0.20 +
                            (group_overlap >= 4) * 0.15 +
                            (group_overlap >= 5) * 0.10), 2)


def boost_ranks(mhc, all_data, dataframe, boost):
    """
    Boost the ranks in the input data frame based on biologically relevant criteria.

    :param str mhc: MHCI or MHCII?
    :param pandas.core.groupby.DataFrameGroupBy all_data: All the data used to generate dataframe
    :param pandas.DatafFrame dataframe: The dataframe produced by :module: `prepare_dataframe`
    :param Boost boost: The relative importance of different boosting criteria.
    :return: A sorted data frame
    :rtype: pandas.DataFrame
    """
    _log.info('Boosting ranks using %s.', repr(boost))
    assert mhc in ('mhci', 'mhcii')
    dataframe['epitope_name'] = dataframe.index
    dataframe['naive_rank'] = range(1, len(dataframe) + 1)
    dataframe['boosted_rank'] = range(1, len(dataframe) + 1)
    dataframe.index = list(dataframe.boosted_rank)
    dataframe.index.name = 'index'

    # Sanity check variables
    max_rank = dataframe.boosted_rank.max()
    expected_sum_of_ranks = (max_rank * (max_rank + 1)) / 2

    # 3-pass approach
    for p in range(0, 3):
        dataframe['temp_rank'] = dataframe.boosted_rank.copy(deep=True)
        for index in dataframe.index:
            if index == 1:
                # We don't process the first index
                continue
            total_boost = 0.0
            # boost by all high peptides
            group = all_data.get_group(dataframe.epitope_name[index])
            total_boost += boost_npa(group, boost.npa)
            total_boost += boost_nph(group, boost.nph)
            total_boost += boost_nmhc(group, boost.nMHC, boost.med_max_MHC)
            total_boost += boost_tpm(dataframe['TPM'][index], boost.TPM, boost.medianTPM)
            total_boost += boost_overlap(dataframe['overlap'][index], boost.overlap)

            total_boost = int(math.floor(index * total_boost))
            if total_boost == 0:
                continue
            new_rank = index - total_boost
            # Swap ranks
            dd = dataframe.boosted_rank[(dataframe.boosted_rank >= new_rank) &
                                        (dataframe.boosted_rank < index)]
            dataframe.loc[(dataframe.boosted_rank >= new_rank) &
                          (dataframe.boosted_rank < index), 'boosted_rank'] = [x + 1 for x in
                                                                               list(dd)]
            dataframe.loc[index, 'boosted_rank'] = new_rank
        assert dataframe.boosted_rank.sum() == expected_sum_of_ranks
        dataframe = dataframe.sort_values(by='boosted_rank', ascending=True)
        dataframe.index = list(dataframe.boosted_rank)
        dataframe.index.name = 'index'
    # temp_rank has no more use
    dataframe.pop('temp_rank')
    _log.info('Successfully boosted ranks in the input data.')
    return dataframe


def write_results(mhc, all_data=None, boosted_ranks=None):
    """
    Write the results of runnin rank boost on the sample.

    :param mhc: MHCI or MHCII?
    :param pandas.core.groupby.DataFrameGroupBy all_data: All the data used to generate
          boosted_ranks
    :param boosted_ranks: The data frame created by :meth:`boost_ranks`.
    :returns: None
    """
    _log.info('Writing output files.')
    columns = ['HUGO_gene', 'mutations', 'epitope_name(s)', 'IAR_sequence', 'binding_score',
               'num_pept', 'num_MHC']
    columns.extend(['overlap'] if mhc == 'mhci' else [])
    columns.extend(['T/N', 'TPM', 'naive_rank', 'boosted_rank', 'binding_MHCs'])

    concise_file_name = os.path.join(os.getcwd(), ''.join([mhc, '_rankboost_concise_results.tsv']))
    detailed_file_name = os.path.join(os.getcwd(), ''.join([mhc,
                                                            '_rankboost_detailed_results.txt']))
    if boosted_ranks is None:
        boosted_ranks = pandas.DataFrame({}, columns=columns)
        open(detailed_file_name, 'w').close()  # Emulates system touch
    else:
        # Only print a detailed file if there are results
        with open(detailed_file_name, 'w') as f:
            for row in boosted_ranks.itertuples():
                print('#', '\t'.join([row.HUGO_gene, row.mutations, row.epitope_name,
                                      row.ENSEMBL_gene]), sep='',
                      file=f)
                epitope_data = all_data.get_group(row.epitope_name)
                epitope_data = epitope_data.sort_values(by=['binding_score', 'peplen'],
                                                        ascending=[1, 0])
                if 'normal_binding_score' not in epitope_data:
                    epitope_data.normal_binding_score = 'NA'
                print(epitope_data.to_string(index=False, header=True,
                                             columns=['allele', 'peptide_seq', 'binding_score',
                                                      'normal_binding_score', 'peplen', 'mutations'
                                                      ]), file=f)
                print('', file=f)
    # Always print a concise file.  Even if there are no results because the input was empty.
    boosted_ranks.to_csv(concise_file_name, sep='\t', header=True, columns=columns, index=False)


def main():
    """
    Rankboost is a program to intelligently rank an input list of peptides predicted to elicit a
    T-Cell response in tumors based on biological features.
    """
    start = time.time()
    parser = argparse.ArgumentParser(prog='Rankboost',
                                     description=main.__doc__,
                                     epilog='Contact Arjun Rao (aarao@ucsc.edu) if you encounter '
                                     'any problems while running Rankboost')
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--mhci', dest='mhci', help='Run analysis on MHCI predictions.',
                        action='store_true')
    inputs.add_argument('--mhcii', dest='mhcii', help='Run analysis on MHCII predictions.',
                        action='store_true')
    parser.add_argument('--predictions', '-P', dest='prediction_file', help='The input predictions '
                        'file that is created by parsing the results of running the IEDB '
                        'MHC-prediction suite on all alleles in the patient. The input format is '
                        '"allele", "peptide_seq", "peptide_name", "core", "ic50", "binding_score", '
                        '["normal_binding_score"], "ensembl_gene", "HUGO_gene", "mutations" -- '
                        'where the square brackets describe and optional field.',
                        type=str,
                        required=True)
    parser.add_argument('--expression', '-E', dest='expression_file', help='The isoform-level '
                        'expression obtained through RSEM (or in the same format as RSEM).',
                        required=True)
    parser.add_argument('--peptides', '-F', dest='peptide_file', help='The 10- or 15-mer peptides '
                        'file generated using Transgene (9- and 10-mer predictions use the 10-mer '
                        'fasta and the 15-mer uses the 15-mer fasta.')
    parser.add_argument('--expression_filter_threshold', '-f', dest='efilt_threshold',
                        help='The percentage of expressed below which we must remove candidate '
                             'genes.', type=float, default=10, required=False)
    parser.add_argument('--ratios', '-r', dest='ratios', help='The relative ratios for the '
                        'different ranking criteria to use during boosting.  The values will be '
                        'converted to fractions. The format is a json-like string of elements '
                        'having keys "npa" (number of peptides that make up an IAR), "nph" (number '
                        'of high-binders that make up an IAR), "nMHC" (the number of MHCs binding '
                        'to peptides from the IAR), "TPM" (Expression in TPM), "overlap" (The '
                        'number of 9/10-mer overlap events -- MHCI only), "tndelta" (The number of '
                        'tumor score >= wt score + delta events).', type=str, default='{"npa": 1, '
                        '"nph": 1, "nMHC": 1, "TPM": 1, "overlap":1, "tndelta": 1}')
    parser.add_argument('--log_level', dest='log_level', help='The level of logging above which '
                        'messages should be printed.', required=False,
                        choices={'DEBUG', 'INFO', 'WARNING', 'ERROR'}, default='INFO')
    params = parser.parse_args()

    logging.basicConfig(level=getattr(logging, params.log_level), format='%(levelname)s: '
                                                                         '%(message)s')
    assert os.path.exists(params.prediction_file)
    params.prediction_file = os.path.abspath(params.prediction_file)
    assert os.path.exists(params.expression_file)
    params.expression_file = os.path.abspath(params.expression_file)
    assert os.path.exists(params.peptide_file)
    params.peptide_file = os.path.abspath(params.peptide_file)
    assert params.efilt_threshold < 100.0
    params.efilt_threshold = round(params.efilt_threshold/100.0, 2)
    try:
        params.ratios = yaml.load(params.ratios)
    except ValueError:
        _log.error('Is your ratios argument properly json-formatted?')
        raise
    assert not(set(params.ratios.keys()) - set(base_boost.keys()))
    base_boost.update(params.ratios)
    params.ratios = base_boost
    for x, y in params.ratios.iteritems():
        assert isinstance(y, (float, int)), 'problem parsing %s:%s in ratios' % (x, y)
    sum_of_ratios = sum(params.ratios.values())
    assert sum_of_ratios != 0.0, 'Cannot boost without ratios.'
    for i, v in params.ratios.items():
        v = round(0.55 * v / sum_of_ratios, 2)
        params.ratios[i] = v
    boost = Boost(**params.ratios)
    mhc = 'mhci' if params.mhci else 'mhcii'
    boosted_ranks = all_data = None
    try:
        all_data, dataframe, boost = prepare_dataframe(mhc,
                                                       params.prediction_file,
                                                       params.expression_file,
                                                       params.peptide_file,
                                                       params.efilt_threshold,
                                                       boost)
    except EmptyInputException as e:
        _log.warning(e.message)
        if e.filetype == 'predictions':
            _log.warning('No input binding predictions received for boosting.')
        else:
            raise
    except NoExpressedActionableGenesException as e:
        _log.warning(e.message)
    else:
        boosted_ranks = boost_ranks(mhc,
                                    all_data,
                                    dataframe,
                                    boost)
    finally:
        write_results(mhc, all_data, boosted_ranks)
        time_taken = time.time() - start
        _log.info('Rankboost completed. Took %s minutes and %s seconds to run.',
                  int(time_taken / 60), round(time_taken % 60, 2))


if __name__ == '__main__':
    main()
