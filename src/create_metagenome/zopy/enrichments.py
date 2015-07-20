#! /usr/bin/env python

import os, sys, re, glob, argparse
from numpy import median
from scipy.stats import fisher_exact, mannwhitneyu
from zopy.fdr import pvalues2qvalues

def generate_background ( annotations ):
    """ extract the list of members from an annotations [term][member] nested dictionary """
    background = {}
    for term, dictMembers in annotations.items():
        for m in dictMembers:
            background[m] = 1
    return background

def fisher_enrich ( sample, annotations, depletions=True,
                    background=None, restrict=False, min_fold=None, min_overlap=None, fdr=None ):
    """
    Given a sample of elements ( members ), and an annotation mapping [term][member] ( nested dictionary )
    Compute significant enrichments among the members
    : If not specified, background is all annotated terms
    : If not specified, analysis is "open" ( terms can appear in sample/annotations not in background )
    : Default is not FDR corrected for multiple tests
    : Default is positive enrichment of 1.1+, pvalue < 0.05, and overlap must be at least two members
    """
    if not background:
        background = generate_background( annotations )
    if restrict:
        # restrict sample and annotation space to background; may result in empty annotations ( which are removed )
        sample = [k for k in sample if k in background]
        annotations = {term:{k:1 for k in annotations[term] if k in background} for term in annotations}
        annotations = {term:dictMembers for term, dictMembers in annotations.items() if len( dictMembers ) > 0}
    # calculate results ( enrichment stats for each term in the annotations map )
    results = []
    for term, members in annotations.items():
        # overlap between members with term and members of sample
        overlap = list( set( sample ).__and__( set( members.keys() ) ) )
        # counts
        count_overlap         = len( overlap )
        count_background      = len( background )
        count_sample          = len( sample )
        count_term            = len( members )
        count_sample_not_term = count_sample - count_overlap
        count_term_not_sample = count_term - count_overlap
        count_remainder       = count_background - count_overlap - count_term_not_sample - count_sample_not_term
        # frequencies
        freq_sample      = count_overlap / float( count_sample )
        freq_background  = count_term / float( count_background )
        # fold enrichment
        fold_enrichment  = freq_sample / freq_background
        # contingency table for fisher exact
        table = [
            [ count_overlap,         count_sample_not_term ],
            [ count_term_not_sample, count_remainder       ]
            ]
        # calculate pvalue and store results
        pvalue = fisher_exact( table )[1] # [0] is an odds ratio; do not want
        results.append( [term, len( overlap ), fold_enrichment, pvalue] )
    # convert pvalues to qvalues
    pvalues = [r[-1] for r in results]
    qvalues = pvalues2qvalues( pvalues )
    # attach qvalues
    for i, result in enumerate( results ):
        results[i].append( qvalues[i] )
    # filter results
    results2 = []
    for result in sorted( results, key=lambda x: x[-1] ):
        include = True
        term, overlap, fe, pvalue, qvalue = result
        if min_overlap is not None and overlap < min_overlap:
            include = False
        if min_fold is not None and 1 / float( min_fold ) < fe < min_fold:
            include = False
        if not depletions and fe < 1:
            include = False
        if fdr is not None and qvalue > fdr:
            include = False
        if include:
            results2.append( result )
    return results2

def rank_enrich ( values, annotations, min_overlap=1, depletions=True, fdr=None ):
    """
    Find terms where items in term have different value than items outside of the term
    000011112233334445566777888999 <= values
    ..........|.........|.|||..||| <= items in term
    =====
    values      is a mapping [item]=value
    annotations is a mapping [term]=[item]
    """
    items = set( values.keys() )
    results = []
    for term, members in annotations.items():
        overlap = items.__and__( set( members.keys() ) )
        x = [value for key, value in values.items() if key in overlap ]
        y = [value for key, value in values.items() if key not in overlap ]
        xmed = median( x )
        ymed = median( y )
        pvalue = mannwhitneyu( x, y )[1] # [0] is the test stat; do not want
        pvalue *= 2 # **** python's mannwhitneyu is 1-sided by default ****
        results.append( [term, len( overlap ), xmed, ymed, pvalue] )
    # convert pvalues to qvalues
    pvalues = [r[-1] for r in results]
    qvalues = pvalues2qvalues( pvalues )
    # attach qvalues
    for i, result in enumerate( results ):
        results[i].append( qvalues[i] )
    # filter results
    results2 = []
    for result in sorted( results, key=lambda x: x[-1] ):
        include = True
        term, overlap, xmed, ymed, pvalue, qvalue = result
        if overlap < min_overlap:
            include = False
        if fdr is not None and qvalue > fdr:
            include = False
        if not depletions and xmed < ymed:
            include = False
        if include:
            results2.append( result )
    return results2
