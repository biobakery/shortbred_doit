#! /usr/bin/env python

"""
"""

def steps ( n, alpha=0.05 ):
    return [alpha * i / float( n ) for i in range( 1, n+1 )]

def fdr_cutoff ( pvalues, alpha = 0.05 ):
    n = len( pvalues )
    pvalues = sorted( pvalues )
    max_step = alpha / float( n )
    for pvalue, step in zip( pvalues, steps( n, alpha ) ):
        if pvalue <= step:
            max_step = pvalue
    return max_step

def pvalues2qvalues ( pvalues, adjusted=False ):
    n = len( pvalues )
    # after sorting, index[i] is the original index of the ith-ranked value
    index = range( n )
    index = sorted( index, key=lambda i: pvalues[i] )
    pvalues = sorted( pvalues )
    qvalues = [pvalues[i-1] * n / i for i in range( 1, n+1 )]
    # adjust qvalues to enforce monotonic behavior?
    # q( i ) = min( q( i..n ) )
    if adjusted:
        qvalues.reverse()
        for i in range( 1, n ):
            if qvalues[i] > qvalues[i-1]:
                qvalues[i] = qvalues[i-1]
        qvalues.reverse()
    # rebuild qvalues in the original order
    ordered_qvalues = [None for q in qvalues]
    for i, q in enumerate( qvalues ):
        ordered_qvalues[index[i]] = q
    return ordered_qvalues

def pdict2qdict ( pdict, adjusted=False ):
    keys, pvalues = [], []
    for key, p in pdict.items():
        keys.append( key )
        pvalues.append( p )
    qvalues = pvalues2qvalues( pvalues, adjusted=adjusted )
    return {key:q for key, q in zip( keys, qvalues )}
