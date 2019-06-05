#!/usr/bin/env python2

"""
    SAMtools pileup to .profile converter.

    Note: 
    The corrected coverage, counting only reads, touching (w.r.t. distance tolerance) the
    3' end is NOT included in the output of this procedure.
"""

from __future__ import division # enable true integer division (<int> / <int> = <float>)

__author__ =     'Thomas Kemmer, Ralf Hauenschild'
__maintainer__ = 'Thomas Kemmer'
__email__ =      'thkemmer@uni-mainz.de'
__version__ =    '2.0'

import sys

if len(sys.argv) != 3:
	print 'Usage: pileup2profile.py <input file> <output file>'
	sys.exit(2)
infile = sys.argv[1]
outfile = sys.argv[2]

import re
r1 = re.compile(r"[+-][0-9]+")
r2 = re.compile(r"(?:(?<=\^).|[^^$atgcn,.])", re.IGNORECASE)

# remove indels and quality scores
def strip_pileup(s):
    for start, group in sorted([(e.start(), e.group()) for e in r1.finditer(s)], reverse=True):
        s = s[:start] + s[start + len(group) + abs(int(group)):]
    return r2.sub('', s)

# compute base distribution
def cnt_bases(s):
    cnt = dict()
    for c in s:
        cnt[c] = cnt.get(c, 0) + 1
    return cnt

with open(infile) as fin, open(outfile, 'w') as fout:
    fout.write('# reference segment\tposition\treference base\tcoverage\trelative coverage ' +
               'change\tnorm of relative coverage change\tmatches\tA mismatches\tG mismatches' +
               '\tT mismatches\tC mismatches\tN mismatches\ta mismatches\tg mismatches\t' +
               't mismatches\tc mismatches\tn mismatches\tarrest rate\n')
    oldarate = None
    oldcov   = 0
    for line in fin:
        if len(line) == 0: continue
        refseg, pos, refbase, cov, pile, _ = line.split('\t')
        cov = int(cov)
        dcov = cov/oldcov-1 if oldcov != 0 else 0
        cnt = cnt_bases(strip_pileup(pile))
        ends  = cov - cnt.get("$", 0)
        arate = cnt.get('^', 0)/ends if ends != 0 else 0
        
        # append arrest rate to previous line
        if oldarate is not None:
            fout.write('%.3f\n' % arate)
        
        fout.write('%s\t%s\t%s\t%d\t%f\t' % (refseg, pos, refbase, cov, dcov))
        fout.write('%f\t' % (dcov/oldcov if oldcov != 0 else 0))
        fout.write('%d\t' % (cnt.get('.', 0) + cnt.get(',', 0)))
        [fout.write('%d\t' % cnt.get(base, 0)) for base in "AGTCNagtcn"]
        
        oldcov   = cov
        oldarate = arate
    # write last arrest rate
    if oldarate is not None:
        fout.write('%.3f\n' % oldarate)
