#!/usr/bin/env python2

"""
    Candidate filtering from .profile files.
"""

from __future__ import division # enable true integer division (<int> / <int> = <float>)

__author__ =     'Thomas Kemmer, Ralf Hauenschild'
__maintainer__ = 'Thomas Kemmer'
__email__ =      'thkemmer@uni-mainz.de'
__version__ =    '2.0'

import sys

class Candidate(object):
    """Potential candidate"""
    refseg  = None
    pos     = None
    arate   = None
    relmism = None
    amism   = 0.
    gmism   = 0.
    tmism   = 0.
    cmism   = 0.
    cov     = None
    base    = None
    rh5     = None
    rh6     = None # to be set later on!
    div     = None # to be set via compute_diversity()
    
    def __init__(self, line):
        """Parses the given line of a .profile file and sets the corresponding attributes."""
        (self.refseg, self.pos, self.base, self.cov, _, __, matches, Amism, Gmism, Tmism, Cmism, ___, amism,
         gmism, tmism, cmism, ____, self.arate) = line.strip().split('\t')
        (self.pos, self.cov, matches, Amism, Gmism, Tmism, Cmism, amism, gmism, tmism, cmism
         )= map(int, [self.pos, self.cov, matches, Amism, Gmism, Tmism, Cmism, amism, gmism, tmism, cmism])
        self.arate = float(self.arate)
        self.relmism = 1 - matches / self.cov if self.cov != 0 else 0
        mism = Amism + Gmism + Tmism + Cmism + amism + gmism + tmism + cmism
        if mism != 0:
            self.amism = (Amism + amism) / mism
            self.gmism = (Gmism + gmism) / mism
            self.tmism = (Tmism + tmism) / mism
            self.cmism = (Cmism + cmism) / mism
        self.rh5 = self.relmism / self.arate if self.arate != 0 else 0

    def compute_diversity(self):
        """Computes the diversity score for this candidate. Requires rh5 and rh6 to be set!"""
        self.div = (1     * int(.1 <= self.rh5 < 10) +
                    10    * int(self.rh6 >= 2) +
                    100   * int(.2 <= self.arate < .9) +
                    1000  * int(.2 <= self.relmism < .9) +
                    10000 * int(self.relmism >= .1 and max(self.amism, self.gmism, self.cmism, self.tmism) <= .9))

if len(sys.argv) < 13:
    print 'Usage: profile2candidates.py <input file> <output file> <reference base> <visual range> <min coverage>',
    print '<min coverage (3p)> <min relative mismatches> <min arrest rate> <min diversity> <min rh6> <min position>',
    print '<max position> [<exlude positions (comma-separated list)>]'
    sys.exit(2)

refbase    = sys.argv[3].upper()
visrange   = int(sys.argv[4])
mincov     = int(sys.argv[5])
mincov3p   = int(sys.argv[6])
minrelmism = float(sys.argv[7])
minarate   = float(sys.argv[8])
mindiv     = int(sys.argv[9])
minrh6     = float(sys.argv[10])
minpos     = int(sys.argv[11])
maxpos     = int(sys.argv[12])
excludepos = []
if len(sys.argv) > 13 and len(sys.argv[13].strip()) != 0:
    excludepos = map(int, sys.argv[13].strip().split(','))

mismcols = ['A', 'G', 'T', 'C']
mismcols.remove(refbase)
    
with open(sys.argv[1]) as fin, open(sys.argv[2], 'w') as fout:
    fout.write('# Reference base: %s\n' % refbase)
    fout.write('# Visual range: %d\n' % visrange)
    fout.write('# Minimum coverage: %d\n' % mincov)
    fout.write('# Minimum coverage (next position): %d\n' % mincov3p)
    fout.write('# Minimum relative mismatch rate: %f\n' % minrelmism)
    fout.write('# Minimum arrest rate: %f\n' % minarate)
    fout.write('# Minimum diversity: %d\n' % mindiv)
    fout.write('# Minimum RH6 score: %f\n' % minrh6)
    fout.write('# Position search interval: [%d, %d]\n' % (minpos,maxpos))
    fout.write('# Excluded positions: [%s]\n' % ', '.join(map(str, excludepos)))
    fout.write('# \n')
    fout.write('# reference segment\tposition\tarrest rate\trelative mismatch rate\trh5 (relative mismatches per' +
               'arrest rate)\trh6 (fold change arrest rate)\t%s' % '\t'.join([b + ' mismatches' for b in mismcols]) + 
               '\tprevious base\tnext base\tnext but one base\tdiversity\tcoverage\n')
    cwindow = []
    currentrefseg = None
    for line in fin:
        if line.startswith('#') or len(line) == 0:
            continue   # skip headers and empty lines
        c = Candidate(line)
        if currentrefseg != c.refseg:
            cwindow = []
            currentrefseg = c.refseg
        cwindow.append(c)
        if len(cwindow) != 2 * visrange + 1:
            continue   # candidate window is not fully filled
        c = cwindow[visrange]
        meanarate = 0.5 * sum([e.arate for e in cwindow[:visrange] + cwindow[visrange+1:]]) / visrange
        c.rh6 = c.arate / meanarate if meanarate != 0 else 0.
        c.compute_diversity()
        # check filter criteria and write candidate
        #if (c.base == refbase and c.cov >= mincov and cwindow[visrange+1].cov >= mincov3p and c.relmism >= minrelmism
            #and c.arate >= minarate and c.div >= mindiv and c.rh6 >= minrh6 and minpos <= c.pos <= maxpos and 
            #c.pos not in excludepos):
	if (c.base == refbase and c.cov >= mincov and cwindow[visrange+1].cov >= mincov3p and c.relmism >= minrelmism
	    and c.arate >= minarate and c.div >= mindiv and c.rh6 >= minrh6 and minpos <= c.pos <= maxpos and c.pos not in excludepos):
            fout.write('%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t' % (c.refseg, c.pos, c.arate, c.relmism, c.rh5, c.rh6))
            fout.write('%s\t' % '\t'.join(['%.3f' % getattr(c, '%smism' % b.lower()) for b in mismcols]))
            fout.write('%s\t%s\t%s\t' % (cwindow[visrange-1].base, cwindow[visrange+1].base, cwindow[visrange+2].base))
            fout.write('%05d\t%d\n' % (c.div, c.cov))
        del cwindow[0]
