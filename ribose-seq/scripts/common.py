#! /usr/bin/env python

''' modmap.common: common methods for modmap pipeline
'''

import sys
from pybedtools import BedTool

# constants
STRANDS = ('pos','neg')

def load_coverage(bam_filename, strand, verbose):
    ''' load coverage from bam_filename.
        strand: `pos` or `neg`
        verbose: bool
        '''

    if verbose:
        print >>sys.stderr, ">> loading coverage for %s on %s strand" % \
                            (bam_filename, strand)

    # XXX calc 5' positions in bedGraph format
    kwargs = {'5':True, 'bg':True}

    if strand == 'pos':
        kwargs['strand'] = '+'
    elif strand == 'neg':
        kwargs['strand'] = '-'

    bedtool = BedTool(bam_filename).genome_coverage(**kwargs)

    return bedtool
