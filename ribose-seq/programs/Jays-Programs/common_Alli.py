#! /usr/bin/env python

import sys

from pybedtools import BedTool

def load_coverage(bam_filename, strand, verbose):

	if verbose: print >>sys.stderr, ">> loading coverage for %s on %s strand" % (bam_filename, strand)

	#"bg": Report genome coverage in BedGraph format
	#"5": Calculate genome coverage of 5' positions
	kwargs = {'5':True, 'bg':True}

	#Positive and negative strands	
	if strand == 'pos':
        	kwargs['strand'] = '+'
    	elif strand == 'neg':
		kwargs['strand'] = '-'

	#Calculate coverage of 5' positions in BedGraph format
	#5' positions are where the ribonucleotides are located
	bedtool = BedTool(bam_filename).genome_coverage(**kwargs)

	return bedtool
