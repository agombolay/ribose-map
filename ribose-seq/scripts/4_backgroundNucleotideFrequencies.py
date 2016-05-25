#! /usr/bin/env python

'''
modmap.genome_nuc_freqs
-----------------------
calculate background nucleotide frequencies report a table of raw and
normalized frequencies. 
$ python -m modmap.genome_nuc_freqs --help
'''

#import ipdb
import sys

from collections import Counter, defaultdict
from pyfaidx import Fasta

__author__ = 'Jay Hesselberth'
__contact__ = 'jay.hesselberth@gmail.com'
__version__ = '0.1'

def genome_nuc_freqs(fasta_filename, region_size_min, region_size_max,
                     ignore_chroms, only_chroms, verbose):

    header = ('#region.size','nuc','count','freq')
    print '\t'.join(header)

    nuc_counts = calc_bkgd_counts(fasta_filename, region_size_min,
                                 region_size_max, ignore_chroms,
                                 only_chroms, verbose)

    for region_size in nuc_counts:

        nuc_sum = float(sum(nuc_counts[region_size].values()))

        for nuc, count in sorted(nuc_counts[region_size].items()):
            nuc_freq = count / nuc_sum
            print '%d\t%s\t%s\t%s' % (region_size, nuc, str(count), str(nuc_freq))

def calc_bkgd_counts(fasta_filename, region_size_min,
                    region_size_max, ignore_chroms,
                    only_chroms, verbose):
    ''' calculate nuc frequencies for normalization.
        Returns: dict of nucleotide frequencies.
    '''

    nuc_counts = defaultdict(Counter)

    fasta = Fasta(fasta_filename, as_raw = True)

    for chrom in fasta.keys():

        # skip data based on specified chromosomes
        if chrom in ignore_chroms: continue

        if only_chroms and chrom not in only_chroms: continue

        seq_len = len(fasta[chrom])
        for idx in range(seq_len + 1):

            for region_size in range(region_size_min,
                                     region_size_max + 1):

                nucs = fasta[chrom][idx:idx+region_size]

                nuc_counts[region_size][nucs] += 1

    # remove entries that are not equal to region_size
    for region_size, nuc_dict in nuc_counts.items():
        for nuc, count in nuc_dict.items():
            if len(nuc) != region_size:
                nuc_dict.pop(nuc)

    return nuc_counts

def parse_options(args):
    from optparse import OptionParser, OptionGroup

    usage = "%prog [OPTION]... FASTA_FILENAME"
    version = "%%prog %s" % __version__
    description = ("pre-calculate nucleotide frequences for genome")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    group = OptionGroup(parser, "Variables")

    group.add_option("--region-size-minimum", action="store", type='int',
        default=1,
        help="minimum region size "
        " (default: %default)")

    group.add_option("--region-size-maximum", action="store", type='int',
        default=1,
        help="maximum region size "
        " (default: %default)")

    group.add_option("--ignore-chrom", action="append",
        metavar="CHROM", default=[],
        help="list of chroms to ignore"
        " (default: %default)")

    group.add_option("--only-chrom", action="append",
        metavar="CHROM", default=[],
        help="list of chroms to include"
        " (default: %default)")

    group.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    parser.add_option_group(group)

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("specify FASTA")

    if len(options.ignore_chrom) != 0 and len(options.only_chrom) != 0:
        parser.error("--ignore-chrom and --only-chrom are mutually exclusive ")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    ignore_chroms = set(options.ignore_chrom)
    only_chroms = set(options.only_chrom)

    kwargs = {'region_size_min':options.region_size_minimum,
              'region_size_max':options.region_size_maximum,
              'ignore_chroms':ignore_chroms,
              'only_chroms':only_chroms,
              'verbose':options.verbose}

    fasta_filename = args[0]

    return genome_nuc_freqs(fasta_filename, **kwargs)

if __name__ == '__main__':
    sys.exit(main()) 
