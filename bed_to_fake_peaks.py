#!/usr/bin/python

import argparse
import sys


def main():
    '''

    '''
    # Arg parsing
    desc = "Converts plain BED to peaks file with placeholder values"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("BED", metavar="BED", help="BED file")
    args = parser.parse_args()
    # Parse BED and add placeholder values.
    id_no = 1
    with open(args.BED, 'rU') as handle:
        pvalue_default = 0.01
        for line in handle:
            arow = line.strip('\n').split('\t')
            chrom = arow[0]
            start = arow[1]
            end = arow[2]
            #strand = arow[3]
            strand = "0"
            reads = arow[3]
            peak_id = '_'.join([args.BED[:-4], str(id_no)])
            #pvalue = str(pvalue_default / float(reads))
            out = [chrom, start, end, peak_id, reads, strand]
            #out = [chrom, start, end, pvalue, reads, strand]
            sys.stdout.write('\t'.join(out) + '\n')
            id_no += 1


if __name__ == "__main__":
    main()
