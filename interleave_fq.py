#!/usr/bin/python

'''
Interleaves paired fastq reads.
'''

import argparse
import itertools
import sys
from Bio import SeqIO
from StringIO import StringIO


def interleave(iter1, iter2):
    '''
    Interleaves fastq
    '''
    for (forward, reverse) in itertools.izip(iter1, iter2):
        assert forward.id == reverse.id
        forward.id += "/1"
        reverse.id += "/2"
        yield forward
        yield reverse


def main():
    '''
    Arg parsing and central dispatch.
    '''
    # Arg parsing
    parser = argparse.ArgumentParser(description="Interleave paired fastq.")
    parser.add_argument("first", metavar="FIRST_FASTQ",
                        help="first in fastq pair")
    parser.add_argument("second", metavar="SECOND_FASTQ",
                        help="second in fastq pair")
    parser.add_argument("-n", "--name", metavar="OUTPUT_NAME",
                        default=None,
                        help="name of output fastq [default: stdout]")
    parser.add_argument("-f", "--format",
                        default="fastq",
                        help="format of files [default: fastq]")
    args = parser.parse_args()
    # function dispatching
    records_f = SeqIO.parse(open(args.first, 'rU'), args.format)
    records_r = SeqIO.parse(open(args.second, 'rU'), args.format)
    if args.name:
        with open(args.name, 'w') as handle:
            SeqIO.write(interleave(records_f, records_r), handle, args.format)
    else:
        for record in interleave(records_f, records_r):
            out_handle = StringIO()
            SeqIO.write(record, out_handle, args.format)
            data = out_handle.getvalue()
            sys.stdout.write(data)


if __name__ == "__main__":
    main()
