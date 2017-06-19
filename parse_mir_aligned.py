#!/usr/bin/python

import argparse
from itertools import tee
from itertools import islice
from itertools import izip
from StringIO import StringIO
import sys
import pysam
from Bio import SeqIO


def has_mir(read_name, mir_name, samfile):
    '''
    Checks samfile for whether read has miR aligned to it.
    '''
    try:
        for mir in samfile.fetch(reference=read_name):
            if mir.query_name == mir_name:
                return True
    except ValueError:
        return False

def pairwise(iterable):
    '''
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    '''
    a, b = tee(iterable)
    return izip(islice(a, 0, None, 2), islice(b, 1, None, 2))


def main():
    '''
    Arg parsing and central dispatch.
    '''
    # Arg parsing
    parser = argparse.ArgumentParser(description="Gets reads aligned with miR")
    parser.add_argument("mir", metavar="miR", help="miR name")
    parser.add_argument("bam", metavar="BAM",
                        help="BAM of reverse mapped miRs")
    parser.add_argument("fastq", metavar="FASTQ",
                        help="FASTQ of the chimeric reads")
    args = parser.parse_args()
    # Central Dispatch
    fastq_iter = pairwise(SeqIO.parse(open(args.fastq, 'rU'), "fastq"))
    with pysam.AlignmentFile(args.bam, 'rb') as samfile:
        for first, second in fastq_iter:
            if has_mir(first.id, args.mir, samfile) or \
            has_mir(second.id, args.mir, samfile):
                out_handle = StringIO()
                SeqIO.write((first, second), out_handle, "fastq")
                data = out_handle.getvalue()
                sys.stdout.write(data)


if __name__ == "__main__":
    main()
