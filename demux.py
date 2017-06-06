#! /usr/bin/python

'''
Demuxes miR from chimeric CLEAR-CLIP reads.
'''

import argparse
import operator
import pysam
import sys
from Bio import SeqIO
from collections import defaultdict
from itertools import tee
from itertools import islice
from itertools import izip
from StringIO import StringIO


def determine_chimeric(first, second, samfile):
    '''
    Determines which read is chimeric and outputs which and pos.
    '''
    try:
        first_map = get_read_map(samfile.fetch(reference=first.id))
    except ValueError:
        first_map = None
    try:
        second_map = get_read_map(samfile.fetch(reference=second.id))
    except ValueError:
        second_map = None
    if not first_map and not second_map:
        return None
    first_best, second_best = None, None
    if first_map:
        first_best = filter_aln_to_best(first_map)
        if first_best[0] > 1: first_best = None
    if second_map:
        second_best = filter_aln_to_best(second_map)
        if second_best[0] > 1: second_best = None
    if first_best and second_best:
        return first_best <= second_best, min((first_best, second_best),
                                              key=operator.itemgetter(0))[1]
    elif first_best:
        return True, first_best[1]
    elif second_best:
        return False, second_best[1]
    else:
        return None


def demux_read(read, pos):
    '''
    Slices read from pos to end.
    '''
    return read[pos:]


def filter_aln_to_best(read_map):
    '''
    Filters read alignments to best likely chimera.
    '''
    best_per_mir = {k : min(v, key=operator.itemgetter(0))
                    for k, v in read_map.items()}
    best_mir = min(best_per_mir.items(), key=lambda x: x[1][0])[0]
    return best_per_mir[best_mir]


def get_min_aln(aln_map):
    '''
    Finds the miR with the min pos.
    '''
    if aln_map:
        best_k = min(aln_map, key=aln_map.get)
        best_v = aln_map[best_k]
        return best_k, best_v
    else:
        return None, None


def get_read_map(read_iter):
    '''
    Returns dict of read start and end pos.
    '''
    read_pos_map = defaultdict(list)
    for read in read_iter:
        if not read.is_unmapped:
            read_pos_map[read.query_name].append((read.reference_start,
                                                  read.reference_end))
    return read_pos_map


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
    desc = "Demuxes miR from chimeric CLEAR-CLIP reads."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("BAM", metavar="BAM",
                        help="BAM of reverse mapped miRs")
    parser.add_argument("FASTQ", metavar="FASTQ",
                        help="FASTQ of the chimeric reads")
    parser.add_argument("-s", "--sam", action="store_true",
                        help="Input is in SAM format")
    #parser.add_argument("-n", "--name", metavar="OUTPUT_NAME",
    #                    default=None,
    #                    help="name of output fastq [default: stdout]")
    args = parser.parse_args()
    if args.sam:
        input_format = 'r'
    else:
        input_format = 'rb'
    # Central dispatch
    fastq_iter = pairwise(SeqIO.parse(open(args.FASTQ, 'rU'), "fastq"))
    with pysam.AlignmentFile(args.BAM, input_format) as samfile:
        for first, second in fastq_iter:
            try:
                is_first, pos = determine_chimeric(first, second, samfile)
            except TypeError:
                continue
            if is_first:
                first_out = demux_read(first, pos)
                second_out = second
            else:
                first_out = first
                second_out = demux_read(second, pos)
            out_handle = StringIO()
            SeqIO.write((first_out, second_out), out_handle, "fastq")
            data = out_handle.getvalue()
            sys.stdout.write(data)


if __name__ == "__main__":
    main()
