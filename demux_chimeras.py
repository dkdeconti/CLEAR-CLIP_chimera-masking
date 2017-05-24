#! /usr/bin/python

'''
Demuxes miR from chimeric CLEAR-CLIP reads.
'''

import argparse
import pysam
from Bio import SeqIO
from collections import defaultdict


def demux_reads(fastq, inverted_map):
    '''
    Removes chimeric miR sequence from read.
    '''
    for record in seqIO.parse(open(fastq, 'rU'), "fastq"):
        if record.id not in inverted_map:
            continue
        demuxed_read = None
        yield demuxed_read


def reverse_map_bam(bam):
    '''
    Reverse maps the miR in BAM to FASTQ headers.
    '''
    inverted_map = defaultdict(list)
    with pysam.AlignmentFile(bam, 'r') as samfile:
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            inverted_map[read.reference_name].append((read.query_name,
                                                      read.pos))
            #inverted_map[read.reference_id].append((read.refereancename, read.pos))
    return inverted_map


def main():
    '''
    Arg parsing and central dispatch.
    '''
    desc = "Demuxes miR from chimeric CLEAR-CLIP reads."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("BAM", metavar="BAM",
                        help="BAM of reverse mapped miRs")
    parser.add_argument("FASTQ", metavar="FASTQ",
                        help="FASTQ of the chimeric reads")
    args = parser.parse_args()
    bam2fq_map = reverse_map_bam(args.BAM)
    for k in bam2fq_map:
        print k
        for v in bam2fq_map[k]:
            print '\t', v
    demux_reads(bam2fq_map, args.FASTQ)

if __name__ == "__main__":
    main()
