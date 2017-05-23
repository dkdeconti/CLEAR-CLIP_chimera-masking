#! /usr/bin/python

'''
Demuxes miR from chimeric CLEAR-CLIP reads.
'''

import argparse
import pysam


def reverse_map_bam(bam):
    '''
    Reverse maps the miR in BAM to FASTQ headers.
    ''' 
    samfile = pysam.AlignmentFile(bam, 'rb')
    for read in samfile.fetch():
        print read.reference_id
        pass
    samfile.close()
    return
    pass


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

if __name__ == "__main__":
    main()
