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
from StringIO import StringIO


def consolidate_mirs(mir_map):
    '''
    Chooses whether the miR is first or last.
    '''
    consol_mir_map = {}
    consol_mir_map = defaultdict(dict)
    for basename, mirs in mir_map.items():
        if mirs:
            for mir, stats in mirs.items():
                best = max({"first": stats["first"],
                            "last": stats["last"]}.iteritems(),
                           key=operator.itemgetter(1))
                if best[1] == 0:
                    continue
                consol_mir_map[basename][mir] = (best[0],
                                                 stats[best[0] + "_end"])
    return consol_mir_map


def consolidate_first_last(mir_map):
    '''
    Chooses which miR align to use for demux.
    '''
    consol_first_last = defaultdict(dict)
    for basename, mirs in mir_map.items():
        first_cnt = len([mir for mir, stats in mirs.items()
                         if stats[0] == "first"])
        last_cnt = len([mir for mir, stats in mirs.items()
                        if stats[0] == "last"])
        if first_cnt >= last_cnt:
            consol_first_last[basename]["first"] = max([s[1]
                                                        for s in mirs.values()
                                                        if s[0] == "first"])
        else:
            consol_first_last[basename]["last"] = max([s[1]
                                                       for s in mirs.values()
                                                       if s[0] == "last"])
    return consol_first_last


def demux_reads(fastq, mirlen_map):
    '''
    Removes chimeric miR sequence from read.
    '''
    for record in SeqIO.parse(open(fastq, 'rU'), "fastq"):
        basename = record.id[:-2]
        if basename not in mirlen_map:
            continue
        if record.id[-1] == "1":
            record_order = "first"
        else:
            record_order = "last"
        try:
            yield record[mirlen_map[basename][record_order]:]
        except KeyError:
            yield record


def determine_demux(inverted_map):
    '''
    Determines length of miR and where to excise.
    '''
    first_last_map = find_first_last(inverted_map)
    consol_mirs = consolidate_first_last(consolidate_mirs(first_last_map))
    return consol_mirs


def find_first_last(inverted_map):
    '''
    Finds miRs that are mapped to beginning of read.
    '''
    mir_map = {}
    for read, alns in inverted_map.items():
        basename = read[:-2]
        mir_map[basename] = {}
        for aln in alns:
            #print read, '\t', aln
            mir, start, end = aln
            if mir not in mir_map[basename]:
                mir_map[basename][mir] = {"first": 0, "last": 0,
                                          "first_end": 0, "last_end": 0}
            read_order = read[-1]
            if read_order == "1" and start <= 1:
                current_end = mir_map[basename][mir]["first_end"]
                mir_map[basename][mir]["first"] += 1
                mir_map[basename][mir]["first_end"] = max(end, current_end)
            elif read_order == "2" and start <= 1:
                current_end = mir_map[basename][mir]["last_end"]
                mir_map[basename][mir]["last"] += 1
                mir_map[basename][mir]["last_end"] = max(end, current_end)
    return mir_map


def invert_map_bam(samfile, reference_name):
    '''
    Reverse maps the miR in BAM to FASTQ headers.
    '''
    inverted_map = defaultdict(list)
    for read in samfile.fetch(reference=reference_name):
        if read.is_unmapped:
            continue
        inverted_map[read.reference_name].append((read.query_name,
                                                  read.reference_start,
                                                  read.reference_end))
    return inverted_map


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
    parser.add_argument("-n", "--name", metavar="OUTPUT_NAME",
                        default=None,
                        help="name of output fastq [default: stdout]")
    args = parser.parse_args()
    # Central dispatch
    if args.sam:
        input_format = 'r'
    else:
        input_format = 'rb'
    with pysam.AlignmentFile(args.BAM, input_format) as samfile:
        for reference_name in samfile.references:
            bam2fq_map = invert_map_bam(samfile, reference_name)
            mirlen_map = determine_demux(bam2fq_map)
            demuxed_reads = demux_reads(args.FASTQ, mirlen_map)
    if args.name:
        with open(args.name, 'w') as handle:
            SeqIO.write(demuxed_reads, handle, "fastq")
    else:
        out_handle = StringIO()
        SeqIO.write(demuxed_reads, out_handle, "fastq")
        data = out_handle.getvalue()
        sys.stdout.write(data)

if __name__ == "__main__":
    main()
