#!/usr/bin/python

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO


def filter_fasta(fasta, mirs):
    '''

    '''
    fasta_iter = SeqIO.parse(open(fasta, 'rU'), "fasta")
    out_handle = StringIO()
    for record in fasta_iter:
        if record.id in mirs:
            transcribed_record = SeqRecord(record.seq.back_transcribe())
            transcribed_record.id = record.id
            transcribed_record.name = record.name
            transcribed_record.description = record.description
            transcribed_record.features = record.features
            SeqIO.write(transcribed_record, out_handle, "fasta")
            data = out_handle.getvalue()
            sys.stdout.write(data)


def get_mir_ids(mir_gff):
    '''
    Pulls miR IDs from GFF.
    '''
    mirs = []
    with open(mir_gff, 'rU') as handle:
        for line in handle:
            if line[0] == "#":
                continue
            arow = line.strip('\n').split('\t')
            mir_type = arow[2]
            if mir_type != "miRNA":
                continue
            desc = arow[-1]
            mirs.append(desc.split(';')[2][5:])
    return tuple(mirs)


def main():
    '''
    Arg parsing and central dispatch.
    '''
    # Arg parsing
    desc = "Filters miR for only Mir in demuxed read alignments"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("mir_gff", metavar="miR_GFF",
                        help="GFF of miR")
    parser.add_argument("mir_fasta", metavar="miR_FASTA",
                        help="FASTA of miRs")
    args = parser.parse_args()
    # Central dispatch
    mirs = get_mir_ids(args.mir_gff)
    filter_fasta(args.mir_fasta, mirs)


if __name__ == "__main__":
    main()
