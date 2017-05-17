#! /usr/bin/python

import argparse
import re
import sys


def parse_ref_fasta(fasta, species):
    '''
    Parses fasta and filters name by species; prints to stdout.
    '''
    with open(fasta, 'rU') as handle:
        while True:
            fasta_header = handle.readline()
            fasta_seq = handle.readline()
            if not fasta_seq:
                break
            if re.search(species, fasta_header):
                sys.stdout.write(fasta_header)
                sys.stdout.write(fasta_seq)


def main():
    '''
    Arg parsing and central dispatch.
    '''
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("MIRBASE", metavar="miRbase_fasta",
                        help="miRbase fasta reference")
    parser.add_argument("-s", "--species", metavar="SPECIES",
                        help="Species to filter for")
    parser.set_defaults(SPECIES="Homo sapiens")
    args = parser.parse_args()
    # Dispatch
    parse_ref_fasta(args.MIRBASE, args.SPECIES)


if __name__ == "__main__":
    main()
