#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse

from KRATER.Routines import AIndexRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sequences", action="store", dest="seqs", required=True,
                    help="Input fasta file with representative contaminating sequences, for example 18S rRNA")
parser.add_argument("-i", "--index_prefix", action="store", dest="index_prefix", required=True,
                    help="Index prefix")
parser.add_argument("-k", "--kmer_length", action="store", dest="kmer_len", required=True,
                    type=int,
                    help="Kmer length used for index creation")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use. Default: 1")

args = parser.parse_args()

AIndexRoutines.threads = args.threads
AIndexRoutines.scan_for_contamination(args.seqs,
                                      args.index_prefix,
                                      args.kmer_len,
                                      args.output,
                                      aindex_prefix=None,
                                      reads_file=None,
                                      parsing_mode="parse",
                                      external_process_pool=None,
                                      threads=None)
