#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from KRATER.Tools.Kmers import Jellyfish


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with representative contaminating sequences, for example 18S rRNA")
parser.add_argument("-j", "--jellyfish_db", action="store", dest="jf_db", required=True,
                    help="Jellyfish database")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use. Default: 1")
parser.add_argument("-k", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="K-mer length in database. Default: 23")
parser.add_argument("-f", "--min_fraction_of_covered_kmers", action="store", dest="min_fraction_of_covered_kmer",
                    type=float, default=0.9,
                    help="Minimum fraction of covered k-mers to pass filtration. Default: 0.9")
parser.add_argument("-m", "--min_median_kmer_coverage", action="store", dest="min_median_kmer_coverage",
                    type=int, default=20,
                    help="Minimum median k-mer coverage to pass filtration. Default: 20")

args = parser.parse_args()

Jellyfish.threads = args.threads
Jellyfish.scan_for_contamination_jf_pb(args.input, args.jf_db, args.output_prefix,
                                       parsing_mode="parse", splited_full_output_dir="splited_full_output",
                                       threads=None, kmer_length=args.kmer_length,
                                       min_covered_kmers_fraction=args.min_fraction_of_covered_kmers,
                                       min_median_kmer_coverage=args.min_median_kmer_coverage)
