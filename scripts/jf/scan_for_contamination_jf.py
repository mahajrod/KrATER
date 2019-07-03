#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse

from KRATER.Routines import JellyfishRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with representative contaminating sequences, for example 18S rRNA")
parser.add_argument("-j", "--jellyfish_db", action="store", dest="jf_db", required=True,
                    help="Jellyfish database")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-r", "--retain_intermediate_files", action="store_true", dest="retain_intermediate_files",
                    default=False, help="Store intermediate files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use. Default: 1")

args = parser.parse_args()

JellyfishRoutines.threads = args.threads
JellyfishRoutines.scan_for_contamination(args.input, args.jf_db, args.output, splited_input_dir="splited_fasta",
                                         parsing_mode="parse", splited_output_dir="splited_output",
                                         splited_sorted_unique_output="splited_sorted_unique_output",
                                         retain_intermediate_file=args.retain_intermediate_files)
