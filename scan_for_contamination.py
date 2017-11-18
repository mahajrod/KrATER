#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse

from KRATER.Tools.Kmers import Jellyfish


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with representative contaminating sequences, for example 18S rRNA")
parser.add_argument("-j", "--jellyfish_db", action="store", dest="jf_db", required=True,
                    help="Jellyfish database")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-r", "--retain_intermediate_files", action="store_true", dest="retain_intermediate_files",
                    default=False, help="Store intermediate files")

args = parser.parse_args()

Jellyfish.scan_for_contamination(args.input, args.jf_db, args.output, splited_input_dir="splited_fasta",
                                 parsing_mode="parse", splited_output_dir="splited_output",
                                 splited_sorted_unique_output="splited_sorted_unique_output",
                                 retain_intermediate_file=args.retain_intermediate_files)
