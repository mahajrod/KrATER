#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from KRATER.Routines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input 2-column (kmer and count) tab-separated file with kmer counts. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output 2-column (gc content and count) tab-separated file with. Default: stdout")

args = parser.parse_args()

with FileRoutines.metaopen(args.input, "r", buffering=100000000) as in_fd, \
     FileRoutines.metaopen(args.output, "w", buffering=100000000) as out_fd:
    for line in in_fd:
        line_list = line.split() # splits on any space-character and trims \n from thew end
        out_fd.write("{0}\t{1}\n".format(line_list[0].count("G") + line_list[0].count("C"), line_list[1]))

