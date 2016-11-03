#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse

from KRATER.Tools.Kmers import Jellyfish


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with data for histogram")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default - 10")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=5,
                    help="Low limit of histogram without logscale")
parser.add_argument("-g", "-high_limit", action="store", dest="high_limit", type=int, default=100,
                    help="High limit of histogram without logscale")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, required=True,
                    help="Length of kmer used to construct histo file.")
parser.add_argument("-d", "--draw_separated_pictures", action="store_true", dest="draw_separated_pictures",
                    default=False,
                    help="Draw additional separated pictures for double logarithmic and linear scales."
                         "Default - False")
#parser.add_argument("-d", "--draw_peaks_and_gaps", action="store_true", dest="draw_peaks_and_gaps",
#                    help="Draw peaks and gaps")

args = parser.parse_args()

Jellyfish.draw_kmer_distribution(args.input, args.kmer_length, args.output_prefix, output_formats=args.output_formats,
                                 logbase=args.logbase, non_log_low_limit=args.low_limit,
                                 non_log_high_limit=args.high_limit,
                                 draw_separated_pictures=args.draw_separated_pictures) #, draw_peaks_and_gaps=args.draw_peaks_and_gaps)
