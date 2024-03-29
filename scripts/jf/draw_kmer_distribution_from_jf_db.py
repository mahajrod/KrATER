#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse

from KRATER.Tools.Kmers import Jellyfish
from KRATER.Routines import JellyfishRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input jellyfish database")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure. Default: png")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default - 10")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")
parser.add_argument("-j", "--jellyfish_path", action="store", dest="jellyfish_path",
                    help="Path to jellyfish")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=5,
                    help="Low limit of histogram without logscale. Default - 5")
parser.add_argument("-g", "--high_limit", action="store", dest="high_limit", type=int, default=100,
                    help="High limit of histogram without logscale. Default - 100")
parser.add_argument("-d", "--draw_separated_pictures", action="store_true", dest="draw_separated_pictures",
                    default=False,
                    help="Draw additional separated pictures for double logarithmic and linear scales."
                         "Default - False")

parser.add_argument("-n", "--naive", action="store_true", dest="naive", default=False,
                    help="Use naive estimation of genome size. By default, genomescope2 is used.")

parser.add_argument("-u", "--point_number", action="store", dest="point_number", default=3, type=int,
                    help="Number of values on both sides from each bin in histogram to count for detection "
                         "of maximums and minimums. Rise if your histogram is noisy."
                         "Affects only naive genomesize estimation. Default: 3")
parser.add_argument("-z", "--use_second_peak", action="store_true", dest="use_second_peak", default=False,
                    help="Use second peak from histogram for naive genome size estimation."
                         "Affects only naive genomesize estimation.")

parser.add_argument("--ploidy", action="store", dest="ploidy", default=2, type=int,
                    help="Ploidy of the sample. Affects only genomesize estimation by GenomeScope2."
                         " Default: 2")
parser.add_argument("--initial_haploid_coverage", action="store", dest="initial_haploid_coverage", default=None,
                    type=float,
                    help="Initial estimation of the haploid coverage. "
                         "Affects only genomesize estimation by GenomeScope2."
                         " Default: not set")
parser.add_argument("--genomescope_cmd", action="store", dest="genomescope_cmd", default="genomescope.R",
                    help="Script to call GenomeScope2. Affects only genomesize estimation by GenomeScope2."
                         " Default: genomescope.R")

parser.add_argument("--dont_show_genome_size_on_plot", action="store_true", dest="dont_show_genome_size_on_plot", default=False,
                    help="Dont show genome size on plot. Default: False")
parser.add_argument("-x", "--dont_show_genome_size_ci_on_plot", action="store_true",
                    dest="dont_show_genome_size_ci_on_plot", default=False,
                    help="Dont show confidence interval for genome size on plot. Default: False")

args = parser.parse_args()

kmer_table_file = "%s_%i_mer.counts" % (args.output_prefix, args.kmer_length)
kmer_file = "%s_%i_mer.kmer" % (args.output_prefix, args.kmer_length)

histo_file = "%s_%i_mer.histo" % (args.output_prefix, args.kmer_length)
picture_prefix = "%s_%i_mer_histogram" % (args.output_prefix, args.kmer_length)

Jellyfish.histo(args.input, histo_file, upper_count=100000000)

JellyfishRoutines.draw_kmer_distribution(histo_file,
                                         args.kmer_length,
                                         picture_prefix,
                                         output_formats=args.output_formats,
                                         logbase=args.logbase,
                                         non_log_low_limit=args.low_limit,
                                         non_log_high_limit=args.high_limit,
                                         order=args.point_number,
                                         use_second_peak_for_genome_size_estimation=args.use_second_peak,
                                         draw_separated_pictures=args.draw_separated_pictures,
                                         dont_show_genome_size_on_plot=args.dont_show_genome_size_on_plot,
                                         genomescope2=not args.naive, ploidy=args.ploidy,
                                         initial_haploid_coverage=args.initial_haploid_coverage,
                                         genomescope_cmd=args.genomescope_cmd,
                                         show_confidence_interval=not args.dont_show_genome_size_ci_on_plot)
