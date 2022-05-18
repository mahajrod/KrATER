#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse

from KRATER.Routines import JellyfishRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of input files with data for histogram")
parser.add_argument("-a", "--labels", action="store", dest="labels", type=lambda s: s.split(","),
                    help="Comma-separated list of labels for input files")
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
parser.add_argument("-n", "--naive", action="store_true", dest="naive", default=False,
                    help="Use naive estimation of genome size. By default, genomescope2 is used.")

parser.add_argument("-u", "--point_number", action="store", dest="point_number", default=3, type=int,
                    help="Number of values on both sides from each bin in histogram to count for detection "
                         "of maximums and minimums. Rise if your histogram is noisy."
                         "Affects only naive genomesize estimation. Default: 3")
parser.add_argument("-z", "--use_second_peak", action="store_true", dest="use_second_peak", default=False,
                    help="Use second peak from histogram for naive genome size estimation."
                         "Affects only naive genomesize estimation.")

#parser.add_argument("-d", "--draw_peaks_and_gaps", action="store_true", dest="draw_peaks_and_gaps",
#                    help="Draw peaks and gaps")

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
parser.add_argument("--dont_show_genome_size_on_plot", action="store_true", dest="dont_show_genome_size_on_plot",
                    default=False, help="Dont show genome size on plot. Default: False")
parser.add_argument("-x", "--dont_show_genome_size_ci_on_plot", action="store_true",
                    dest="dont_show_genome_size_ci_on_plot", default=False,
                    help="Dont show confidence interval for genome size on plot. Default: False")
args = parser.parse_args()

JellyfishRoutines.draw_kmer_distribution(args.input,
                                         args.kmer_length,
                                         args.output_prefix,
                                         output_formats=args.output_formats,
                                         logbase=args.logbase,
                                         non_log_low_limit=args.low_limit,
                                         label_list=args.labels,
                                         non_log_high_limit=args.high_limit,
                                         order=args.point_number,
                                         draw_separated_pictures=args.draw_separated_pictures,
                                         use_second_peak_for_genome_size_estimation=args.use_second_peak,
                                         dont_show_genome_size_on_plot=args.dont_show_genome_size_on_plot,
                                         genomescope2=not args.naive, ploidy=args.ploidy,
                                         initial_haploid_coverage=args.initial_haploid_coverage,
                                         genomescope_cmd=args.genomescope_cmd,
                                         show_confidence_interval=not args.dont_show_genome_size_ci_on_plot) #, draw_peaks_and_gaps=args.draw_peaks_and_gaps)
