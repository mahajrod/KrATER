#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from copy import deepcopy
from Bio import SeqIO


from KRATER.Tools.Kmers import Jellyfish

from KRATER.Routines import JellyfishRoutines


def rev_com_generator(rec_dict, yield_original_record=False, rc_suffix="_rc", add_rc_suffix="False"):
    for record_id in rec_dict:
        reverse_seq = rec_dict[record_id].seq.reverse_complement()
        record = deepcopy(rec_dict[record_id])
        record.seq = reverse_seq
        record.id = record_id + rc_suffix if yield_original_record or add_rc_suffix else record.id
        if yield_original_record:
            yield rec_dict[record_id]
        yield record


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of fasta or fastq files or directories containing them.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default - 10")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=str, default="1000000",
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L). Default - 1000000")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default - 1")
parser.add_argument("-b", "--count_both_strands", action="store_true", dest="count_both_strands",
                    help="Count kmers in both strands. NOTICE: only mer or its reverse-complement, whichever "
                         "comes first lexicographically, is stored and the count value is the number of "
                         "occurrences of both. So this option is not suitable for generating sets of forward "
                         "and reverse-complement kmers. For this case use -r/--add_reverse_complement option. "
                         "Not compatible with -r/--add_reverse_complement option.")
parser.add_argument("-r", "--add_reverse_complement", action="store_true", dest="add_rev_com",
                    help="Add reverse-complement sequences before counting kmers. "
                         "Works only for fasta sequences. "
                         "Not compatible with -b/--count_both_strands option")
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
parser.add_argument("--turn_on_timelog", action="store_true", dest="turn_on_timelog", default=False,
                    help="Turn on timelog. Default: False")
parser.add_argument("--generators", action="store", dest="generators",
                    help="File with commands for fastq generators. Default: None")
parser.add_argument("-x", "--dont_show_genome_size_ci_on_plot", action="store_true",
                    dest="dont_show_genome_size_ci_on_plot", default=False,
                    help="Dont show confidence interval for genome size on plot. Default: False")


args = parser.parse_args()

args.input = Jellyfish.make_list_of_path_to_files(args.input)
if args.count_both_strands and args.add_rev_com:
    raise ValueError("Options -b/--count_both_strands and -r/--add_reverse_complement are not compatible")

if args.add_rev_com:
    file_with_rev_com = args.output_prefix + "_with_rev_com.fasta"
    record_dict = SeqIO.index_db("temp_index.idx", args.input, format="fasta")
    SeqIO.write(rev_com_generator(record_dict, yield_original_record=True), file_with_rev_com, "fasta")
    args.output_prefix += "_with_rev_com"
    os.remove("temp_index.idx")

output_prefix = "%s.k%i" % (args.output_prefix, args.kmer_length)
base_file = "%s.jf" % output_prefix
kmer_table_file = "%s.counts" % output_prefix
kmer_file = "%s.kmer" % output_prefix

histo_file = "%s.histo" % output_prefix
picture_prefix = "%s.histogram" % output_prefix

Jellyfish.threads = args.threads
Jellyfish.timelog = "%s.jellyfish.time.log" % output_prefix if args.turn_on_timelog else None
Jellyfish.path = args.jellyfish_path if args.jellyfish_path else ""
Jellyfish.count(args.input if not args.add_rev_com else file_with_rev_com, base_file,
                kmer_length=args.kmer_length, hash_size=args.hash_size,
                count_both_strands=args.count_both_strands,
                generators=args.generators)
Jellyfish.histo(base_file, histo_file, upper_count=100000000)
JellyfishRoutines.draw_kmer_distribution(histo_file, args.kmer_length, picture_prefix,
                                         output_formats=args.output_formats,
                                         logbase=args.logbase, non_log_low_limit=args.low_limit,
                                         non_log_high_limit=args.high_limit, order=args.point_number,
                                         use_second_peak_for_genome_size_estimation=args.use_second_peak,
                                         draw_separated_pictures=args.draw_separated_pictures,
                                         dont_show_genome_size_on_plot=args.dont_show_genome_size_on_plot,
                                         genomescope2=not args.naive, ploidy=args.ploidy,
                                         initial_haploid_coverage=args.initial_haploid_coverage,
                                         genomescope_cmd=args.genomescope_cmd,
                                         show_confidence_interval=not args.dont_show_genome_size_ci_on_plot)
