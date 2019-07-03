#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import numpy as np
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()

from aindex import *

from RouToolPa.Tools.Abstract import Tool
from RouToolPa.Parsers.Sequence import CollectionSequence


class AIndexRoutines(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "aindex", path=path, max_threads=max_threads)

    def scan_for_contamination(self,
                               sequence_file,
                               index_prefix,
                               kmer_length,
                               output_file,
                               aindex_prefix=None,
                               reads_file=None,
                               parsing_mode="parse",
                               external_process_pool=None,
                               threads=None):

        index_settings = {
                          "index_prefix": index_prefix,
                          "aindex_prefix": aindex_prefix,
                          "reads_file": reads_file,
                          }

        print("Parsing fasta file...")

        sequence_collection = CollectionSequence(in_file=sequence_file, parsing_mode=parsing_mode)
        sequence_collection.get_stats_and_features(count_gaps=False, sort=False, min_gap_length=1)

        index = load_aindex(index_settings,
                            skip_reads=False if reads_file else True,
                            skip_aindex=False if reads_file else True)

        def get_kmer_coverage(seq_id):
            return (seq_id, [index[sequence_collection.records[seq_id][i:i + kmer_length]] \
                             for i in xrange(sequence_collection.lengths.at[seq_id, "length"] - kmer_length + 1)])

        process_pool = external_process_pool if external_process_pool else mp.Pool(threads if threads else self.threads)

        results = process_pool.map(get_kmer_coverage, sequence_collection.scaffolds, chunksize=1)

        print("Analyzing results...")
        with open(output_file, "w") as out_fd:
            out_fd.write("#record_id\tkmer_number\tcovered_positions\tcovered_positions,%\t"
                         "kmer_mean_coverage\tkmer_median_coverage\tdescription\n")
            for seq_id, kmer_coverage_list in results:
                kmer_number = sequence_collection.lengths.at[seq_id, "length"]-kmer_length+1
                covered_positions = kmer_number - kmer_coverage_list.count(0)

                mean_kmer_coverage = np.mean(kmer_coverage_list)
                median_kmer_coverage = np.median(kmer_coverage_list)

                out_fd.write("%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%s\n" % (seq_id,
                                                                     kmer_number,
                                                                     covered_positions,
                                                                     100*float(covered_positions)/float(kmer_number),
                                                                     mean_kmer_coverage,
                                                                     median_kmer_coverage,
                                                                     sequence_collection.description[seq_id]))
