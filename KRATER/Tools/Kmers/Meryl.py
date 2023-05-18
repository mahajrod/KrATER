#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from KRATER.Tools.Abstract import Tool


class Meryl(Tool): # couple of jellyfish specific arguments were retained in functions for compartibility
    """
    Several subcommands are not implemented: qhisto, qdump, qmerge, cite
    Not all options were implemented for count subcommand

    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "meryl", path=path, max_threads=max_threads)

    def count(self, in_file, out_file, kmer_length=23, memory="100g", hash_size=1000000, count_both_strands=False,
              lower_count=None, upper_count=None, generators=None):
        # IMPORTANT! Not all options were implemented

        input_files = [in_file] if isinstance(in_file, str) else in_file

        options =  " k=%i " % kmer_length
        options += " memory=%s " % memory
        options += " threads=%i " % self.threads
        options += " count "
        options += " output %s " % out_file
        options += " %s " % " ".join(input_files)

        self.execute(options)

    def histo(self, in_file, out_file, bin_width=1, lower_count=1, upper_count=100000000,
              include_absent_kmers=False):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options =  " memory=%s " % memory
        options += " threads=%i " % self.threads
        options += " histogram "
        options += " %s " % in_file
        options += " > %s" % out_file

        self.execute(options)

    def dump(self, in_file, out_file, lower_count=None, upper_count=None):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options =  " memory=%s " % memory
        options += " threads=%i " % self.threads
        options += " print "
        options += " less-than %i " % upper_count
        options += " greater-than %i " % lower_count
        options += " %s " % in_file
        if out_file != "stdout":
            options += " > %s " % out_file

        self.execute(options)
