#!/usr/bin/env python
__author__ = 'mahajrod'

import re
import os
import sys
import bz2
import gzip
import shutil
from collections import Iterable, OrderedDict

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

class FileRoutines:
    def __init__(self):
        self.filetypes_dict = {"fasta": [".fa", ".fasta", ".fa", ".pep", ".cds"],
                               "fastq": [".fastq", ".fq"],
                               "genbank": [".gb", ".genbank"],
                               "newick": [".nwk"],
                               "gz": [".gz"],
                               "bzip": [".bz2"],
                               "bam": [".bam"],
                               "sam": [".sam"],
                               "cram": [".cram"]}

    @staticmethod
    def metaopen(filename, flags, buffering=None, compresslevel=5):
        if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
            if isinstance(filename, file):
                return filename
            else:
                raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
        elif filename[-3:] == ".gz":
            return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
        elif filename[-4:] == ".bz2":
            return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
        else:
            if buffering is not None:
                return open(filename, flags, buffering=buffering)
            else:
                return open(filename, flags)

    @staticmethod
    def add_external_extraction_to_filename(filename):
        if filename[-3:] == ".gz":
            return "<(gunzip -c %s)" % filename
        elif filename[-4:] == ".bz2":
            return "<(bunzip2 -c %s)" % filename
        else:
            return filename

    @staticmethod
    def add_external_extraction_to_filelist(filelist):
        new_filelist = []

        for filename in filelist:
            if filename[-3:] == ".gz":
                new_filelist.append("<(gunzip -c %s)" % filename)
            elif filename[-4:] == ".bz2":
                new_filelist.append("<(bunzip2 -c %s)" % filename)
            else:
                new_filelist.append(filename)

        return new_filelist

    @staticmethod
    def safe_mkdir(dirname,
                   description_filename=None, description_text=None,
                   readme_filename=None, readme_text=None):
        try:
            os.mkdir(dirname)
        except OSError:
            pass

        if not(description_filename is None):
            description_filename = "%s/%s" % (dirname, description_filename)

            if not os.path.isfile(description_filename):
                with open(description_filename, "w") as descr_fd:
                    if not (description_text is None):
                        descr_fd.write(description_text)

        if not(readme_filename is None):
            readme_filename = "%s/%s" % (dirname, readme_filename)

            if not os.path.isfile(readme_filename):
                with open(readme_filename, "w") as descr_fd:
                    if not (readme_text is None):
                        descr_fd.write(readme_text)

    def recursive_mkdir(self, dir_dict=None, out_dir=None,
                        description_filename=None, description_text=None,
                        readme_filename=None, readme_text=None):
        if not(out_dir is None):
            self.safe_mkdir(out_dir)

        if dir_dict:
            for directory in dir_dict:
                dirname = directory if out_dir is None else "%s/%s" % (out_dir, directory)
                self.safe_mkdir(dirname,
                                description_filename=description_filename,
                                description_text=description_text,
                                readme_filename=readme_filename,
                                readme_text=readme_text)

                if isinstance(dir_dict[directory], dict):
                    self.recursive_mkdir(dir_dict[directory],
                                         out_dir=dirname,
                                         description_filename=description_filename,
                                         description_text=description_text,
                                         readme_filename=readme_filename,
                                         readme_text=readme_text)

    def detect_filetype_by_extension(self, filename, filetypes_dict=None):
        filetypes = filetypes_dict if filetypes_dict else self.filetypes_dict
        directory, prefix, extension = self.split_filename(filename)
        for filetype in filetypes:
            if extension in filetypes[filetype]:
                return filetype
        return None