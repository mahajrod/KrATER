#!/usr/bin/env python
import os
import sys
import multiprocessing as mp
from subprocess import PIPE, Popen
import subprocess

print_mutex = mp.Lock()
SHELLPATH = "/bin/bash"


def execute(exe_string):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use mutex to safe write to stdout from multiple threads
    print_mutex.acquire()
    sys.stdout.write("Executing:\n\t%s\n" % exe_string)
    print_mutex.release()

    # os.system(exe_string)
    subprocess.call(exe_string, shell=True, executable=SHELLPATH)


class Tool:  # AlignmentRoutines

    def __init__(self, cmd, path="", max_threads=4, jar_path="", jar=None,
                 max_memory="500m", max_per_thread_memory="500m", timelog=None, tmp_dir=None):
        # SequenceRoutines.__init__(self)
        self.path = self.check_path(path)
        self.cmd = cmd
        self.threads = max_threads

        self.jar_path = self.check_path(jar_path) if jar_path else ""
        self.jar = jar
        self.max_memory = max_memory
        self.timelog = timelog
        self.max_per_thread_memory = max_per_thread_memory
        self.tmp_dir = tmp_dir

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
    def check_path(path_to_check):
        # print (path_to_check)
        # returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    @staticmethod
    def check_dir_path(path_to_check):
        # print (path_to_check)
        # returns path with / at end or blank path
        if path_to_check != "":
            if path_to_check[-1] != "/":
                return path_to_check + "/"
        return path_to_check

    def detect_filetype_by_extension(self, filename, filetypes_dict=None):
        filetypes = filetypes_dict if filetypes_dict else self.filetypes_dict
        directory, prefix, extension = self.split_filename(filename)
        for filetype in filetypes:
            if extension in filetypes[filetype]:
                return filetype
        return None

    @staticmethod
    def split_filename(filepath):
        directory, basename = os.path.split(filepath)
        if directory:
            if directory[-1] != "/":
                directory += "/"
        prefix, extension = os.path.splitext(basename)
        return directory, prefix, extension if filepath else None

    def make_list_of_path_to_files(self, list_of_dirs_and_files, expression=None, recursive=False,
                                   return_absolute_paths=True):
        file_list = []
        for entry in [list_of_dirs_and_files] if isinstance(list_of_dirs_and_files, str) else list_of_dirs_and_files:
            if os.path.isdir(entry):
                files_in_dir = ["%s%s" % (self.check_path(entry), filename)
                                for filename in sorted(filter(expression, os.listdir(entry))
                                                       if expression else os.listdir(entry))]
                if recursive:
                    for filename in files_in_dir:
                        if os.path.isdir(filename):
                            file_list += self.make_list_of_path_to_files([filename],
                                                                         expression=expression,
                                                                         recursive=recursive)
                        else:
                            file_list.append(filename)
                else:
                    file_list += files_in_dir
            elif os.path.exists(entry):
                if expression:
                    if expression(os.path.abspath(entry)):
                        file_list.append(os.path.abspath(entry))

                else:
                    file_list.append(os.path.abspath(entry))
            else:
                print("%s does not exist" % entry)
        # direct conversion to list was added for compatibility with python3
        # in which map function returns map object instead of list
        return list(map(os.path.abspath, file_list)) if return_absolute_paths else file_list

    def execute(self, options="", cmd=None, directory=None, capture_output=False,
                generate_cmd_string_only=False, intepreter=None):
        command = cmd if cmd is not None else self.cmd

        exe_string = ""
        exe_string += " cd %s && " % self.check_dir_path(directory) if directory else ""
        exe_string += ("%s " % intepreter if intepreter else "") + (self.check_path(self.path) if self.path else "") + command + " " + options

        sys.stdout.write("Executing:\n\t%s\n" % exe_string)
        if self.timelog:
            os.system("date >> %s" % self.timelog)
            with open(self.timelog, "a") as time_fd:
                time_fd.write("Command\t%s\n" % exe_string)

        if generate_cmd_string_only:
            return exe_string

        if capture_output:
            return Popen([exe_string], shell=True, stdout=PIPE).stdout  # returns file object
        else:
            subprocess.call(exe_string, shell=True, executable=SHELLPATH)
            if self.timelog:
                os.system("date >> %s" % self.timelog)
            return None
