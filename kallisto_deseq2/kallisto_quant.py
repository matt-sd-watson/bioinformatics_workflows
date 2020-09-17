# This script generate the shell command script for bulk quantification of rna seq libraries
# using kallisto quant with the single read option

from os import listdir
from os.path import isfile, join
import re
import os


# full file path names or fastqs were generated using
# find $(pwd) -maxdepth 1 -type f -not -path '*/\.*' > files_fullpath.txt

# partial file names (without the complete directory path) for the fastqs were generated using
# ls $search_path > filename.txt


def absoluteFilePaths(directory):
    absolute_paths = []
    for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           if ".fastq.gz" in f:
            absolute_paths.append(os.path.abspath(os.path.join(dirpath, f)))
    return absolute_paths


def kallisto_script(directory, output_file):

    partials = [f for f in listdir(directory) if isfile(join(directory, f))]
    partial_fastqs = [f for f in partials if ".fastq.gz" in f]

    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

    partial_name_read = sorted(partial_fastqs, key=alphanum_key)

    partial_split = [i.split('.fastq.gz', 1)[0] for i in partial_name_read]
    partial_no_n = [i.strip('\n') for i in partial_split]
    partial_clean = filter(None, partial_no_n)

    full_no_n = [i.strip('\n') for i in absoluteFilePaths(directory)]
    full_clean = filter(None, full_no_n)
    full_name_sorted = sorted(full_clean, key=alphanum_key)

    kallisto_commands = []

    for i, k in zip(full_name_sorted, partial_clean):
        # shell command line for kallisto rna seq quantification of single read libraries

        if re.search('ERR', i):
            kallisto_commands.append("kallisto quant -i /Users/mattsdwatson/kallisto/mus_musculus/transcriptome.idx -o {}_kallisto --single -l 125 -s 20 {}".format(k, i))
        else:
            kallisto_commands.append(
                "kallisto quant -i /Users/mattsdwatson/kallisto/mus_musculus/transcriptome.idx -o {}_kallisto --single -l 50 -s 10 {}".format(
                    k, i))

    with open(output_file, "w") as handle:
        # include the bash shebang line that is preferred for portability
        handle.write("#!/usr/bin/env bash")
        handle.write("\n")
        for element in kallisto_commands:
            handle.write(element)
            handle.write("\n")


# sample execution of the file function


kallisto_script("/Users/mattsdwatson/mouse_de/fastq/", "kallisto_quant.sh")

absoluteFilePaths("/Users/mattsdwatson/mouse_de/fastq/")

