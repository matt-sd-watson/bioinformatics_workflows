# This script generate the shell command script for bulk quantification of rna seq libraries
# using kallisto quant with the single read option

from os import listdir
from os.path import isfile, join


# full file path names or fastqs were generated using
# find $(pwd) -maxdepth 1 -type f -not -path '*/\.*' > files_fullpath.txt

# partial file names (without the complete directory path) for the fastqs were generated using
# ls $search_path > filename.txt

def kallisto_directory(partial_names, output_file):

    partial_name_read = [f for f in open(partial_names, "r")]
    partial_split = [i.split('.fastq.gz', 1)[0] for i in partial_name_read]

    kallisto_commands = []

    for i in partial_split:
        kallisto_commands.append("mkdir {}_kallisto".format(i))

    with open(output_file, "w") as handle:
        # include the bash shebang line that is preferred for portability
        handle.write("#!/usr/bin/env bash")
        handle.write("\n")
        for element in kallisto_commands:
            handle.write(element)
            handle.write("\n")

    with open("samples.txt", "w") as handle:
        # include the bash shebang line that is preferred for portability
        for element in partial_split:
            handle.write(element)
            handle.write("\n")


# sample execution of the file function


kallisto_directory("/Users/mattsdwatson/mouse_de/fastq/filename.txt",
                "kallisto_directory.sh")

