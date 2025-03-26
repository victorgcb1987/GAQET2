import argparse
import sys

from pathlib import Path

from src.stringtie import run_stringtie, run_gffcompare, calculate_annotation_scores


#Function to create arguments and help
def parse_arguments():
    description = "Module to ....................."
    parser = argparse.ArgumentParser(description=description)

    help_input = "(Required) Input GTF/GFF file"
    parser.add_argument("--gff", "-i", type=str,
                        help=help_input, required=True)

    help_bam = "(Required) Input BAM file"
    parser.add_argument("--bam", "-b", type=str,
                        help=help_bam, required=True)
    
    help_threads = "(Optional) Threads to use"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)

    help_output = "(Required) Output path"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output, required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()

#Function to compile de arguments into a dictionary
def get_arguments():
    parser = parse_arguments()
    return {"gff": Path(parser.gff),
            "bam": Path(parser.bam),
            "threads": parser.threads,
            "output": Path(parser.output)}


def main():
    arguments = get_arguments()
    stringtie = run_stringtie(arguments)
    print(stringtie)
    gffcompare = run_gffcompare(arguments)
    print(gffcompare)
    annotation_scores = calculate_annotation_scores(arguments)
    print(annotation_scores)
    ###


if __name__ == "__main__":
    main()