import argparse
import sys

from pathlib import Path

from src.agat import run_agat

#Function to create arguments and help
def parse_arguments():
    description = "Module to summarise de GFF metrics using AGAT toolkit"
    parser = argparse.ArgumentParser(description=description)

    help_input = "(Required) Input GTF/GFF file"
    parser.add_argument("--gff", "-i", type=str,
                        help=help_input, required=True)
    
    help_gsize = "(Optional) Inform about the genome size in NT in order to compute more statics"
    parser.add_argument("--gs", "-g", type=int,
                        help=help_gsize, default=0)

    help_distribution = "Use this option to plot distribtions in pdf"
    parser.add_argument("--distribution", "-d", help=help_distribution,
                        action="store_true")

    help_plot = "Use this option to plot distribtions in pdf"
    parser.add_argument("--plot", "-p", help=help_plot,
                        action="store_true")

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
            "gs": parser.gs,
            "distribution": parser.distribution,
            "plot": parser.plot,
            "output": Path(parser.output)}


def main():
    arguments = get_arguments()
    agat_results = run_agat(arguments)
    print(agat_results)


if __name__ == "__main__":
    main()