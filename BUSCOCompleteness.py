import argparse
import sys

from pathlib import Path

from src.busco import run_busco

#Function to create arguments and help
def parse_arguments():
    description = "Assess the completeness of a genome annotation using BUSCO"
    parser = argparse.ArgumentParser(description=description)

    help_input = "(Required) Input sequence"
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)
    
    help_lineage = "(Optional) BUSCO lineage to be used. Viridiplantae_odb10 by default"
    parser.add_argument("--lineage", "-l", type=str,
                        help=help_lineage, default="Viridiplantae_odb10")
    
    help_threads = "(Optional) Threads to use"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)
    
    help_mode = "(Required) Specify which BUSCO analysis mode to run"
    parser.add_argument("--mode", "-m", type=str,
                        help=help_mode, required=True)
    
    help_output = "(Required) Output path"
    parser.add_argument("--out", "-o", type=str,
                        help=help_output, required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()

#Function to compile de arguments into a dictionary
def get_arguments():
    parser = parse_arguments()
    return {"input_file": Path(parser.input_file),
            "lineage": parser.lineage,
            "threads": parser.threads,
            "mode": parser.mode,
            "output": Path(parser.out)}

#Main function
def main():
    arguments = get_arguments()
    busco_results = run_busco(arguments)
    print(busco_results)
    #write_busco_results(busco_results, arguments["out"])

if __name__ == "__main__":
    main()