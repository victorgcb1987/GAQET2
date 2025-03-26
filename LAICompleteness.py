import argparse
import sys

from pathlib import Path


import src.LTR_retriever as LTR


#parses command-line arguments
def parse_arguments():
    description = "Module to assess the completeness of the genome assembly in term of LTR elements (LTR Assembly Index, LAI)"
    parser = argparse.ArgumentParser(description=description)

    #inputs
    help_fasta = "(Required) Specify the genome sequence file (FASTA)"
    parser.add_argument("--fasta", "-f", type=str, 
                        help=help_fasta, required=True)
    
    help_gff = "(Required) "
    parser.add_argument("--gff", "-g", type=str,
                        help=help_gff)

    #threads
    help_threads = "Number of threads to run the program"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)

    #output
    help_output = "(Required) Output path"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output) #required=True

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


#stores command-line arguments in a dictionary with proper data types
def get_arguments():
    parser = parse_arguments()
    return {"fasta": Path(parser.fasta),
            "gff": Path(parser.gff),
            "threads": parser.threads,
            "output": Path(parser.output)} 
            ###


#main function
def main():
    arguments = get_arguments()
    #function to create the output directory if this does not exist
    outdir =  LTR.create_outdir(arguments)
    print(outdir)

    suffixerator =  LTR.run_suffixerator(arguments)
    if "returncode" in suffixerator:
        if suffixerator["returncode"] == 1:
            raise RuntimeError("Suffixerator has failed")
    print(suffixerator)
    harvest = LTR.run_harvest(arguments)
    print(harvest)
    
    finder = LTR.run_finder(arguments)
    print(finder)

    cat = LTR.concatenate_outputs(arguments)
    print(concatenate_outputs)

    LTR = LTR.run_LTR_retriever(arguments)
    print(LTR)

    LAI = LTR.LTRrun_LAI(arguments)
    print(LAI)


if __name__ == "__main__":
    main()