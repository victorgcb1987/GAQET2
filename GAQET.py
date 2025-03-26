import argparse
import sys

from csv import DictReader
from pathlib import Path

from src.agat import run_agat
from src.busco import run_busco, run_gffread
from src.LTR_retriever import create_outdir, run_suffixerator, run_harvest, run_finder, concatenate_outputs, run_LTR_retriever, run_LAI
from src.stringtie import run_stringtie, run_gffcompare, calculate_annotation_scores


#Function to create arguments and help
def parse_arguments():
    description = 'Tool with four modules to summarise the metrics and evaluate the quality of a genome annotation'
    parser = argparse.ArgumentParser(description=description)

    help_input = "(Required) Input File of Files"
    parser.add_argument("--input", "-i", type=str,
                        help=help_input, required=True)
    
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
    fof_fpath = Path(parser.input)
    if not fof_fpath.exists():
        raise RuntimeError("Fof does not exist")
    else:
        samples = {}
        with open(fof_fpath) as fof_fhand:
            samples = {line["name"]: line for line in DictReader(fof_fhand, delimiter="\t")}
                       
    return {"input": samples,
            "threads": parser.threads,
            "output": Path(parser.output)}


def main():
    arguments = get_arguments()
    out_dir =  arguments["output"]
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
    stats = {}                                          #diccionario para guardar los resultados de cada modulo por cada especie
    for name, values in arguments["input"].items():
        stats[name] = {}
        name_dir = out_dir / name
        values["output"] = name_dir
        values["threads"] = arguments["threads"]
        if not name_dir.exists():
            name_dir.mkdir(parents=True, exist_ok=True)

        agat_statistics = run_agat(values)
        print(agat_statistics)
        stats[name]["agat_statistics"] = agat_statistics

        gffread_results = run_gffread(values)
        print(gffread_results)
        busco_results = run_busco(values)
        print(busco_results)
        stats[name]["busco_results"] = get_busco_results(busco_results)

        LAI_out_dir =  create_outdir(values)
        print(LAI_out_dir)
        values["LAI_dir"] = LAI_out_dir["out_fpath"]
        suffixerator =  run_suffixerator(values)
        if "returncode" in suffixerator:
            if suffixerator["returncode"] == 1:
                raise RuntimeError("Suffixerator has failed")
        print(suffixerator)
        harvest = run_harvest(values)
        print(harvest)
        finder = run_finder(values)
        print(finder)
        cat = concatenate_outputs(values)
        print(cat)
        LTR = run_LTR_retriever(values)
        print(LTR)
        LAI = run_LAI(values)
        print(LAI)
        stats[name]["LAI"] = LAI

        stringtie = run_stringtie(values)
        print(stringtie)
        gffcompare = run_gffcompare(values)
        print(gffcompare)
        annotation_scores = calculate_annotation_scores(values)
        print(annotation_scores)
        stats[name]["annotation_scores"] = annotation_scores

        



if __name__ == "__main__":
    main()