#!/usr/bin/env python

import argparse
import os
import sys

from pathlib import Path


from src.parsers import (parse_fof, get_pfams_from_db, get_pfams_from_interpro_query, 
                         parse_TEsort_output, classify_pfams, create_summary, write_summary,
                         get_stats)
from src.run import run_gffread, run_TEsorter, remove_stop_codons, run_interpro, run_agat


import subprocess

REXDB_PFAMS = Path(os.path.dirname(os.path.realpath(__file__))).parent / "docs" / "Viridiplantae_2.0_pfams.txt"


CATEGORIES = {"No_TE(PcpM0)": "PcpM0", "Protein_TE_only(PteM0)": "PteM0",
              "Chimeric_Protein_Only(PchM0)": "PchM0", "mRNA_TE_Only(PcpMte)": "PcpMte",
              "Protein_and_mRNA_TE(PteMte)": "PteMte", "Chimeric_Protein_and_mRNA_TE(PchMte)": "PchMte",
              "No_Protein_Domains_mRNA_TE(P0Mte)": "P0Mte"}
EXCLUDE = ["AntiFam", "CDD", "Coils", "FunFam",
               "Gene3D", "Hamap", "MobiDBLite",
               "NCBIfam", "PANTHER", "PIRSF", 
               "PIRSR", "PRINTS", "ProSitePatterns",
               "ProSiteProfiles", "SFLD", "SMART", 
               "SUPERFAMILY"]

#Generating program options
def parse_arguments():
    desc = "Pipeline to identify Transposable Elments (TE) in annotated genes"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) File of Files with the following format:
                    "NAME   FASTA   GFF'''
    parser.add_argument("--input", "-i", type=str,
                        help=help_input,
                        required=True)
    
    help_output_dir = '''(Required) Output dir'''
    parser.add_argument("--output", "-out", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_threads = "(Optional) number of threads. 1 by default"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)
    
    help_database = "(Optional) database for TEsorter. rexdb-plant by default"
    parser.add_argument("--tesorter_database", "-d", type=str,
                        help=help_database, default="rexdb-plant")
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()

def get_arguments():
    parser = parse_arguments()
    output = Path(parser.output)
    if not output.exists():
        output.mkdir(parents=True)
    return {"input": parser.input,
            "out": output,
            "threads": parser.threads,
            "tesorter_database": parser.tesorter_database}


def create_header():
    header = ["Run", "Genome", "Annotation", "Annotated_transcripts"]
    for key in CATEGORIES:
        header.append(f"{key}_N")
    for key in CATEGORIES:
        header.append(f"{key}_%")
    header += ["Summary_N", "Summary_%"]
    return "\t".join(header)+"\n"


def get_row(label, genome, annotation, stats):
    inverse_categories = {value: key for key, value in CATEGORIES.items()}
    categories = ["T"] + [key for key in inverse_categories]
    values = [str(stats["num_transcripts"])] + [str(stats[key]) for key in inverse_categories]
    per_values = [str(stats["num_transcripts"])] + [str(round(float(stats[key]/stats["num_transcripts"])*100, 2)) for key in inverse_categories]
    summary = "{0}: {1};{2}: {3};{4}: {5};{6}: {7};{8}: {9};{10}: {11};{12}: {13}"
    row = [label, genome, annotation]
    row += values
    row += per_values[1:]
    row += [summary.format(*[item for pair in zip(categories, values) for item in pair])]
    row += [summary.format(*[item for pair in zip(categories, per_values) for item in pair])]    
    return "\t".join(row)+"\n"


def run_detenga(config, protein_sequences, mrna_sequences):
    report = {"TEsorter": {}, "Stop codons removed": {},
              "InterproScan": {}}
    outdir = Path(config["Basedir"]) / "DETENGA_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    #Run TEsorter
    base_dir = Path(os.getcwd())
    tesorter_outfile = Path("{}.{}.cls.tsv".format(mrna_sequences, config["DETENGA_db"]))
    os.chdir(outdir)
    cmd = "TEsorter {} -db {} -p {}".format(mrna_sequences, config["DETENGA_db"], str(config["Threads"]))
    if tesorter_outfile.is_file():
        msg = "DeTEnGA TEsorter step already done"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "DeTEnGA TEsorter step run successfully"
        else:
            msg = "DeTEnGA TEsorter step Failed: \n {}".format(run_.stderr)
    report["tesorter"] = {"command": cmd,
                          "status": msg,
                          "outfile": tesorter_outfile}
    os.chdir(base_dir)
    
    #REMOVE stop codons
    stop_codons_outfile = outdir / "{}.pep.nostop.fasta".format(Path(config["Assembly"]).stem)
    if stop_codons_outfile.is_file():
        msg = "DeTEnGA Removing stop codons step already done"
    else:
        try:
            id = ""
            sequences_log = []
            stop = False
            original_len = 0
            new_len = 0
            with open(stop_codons_outfile, "w") as out_fhand:
                with open(protein_sequences) as seqs_fhand:
                    for line in seqs_fhand:
                        if line.startswith(">"):
                            if id:
                                sequences_log.append("{}\t{}\t{}\n".format(id, original_len, new_len))
                                original_len = 0
                                new_len = 0
                            id = line.rstrip()[1:]
                            out_fhand.write(line)
                            stop = False
                        else:
                            original_len += len(line.rstrip())
                            if stop:
                                continue
                            else:
                                stop_codons = [".", "*"]
                                for symbol in stop_codons:
                                    if symbol in line:
                                        stop = True
                                        seq = line.split(symbol)[0]
                                        new_len += len(seq)
                                        out_fhand.write(seq+"\n")
                                if not stop:
                                    out_fhand.write(line)
                                    new_len += len(line.rstrip())
            msg = "DeTEnGA Removing stop codons step run successfully"
        except Exception as error:
            msg = "DeTEnGA Removing stop codons step Failed: \n {}".format(error)
    report["stop_codons"] = {"command": "",
                            "status": msg,
                            "outfile": stop_codons_outfile}
    
   
    #Run interproscan
    interpro_outfile = outdir / "{}.pep.nostop.fasta.tsv".format(Path(config["Assembly"]).stem)
    base_dir = Path(os.getcwd())
    os.chdir(outdir)
    cmd = "interproscan.sh -i {} -cpu {} -exclappl {} --disable-precalc".format(interpro_outfile, 
                                                                                config["Threads"], 
                                                                                ",".join(EXCLUDE))
    if interpro_outfile.is_file():
        msg = "DeTEnGA InteproScan analysis step already done"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "DeTEnGA InteproScan analysis step run successfully"
        else:
            msg = "DeTEnGA InteproScan analysis step Failed: \n {}".format(run_.stderr)
    report["InterproScan"] = {"command": cmd,
                              "status": msg,
                              "outfile": interpro_outfile}
        
    os.chdir(base_dir)
    return report
    

#     msg = "##STEP 5: merging evidences from interpro and TEsorter\n"
#     print(msg)
#     log_fhand.write(msg)
#     log_fhand.flush()
#     TE_pfams = get_pfams_from_db(REXDB_PFAMS)
#     summaries = {}
#     for label in sequences:
#         if label in interpro_results and label in TEsorter_results:
#             with open(TEsorter_results[label]["out_fpath"]) as TEsorter_fhand:
#                 te_sorter_output = parse_TEsort_output(TEsorter_fhand)
        
#             with open(interpro_results[label]["out_fpath"]) as interpro_fhand:
#                 interpro = get_pfams_from_interpro_query(interpro_fhand)
#                 classified_pfams = classify_pfams(interpro, TE_pfams)
    
#             te_summary = create_summary(classified_pfams, te_sorter_output)
    
#             out_fpath = Path(out_dir / label / "{}_TE_summary.csv".format(label))
#             with open(out_fpath, "w") as out_fhand:
#                 write_summary(te_summary, out_fhand)
#                 summaries[label] = out_fpath
#                 msg = "TE Summary for {} written in {}\n".format(label, out_fpath)
#                 log_fhand.write(msg)
#                 log_fhand.flush()
    
#     msg = "##STEP 6: Running stats on annotation files\n"
#     print(msg)
#     log_fhand.write(msg)
#     log_fhand.flush()
#     agat_results = run_agat(summaries, files)
#     with open(args["out"]/ "combined_summaries.tsv", "w") as combined_summaries_fhand:
#         header = create_header()
#         combined_summaries_fhand.write(header)
#         for label, results in agat_results.items():
#             stats = get_stats(results["out_fpath"], summaries[label])
#             genome = files[label]["assembly"].stem
#             annotation = files[label]["annotation"].stem
#             row = get_row(label, genome, annotation, stats)
#             combined_summaries_fhand.write(row)


# if __name__ == "__main__":
#     main()

#     base_dir = Path(os.getcwd())
#     for label, values in sequences_input.items():
#         results = {}
#         input_mrna = values["out_fpath"]["mrna"]
#         out_mrna = Path("{}.{}.cls.tsv".format(input_mrna, database))
#         os.chdir(out_mrna.parents[0].absolute())
#         cmd = "TEsorter {} -db {} -p {}".format(input_mrna, database, str(threads))
#         if out_mrna.exists():
#             returncode = 99
#             msg = "File {} already exists\n".format(str(out_mrna))
#         else:
#             run_ = run(cmd, capture_output=True, shell=True)
#             returncode = run_.returncode
#             log = out_mrna.parents[0] / "{}_TEsorter.log.txt".format(label)
#             if returncode == 0:
#                 msg = "Done, check {} for details \n".format(str(log))
#             else:
#                 msg = run_.stderr.decode()
#             with open(log, "w") as log_fhand:
#                 log_fhand.write(run_.stderr.decode())
#         results = {"command": cmd, "returncode": returncode,
#                    "msg": msg, "out_fpath": out_mrna}
#         tesorter_results[label] = results
#         os.chdir(base_dir)
#     return tesorter_results