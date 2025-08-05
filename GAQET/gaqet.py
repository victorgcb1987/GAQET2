#!/usr/bin/env python

import argparse
import sys
import time

from argparse import RawTextHelpFormatter
from datetime import datetime
from pathlib import Path


from pathlib import Path
from yaml import safe_load as load_yaml

from src.agat import run_agat, get_longest_isoform, split_annotation
from src.busco import run_busco
from src.error_check import correct_fasta_length
from src.detenga import run_detenga
from src.dependencies import check_dependencies
from src.gffread import run_gffread
from src.omark import run_omark
from src.psauron import run_psauron
from src.seqtk import reformat_fasta_file
from src.YAML import report_yaml_file
from src.homology import run_protein_homology
from src.agat_parsers import parse_agat_stats, parse_agat_incomplete, parse_agat_premature
from src.busco_parsers import busco_stats
from src.detenga_parsers import detenga_stats
from src.homology_parsers import protein_homology_stats
from src.omark_parsers import omark_stats
from src.psauron_parsers import psauron_stats


BULLET_OK = "\t✓\t"
BULLET_FIX = "\tERROR!\t"
HEADER = "-"*5
AVAILABLE_ANALYSIS = ["AGAT", "BUSCO", "PSAURON",
                      "DETENGA", "OMARK", "PROTHOMOLOGY"]
VERSION = "v1.11.16"


def parse_arguments():
    description = '''\t\t\t###############\n\t\t\t##   GAQET   ##\n\t\t\t###############\n
            Genome Annotation Quality Evaluation Tools\n
            
            Needs a YAML configuration file to work.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

    help_config = "(Required) YAML configuration file"
    parser.add_argument("--yaml", "-i",
                        help=help_config, required=True)
    
    help_genome = "(Optional) Override YAML species"
    parser.add_argument("--species", "-s", type=str,
                        help=help_genome, default="")
     
    help_genome = "(Optional) Override YAML assembly"
    parser.add_argument("--genome", "-g", type=str,
                        help=help_genome, default="")

    help_annot = "(Optional) Override YAML Genome Annotation"
    parser.add_argument("--annotation", "-a", type=str,
                        help=help_annot, default="")
    
    help_taxid = "(Optional) Override NCBI taxid"
    parser.add_argument("--taxid", "-t", type=str,
                        help=help_taxid, default="")
    
    help_outbase = "(Optional) Override YAML outbase"
    parser.add_argument("--outbase", "-o", type=str,
                        help=help_outbase, default="")
    
    help_busco_switch =  "(Optional) disable Filter BUSCO to work with only complete proteins. False by default"
    parser.add_argument("--disable_busco_filter", "-f", action='store_true',
                        help=help_busco_switch)

    
    help_version = "Print version and exit"
    parser.add_argument("--version", "-v", type=str,
                        help=help_version, default="")
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def create_basedir_fpath():
    suffix = datetime.now().isoformat().split(".")[0].replace("-","").replace(":","")
    basedir = Path("./").absolute() / "AnnotationQC_{}".format(suffix)
    return basedir


def get_arguments():
    parser = parse_arguments()
    yaml = parser.yaml
    with open(Path(parser.yaml)) as yaml_fhand:
        yaml = load_yaml(yaml_fhand)

    if parser.species:
        yaml["ID"] = parser.species
    if parser.genome:
        yaml["Assembly"] = parser.genome
    if parser.annotation:
        yaml["Annotation"] = parser.annotation
    if parser.taxid:
        yaml["OMARK_taxid"] = parser.taxid
    if parser.outbase:
        yaml["Basedir"] = parser.outbase
    else:
        basedir = create_basedir_fpath()
        yaml["Basedir"] = basedir
    yaml["disable_busco_filter"] = parser.disable_busco_filter
    config_report = report_yaml_file(yaml)
    #subs whitespaces with _ in ID∫
    yaml["ID"] = "_".join(yaml["ID"].split())
    binary = sys.argv[0]
    command_used = f"{binary} -i {parser.yaml} -s {yaml['ID']} "
    command_used += f"-t {yaml['OMARK_taxid']} "
    command_used += f"-g {yaml['Assembly']} "
    command_used += f"-a {yaml['Annotation']} "
    command_used += f"-o {yaml['Basedir']}"
    return yaml, config_report, command_used


def emit_msg(string, log_fhand):
    for line in string.split("\n"):
        if BULLET_OK in line:
            print('\033[92m'+line+'\033[0m')
        elif BULLET_FIX in line:
            print('\033[31m'+line+'\033[0m')
        else:
            print(line)
    log_fhand.write("\n"+ string)
    log_fhand.flush()


def main():
    sys.tracebacklimit = 0
    if '--version' in sys.argv or "-v" in sys.argv:
        print(VERSION)
        sys.exit(0)
    arguments, config_report, command_used = get_arguments()
    basedir = Path(arguments["Basedir"])
    if not basedir.exists():
        basedir.mkdir(parents=True, exist_ok=True)

    log_fpath = basedir / "GAQET.log.txt" 
    log_fhand = open(log_fpath, "w")
    error_msg = "GAQET has failed, {} for details".format(log_fpath.resolve())
    
    start = time.ctime()
    header = "\t\t\t###############\n\t\t\t##   GAQET   ##\n\t\t\t###############\n\n" + VERSION + "\n" + start + "\n"
    emit_msg(header, log_fhand)
    emit_msg(config_report + "\n", log_fhand)
    emit_msg("You can redo this analysis using the following command:\n", log_fhand)
    emit_msg(f"CMMD: {command_used}\n", log_fhand)

    if "ERROR!" in config_report:
        raise RuntimeError(error_msg)
    
    dependencies = check_dependencies(arguments)
    for analysis, msg in dependencies.items():
        if analysis == "ok":
            continue
        emit_msg(msg, log_fhand)
    if not dependencies["ok"]:
        msg = "Some dependencies are missing. GAQET has stopped working"
        emit_msg(msg, log_fhand)
        raise RuntimeError(msg)
    emit_msg("#Results will be stored at {}\n".format(basedir.resolve()), log_fhand)
   
    #Run analysis
    overall_start = time.time()
    start_time = time.time()
    emit_msg(HEADER + "Splitting annotation by features"+ HEADER + "\n", log_fhand)
    mrna_features = split_annotation(arguments)
    end_time = time.time()
    status = mrna_features["status"]
    emit_msg("#Separate annotation by type, command used: \n\t{}\n".format(mrna_features["command"]), log_fhand)
    if "Failed" in status:
        emit_msg(BULLET_FIX + status + "\n", log_fhand)
    else:
        emit_msg(BULLET_OK + status + "\n", log_fhand)
    emit_msg("Time consumed getting  splitting annotation: {}s\n".format(round(end_time-start_time,2)), log_fhand) 
    #From now, we are using only mRNA features
    arguments["Annotation"] = mrna_features["outfile"]

    start_time = time.time()
    emit_msg(HEADER + "Checking if Assembly file has a correct length format"+ HEADER + "\n", log_fhand)
    correct_seq_length = correct_fasta_length(arguments)
    if not correct_seq_length:
            emit_msg(HEADER + "Reformatting Assembly file with seqtk"+ HEADER + "\n", log_fhand)
            reformatted_assembly = reformat_fasta_file(arguments)
            end_time = time.time()
            status = reformatted_assembly["status"]
            emit_msg("Reformat Assembly file, command used: \n\t{}\n".format(reformatted_assembly["command"]), log_fhand)
            if "Failed" in status:
                emit_msg(BULLET_FIX + status + "\n", log_fhand)
            else:
                emit_msg(BULLET_OK + status + "\n", log_fhand)
                arguments["Assembly"] = reformatted_assembly["outfile"]
            emit_msg("Time consumed reformatting Assembly file : {}s\n".format(round(end_time-start_time,2)), log_fhand) 

    start_time = time.time()
    emit_msg(HEADER + "Getting longest isoforms"+ HEADER + "\n", log_fhand)
    longest_isoform = get_longest_isoform(arguments)
    end_time = time.time()
    status = longest_isoform["status"]
    emit_msg("#Getting longest isoforms, command used: \n\t{}\n".format(longest_isoform["command"]), log_fhand)
    if "Failed" in status:
        emit_msg(BULLET_FIX + status + "\n", log_fhand)
    else:
        emit_msg(BULLET_OK + status + "\n", log_fhand)
    emit_msg("Time consumed getting longest isoforms: {}s\n".format(round(end_time-start_time,2)), log_fhand) 


    #change annotation to work with the longest isoform
    arguments["Annotation_Longest"] = longest_isoform["outfile"]

    start_time = time.time()
    emit_msg(HEADER + "Extracting CDS and protein sequences" + HEADER + "\n", log_fhand)
    gffread = run_gffread(arguments)
    end_time = time.time()
    for kind, values in gffread.items():
        status =  values["status"]
        emit_msg("#{} extraction, command used: \n\t{}\n".format(kind, values["command"]), log_fhand)
        if "Failed" in status:
            emit_msg(BULLET_FIX + status + "\n", log_fhand)
        else:
            emit_msg(BULLET_OK + status + "\n", log_fhand)
    emit_msg("Time consumed extracting CDS and proteins: {}s\n".format(round(end_time-start_time, 2)), log_fhand)


    for analysis in arguments["Analysis"]:
        if analysis == "AGAT":
            emit_msg(HEADER + "Running AGAT on the GFF file"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            agat = run_agat(arguments)
            end_time = time.time()
            for mode, values in agat.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(mode, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running AGAT: {}s\n".format(round(end_time-start_time, 2)), log_fhand)       


        if analysis == "BUSCO":
            emit_msg(HEADER + "Running BUSCO"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            if not arguments["disable_busco_filter"]:
                busco_sequences = gffread["proteins_longest_busco"]["outfile"]
            else:
                busco_sequences = gffread["proteins_longest_isoform"]["outfile"]
            busco = run_busco(arguments, busco_sequences)
            end_time = time.time()
            for lineage, values in busco.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(lineage, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running BUSCO: {}s\n".format(round(end_time-start_time, 2)), log_fhand)  


        if analysis == "PSAURON":
            emit_msg(HEADER + "Running PSAURON"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            psauron = run_psauron(arguments, gffread["cds"]["outfile"])
            end_time = time.time()
            status = psauron["status"]
            emit_msg("#PSAURON command used: \n\t{}\n".format(psauron["command"]), log_fhand)
            if "Failed" in status:
                emit_msg(BULLET_FIX + status + "\n", log_fhand)
            else:
                emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running PSAURON: {}s\n".format(round(end_time-start_time,2)), log_fhand) 

        if analysis == "OMARK":
            emit_msg(HEADER + "Running OMARK"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            omark = run_omark(arguments, gffread["proteins_longest_isoform"]["outfile"])
            end_time = time.time()
            for analysis, values in omark.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(analysis, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running OMARK: {}s\n".format(round(end_time-start_time,2)), log_fhand) 

        if analysis == "DETENGA":
            emit_msg(HEADER + "Running DeTEnGA"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            detenga = run_detenga(arguments, gffread["proteins"]["outfile"], gffread["mrna"]["outfile"])
            end_time = time.time()
            for analysis, values in detenga.items():
                status = values["status"]
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running DeTEnGA: {}s\n".format(round(end_time-start_time, 2)), log_fhand)
        

        if analysis == "PROTHOMOLOGY":
            emit_msg(HEADER + "Running Protein homology"+ HEADER + "\n", log_fhand)
            start_time = time.time()
            protein_homology = run_protein_homology(arguments, gffread["proteins"]["outfile"])
            end_time = time.time()
            for db, values in protein_homology.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(db, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
            emit_msg("Time consumed Running Protein homology: {}s\n".format(round(end_time-start_time, 2)), log_fhand)


    #Get results from analysis
    results = {}
    for analysis in arguments["Analysis"]:
        if analysis == "AGAT":
            results.update(parse_agat_stats(agat))
            results.update(parse_agat_premature(agat))
            results.update(parse_agat_incomplete(agat))
        if analysis == "BUSCO":
            busco_results = busco_stats(busco)
            for lineage, stats in busco_results.items():
                results["Annotation_BUSCO_{}".format(lineage)] = stats
        if analysis == "OMARK":
            results.update(omark_stats(omark))
        if analysis == "PSAURON":
            results.update(psauron_stats(psauron))
        if analysis == "DETENGA":
            results.update(detenga_stats(results["Transcript_Models (N)"], 
                                         detenga["create_summary"]["outfile"]))
        if analysis == "PROTHOMOLOGY":
            results.update(protein_homology_stats(protein_homology, results["Transcript_Models (N)"]))


    outfile = Path(arguments["Basedir"]) / "{}_GAQET.stats.tsv".format(arguments["ID"])
    with open(outfile, "w") as out_fhand:
        header = ["Species", "NCBI_TaxID", "Assembly_Version", "Annotation_Version"] + [stats for stats in results]
        out_fhand.write("{}\n".format("\t".join(header)))
        row = [arguments["ID"], str(arguments["OMARK_taxid"]), Path(arguments["Assembly"]).name, Path(arguments["Annotation"]).name]
        row += [str(value) for stats, value in results.items()]
        out_fhand.write("{}\n".format("\t".join(row)))
    overall_end = time.time()
    emit_msg("GAQET finished successfully at {}, runtime: {} minutes".format(time.ctime(), str((overall_end-overall_start)/60)), log_fhand)


if __name__ == "__main__":
    main()