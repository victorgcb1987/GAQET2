import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import sys

from pathlib import Path
from yaml import safe_load as load_yaml



from src.agat import run_agat
from src.busco import run_busco
from src.detenga import run_detenga
from src.gffread import run_gffread
from src.omark import run_omark
from src.psauron import run_psauron
from src.YAML import report_yaml_file
from src.agat_parsers import parse_agat_stats, parse_agat_incomplete, parse_agat_premature
from src.busco_parsers import busco_stats
from src.detenga_parsers import detenga_stats
from src.omark_parsers import omark_stats
from src.psauron_parsers import psauron_stats


from pathlib import Path


BULLET_OK = "\t✓\t"
BULLET_FIX = "\tERROR!\t"
HEADER = "-"*5
AVAILABLE_ANALYSIS = ["AGAT", "BUSCO", "PSAURON",
                      "DETENGA", "OMARK", "PROTHOMOLOGY"]

def parse_arguments():
    description = '''\t\t\t###############\n\t\t\t##   GAQET   ##\n\t\t\t###############\n
            Genome Annotation Quality Evaluation Tools\n
            
            Needs a YAML configuration file to work.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

    help_config = "(Required) YAML configuration file"
    parser.add_argument("--yaml", "-i",
                        help=help_config, required=True)
    
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
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    yaml = parser.yaml
    with open(Path(parser.yaml)) as yaml_fhand:
        yaml = load_yaml(yaml_fhand)
        
    if parser.genome:
        yaml["Assembly"] = parser.genome
    if parser.annotation:
        yaml["Annotation"] = parser.annotation
    if parser.taxid:
        yaml["OMARK_taxid"] = parser.taxid
    if parser.outbase:
        yaml["Basedir"] = parser.outbase
    config_report = report_yaml_file(yaml)
    #subs whitespaces with _ in ID∫
    yaml["ID"] = "_".join(yaml["ID"].split())
    return yaml, config_report

def emit_msg(string, log_fhand):
    print(string)
    log_fhand.write("\n"+ string)

def main():
    sys.tracebacklimit = 1
    arguments, config_report = get_arguments()
    basedir = Path(arguments["Basedir"])
    if not basedir.exists():
        basedir.mkdir(parents=True, exist_ok=True)

    log_fpath = basedir / "GAQET.log.txt" 
    log_fhand = open(log_fpath, "w")
    error_msg = "GAQET has failed, {} for details".format(log_fpath.resolve())
    
    header = "\t\t\t###############\n\t\t\t##   GAQET   ##\n\t\t\t###############\n\n" + config_report + "\n"
    emit_msg(header, log_fhand)

    if "ERROR!" in config_report:
        raise RuntimeError(error_msg)
    
    emit_msg("#Results will be stored at {}\n".format(basedir.resolve()), log_fhand)
   
    emit_msg(HEADER + "Extracting CDS and protein sequences" + HEADER + "\n", log_fhand)
    
    #Run analysis
    gffread = run_gffread(arguments)
    for kind, values in gffread.items():
        status =  values["status"]
        emit_msg("#{} extraction, command used: \n\t{}\n".format(kind, values["command"]), log_fhand)
        if "Failed" in status:
            emit_msg(BULLET_FIX + status + "\n", log_fhand)
            emit_msg(HEADER + "GAQET has stopped working", log_fhand)
            raise RuntimeError(error_msg)
        else:
            emit_msg(BULLET_OK + status + "\n", log_fhand)

    for analysis in arguments["Analysis"]:
        if analysis == "AGAT":
            emit_msg(HEADER + "Running AGAT on the GFF file"+ HEADER + "\n", log_fhand)
            agat = run_agat(arguments)
            for mode, values in agat.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(mode, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                    emit_msg(HEADER + "GAQET has stopped working", log_fhand)
                    raise RuntimeError(error_msg)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)

        if analysis == "BUSCO":
            emit_msg(HEADER + "Running BUSCO"+ HEADER + "\n", log_fhand)
            busco = run_busco(arguments, gffread["proteins"]["outfile"])
            for lineage, values in busco.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(lineage, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                    emit_msg(HEADER + "GAQET has stopped working", log_fhand)
                    raise RuntimeError(error_msg)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)

        if analysis == "PSAURON":
            emit_msg(HEADER + "Running PSAURON"+ HEADER + "\n", log_fhand)
            psauron = run_psauron(arguments, gffread["cds"]["outfile"])
            status = psauron["status"]
            emit_msg("#{} command used: \n\t{}\n".format(lineage, values["command"]), log_fhand)
            if "Failed" in status:
                emit_msg(BULLET_FIX + status + "\n", log_fhand)
                emit_msg(HEADER + "GAQET has stopped working", log_fhand)
                raise RuntimeError(error_msg)
            else:
                emit_msg(BULLET_OK + status + "\n", log_fhand)

        if analysis == "OMARK":
            emit_msg(HEADER + "Running OMARK"+ HEADER + "\n", log_fhand)
            omark = run_omark(arguments, gffread["proteins"]["outfile"])
            for analysis, values in omark.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(analysis, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                    emit_msg(HEADER + "GAQET has stopped working", log_fhand)
                    raise RuntimeError(error_msg)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)

        if analysis == "DETENGA":
            emit_msg(HEADER + "Running DeTEnGA"+ HEADER + "\n", log_fhand)
            detenga = run_detenga(arguments, gffread["proteins"]["outfile"], gffread["mrna"]["outfile"])
            for analysis, values in detenga.items():
                status = values["status"]
                emit_msg("#{} command used: \n\t{}\n".format(analysis, values["command"]), log_fhand)
                if "Failed" in status:
                    emit_msg(BULLET_FIX + status + "\n", log_fhand)
                    emit_msg(HEADER + "GAQET has stopped working", log_fhand)
                    raise RuntimeError(error_msg)
                else:
                    emit_msg(BULLET_OK + status + "\n", log_fhand)
        
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
    outfile = Path(arguments["Basedir"]) / "{}_GAQET.stats.txt"
    with open(outfile, "w") as out_fhand:
        header = ["Species", "NCBI_TaxID", "Assembly_Version", "Annotation_Version"] + [stats for stats in results]
        outfile.write("{}\n".format("\t".join(header)))
        row = [arguments["ID"], arguments["OMARK_taxid"], Path(arguments["Assembly"]).name, Path(arguments["Annotation"]).name]
        row += [value for stats, value in results.items()]
        outfile.write("{}\n".format("\t".join(row)))

            
if __name__ == "__main__":
    main()