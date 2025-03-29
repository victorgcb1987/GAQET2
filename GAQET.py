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


    




    
    

    # for name, values in arguments["input"].items():
    #     stats[name] = {}
    #     name_dir = out_dir / name
    #     values["output"] = name_dir
    #     values["threads"] = arguments["threads"]
    #     if not name_dir.exists():
    #         name_dir.mkdir(parents=True, exist_ok=True)

    #     agat_statistics = run_agat(values)
    #     print(agat_statistics)
    #     stats[name]["agat_statistics"] = agat_statistics

    #     gffread_results = run_gffread(values)
    #     print(gffread_results)
    #     busco_results = run_busco(values)
    #     print(busco_results)
    #     stats[name]["busco_results"] = get_busco_results(busco_results)

    #     LAI_out_dir =  create_outdir(values)
    #     print(LAI_out_dir)
    #     values["LAI_dir"] = LAI_out_dir["out_fpath"]
    #     suffixerator =  run_suffixerator(values)
    #     if "returncode" in suffixerator:
    #         if suffixerator["returncode"] == 1:
    #             raise RuntimeError("Suffixerator has failed")
    #     print(suffixerator)
    #     harvest = run_harvest(values)
    #     print(harvest)
    #     finder = run_finder(values)
    #     print(finder)
    #     cat = concatenate_outputs(values)
    #     print(cat)
    #     LTR = run_LTR_retriever(values)
    #     print(LTR)
    #     LAI = run_LAI(values)
    #     print(LAI)
    #     stats[name]["LAI"] = LAI

    #     stringtie = run_stringtie(values)
    #     print(stringtie)
    #     gffcompare = run_gffcompare(values)
    #     print(gffcompare)
    #     annotation_scores = calculate_annotation_scores(values)
    #     print(annotation_scores)
    #     stats[name]["annotation_scores"] = annotation_scores


if __name__ == "__main__":
    main()