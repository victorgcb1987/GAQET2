import hashlib
import sys 

from argparse import ArgumentParser, RawTextHelpFormatter
from glob import glob
from pathlib import Path
from yaml import safe_load as load_yaml

from src.agat import run_agat_reviewer
from src.agat_parsers import generate_additional_features_reports
from src.YAML import report_yaml_reviewer_file


BULLET_OK = "\t✓\t"
BULLET_FIX = "\tERROR!\t"
HEADER = "-"*5

VERSION = "v1.0.0"


def parse_arguments():
    description = '''\t\t\t###############\n\t\t\t##   GAQET REVIEW   ##\n\t\t\t###############\n
            Parse GAQET results for reviewing\n
            
            Needs a YAML configuration file to work.'''
    parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

    help_config = "(Required) YAML configuration file"
    parser.add_argument("--yaml", "-i",
                        help=help_config, required=True)
    
    help_species = "(Optional) Override YAML species"
    parser.add_argument("--species", "-s", type=str,
                        help=help_species, default="")

    help_class = "(Optional) Override YAML Class"
    parser.add_argument("--clas", "-c", type=str,
                        help=help_class, default="")
   
    help_order = "(Optional) Override YAML Order"
    parser.add_argument("--order", "-d", type=str,
                        help=help_order, default="")
   
    help_genome = "(Optional) Override YAML assembly"
    parser.add_argument("--genome", "-g", type=str,
                        help=help_genome, default="")
    
    help_annot = "(Optional) Override YAML Genome Annotation"
    parser.add_argument("--annotation", "-a", type=str,
                        help=help_annot, default="")
    
    help_repeats = "(Optional) Override YAML repeats"
    parser.add_argument("--repeats", "-r", type=str,
                        help=help_repeats, default="")
    
    help_taxid = "(Optional) Override YAML NCBI taxid"
    parser.add_argument("--taxid", "-t", type=str,
                        help=help_taxid, default="")
    
    help_tolid = "(Optional) Override YAML ToLID taxid"
    parser.add_argument("--tolid", "-l", type=str,
                        help=help_tolid, default="")
    
    help_outbase = "(Optional) Override YAML outbase"
    parser.add_argument("--outbase", "-o", type=str,
                        help=help_outbase, default="")
    
    help_version = "Print version and exit"
    parser.add_argument("--version", "-v", type=str,
                        help=help_version, default="")
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    yaml = parser.yaml
    with open(Path(parser.yaml)) as yaml_fhand:
        yaml = load_yaml(yaml_fhand)

    if parser.species:
        yaml["Species"] = parser.species
    if parser.clas:
        yaml["Class"] = parser.clas
    if parser.order:
        yaml["Order"] = parser.order
    if parser.genome:
        yaml["Assembly"] = parser.genome
    if parser.annotation:
        yaml["Annotation"] = parser.annotation
    if parser.taxid:
        yaml["NCBI_taxid"] = parser.taxid
    if parser.tolid:
        yaml["ToLID"] = parser.tolid
    if parser.outbase:
        yaml["Basedir"] = parser.outbase
    if parser.repeats:
        yaml["Repeats"] = parser.repeats
    
    print(yaml)

    config_report = report_yaml_reviewer_file(yaml)
    #subs whitespaces with _ in ID∫
    yaml["Species"] = "_".join(yaml["Species"].split())
    binary = sys.argv[0]
    command_used = f"{binary} -i {parser.yaml} -s {yaml['Species']} "
    command_used += f"-c {yaml['Class']} "
    command_used += f"-d {yaml['Order']} "
    command_used += f"-t {yaml['NCBI_taxid']} "
    command_used += f"-l {yaml['ToLID']} "
    command_used += f"-g {yaml['Assembly']} "
    command_used += f"-a {yaml['Annotation']} "
    command_used += f"-o {yaml['Basedir']}"
    command_used += f"-r {yaml['Repeats']}"
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

def get_additional_features(basedir):
    additional_features = {"trna": "trna.gff", 
                           "rRNA": "rRNA.gff",
                           "miRNA": "miRNA.gff",
                           "lncRNA": "lnc_rna.gff",
                           }
    features_to_analyze = {}
    for feature, filename in additional_features.items():
        filepath = basedir / "input_sequences" / filename
        if filepath.is_file():
            features_to_analyze[feature] = filepath
    return features_to_analyze

def generate_reviewer_metrics(additional_metrics, arguments, outdir):
    with open(outdir / "Metadata.tsv") as out_fhand:
        with open(arguments["Annotation"], "rb") as annotation_fhand:
            annotation_md5 = hashlib.file_digest(annotation_fhand, "sha256").hexdigest()
        with open(arguments["Assembly"], "rb") as assembly_fhand:
            assembly_md5 = hashlib.file_digest(assembly_fhand, "sha256").hexdigest()
        out_fhand.write(f"Species\t{arguments['Species']}\n")
        out_fhand.write(f"Class\t{arguments['Class']}\n")
        out_fhand.write(f"Order\t{arguments['Order']}\n")
        out_fhand.write(f"NCBI_taxid\t{arguments['NCBI_taxid']}\n")
        out_fhand.write(f"ToLID\t{arguments['ToLID']}\n")
        out_fhand.write(f"Annotation_file\t{arguments['ToLID']}\n")



def main():
    sys.tracebacklimit = 1
    if '--version' in sys.argv or "-v" in sys.argv:
        print(VERSION)
        sys.exit(0)
    arguments, config_report, command_used = get_arguments()
    basedir = Path(arguments["Basedir"])
    outdir = basedir / "REVIEWER"
    general_metrics = glob(basedir/"*.GAQET.stats.tsv")[0]
    generate_reviewer_metrics(general_metrics, arguments, outdir)
    additional_features = get_additional_features(basedir)
    additional_agat_reports = run_agat_reviewer(additional_features, basedir)
    generate_additional_features_reports(additional_agat_reports, outdir)


if __name__ == "__main__":
    main()