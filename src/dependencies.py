from os import environ as env
from os.path import (exists, join)
from shutil import which


BINARIES = {"AGAT": ["agat_sp_statistics.pl", 
                     "agat_sp_flag_premature_stop_codons.pl",
                     "agat_sp_filter_incomplete_gene_coding_models.pl"], 
            "BUSCO": ["busco"], 
            "PSAURON": ["psauron"], 
            "DETENGA": ["TEsorter", "interproscan.sh"], 
            "OMARK": ["omamer", "omark"], 
            "PROTHOMOLOGY": ["diamond"]}

HEADER = "-"*5
BULLET_OK = "\t✓\t"
BULLET_FIX = "\tERROR!\t"


def check_dependencies(config):
    report = {"ok": True}
    msg = HEADER+"Checking binaries for {}".format(analysis) + HEADER+ "\n"
    if which["gffread"]:
        msg += BULLET_OK + "Binary gffread found" + "\n"
    else:
        msg += BULLET_FIX + "Binary gffread not found" + "\n"
        report["ok"] = False
    report["gffread"] = msg       
    for analysis in config["Analysis"]:
        msg = HEADER+"Checking binaries for {}".format(analysis) + HEADER = "-"*5 + "\n"
        for binary in analysis:
            if which(binary):
                msg += BULLET_OK + "Binary {} found" + "\n"
            else:
                msg += BULLET_FIX + "Binary {} not found" + "\n"
                report["ok"] = False
        report[analysis] = msg    
    return report