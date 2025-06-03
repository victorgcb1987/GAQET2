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
    msg = HEADER+"Checking binaries for gffread" + HEADER+ "\n"
    if which("gffread"):
        msg += BULLET_OK + "Binary gffread found" + "\n"
    else:
        msg += BULLET_FIX + "Binary gffread not found" + "\n"
        report["ok"] = False
    report["gffread"] = msg       
    msg = HEADER+"Checking binaries for seqtk" + HEADER+ "\n"
    if which("seqtk"):
        msg += BULLET_OK + "Binary seqtk found" + "\n"
    else:
        msg += BULLET_FIX + "Binary seqtk not found" + "\n"
        report["ok"] = False
    report["seqtk"] = msg
    for analysis in config["Analysis"]:
        msg = HEADER+"Checking binaries for {}".format(analysis) + HEADER  + "\n"
        for binary in BINARIES[analysis]:
            if which(binary):
                msg += BULLET_OK + "Binary {} found".format(binary) + "\n"
            else:
                msg += BULLET_FIX + "Binary {} not found".format(binary) + "\n"
                report["ok"] = False
        report[analysis] = msg    
    return report