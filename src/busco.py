import os
import subprocess
from pathlib import Path


def run_busco(arguments, protein_sequences):
    execution_path = os.getcwd()
    proteins_path = protein_sequences.resolve()
    outdir = Path(arguments["Basedir"]) / "BUSCOCompleteness_run"
    if not outdir.exists():
        outdir.mkdir()
    os.chdir(outdir)
    report = {}
    for lineage in arguments["BUSCO_lineages"]:
        #lineage_outdir = outdir / lineage
        #Busco have problems with fullpaths
        #if not lineage_outdir.exists():
        #    lineage_outdir.mkdir(parents=True, exist_ok=True)

        outfile = Path(lineage) / "run_{}".format(lineage) / "short_summary.txt"
        
        cmd = "busco --cpu {} -i {} -o {} -m prot -l {} --tar".format(arguments["Threads"],
                                                                      proteins_path,
                                                                      lineage,
                                                                      lineage)
        if outfile.exists():
            msg = "Busco on lineage {} done already".format(lineage)
        else:
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            if run_.returncode == 0:
                msg = "BUSCO analysis with lineage {} run successfully".format(lineage)
            else:
                msg = "BUSCO analysis with lineage {} Failed: \n {}".format(lineage, run_.stderr)
        report[lineage] = {"command": cmd,
                           "status": msg,
                           "outfile": outfile.resolve()}
    os.chdir(execution_path)
    return report