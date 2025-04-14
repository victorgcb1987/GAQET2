import os
import subprocess
from pathlib import Path


def run_busco(arguments, protein_sequences):
    report = {}
    outdir = Path(arguments["Basedir"]) / "BUSCOCompleteness_run"
    for lineage in arguments["BUSCO_lineages"]:
        lineage_outdir = outdir / lineage
        #Busco have problems with fullpaths
        busco_dir = "{}/BUSCOCompleteness_run/{}".format(str(arguments["Basedir"]).split("/")[-1], lineage)
        if not lineage_outdir.exists():
            lineage_outdir.mkdir(parents=True, exist_ok=True)

        outfile = lineage_outdir / "run_{}".format(lineage) / "short_summary.txt"
        
        cmd = "busco --cpu {} -i {} -o {} -m prot -l {} --force".format(arguments["Threads"],
                                                                            protein_sequences,
                                                                            busco_dir,
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
                           "outfile": outfile}
    return report