import os
import subprocess
from pathlib import Path


def run_busco(arguments, protein_sequences):
    execution_path = os.getcwd()
    print(protein_sequences)
    outdir = Path(arguments["Basedir"]) / "BUSCOCompleteness_run"
    if not outdir.exists():
        outdir.mkdir()
    os.chdir(outdir)
    report = {}
    for lineage in arguments["BUSCO_lineages"]:
        lineage_outdir = outdir / lineage
        #Busco have problems with fullpaths
        if not lineage_outdir.exists():
            lineage_outdir.mkdir(parents=True, exist_ok=True)

        outfile = lineage_outdir / "run_{}".format(lineage) / "short_summary.txt"
        
        cmd = "busco --cpu {} -i {} -o run_{} -m prot -l {} --force --tar".format(arguments["Threads"],
                                                                            protein_sequences.resolve(),
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
                           "outfile": outfile}
    os.chdir(execution_path)
    return report