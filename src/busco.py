import subprocess
from pathlib import Path


def run_busco(arguments, protein_sequences):
    report = {}
    outdir = Path(arguments["Basedir"]) / "BUSCOCompleteness"
    for lineage in arguments["BUSCO_lineages"]:
        lineage_outdir = outdir / lineage
        if not lineage_outdir.exists():
            lineage_outdir.mkdir(parents=True, exist_ok=True)
            outfile = lineage_outdir / "run_{}".format(lineage) / "short_summary.txt"
            cmd = "busco --cpu {} -i {} -o {} -m prot -l {}".format(arguments["Threads"],
                                                                    protein_sequences,
                                                                    lineage_outdir,
                                                                    lineage)
            if outfile.exists():
                msg = "Busco on lineage {} done already".format(lineage)
            else:
                run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
                if run_.returncode == 0:
                    msg = "BUSCO analysis with lineage {} run successfully".format(lineage)
                else:
                    msg = "BUSCO analysis with lineage {} Failed: \n {}".format(run_.stderr)
                report[lineage]["command"] = cmd
                report[lineage]["status"] = msg
                report[lineage]["outfile"] = outfile
    return report