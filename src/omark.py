import subprocess

from pathlib import Path


def run_omark(arguments, protein_sequences):
    report = {"OMAMER": {}, "OMARK": {}}
    outdir = Path(arguments["Basedir"]) / "OMARK"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    #Run OMAMER
    omamer_outfile = outdir / "{}_proteins.omamer".format(arguments["ID"])

    cmd = "omamer search --db {} --query {} --out {}".format(arguments["OMARK_db"],
                                                             protein_sequences,
                                                             omamer_outfile)
    if omamer_outfile.exists():
        msg = "OMAMER search analysis done already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "OMAMER search analysis run successfully"
        else:
            msg = "OMAMER search analysis Failed: \n {}".format(run_.stderr)
    report["OMAMER"] = {"command": cmd,
                        "status": msg,
                        "outfile": omamer_outfile}
    
    #Run OMARK
    omark_outfile = outdir / "{}_proteins.omark".format(arguments["ID"]) / "{}_proteins_detailed_summary.txt".format(arguments["ID"])

    cmd = "omark  -f {} -d {} -t {} -o {}".format(omamer_outfile,
                                                  arguments["OMARK_db"],
                                                  arguments["OMARK_taxid"],
                                                  omark_outfile)
    if omark_outfile.is_file():
        msg = "OMARK analysis done already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "OMARK analysis run successfully"
        else:
            msg = "OMARK analysis Failed: \n {}".format(run_.stderr)
    report["OMARK"] = {"command": cmd,
                        "status": msg,
                        "outfile": omark_outfile}
   

    return report