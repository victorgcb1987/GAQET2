import subprocess

from pathlib import Path


def succes(outfile):
    if outfile.is_file():
        with open(outfile) as fhand:
             for line in fhand:
                if "psauron score" in line:
                    return True
    return False



def run_psauron(arguments, cds_sequences):
    report = {}
    outdir = Path(arguments["Basedir"]) / "PSAURON_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "{}.cds.psauron.csv".format(arguments["ID"])
    cmd = "psauron -i {} -o {}".format(cds_sequences, outfile)

    if outfile.is_file():
        msg = "PSAURON analysis done already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if succes(outfile):
            msg = "PSAURON analysis run successfully"
        else:
            msg = "PSAURON analysis Failed: \n {}".format(run_.stderr)
    report = {"command": cmd,
              "status": msg,
              "outfile": outfile}
    return report