import subprocess

from pathlib import Path


def run_psauron(arguments, cds_sequences):
    report = {}
    outdir = Path(arguments["Basedir"]) / "PSAURON_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "{}.cds.psauron.csv".format(arguments["ID"])
    cmd = "psauron -i {} -o {}".format(cds_sequences, outfile)

    if outfile.exists():
        msg = "PSAURON analysis done already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "PSAURON analysis run successfully"
        else:
            msg = "PSAURON analysis Failed: \n {}".format(run_.stderr)
    report = {"command": cmd,
              "status": msg,
              "outfile": outfile}
    return report