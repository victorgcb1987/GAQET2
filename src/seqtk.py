import subprocess
from pathlib import Path


def reformat_fasta_file(config):
    report = {}
    outdir = Path(config["Basedir"]) / "input_sequences"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "{}.reformatted.fasta".format(Path(config["Assembly"]).stem)
    cmd = "seqtk seq -l 80 {} > {}".format(Path(config["Assembly"]), outfile)
    if outfile.exists():
        msg = "Assembly file reformatted already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)   
        if run_.returncode == 0:
            msg = "Assembly file reformatted successfully"
        else:
            msg = "Assembly file reformating Failed: \n {}".format(run_.stdout)
    report = {"command": cmd, "status": msg, 
              "outfile": outfile}
    return report