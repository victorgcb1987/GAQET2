import subprocess

from pathlib import Path


def run_gffread(config):
    report = {"cds": {"mode": "x", "command": "", "status": "", "outfile": ""}, 
              "proteins": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna": {"mode": "w", "command": "", "status": "", "outfile": ""}}
    outdir = Path(config["Basedir"]) / "input_sequences"

    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    for kind, values in report.items():
        outfile = outdir / "{}.{}.fasta".format(Path(config["Assembly"]).stem, kind)
        cmd = "gffread -{} {} -J -g {} {}".format(values["mode"],
                                               outfile,
                                               config["Assembly"],
                                               config["Annotation"])
        if outfile.exists():
            msg = "{} sequences already extracted".format(kind)
        else:
        #Run BUSCO with command
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #Is process has gone well
            if run_.returncode == 0:
                msg = "GFFread, mode {} run successfully".format(kind)
        #But if not
            else:
                msg = "GFFread, mode {} Failed: \n {}".format(run_.stderr)
        report[kind]["command"] = cmd
        report[kind]["status"] = msg
        report[kind]["outfile"] = outfile

        #Return command, final message and output dir path
    return report
