import subprocess
from pathlib import Path


def run_gffread(config):
    report = {"cds": {"mode": "x", "command": "", "status": "", "outfile": ""}, 
              "proteins": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna": {"mode": "w", "command": "", "status": "", "outfile": ""},
              "cds_longest_isoform": {"mode": "x", "command": "", "status": "", "outfile": ""}, 
              "proteins_longest_isoform": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna_longest_isoform": {"mode": "w", "command": "", "status": "", "outfile": ""}}
    outdir = Path(config["Basedir"]) / "input_sequences"


    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)


    for kind, values in report.items():
        outfile = outdir / "{}.{}.fasta".format(Path(config["Assembly"]).stem, kind)
        if "longest" in kind:
            annotation = config["Annotation_Longest"]
        else:
            annotation = config["Annotation"]
        cmd = "gffread -{} {} -J -g {} {}".format(values["mode"],
                                               outfile,
                                               config["Assembly"],
                                               annotation)
        if outfile.exists():
            msg = "{} sequences already extracted".format(kind)
        else:
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            if run_.returncode == 0:
                msg = "GFFread, mode {} run successfully".format(kind)
            else:
                msg = "GFFread, mode {} Failed: \n {}".format(kind, run_.stderr)
        report[kind]["command"] = cmd
        report[kind]["status"] = msg
        report[kind]["outfile"] = outfile
    return report
