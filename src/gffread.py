import subprocess

from collections import defaultdict
from pathlib import Path



def reformat_annotation(config):
    outdir = Path(config["Basedir"]) / "input_sequences"
    annotation = config["Annotation"]
    outfile = Path(config["Basedir"]) / "input_sequences" / "reformatted_annotation.gff3"
    report = {"transcripts_to_mRNA": [], "outfile": outfile}
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    with open(outfile, "w") as out_fhand:
        with open(annotation) as annot_fhand:
            for line in annot_fhand:
                if line.startswith("#"):
                    out_fhand.write(line)
                else:
                    line = line.split()
                    if line[2] == "transcript":
                        line[2] = "mRNA"
                        report["transcripts_to_mRNA"].append(line[-1].rstrip())
                    out_fhand.write("\t".join(line)+"\n")
    return report

def run_gffread(config):
    report = {"cds": {"mode": "x", "command": "", "status": "", "outfile": ""}, 
              "proteins": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna": {"mode": "w", "command": "", "status": "", "outfile": ""},
              "cds_longest_isoform": {"mode": "x", "command": "", "status": "", "outfile": ""}, 
              "proteins_longest_isoform": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna_longest_isoform": {"mode": "w", "command": "", "status": "", "outfile": ""},
              "proteins_longest_busco": {"mode": "y", "command": "", "status": "", "outfile": ""}}
    outdir = Path(config["Basedir"]) / "input_sequences"


    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)


    for kind, values in report.items():
        outfile = outdir / "{}.{}.fasta".format(Path(config["Assembly"]).stem, kind)
        if "longest" in kind:
            annotation = config["Annotation_Longest"]
        else:
            annotation = config["Annotation"]
        if "busco" not in kind:
            cmd = "gffread -{} {} -J -g {} {}".format(values["mode"],
                                                      outfile,
                                                      config["Assembly"],
                                                      annotation)
        else:
            cmd = "gffread -{} {} -g {} {}".format(values["mode"],
                                                      outfile,
                                                      config["Assembly"],
                                                      annotation)
        if outfile.exists():
            msg = "{} sequences already extracted".format(kind)
        else:
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            if run_.returncode == 0:
                if "busco" in kind:
                    outfile_renamed = outdir / "{}.{}.renamed.fasta".format(Path(config["Assembly"]).stem, kind)
                    with open(outfile) as fhand:
                        with open(outfile_renamed, "w") as out_fhand:
                            seen = defaultdict(int)
                            for line in fhand:
                                if line.startswith(">"):
                                    header = line.strip()
                                    base_id = header[1:].split()[0]
                                    seen[base_id] += 1
                                    out_fhand.write(f">{base_id}_{seen[base_id]}\n")
                                else:
                                    out_fhand.write(line)
                    outfile = outfile_renamed
                msg = "GFFread, mode {} run successfully".format(kind)
            else:
                msg = "GFFread, mode {} Failed: \n {}".format(kind, run_.stderr)
        report[kind]["command"] = cmd
        report[kind]["status"] = msg
        report[kind]["outfile"] = outfile
    return report
