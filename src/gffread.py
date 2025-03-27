import subprocess

from pathlib import Path



def run_gffread(config):
    report = {"commands": [], "status": [], "outfiles": []}
    modes = {"cds": "x", "proteins": "y"}
    outdir = config["Basedir"] / "input_sequences"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
        
    for mode, arg in modes.items():
        outfile = outdir / "{}.{}.fasta".format(Path(config["Assembly"]).stem, mode)
        cmd = "gffread -{} {} -g {} {}".format(arg,
                                               outfile,
                                               config["Assembly"],
                                               config["Annotation"])
        if outfile.exists():
            msg = "{} sequences already extracted".format(mode)
        else:
        #Run BUSCO with command
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #Is process has gone well
            if run_.returncode == 0:
                msg = "GFFread, mode {} run successfully".format(mode)
        #But if not
            else:
                msg = "GFFread, mode {} Failed: \n {}".format(run_.stderr)
        report["commands"].append(cmd)
        report["status"].append(msg)
        report["outfiles"].append(outfile)

        #Return command, final message and output dir path
        return report
