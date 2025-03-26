import subprocess
from pathlib import Path

def run_gffread(arguments):
    outdir = arguments["output"] / "input_sequences"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "{}.proteins.fasta".format(Path(arguments["ref_assembly"]).stem)
    cmd = "gffread -y {} -g {} {}".format(outfile, 
                                            arguments["ref_assembly"],
                                            arguments["annotation"])

    #Add input sequences file for BUSCO
    arguments["input"] = outfile
    if outfile.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "Extract sequences already done",
                "out_fpath": outfile}
    #But if is not done
    else:
        #Run BUSCO with command
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "GFFread run successfully"
        #But if not
        else:
            msg = "GFFread Failed: \n {}".format(run_.stderr)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "out_fpath": outfile, "returncode": run_.returncode}


def run_busco(arguments):
    #Creating output dir
    outdir = arguments["output"] / "BUSCOCompleteness"
    cmd = "busco -i {} -c {} -o {} --mode prot -l {}".format(arguments["input"],
                                                             arguments["threads"],
                                                              outdir,
                                                            arguments["lineage"])
    if outdir.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "BUSCO already done",
                "out_fpath": outdir}
    else: 
    #Command to run BUSCO
        
        #Run BUSCO with command
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "BUSCO run successfully"
        #But if not
        else:
            msg = "BUSCO Failed: \n {}".format(run_.stderr)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "out_fpath": outdir, "returncode": run_.returncode}


def get_busco_results(busco_results):
    with open(busco_results) as input:
        for line in input:
            if "%" in line:
                return line.strip()