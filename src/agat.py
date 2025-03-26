import subprocess

def run_agat(arguments):
    #Creating output dir
    outdir = arguments["output"] / "RunAgat" 
    outdir.mkdir(parents=True, exist_ok=True)
    out_fpath = outdir / "ResultAgat.txt"

    #Creating command to run AGAT as a list
    cmd = ["agat_sp_statistics.pl", "--gff", "{}".format(arguments["annotation"]), "-o", "{}".format(out_fpath)]

    # #Adding "distribution" to command if it has been selected
    # if arguments["distribution"]:
    #     dist_arg = ["-d"]
    #     cmd += dist_arg
    # #Adding "plot" to command if it has been selected
    # if arguments["plot"]:
    #     plot_arg = ["-p"]
    #     cmd += plot_arg

    # #Adding "genome size" to command if it has a different value than default
    # if arguments["gs"] > 0:
    #     size_arg = ["-g {}".format(arguments["gs"])]
    #     cmd += size_arg

#Check if AGAT is already done
    if out_fpath.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "AGAT already done",
                "out_fpath": out_fpath, "returncode": 99}
    #But if is not done
    else:
        #Run AGAT with command
        command = ' '.join(cmd)
        print(command)
        run_ = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
        #Is process has gone well
        print(run_.returncode)
        if run_.returncode == 0:
            msg = "AGAT run successfully"
        #But if not
        else:
            msg = "AGAT Failed: \n {}".format(run_.stdout)
        #Return command, final message and output dir path
        return {"command": command, "msg": msg,
                "out_fpath": out_fpath}

#def get_agat_results(agat_results):


#    return = "NG:1212112;NT:34345334"
