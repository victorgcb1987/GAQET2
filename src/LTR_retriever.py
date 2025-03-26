import os
import subprocess
import shutil

from pathlib import Path


def create_outdir(arguments):
    #Output directory (to save LTR_retriever input and output files) path
    outdir = arguments["output"] / "LAICompleteness"
    outfile = outdir / Path(arguments["ref_assembly"]).name

    if not outdir.exists():
        outdir.mkdir(parents=True)
    if not outfile.exists():
        cmd = f"ln -s {str(arguments['ref_assembly'])} {str(outfile)}"
        run_ = subprocess.run(cmd, shell=True)
    msg = "The output directory for LAICompleteness has been created"

    return {"msg": msg, "out_fpath": outdir}




def run_suffixerator(arguments):
    #output dir path + name for output files
    index = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name
    #suffixerator command
    cmd = "gt suffixerator -db {} -indexname {} -tis -suf -lcp -des -ssp -sds -dna".format(arguments["ref_assembly"], index)

    #Check if suffixerator is already done
    md5 = arguments["LAI_dir"] / "{}.md5".format(Path(arguments["ref_assembly"]).name)
    print(md5)
    if md5.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "suffixerator already done",
                "LAI_dir": index}
    #But if is not done
    else:
        #Run suffixerator
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #If process has gone well, send this message
        if run_.returncode == 0:
            msg = "suffixerator ran successfully"
        #Otherwise, send this error message
        else:
            msg = "suffixerator Failed: \n {}".format(run_.stderr)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "LAI_dir": index, "returncode": run_.returncode}




def run_harvest(arguments):
    #taking suffixerator output dir to use it as input
    index = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name
    #creating output: output path + file .harvest.scn /w ref_assembly file name
    out = arguments["LAI_dir"] / "{}.harvest.scn".format(Path(arguments["ref_assembly"]).name)
    #harvest command
    cmd = "gt ltrharvest -index {} -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > {}".format(index, out)

    #Check if HARVEST is already done
    if out.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "harvest already done",
                "LAI_dir": out}
    #Otherwise, send this error message
    else:
        #Run harvest
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #If process has gone well, send this message
        if run_.returncode == 0:
            msg = "HARVEST ran successfully"
        #But if not
        else:
            msg = "HARVEST Failed: \n {}".format(run_.stderr)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "LAI_dir": out, "returncode": run_.returncode}




def run_finder(arguments):
    #finder command
    cwd = Path(os.getcwd())
    cmd = "LTR_FINDER_parallel -seq {} -threads {} -harvest_out -size 1000000 -time 300".format(arguments["ref_assembly"],
                                                                                                arguments["threads"])

    #Check if FINDER is already done
    out_file = arguments["LAI_dir"] / "{}.finder.combine.scn".format(Path(arguments["ref_assembly"]).name)
    if out_file.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "harvest already done",
                "LAI_dir": arguments["LAI_dir"]}
    #But if is not done
    else:
        #Change the working directory to the "output" path
        os.chdir(arguments["LAI_dir"])
        #Run finder
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #If process has gone well, send this message
        if run_.returncode == 0:
            msg = "FINDER ran successfully"
        #Otherwise, send this error message
        else:
            msg = " FINDER Failed: \n {}".format(run_.stderr)
        #Restore the original working directory
        os.chdir(cwd)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "LAI_dir": arguments["LAI_dir"], "returncode": run_.returncode}




def concatenate_outputs(arguments):
    outpath = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name
    cmd = "cat {}.harvest.scn {}.finder.combine.scn > {}.rawLTR.scn".format(outpath, 
                                                                            outpath, 
                                                                            outpath)

    #Check if "cat" is already done
    out_file = arguments["LAI_dir"] / "{}.rawLTR.scn".format(Path(arguments["ref_assembly"]).name)
    if out_file.exists():
        #Show a message if it is
        return {"command": cmd, "msg": "Concatenation of the output files from Harvest and Finder is already done.",
                "LAI_dir": arguments["LAI_dir"]}
    #But if is not done
    else:
        #Run command
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #If process has gone well, send this message
        if run_.returncode == 0:
            msg = "Concatenation of the output files from Harvest and Finder successfully completed."
        #Otherwise, send this error message
        else:
            msg = "Failed: \n {}".format(run_.stderr)
        #Return command, final message and output dir path
        return {"command": cmd, "msg": msg,
                "LAI_dir": arguments["LAI_dir"], "returncode": run_.returncode}




def run_LTR_retriever(arguments):
    cwd = Path(os.getcwd())
    cmd = "LTR_retriever -genome {} -inharvest {}.rawLTR.scn -threads {}".format(Path(arguments["ref_assembly"]).name,
                                                                                Path(arguments["ref_assembly"]).name,
                                                                                arguments["threads"])
    outfile = arguments["LAI_dir"] / "{}.mod.pass.list".format(Path(arguments["ref_assembly"]).name)
    print(outfile)
    if outfile.exists():
        return {"command": cmd, "msg": "LTR_retriever already done",
                "LAI_dir": arguments["LAI_dir"]}
    else:
        os.chdir(arguments["LAI_dir"])
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "LTR_retriever ran successfully"
        else:
            msg = "LTR_retriever Failed: \n {}".format(run_.stderr)
        os.chdir(cwd)
        return {"command": cmd, "msg": msg,
                "LAI_dir": arguments["LAI_dir"], "returncode": run_.returncode}




def run_LAI(arguments):
    cwd = Path(os.getcwd())
    cmd = "LAI -genome {} -intact {}.mod.pass.list -all {}.mod.out".format(Path(arguments["ref_assembly"]).name,
                                                                            Path(arguments["ref_assembly"]).name,
                                                                            Path(arguments["ref_assembly"]).name)
    outfile = arguments["LAI_dir"] / "{}.mod.out.LAI".format(Path(arguments["ref_assembly"]).name)
    if outfile.exists():
        return {"command": cmd, "msg": "LAI already done",
                "LAI_dir": arguments["LAI_dir"]}
    else:
        os.chdir(arguments["LAI_dir"])
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "LAI ran successfully"
        else:
            msg = "LAI Failed: \n {}".format(run_.stderr)
        os.chdir(cwd)
        return {"command": cmd, "msg": msg,
                "out_fpath": outfile, "returncode": run_.returncode,
                }