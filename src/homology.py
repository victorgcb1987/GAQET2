import subprocess

from pathlib import Path

def run_protein_homology(config, protein_sequences):
    outdir = Path(config["Basedir"]) / "DIAMOND_run"
    results = {}
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    for db in config["PROTHOMOLOGY_tags"]:
        for tag, db_fpath in db.items():
            outfile = outdir / "{}.proteins.dmd.{}.o6.txt".format(config["ID"], tag)
            cmd = "diamond blastp --threads {} --db {} --query {} --out {}".format(config["Threads"],
                                                                                str(db_fpath),
                                                                                str(protein_sequences),
                                                                                str(outfile))
            if outfile.is_file():
                msg = "Protein homology analysis with {} already done".format(tag)
            else:
                run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            #Is process has gone well
                if run_.returncode == 0:
                    msg = "Protein homology analysis with {} run successfully".format(tag)
            #But if not
                else:
                    msg = "Protein homology analysis with {} Failed: \n {}".format(tag, run_.stderr)
        results[tag] = {"command": cmd, "status": msg, "outfile": outfile}
    return results