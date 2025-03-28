import subprocess

from pathlib import Path


def run_agat(config):
    report = {}
    outdir = Path(config["Basedir"]) / "AGAT_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    annot = config["Annotation"]
    assembly = config["Assembly"]

    #Running AGAT STATS
    stats_outfile = outdir / "{}.01_agat_stats.txt".format(config["ID"])
    cmd = "agat_sp_statistics.pl --gff {} -o {}".format(annot, stats_outfile)

    if stats_outfile.is_file():
        msg = "AGAT stats already done"

    else:
        run_ = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "AGAT stats run successfully"
        #But if not
        else:
            msg = "AGAT stats Failed: \n {}".format(run_.stdout)
    
    report["AGAT stats"] = {"command": cmd, "status": msg, 
                            "outfile": stats_outfile}
    
    #Running AGAT premature stop codons
    premature_stop_outfile = outdir / "{}.01_agat_premature_stop.txt".format(config["ID"])
    cmd = "agat_sp_flag_premature_stop_codons.pl --gff {} --fasta {} -o {}".format(annot, 
                                                                                   assembly,
                                                                                   premature_stop_outfile)
    
    if premature_stop_outfile.is_file():
        msg = "AGAT premature stop codons analysis already done"

    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "AGAT premature stop codons analysis run successfully"
        #But if not
        else:
            msg = "AGAT premature stop codons analysis Failed: \n {}".format(run_.stderr)
    
    report["AGAT stop codons"] = {"command": cmd, "status": msg, 
                                  "outfile": premature_stop_outfile}
    

    #Running AGAT incomplete CDS
    incomplete_cds_outfile = outdir / "{}.01_agat_incomplete.txt".format(config["ID"])
    cmd = "agat_sp_filter_incomplete_gene_coding_models.pl --add_flag --gff {} ".format(annot)
    cmd += "--fasta {} -o {}".format(assembly, incomplete_cds_outfile)
    
    if incomplete_cds_outfile.is_file():
        msg = "AGAT incomplete CDS analysis already done"
 
    else:
        run_ = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "AGAT incomplete CDS analysis run successfully"
        #But if not
        else:
            msg = "AGAT incomplete CDS analysis Failed: \n {}".format(run_.stdout)
    
    report["AGAT incomplete CDS"] = {"command": cmd, "status": msg, 
                                     "outfile": incomplete_cds_outfile}
    return report