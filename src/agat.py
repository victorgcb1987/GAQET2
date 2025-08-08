import subprocess
from pathlib import Path


def get_longest_isoform(config):
    report = {}
    outdir = Path(config["Basedir"]) / "input_sequences"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "{}.longest_isoform.gff3".format(Path(config["Assembly"]).stem)
    cmd = "agat_sp_keep_longest_isoform.pl -gff {} -o {}".format(Path(config["Annotation"]), outfile)
    if outfile.exists():
        msg = "Longest isoform from annotation file selected already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)   
        if run_.returncode == 0:
            msg = "AGAT longest isoform run successfully"
        else:
            msg = "AGAT longest isoform Failed: \n {}".format(run_.stdout)
    report = {"command": cmd, "status": msg, 
              "outfile": outfile}
    return report
    

def split_annotation(config):
    report = {}
    outdir = Path(config["Basedir"]) / "input_sequences"
    outfile = outdir / "mrna.gff"
    cmd = "agat_sp_separate_by_record_type.pl  --gff {} -o {}".format(config["Annotation"],
                                                                      outdir)
    if outfile.exists():
        msg = "Annotation splitting by feature type done already"
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)   
        if run_.returncode == 0:
            msg = "AGAT separate by type run successfully"
        else:
            msg = "AGAT separate by type Failed: \n {}".format(run_.stdout)
    report = {"command": cmd, "status": msg, 
              "outfile": outfile}
    return report


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
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
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
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "AGAT premature stop codons analysis run successfully"
        #But if not
        else:
            msg = "AGAT premature stop codons analysis Failed: \n {}".format(run_.stdout)
    
    report["AGAT stop codons"] = {"command": cmd, "status": msg, 
                                  "outfile": premature_stop_outfile}
    

    #Running AGAT incomplete CDS
    incomplete_cds_outfile = outdir / "{}.01_agat_incomplete.txt".format(config["ID"])
    cmd = "agat_sp_filter_incomplete_gene_coding_models.pl --add_flag --gff {} ".format(annot)
    cmd += "--fasta {} -o {}".format(assembly, incomplete_cds_outfile)
    
    if incomplete_cds_outfile.is_file():
        msg = "AGAT incomplete CDS analysis already done"
 
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        #Is process has gone well
        if run_.returncode == 0:
            msg = "AGAT incomplete CDS analysis run successfully"
        #But if not
        else:
            msg = "AGAT incomplete CDS analysis Failed: \n {}".format(run_.stdout)
    
    report["AGAT incomplete CDS"] = {"command": cmd, "status": msg, 
                                     "outfile": incomplete_cds_outfile}
    return report


def run_agat_reviewer(features, basedir):
    report = {}
    outdir = basedir / "AGAT_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    for feature, filepath in features.items():

    #Running AGAT STATS
        stats_outfile = outdir / "{}.01_agat_stats.txt".format(feature)
        cmd = "agat_sp_statistics.pl --gff {} -o {}".format(filepath, stats_outfile)

        if stats_outfile.is_file():
            msg = f"AGAT stats for {feature} already done"

        else:
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            #Is process has gone well
            if run_.returncode == 0:
                msg = f"AGAT stats for {feature} run successfully"
            #But if not
            else:
                msg = f"AGAT stats for {feature} Failed: \n {run_.stdout}"
    
        report[feature] = {"command": cmd, "status": msg, 
                           "outfile": stats_outfile}
    return report