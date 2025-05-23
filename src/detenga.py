import argparse
import os
import sys

from pathlib import Path

from src.detenga_parsers import (get_pfams_from_interpro_query, parse_TEsort_output, 
                         classify_pfams, create_summary, write_summary, get_pfams_from_db)


import subprocess

REXDB_PFAMS = {"rexdb-plant": Path(os.path.dirname(os.path.realpath(__file__))).parent / "docs" / "Viridiplantae_2.0_pfams.txt",
               "rexdb-metazoa": Path(os.path.dirname(os.path.realpath(__file__))).parent / "docs" / "Metazoa_3.1_pfams.txt",
               "rexdb": Path(os.path.dirname(os.path.realpath(__file__))).parent / "docs" / "Combined_pfams.txt"}


EXCLUDE = ["AntiFam", "CDD", "Coils", "FunFam",
               "Gene3D", "Hamap", "MobiDBLite",
               "NCBIfam", "PANTHER", "PIRSF", 
               "PIRSR", "PRINTS", "ProSitePatterns",
               "ProSiteProfiles", "SFLD", "SMART", 
               "SUPERFAMILY"]


def run_detenga(config, protein_sequences, mrna_sequences):
    report = {"TEsorter": {}, "Stop codons removed": {},
              "InterproScan": {}, "classify_interpro": {},
              "classify_tesorter": {}, "create_summary": {}}
    outdir = Path(config["Basedir"]) / "DETENGA_run"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

    #Run TEsorter
    base_dir = Path(os.getcwd())
    tesorter_outfile = outdir / "{}.{}.cls.tsv".format(mrna_sequences.name, config["DETENGA_db"])
    cmd = "TEsorter {} -db {} -p {}".format(mrna_sequences.absolute(), config["DETENGA_db"], str(config["Threads"]))

    if tesorter_outfile.is_file():
        msg = "DeTEnGA TEsorter step already done"
    else:
        os.chdir(outdir)
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg = "DeTEnGA TEsorter step run successfully"
        else:
            msg = "DeTEnGA TEsorter step Failed: \n {}".format(run_.stderr)
        os.chdir(base_dir)
    report["TEsorter"] = {"command": cmd,
                          "status": msg,
                          "outfile": tesorter_outfile}
    
    
    #REMOVE stop codons
    stop_codons_outfile = outdir / "{}.pep.nostop.fasta".format(Path(config["Assembly"]).stem)
    if stop_codons_outfile.is_file():
        msg = "DeTEnGA Removing stop codons step already done"
    else:
        try:
            id = ""
            sequences_log = []
            stop = False
            original_len = 0
            new_len = 0
            with open(stop_codons_outfile, "w") as out_fhand:
                with open(protein_sequences) as seqs_fhand:
                    for line in seqs_fhand:
                        if line.startswith(">"):
                            if id:
                                sequences_log.append("{}\t{}\t{}\n".format(id, original_len, new_len))
                                original_len = 0
                                new_len = 0
                            id = line.rstrip()[1:]
                            out_fhand.write(line)
                            stop = False
                        else:
                            original_len += len(line.rstrip())
                            if stop:
                                continue
                            else:
                                stop_codons = [".", "*"]
                                for symbol in stop_codons:
                                    if symbol in line:
                                        stop = True
                                        seq = line.split(symbol)[0]
                                        new_len += len(seq)
                                        out_fhand.write(seq+"\n")
                                if not stop:
                                    out_fhand.write(line)
                                    new_len += len(line.rstrip())
            msg = "DeTEnGA Removing stop codons step run successfully"
        except Exception as error:
            msg = "DeTEnGA Removing stop codons step Failed: \n {}".format(error)
    report["Stop codons removed"] = {"command": "",
                            "status": msg,
                            "outfile": stop_codons_outfile}
    
   
    #Run interproscan
    interpro_outfile = outdir / "{}.pep.nostop.fasta.tsv".format(Path(config["Assembly"]).stem)
    base_dir = Path(os.getcwd())

    cmd = "interproscan.sh -i {} -cpu {} -exclappl {} --disable-precalc".format(stop_codons_outfile.absolute(), 
                                                                                config["Threads"], 
                                                                                ",".join(EXCLUDE))
    if interpro_outfile.is_file():
        msg = "DeTEnGA InteproScan analysis step already done"
    else:
        os.chdir(outdir)
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "DeTEnGA InteproScan analysis step run successfully"
        else:
            msg = "DeTEnGA InteproScan analysis step Failed: \n {}".format(run_.stderr)
        os.chdir(base_dir)
    report["InterproScan"] = {"command": cmd,
                              "status": msg,
                              "outfile": interpro_outfile}

    try:
        with open(report["TEsorter"]["outfile"]) as tesorter_fhand:
            te_sorter_output = parse_TEsort_output(tesorter_fhand)
            msg = "DeTEnGA Parse TEsorter step run succesfully"
    except Exception as error:
        msg = "DeTEnGA Parse TEsorter step Failed: \n {}".format(error)
    report["classify_tesorter"] = {"command": "",
                                   "status": msg,
                                   "outfile": ""}
         
    try:
        with open(report["InterproScan"]["outfile"]) as interpro_fhand:
            TE_pfams = get_pfams_from_db(REXDB_PFAMS[config["DETENGA_db"]])
            intepro_pfams = get_pfams_from_interpro_query(interpro_fhand)
            classified_pfams = classify_pfams(intepro_pfams, TE_pfams)
            msg = "DeTEnGA Parse Interpro step run succesfully"
    except Exception as error:
        msg = "DeTEnGA Parse Interpro step Failed: \n {}".format(error)
    report["classify_interpro"] = {"command": "",
                                   "status": msg,
                                   "outfile": ""}
    try:
        te_summary = create_summary(classified_pfams, te_sorter_output)
        outfile = outdir / "{}_TE_summary.csv".format(config["ID"])
        with open(outfile, "w") as out_fhand:
            write_summary(te_summary, out_fhand)
            msg = "DeTEnGA create summary step done"
    except Exception as error:
        msg = "DeTEnGA create summary step done Failed: \n {}".format(error)
    report["create_summary"] = {"command": "",
                                "status": msg,
                                "outfile": outfile}
    return report