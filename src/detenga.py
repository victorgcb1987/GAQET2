import argparse
import os
import sys

from Bio import SeqIO
from pathlib import Path

from src.detenga_parsers import (get_pfams_from_interpro_query, parse_TEsort_output, 
                         classify_pfams, create_summary, write_summary, get_pfams_from_db)

from importlib.resources import files

import subprocess

REXDB_PFAMS = {"rexdb-plant": Path(files('GAQET').joinpath("docs/Viridiplantae_2.0_pfams.txt")),
               "rexdb-metazoa": Path(files('GAQET').joinpath("docs/Metazoa_3.1_pfams.txt")),
               "rexdb": Path(files('GAQET').joinpath("docs/Combined_pfams.txt"))}


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

    filtered_mRNA_outfile = outdir / "{}.mRNA.noNs.no100k.fasta".format(Path(config["Assembly"]).stem)
    sequences_removed_Ns = []
    sequences_removed_length = []
    msg = ""
    
    with open(filtered_mRNA_outfile, "w") as filtered_mRNA_outfhand:
        for record in SeqIO.parse(mrna_sequences.absolute(), "fasta"):
            if "N" in record.seq.upper():
                sequences_removed_Ns.append(record.id)
            elif len(record.seq) >= 100000:
                sequences_removed_length.append(record.id)
            else:
                SeqIO.write(record, filtered_mRNA_outfhand, "fasta")
            

    if sequences_removed_Ns:
        msg = "Warning: the following sequences have Ns in the sequence and has been removed from TEsorter analysis:\n {}\n".format("\n".join(sequences_removed_Ns))
    if sequences_removed_length:
         msg = "Warning: the following sequences have more than 100k residues has been removed from TEsorter analysis:\n {}\n".format("\n".join(sequences_removed_length))

    base_dir = Path(os.getcwd())
    tesorter_outfile = outdir / "{}.{}.cls.tsv".format(filtered_mRNA_outfile.name, config["DETENGA_db"])
    cmd = "TEsorter {} -db {} -p {}".format(filtered_mRNA_outfile.absolute(), config["DETENGA_db"], str(config["Threads"]))

    if tesorter_outfile.is_file():
        msg += "DeTEnGA TEsorter step already done"
    else:
        os.chdir(outdir)
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
        if run_.returncode == 0:
            msg += "DeTEnGA TEsorter step run successfully"
        else:
            msg += "DeTEnGA TEsorter step Failed: \n {}".format(run_.stderr)
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
    outfile = outdir / "{}_TE_summary.csv".format(config["ID"])
    try:
        te_summary = create_summary(classified_pfams, te_sorter_output)
        with open(outfile, "w") as out_fhand:
            write_summary(te_summary, out_fhand)
            msg = "DeTEnGA create summary step done"
    except Exception as error:
        msg = "DeTEnGA create summary step done Failed: \n {}".format(error)
    report["create_summary"] = {"command": "",
                                "status": msg,
                                "outfile": outfile}
    return report