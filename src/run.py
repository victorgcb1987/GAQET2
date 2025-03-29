import os


from pathlib import Path
from subprocess import run

def run_gffread(fof, output):
    results_catalog = {}
    for label, values in fof.items():
        out_dir = output / label
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)
        prefix  = values["assembly"].stem
        mrna_out = out_dir / "{}.mRNA.fasta".format(prefix)
        pep_out = out_dir / "{}.pep.fasta".format(prefix)
        
        
        cmd = "gffread -w {} -g {} {}".format(str(mrna_out), 
                                              str(values["assembly"]),
                                              str(values["annotation"]))
        if mrna_out.exists():
            returncode = 99
            msg = "File {} already exists\n".format(str(mrna_out))
        else:
            run_ = run(cmd, capture_output=True, shell=True)
            msg = run_.stderr.decode()
            returncode = run_.returncode
        
        results = {"command": {"mrna": cmd}, "returncode": {"mrna": returncode},
                   "msg": {"mrna": msg}, "out_fpath": {"mrna": mrna_out.absolute()}}
        
        cmd = "gffread -y {} -g {} {}".format(str(pep_out), 
                                              str(values["assembly"]),
                                              str(values["annotation"]))
        if pep_out.exists():
            returncode = 99
            msg = "File {} already exists\n".format(str(pep_out))
        else:
            run_ = run(cmd, capture_output=True, shell=True)
            msg = run_.stderr.decode()
            returncode = run_.returncode

        results["command"]["protein"] = cmd
        results["returncode"]["protein"] = returncode
        results["msg"]["protein"] = msg
        results["out_fpath"]["protein"] = pep_out.absolute()
        results_catalog[label] = results

    return results_catalog


def run_TEsorter(sequences_input, database, threads):
    base_dir = Path(os.getcwd())
    tesorter_results = {}
    for label, values in sequences_input.items():
        results = {}
        input_mrna = values["out_fpath"]["mrna"]
        out_mrna = Path("{}.{}.cls.tsv".format(input_mrna, database))
        os.chdir(out_mrna.parents[0].absolute())
        cmd = "TEsorter {} -db {} -p {}".format(input_mrna, database, str(threads))
        if out_mrna.exists():
            returncode = 99
            msg = "File {} already exists\n".format(str(out_mrna))
        else:
            run_ = run(cmd, capture_output=True, shell=True)
            returncode = run_.returncode
            log = out_mrna.parents[0] / "{}_TEsorter.log.txt".format(label)
            if returncode == 0:
                msg = "Done, check {} for details \n".format(str(log))
            else:
                msg = run_.stderr.decode()
            with open(log, "w") as log_fhand:
                log_fhand.write(run_.stderr.decode())
        results = {"command": cmd, "returncode": returncode,
                   "msg": msg, "out_fpath": out_mrna}
        tesorter_results[label] = results
        os.chdir(base_dir)
    return tesorter_results


def remove_stop_codons(sequences):
    out_fpath = Path("{}/{}.nostop.fasta".format(sequences.parents[0], sequences.stem))
    log_file = Path("{}/internal_stop_codons.log.txt".format(sequences.parents[0]))
    if out_fpath.exists():
        return {"command":  "Remove internal stop codons from {}".format(str(sequences)),
                "msg": "File exists already, check {} for details".format(log_file),
                "out_fpath": out_fpath}

    else:
        id = ""
        sequences_log = []
        stop = False
        original_len = 0
        new_len = 0
        with open(out_fpath, "w") as out_fhand:
            with open(sequences) as seqs_fhand:
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
    with open(log_file, "w") as log_fhand:
        for line in sequences_log:
            log_fhand.write(line)
    return {"command": "Remove internal stop codons", 
            "msg": "Done, check {} for details".format(log_file),
            "out_fpath": out_fpath}    



def run_interpro(sequences, threads):
    exclude = ["AntiFam", "CDD", "Coils", "FunFam",
               "Gene3D", "Hamap", "MobiDBLite",
               "NCBIfam", "PANTHER", "PIRSF", 
               "PIRSR", "PRINTS", "ProSitePatterns",
               "ProSiteProfiles", "SFLD", "SMART", 
               "SUPERFAMILY"]
    base_dir = Path(os.getcwd())
    
    interpro_results = {}
    for label, values in sequences.items():
        sequences = values["out_fpath"]
        os.chdir(sequences.parents[0].absolute())
        out_fpath = Path("{}.tsv".format(values["out_fpath"]))
        log_fpath = Path("{}/interpro.log.txt".format(out_fpath.parents[0]))
        cmd = "interproscan.sh -i {} -cpu {} -exclappl {} --disable-precalc > {}".format(str(values["out_fpath"]), 
                                                                                         threads, ",".join(exclude),
                                                                                         log_fpath)
        if out_fpath.exists():
            returncode = 99
            msg = "File {} already exists, check log {} for details".format(str(out_fpath),
                                                                           str(log_fpath))
        else:
            run_ = run(cmd, shell=True, capture_output=True)
            returncode = run_.returncode
            if returncode == 0:
                msg = "Done, check {} for details".format(log_fpath)
            else:
                msg = run_.stdout.decode()
        interpro_results[label] = {"command": cmd, "msg": msg,
                                   "out_fpath": out_fpath, "returncode": returncode}
        
        os.chdir(base_dir)
    return interpro_results
    
        
def run_agat(summaries, annotations):
    agat_results = {}
    for label, summary in summaries.items():
        base_dir = summary.parents[0].absolute()
        agat_out = base_dir / "{}.agat.stats.txt".format(label)
        annot_file = annotations[label]["annotation"]
        cmd = "agat_sp_statistics.pl --gff {} -o {}".format(str(annot_file), 
                                                            agat_out)
        if agat_out.exists():
            returncode = 99
            msg = "File {} already exists".format(str(agat_out))
        else:
            run_ = run(cmd, shell=True, capture_output=True)
            returncode = run_.returncode
            if returncode == 0:
                msg = "Done"
            else:
                msg = run_.stdout.decode()
        agat_results[label] = {"command": cmd, "msg": msg,
                               "out_fpath": agat_out, "returncode": returncode}
    return agat_results 