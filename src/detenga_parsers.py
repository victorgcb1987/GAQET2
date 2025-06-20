import os
import re

from collections import defaultdict
from csv import DictReader


CATEGORIES = {"No_TE(PcpM0)": "PcpM0", "Protein_TE_only(PteM0)": "PteM0",
              "Chimeric_Protein_Only(PchM0)": "PchM0", "mRNA_TE_Only(PcpMte)": "PcpMte",
              "Protein_and_mRNA_TE(PteMte)": "PteMte", "Chimeric_Protein_and_mRNA_TE(PchMte)": "PchMte",
              "No_Protein_Domains_mRNA_TE(P0Mte)": "P0Mte"}


def get_pfams_from_db(fpath):
    pfams = {}
    with open(fpath) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            line = line.split()
            pfams[line[0]] = " ".join(line[1:])
    return pfams

def create_header():
    header = ["Run", "Genome", "Annotation", "Annotated_transcripts"]
    for key in CATEGORIES:
        header.append(f"{key}_N")
    for key in CATEGORIES:
        header.append(f"{key}_%")
    header += ["Summary_N", "Summary_%"]
    return "\t".join(header)+"\n"


def get_row(stats):
    results = {}
    inverse_categories = {value: key for key, value in CATEGORIES.items()}
    categories = ["T"] + [key for key in inverse_categories]
    values = [str(stats["num_transcripts"])] + [str(stats[key]) for key in inverse_categories]
    per_values = [str(stats["num_transcripts"])] + [str(round(float(stats[key]/stats["num_transcripts"])*100, 2)) for key in inverse_categories]
    summary = "{0}: {1};{2}: {3};{4}: {5};{6}: {7};{8}: {9};{10}: {11};{12}: {13}"
    results["DETENGA_FPV"] = ";".join([summary.format(*[item for pair in zip(categories, values) for item in pair])])
    results["DETENGA_FP%"] = ";".join([summary.format(*[item for pair in zip(categories, per_values) for item in pair])])  
    return results


def get_pfams_from_interpro_query(fhand):
    genes = defaultdict(list)
    for line in fhand:
        line = line.split("\t")
        if line[3] == "Pfam":
            gen, code, description, start, end, evalue = line[0], line[4], line[5], line[6], line[7], float(line[8])
            if evalue <= 0.005:
                genes[gen].append([code, description, start, end])
    sorted_genes = {
                    key: sorted(value, key=lambda x: int(x[2]))  # Ordena por el tercer valor (convertido a entero)
                    for key, value in genes.items()}
    return sorted_genes


def classify_pfams(interpro, te_pfams):
    for gene, pfams in interpro.items():
        for pfam in pfams:
            if pfam[0] in te_pfams:
                pfam.append("TE")
            else:
                pfam.append("NT")
    return interpro


def parse_TEsort_output(fhand):
    output = defaultdict(list)
    for line in DictReader(fhand,delimiter="\t"):
        output[line["#TE"]] = {"domains": line["Domains"], 
                               "complete": line["Complete"],
                               "classification": "{}|{}|{}".format(line["Order"],
                                                                   line["Superfamily"],
                                                                   line["Clade"]),
                               "strand": line["Strand"]}
    return output


def create_summary(interpro_classified, tesort_output):
    summary = []
    for transcript, values in interpro_classified.items():
        row = {"transcript": transcript}
        transposable = False
        no_transposable = False
        pfams_ids = []
        pfams_descriptions = []
        for value in values:
            pfams_ids.append(value[0])
            pfams_descriptions.append(value[1])
            if "TE" in value:
                transposable = True
            if "NT" in value:
                no_transposable = True
        if transposable and not no_transposable:
            status = "transposable_element"
        if not transposable and no_transposable:
            status = "coding_sequence"
        if transposable and no_transposable:
            status = "mixed"
        row["interpro_status"] = status
        row["pfams_ids"] = "|".join(pfams_ids)
        row["pfams_descriptions"] = "|".join(pfams_descriptions)
        transcript_tesort = tesort_output.get(transcript, None)
        if transcript_tesort is not None:
            row["tesort_domains"] = transcript_tesort["domains"]
            row["tesort_complete"] = transcript_tesort["complete"]
            row["tesort_class"] = transcript_tesort["classification"]
            row["tesort_strand"] = transcript_tesort["strand"]
        else:
            row["tesort_domains"] = "NA"
            row["tesort_complete"] = "NA"
            row["tesort_class"] = "NA"
            row["tesort_strand"] = "NA"
        row["detenga_status"] = detenga_status(row)
        summary.append(row)
    for transcript, values in tesort_output.items():
        if transcript not in interpro_classified:
            row = {"transcript": transcript,
                   "interpro_status": "NA",
                   "pfams_ids": "NA",
                   "pfams_descriptions": "NA",
                   "tesort_domains": values["domains"],
                   "tesort_complete": values["complete"],
                   "tesort_class": values["classification"],
                   "tesort_strand": values["strand"]}
            row["detenga_status"] = detenga_status(row) 
            summary.append(row)
    return summary
        

def detenga_status(row):
    status = "NA"
    if row["interpro_status"] == "coding_sequence" and row["tesort_domains"] == "NA":
        status = "PcpM0" 
    if row["interpro_status"] == "transposable_element" and row["tesort_domains"] == "NA":
        status = "PteM0"
    if row["interpro_status"] == "mixed" and row["tesort_domains"] == "NA":
        status = "PchM0" 
    if row["interpro_status"] == "coding_sequence" and row["tesort_domains"] != "NA":
        status = "PcpMte"
    if row["interpro_status"] == "transposable_element" and row["tesort_domains"] != "NA":
        status = "PteMte"
    if row["interpro_status"] == "mixed" and row["tesort_domains"] != "NA":
        status = "PchMte"
    if row["interpro_status"] == "NA" and row["tesort_domains"] != "NA":
        status = "P0Mte"
    return status


def write_summary(summary, out_fhand):
    out_fhand.write("Transcript_ID;Interpro_status;TEsort_class;PFAM_domains;")
    out_fhand.write("PFAM_descriptions;TEsort_domains;TEsort_completness;")
    out_fhand.write("TEsort_strand;DeTEnGA_status\n")
    out_fhand.flush()
    for row in summary:
        line_total = "" 
        line_total += "{};{};{};{};".format(row["transcript"], row["interpro_status"],
                                            row["tesort_class"], row["pfams_ids"])
        line_total += "{};{};{};{};{}\n".format(row["pfams_descriptions"].replace(";", ","), row["tesort_domains"],
                                                row["tesort_complete"], row["tesort_strand"],
                                                row["detenga_status"])
        out_fhand.write(line_total)
        out_fhand.flush()


def detenga_stats(num_transcripts, summary):
    stats = {"PcpM0": 0, "PteM0": 0, "PchM0": 0, 
             "PcpMte": 0, "PteMte": 0, "PchMte": 0, 
             "P0Mte":0, "num_transcripts": num_transcripts}
    try:
        for row in DictReader(open(summary), delimiter=";"):
            stats[row["DeTEnGA_status"]] += 1
        return get_row(stats)
    except FileNotFoundError:
        return {"DETENGA_FPV":  "FAILED",
                "DETENGA_FP%": "FAILED"}