import subprocess

from collections import defaultdict
from pathlib import Path



def split_gff_line(line):
    """The 9 GFF columns. Columns are tab-separated and may contain spaces (e.g. source "Curated Genomic");
    lines with fewer than 9 tab-separated columns are split on whitespace (space-separated files)."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9:
        fields = line.split()
    if len(fields) < 9:
        return None
    return fields[0:8] + [" ".join(fields[8:])]


def gff_attribute(attributes, key):
    for item in attributes.split(";"):
        name, _, value = item.strip().partition("=")
        if name == key:
            return value
    return ""


def unknown_strand_features(annotation):
    """IDs of the features with strand "?" (e.g. trans-spliced organelle genes in NCBI GFFs) and of all their
    descendants. gffread stops at the first "?" strand and writes no sequence at all, so they are removed."""
    children = defaultdict(list)
    unknown = set()
    with open(annotation) as annot_fhand:
        for line in annot_fhand:
            if line.startswith("#"):
                continue
            fields = split_gff_line(line)
            if fields is None:
                continue
            feature_id = gff_attribute(fields[8], "ID")
            for parent in filter(None, gff_attribute(fields[8], "Parent").split(",")):
                children[parent].append(feature_id)
            if fields[6] == "?" and feature_id:
                unknown.add(feature_id)
    removed, pending = set(), list(unknown)
    while pending:
        feature_id = pending.pop()
        if feature_id in removed:
            continue
        removed.add(feature_id)
        pending.extend(child for child in children.get(feature_id, []) if child)
    return removed


def reformat_annotation(config):
    """Writes the annotation with tab-separated columns, "transcript" features renamed to mRNA and the
    features with unknown strand ("?") removed together with their descendants (listed in
    {ID}_unknown_strand_features.tsv)."""
    outdir = Path(config["Basedir"]) / "reformatted_annotation"
    annotation = config["Annotation"]
    outfile = outdir / "reformatted_annotation.gff3"
    removed_file = Path(config["Basedir"]) / "{}_unknown_strand_features.tsv".format(config["ID"])
    report = {"transcripts_to_mRNA": [], "outfile": outfile, "removed_outfile": removed_file,
              "stats": {"Features removed for unknown strand (N)": 0, "Transcripts removed for unknown strand (N)": 0}}
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    to_remove = unknown_strand_features(annotation)
    removed = []
    with open(outfile, "w") as out_fhand:
        with open(annotation) as annot_fhand:
            for line in annot_fhand:
                if line.startswith("#"):
                    out_fhand.write(line)
                    continue
                line = split_gff_line(line)
                if line is None:
                    continue
                description = line[8]
                feature_id = gff_attribute(description, "ID")
                parents = [parent for parent in gff_attribute(description, "Parent").split(",") if parent]
                if (line[6] == "?" or (feature_id and feature_id in to_remove)
                        or any(parent in to_remove for parent in parents)):
                    removed.append((line[0], line[2], line[3], line[4], line[6], feature_id or ",".join(parents)))
                    continue
                if line[2] == "transcript":
                    line[2] = "mRNA"
                    report["transcripts_to_mRNA"].append(description)
                line = "\t".join(line[0:8]) + "\t" + description + "\n"
                out_fhand.write(line)
    with open(removed_file, "w") as out_fhand:
        out_fhand.write("Seqid\tType\tStart\tEnd\tStrand\tID\n")
        for row in removed:
            out_fhand.write("\t".join(row) + "\n")
    report["removed"] = removed
    report["stats"]["Features removed for unknown strand (N)"] = len(removed)
    report["stats"]["Transcripts removed for unknown strand (N)"] = sum(row[1] in ("mRNA", "transcript") for row in removed)
    return report

STOP_CODONS = {"TAA", "TAG", "TGA"}
# Extracted without -J (all coding transcripts): BUSCO input, and the reference sets used to
# find the transcripts that -J leaves out of the other analyses
NOT_FILTERED = {"proteins_longest_busco", "cds_all", "cds_longest_isoform_all"}
EXCLUSION_REASONS = {"no_start": "no START", "no_stop": "no STOP", "in_frame_stop": "in-frame STOP",
                     "length_not_multiple_of_3": "CDS length not multiple of 3", "other": "other"}


def run_gffread(config):
    report = {"cds": {"mode": "x", "command": "", "status": "", "outfile": ""},
              "proteins": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna": {"mode": "w", "command": "", "status": "", "outfile": ""},
              "cds_longest_isoform": {"mode": "x", "command": "", "status": "", "outfile": ""},
              "proteins_longest_isoform": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "mrna_longest_isoform": {"mode": "w", "command": "", "status": "", "outfile": ""},
              "proteins_longest_busco": {"mode": "y", "command": "", "status": "", "outfile": ""},
              "cds_all": {"mode": "x", "command": "", "status": "", "outfile": ""},
              "cds_longest_isoform_all": {"mode": "x", "command": "", "status": "", "outfile": ""}}
    outdir = Path(config["Basedir"]) / "input_sequences"


    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)


    for kind, values in report.items():
        outfile = outdir / "{}.{}.fasta".format(Path(config["Assembly"]).stem, kind)
        if "longest" in kind:
            annotation = config["Annotation_Longest"]
        else:
            annotation = config["Annotation"]
        if kind not in NOT_FILTERED:
            cmd = "gffread -{} {} -J -g {} {}".format(values["mode"],
                                                      outfile,
                                                      config["Assembly"],
                                                      annotation)
        else:
            cmd = "gffread -{} {} -g {} {}".format(values["mode"],
                                                      outfile,
                                                      config["Assembly"],
                                                      annotation)
        if outfile.exists():
            msg = "{} sequences already extracted".format(kind)
        else:
            run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL)
            if run_.returncode == 0:
                if "busco" in kind:
                    outfile_renamed = outdir / "{}.{}.renamed.fasta".format(Path(config["Assembly"]).stem, kind)
                    with open(outfile) as fhand:
                        with open(outfile_renamed, "w") as out_fhand:
                            seen = defaultdict(int)
                            for line in fhand:
                                if line.startswith(">"):
                                    header = line.strip()
                                    base_id = header[1:].split()[0]
                                    seen[base_id] += 1
                                    out_fhand.write(f">{base_id}_{seen[base_id]}\n")
                                else:
                                    out_fhand.write(line)
                    outfile = outfile_renamed
                msg = "GFFread, mode {} run successfully".format(kind)
            else:
                msg = "GFFread, mode {} Failed: \n {}".format(kind, run_.stderr)
        report[kind]["command"] = cmd
        report[kind]["status"] = msg
        report[kind]["outfile"] = outfile
    return report


def read_fasta(fpath):
    sequences = {}
    seq_id = None
    with open(fpath) as fhand:
        for line in fhand:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                sequences[seq_id] = []
            elif seq_id is not None:
                sequences[seq_id].append(line.strip())
    return {seq_id: "".join(parts) for seq_id, parts in sequences.items()}


def cds_problems(cds):
    """Why gffread -J would leave a CDS out (standard genetic code, as gffread uses)."""
    cds = cds.upper()
    problems = []
    if len(cds) % 3:
        problems.append("length_not_multiple_of_3")
    codons = [cds[index:index + 3] for index in range(0, len(cds) - len(cds) % 3, 3)]
    if not codons or codons[0] != "ATG":
        problems.append("no_start")
    if not codons or codons[-1] not in STOP_CODONS:
        problems.append("no_stop")
    if any(codon in STOP_CODONS for codon in codons[:-1]):
        problems.append("in_frame_stop")
    return problems or ["other"]


def excluded_transcripts(config, gffread):
    """Coding transcripts that gffread -J leaves out of the protein analyses (DeTEnGA, DIAMOND,
    PSAURON, OMArk), with the reasons. BUSCO is not affected (its proteins are extracted without -J).

    Returns the stats for the GAQET table and {transcript: [reasons]}; writes the list of excluded
    transcripts to {ID}_excluded_from_protein_analyses.tsv."""
    all_cds = read_fasta(gffread["cds_all"]["outfile"])
    kept = set(read_fasta(gffread["proteins"]["outfile"]))
    longest_all = set(read_fasta(gffread["cds_longest_isoform_all"]["outfile"]))
    longest_kept = set(read_fasta(gffread["proteins_longest_isoform"]["outfile"]))
    longest_excluded = longest_all - longest_kept

    excluded = {transcript: cds_problems(cds) for transcript, cds in all_cds.items() if transcript not in kept}
    outfile = Path(config["Basedir"]) / "{}_excluded_from_protein_analyses.tsv".format(config["ID"])
    with open(outfile, "w") as out_fhand:
        out_fhand.write("Transcript_ID\tLongest_isoform\tReasons\n")
        for transcript, reasons in excluded.items():
            out_fhand.write("{}\t{}\t{}\n".format(transcript, "yes" if transcript in longest_excluded else "no",
                                                  ",".join(reasons)))

    total = len(all_cds)
    # coding transcripts (mRNAs with a CDS): denominator of the per-transcript protein and CDS metrics
    # (DIAMOND, DeTEnGA, AGAT START/STOP); mRNAs without CDS (e.g. NCBI pseudogene mRNAs) are not counted
    stats = {"Coding transcripts (N)": total,
             "Coding transcripts excluded from protein analyses (N)": len(excluded),
             "% Coding transcripts excluded from protein analyses": round(100 * len(excluded) / total, 2) if total else 0.0}
    for reason, label in EXCLUSION_REASONS.items():
        stats["Excluded transcripts with {} (N)".format(label)] = sum(reason in reasons for reasons in excluded.values())
    stats["Longest isoforms excluded from OMArk (N)"] = len(longest_excluded)
    return {"stats": stats, "transcripts": excluded, "outfile": outfile,
            "coding": set(all_cds), "coding_longest": longest_all}
