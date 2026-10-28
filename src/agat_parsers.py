import re

from pathlib import Path

from src.error_check import operation_failed

def parse_agat_stats(agat_results):
    results = {
        "Gene_Models (N)": 0,
        "Transcript_Models (N)": 0,
        "CDS_Models (N)": 0,
        "Exons (N)": 0,
        "UTR5' (N)": 0,
        "UTR3' (N)": 0,
        "Both sides UTR' (N)": 0,
        "Overlapping_Gene_Models (N)": 0,
        "Single Exon Gene Models (N)": 0,
        "Single Exon Transcripts (N)": 0,
        "Total Gene Space (Mb)": 0,
        "Mean Gene Model Length (bp)": 0,
        "Mean CDS Model Length (bp)": 0,
        "Mean Exon Length (bp)": 0,
        "Mean Intron Length (bp)": 0,
        "Longest Gene Model Length (bp)": 0,
        "Longest CDS Model Length (bp)": 0,
        "Longest Intron Length (bp)": 0,
        "Shortest Gene Model Length (bp)": 0,
        "Shortest CDS Model Length (bp)": 0,
        "Shortest Intron Length (bp)": 0
    }
    mapping = {
        "Number of gene": "Gene_Models (N)",
        "Number of mrna": "Transcript_Models (N)",
        "Number of cds": "CDS_Models (N)",
        "Number of exon": "Exons (N)",
        "Number of five_prime_utr": "UTR5' (N)",
        "Number of three_prime_utr": "UTR3' (N)",
        "Number of mrnas with utr both sides": "Both sides UTR' (N)", 
        "Number gene overlapping": "Overlapping_Gene_Models (N)",
        "Number of single exon gene": "Single Exon Gene Models (N)",
        "Number of single exon mrna": "Single Exon Transcripts (N)",
        "Total gene length (bp)": "Total Gene Space (Mb)",
        "mean gene length (bp)": "Mean Gene Model Length (bp)",
        "mean cds length (bp)": "Mean CDS Model Length (bp)",
        "mean exon length (bp)": "Mean Exon Length (bp)",
        "mean intron in cds length (bp)": "Mean Intron Length (bp)",
        "Longest gene (bp)": "Longest Gene Model Length (bp)",
        "Longest cds (bp)": "Longest CDS Model Length (bp)",
        "Longest intron into cds part (bp)": "Longest Intron Length (bp)",
        "Shortest gene (bp)": "Shortest Gene Model Length (bp)",
        "Shortest cds piece (bp)": "Shortest CDS Model Length (bp)",
        "Shortest intron into cds part (bp)": "Shortest Intron Length (bp)"
    }

    error = operation_failed(agat_results["AGAT stats"])
    if error:
        return error
    with open(agat_results["AGAT stats"]["outfile"], 'r') as stats_fhand:
        start_mrna = False
        for line in stats_fhand:
            if "-" in line:
                if "mrna" in line:
                    start_mrna = True
                else:
                    start_mrna = False
            if not line.rstrip():
                continue
            if ':' in line:
                break
            if start_mrna:
                try:
                    key, val = line.rsplit(maxsplit=1)
                    key = key.strip()
                    val = int(val.strip())
                    if key in mapping:
                        result_key = mapping[key]
                        if result_key == "Total Gene Space (Mb)":
                            results[result_key] = round(val / 1_000_000, 2)
                        else:
                            results[result_key] = val
                except ValueError:
                    continue  # Skip lines that can't be parsed

    return results


def parse_agat_incomplete(agat_results):
    error = operation_failed(agat_results["AGAT incomplete CDS"])
    if error:
        return {"Models START missing": error,
                "Models STOP missing": error,
                "Models START & STOP missing": error}
    results = {"Models START missing": 0,
               "Models STOP missing": 0,
               "Models START & STOP missing":0}
    with open(agat_results["AGAT incomplete CDS"]["outfile"]) as fhand:
        for line in fhand:
            if "incomplete=1" in line:
                results["Models START missing"] += 1
            if "incomplete=2" in line:
                results["Models STOP missing"] += 1
            if "incomplete=3" in line:
                results["Models START & STOP missing"] += 1
    return results


def parse_attributes(text):
    attributes = {}
    for item in text.strip().split(";"):
        key, sep, value = item.partition("=")
        if sep:
            attributes[key.strip()] = value.strip()
    return attributes


def parse_agat_introns(agat_results, artifact_threshold, review_threshold, summary_prefix=None):
    """Intron length metrics for two cutoffs: potential artifacts and the floor.

    artifact_threshold: introns with length <= this value (bp) are in the range of potential artifacts.
    review_threshold: introns with length <= this value (bp) are within the floor (the suspicious
    intron cutoff, the review band), so this group includes the potential artifacts.

    For each cutoff it reports, as numbers and percentages:
      - introns with that length (% of all introns),
      - transcripts with at least one such intron (% of all transcripts),
      - genes with at least one transcript (isoform) with such an intron (% of all genes).
    Transcripts and genes are counted from the mRNA features of the file, so single-exon
    transcripts are in the denominator.

    If summary_prefix is given, two tables list what is inside the floor cutoff (which includes the
    potential artifacts): <prefix>_transcripts.tsv (one row per transcript with an intron <= the
    floor cutoff) and <prefix>_genes.tsv (one row per gene with at least one such transcript).

    Relies on the GFF3 produced by agat_sp_add_introns.pl (report key
    "AGAT introns"), which materializes explicit intron features so their
    individual lengths can be measured directly, rather than only the
    aggregate mean/longest/shortest values AGAT stats reports.
    """
    cutoffs = [("Introns <= {:g}bp (potential artifacts)".format(artifact_threshold), artifact_threshold),
               ("Introns <= {:g}bp (floor)".format(review_threshold), review_threshold)]

    def columns(label):
        lowered = label.lower()
        return [label + " (N)", "% " + label,
                "Transcripts with " + lowered + " (N)", "% Transcripts with " + lowered,
                "Genes with " + lowered + " (N)", "% Genes with " + lowered]

    error = operation_failed(agat_results["AGAT introns"])
    if error:
        results = {"Introns (N)": error}
        for label, threshold in cutoffs:
            results.update({column: error for column in columns(label)})
        return results

    transcript_gene = {}
    introns = []
    with open(agat_results["AGAT introns"]["outfile"]) as fhand:
        for line in fhand:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "mRNA":
                attributes = parse_attributes(fields[8])
                if "ID" in attributes:
                    transcript_gene[attributes["ID"]] = attributes.get("Parent", "").split(",")[0]
            elif fields[2] == "intron":
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError:
                    continue
                parents = parse_attributes(fields[8]).get("Parent", "").split(",")
                introns.append((abs(end - start) + 1, [parent for parent in parents if parent],
                                fields[0], start, end, fields[6]))

    total_introns = len(introns)
    total_transcripts = len(transcript_gene)
    total_genes = len({gene for gene in transcript_gene.values() if gene})
    results = {"Introns (N)": total_introns}
    for label, threshold in cutoffs:
        selected = [intron[1] for intron in introns if intron[0] <= threshold]
        transcripts = {parent for parents in selected for parent in parents}
        genes = {transcript_gene[transcript] for transcript in transcripts if transcript_gene.get(transcript)}
        counts = [len(selected), len(transcripts), len(genes)]
        totals = [total_introns, total_transcripts, total_genes]
        names = columns(label)
        for index, (count, total) in enumerate(zip(counts, totals)):
            results[names[2 * index]] = count
            results[names[2 * index + 1]] = round(100 * count / total, 2) if total else 0.0
    if summary_prefix is not None:
        write_short_intron_summary(summary_prefix, introns, transcript_gene, artifact_threshold, review_threshold)
    return results


def write_short_intron_summary(prefix, introns, transcript_gene, artifact_threshold, review_threshold):
    """Tables of the transcripts and genes with introns within the cutoffs (see parse_agat_introns)."""
    by_transcript = {}
    for length, parents, seqid, start, end, strand in introns:
        for parent in parents:
            by_transcript.setdefault(parent, []).append((length, seqid, min(start, end), max(start, end), strand))
    transcripts_per_gene = {}
    for transcript, gene in transcript_gene.items():
        transcripts_per_gene.setdefault(gene, []).append(transcript)

    transcript_rows, gene_rows = [], {}
    for transcript in sorted(by_transcript):
        short = sorted(i for i in by_transcript[transcript] if i[0] <= review_threshold)
        if not short:
            continue
        artifacts = [i for i in short if i[0] <= artifact_threshold]
        gene = transcript_gene.get(transcript, "")
        transcript_rows.append([transcript, gene, short[0][1], short[0][4], len(by_transcript[transcript]),
                                len(artifacts), len(short), min(i[0] for i in short),
                                ",".join("{}-{}:{}".format(i[2], i[3], i[0]) for i in short)])
        info = gene_rows.setdefault(gene, {"artifact": [], "floor": [], "shortest": short[0][0]})
        info["floor"].append(transcript)
        if artifacts:
            info["artifact"].append(transcript)
        info["shortest"] = min(info["shortest"], short[0][0])

    prefix = str(prefix)
    with open(prefix + "_transcripts.tsv", "w") as fhand:
        fhand.write("\t".join(["Transcript_ID", "Gene_ID", "Seqid", "Strand", "Introns (N)",
                               "Introns <= {:g}bp (potential artifacts) (N)".format(artifact_threshold),
                               "Introns <= {:g}bp (floor) (N)".format(review_threshold),
                               "Shortest_intron (bp)", "Introns_within_floor (start-end:length)"]) + "\n")
        for row in transcript_rows:
            fhand.write("\t".join(str(value) for value in row) + "\n")
    with open(prefix + "_genes.tsv", "w") as fhand:
        fhand.write("\t".join(["Gene_ID", "Transcripts (N)",
                               "Transcripts with introns <= {:g}bp (potential artifacts) (N)".format(artifact_threshold),
                               "Transcripts with introns <= {:g}bp (floor) (N)".format(review_threshold),
                               "Shortest_intron (bp)", "Transcripts_within_floor"]) + "\n")
        for gene in sorted(gene_rows):
            info = gene_rows[gene]
            fhand.write("\t".join([gene, str(len(transcripts_per_gene.get(gene, []))), str(len(info["artifact"])),
                                   str(len(info["floor"])), str(info["shortest"]), ",".join(info["floor"])]) + "\n")


def parse_agat_premature(agat_results):
    """Number of mRNAs with premature (in-frame) stop codons.

    agat_sp_flag_premature_stop_codons.pl writes its counts to <output>_report.txt, not to the GFF
    given with -o ("We found N cases where mRNAs contain premature stop codons"). If that report is
    missing (other AGAT versions), the mRNAs carrying the "pseudo" attribute in the GFF are counted."""
    error = operation_failed(agat_results["AGAT stop codons"])
    if error:
        return {"Models with early STOP": error}
    gff_fpath = Path(agat_results["AGAT stop codons"]["outfile"])
    report_fpath = gff_fpath.with_name(gff_fpath.stem + "_report.txt")
    if report_fpath.is_file():
        with open(report_fpath) as fhand:
            for line in fhand:
                match = re.search(r"We found (\d+) cases where mRNAs contain premature stop codons", line)
                if match:
                    return {"Models with early STOP": int(match.group(1))}
    count = 0
    with open(gff_fpath) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2].lower() in ("mrna", "transcript") and \
                    "pseudo" in parse_attributes(fields[8]):
                count += 1
    return {"Models with early STOP": count}


def generate_additional_features_reports(features, outdir):
    translations = {"lncRNA": "lnc_rna", "miRNA": "mirna",
                    "rRNA": "rrna", "tRNA": "trna"}
    for feature, results in features.items():
        metrics = {"Gene models (N)": "NA", "Average Gene length (bp)": "NA",
                   "Transcript models (N)": "NA", "Exons (N)": "NA", 
                   "Average exons per transcript (N)": "NA", "Single exon gene models (N)": "NA"}
        with open(results["outfile"]) as fhand:
            translated_feature = translations[feature]
            for line in fhand:
                if line.strip():
                    value = line.rstrip().split()[-1]
                if ":" in line:
                    break
                elif "Number of gene" in line:
                    metrics["Gene models (N)"] = value
                elif "mean gene length (bp)" in line:
                    metrics["Average Gene length (bp)"] = value
                elif f"Number of {translated_feature}" in line:
                    metrics["Transcript models (N)"] = value
                elif "Number of exon" in line:
                    metrics["Exons (N)"] = value
                elif f"Number of single exon {translated_feature}" in line:
                    metrics["Single exon gene models (N)"] = value
                elif f"mean exons per {translated_feature}" in line:
                    metrics["Average exons per transcript (N)"] = value 
            with open(outdir / f"{feature}.metrics.tsv", "w") as out_fhand:
                out_fhand.write(f"Feature\t{feature}\n")
                for metric, value in metrics.items():
                    out_fhand.write(f"{metric}\t{value}\n")