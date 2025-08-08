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


def parse_agat_premature(agat_results):
    error = operation_failed(agat_results["AGAT stop codons"])
    if error:
        return {"Models with early STOP": error}
    results = {"Models with early STOP": 0}
    with open(agat_results["AGAT stop codons"]["outfile"]) as fhand:
        for line in fhand:
            if "genes have been flagged as pseudogene" in line:
                results["Models with early STOP"] = int(line.split())
    return results


def generate_additional_features_reports(features, outdir):
    for feature, filepath in features.items():
        metrics = {"Gene models (N)": "NA", "Average Gene length (bp)": "NA",
                   "Transcript models (N)": "NA", "Exons (N)": "NA", 
                   "Average exons per transcript (N)": "NA", "Single exon gene models (N)": "NA"}
        with open(filepath) as fhand:
            for line in fhand:
                value = line.rstrip().split()[-1]
                if ":" in line:
                    break
                elif "Number of gene" in line:
                    metrics["Gene models (N)"] = value
                elif "mean gene length (bp)" in line:
                    metrics["Average Gene length (bp)"] = value
                elif f"Number of {feature}" in line:
                    metrics[ "Transcript models (N)": 0] = value
                elif "Number of exon" in line:
                    metrics["Exons (N)"] = value
                elif f"Number of single exon {feature}" in line:
                    metrics["Single exon gene models (N)"] = value
                elif f"mean exons per {feature}" in line:
                    metrics["Average exons per transcript (N)"] = value 
            with open(outdir / f"{feature}.metrics.tsv") as out_fhand:
                out_fhand.write(f"Feature\t{feature}\n")
                for metric, value in metrics.items():
                    out_fhand.write(f"{metric}\t{value}\n")