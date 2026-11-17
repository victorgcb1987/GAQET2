"""Model-quality metrics re-counted on one transcript per gene (the longest isoform).

BUSCO and OMArk always receive one protein per gene (the longest isoform): with every isoform they
would report multi-isoform genes as duplicated. DIAMOND, DeTEnGA, PSAURON, AGAT START/STOP and the
intron metrics are measured on every transcript. This module repeats those per-transcript counts
on the longest isoform of each gene, from the files the analyses already wrote, so nothing is run
twice. The longest isoforms are the ones selected by agat_sp_keep_longest_isoform.pl, the same
representatives given to BUSCO and OMArk. The denominator is their number (one per gene): coding longest
isoforms for the protein and CDS metrics, every longest isoform for the intron metrics.

Every column is named LONGEST_PREFIX + <name of the per-transcript column>, so gaqet_plot.py can
use them in place of the per-transcript columns for the longest-isoform score.
"""
import csv

from src.agat_parsers import parse_attributes
from src.detenga_parsers import get_row, CATEGORIES
from src.error_check import operation_failed


LONGEST_PREFIX = "Longest isoform: "
MEASURED_COLUMN = LONGEST_PREFIX + "Measured transcripts (N)"
FAILED = "FAILED"


def read_gff_transcripts(gff_fpath):
    """IDs of the mRNA features of a GFF3 file."""
    transcripts = set()
    with open(gff_fpath) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] == "mRNA":
                transcript_id = parse_attributes(fields[8]).get("ID")
                if transcript_id:
                    transcripts.add(transcript_id)
    return transcripts


def percentage(count, total, digits=2):
    return round(100 * count / total, digits) if total else 0.0


def homology_longest(protein_homology, longest):
    """Same rule as protein_homology_stats (any hit with e-value < 1e-20), over the longest isoforms."""
    results = {}
    for tag, values in protein_homology.items():
        column = LONGEST_PREFIX + "ProteinsWith{}Hits (%)".format(tag)
        if operation_failed(values):
            results[column] = FAILED
            continue
        proteins = set()
        with open(values["outfile"]) as fhand:
            for line in fhand:
                parts = line.rstrip().split()
                if len(parts) > 10:
                    try:
                        if float(parts[10]) < 1e-20:
                            proteins.add(parts[0])
                    except ValueError:
                        continue
        results[column] = percentage(len(proteins & longest), len(longest))
    return results


def psauron_longest(psauron, longest):
    """Percentage of the longest-isoform CDS that PSAURON calls coding. As in the per-transcript PSAURON
    score, the denominator is the CDS PSAURON analysed (transcripts with a complete CDS)."""
    column = LONGEST_PREFIX + "PSAURON SCORE"
    if operation_failed(psauron):
        return {column: FAILED}
    analysed, coding = 0, 0
    with open(psauron["outfile"]) as fhand:
        for row in csv.reader(fhand):
            if len(row) < 2 or row[1] not in ("True", "False") or row[0] not in longest:
                continue
            analysed += 1
            coding += row[1] == "True"
    if not analysed:
        return {column: FAILED}
    return {column: percentage(coding, analysed, digits=1)}


def detenga_longest(detenga, longest):
    """DeTEnGA classes of the longest isoforms, as DETENGA_FPV / DETENGA_FP% cells (T = genes)."""
    columns = [LONGEST_PREFIX + "DETENGA_FPV", LONGEST_PREFIX + "DETENGA_FP%"]
    summary = detenga.get("create_summary", {})
    if not summary or operation_failed(summary):
        return {column: FAILED for column in columns}
    stats = {key: 0 for key in CATEGORIES.values()}
    stats["num_transcripts"] = len(longest)
    try:
        with open(summary["outfile"]) as fhand:
            for row in csv.DictReader(fhand, delimiter=";"):
                if row["Transcript_ID"] in longest:
                    stats[row["DeTEnGA_status"]] += 1
    except FileNotFoundError:
        return {column: FAILED for column in columns}
    row = get_row(stats)
    return {LONGEST_PREFIX + key: value for key, value in row.items()}


def agat_incomplete_longest(agat, longest):
    """START/STOP missing counts (AGAT incomplete=1/2/3) among the longest isoforms."""
    names = {"1": "Models START missing", "2": "Models STOP missing", "3": "Models START & STOP missing"}
    columns = {code: LONGEST_PREFIX + name for code, name in names.items()}
    if operation_failed(agat["AGAT incomplete CDS"]):
        return {column: FAILED for column in columns.values()}
    flagged = {code: set() for code in names}
    with open(agat["AGAT incomplete CDS"]["outfile"]) as fhand:
        for line in fhand:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue
            attributes = parse_attributes(fields[8])
            code = attributes.get("incomplete")
            if code in flagged and attributes.get("ID"):
                flagged[code].add(attributes["ID"])
    return {columns[code]: len(ids & longest) for code, ids in flagged.items()}


def introns_longest(agat, longest, artifact_threshold, floor_threshold=None):
    """Longest isoforms with at least one intron <= the potential-artifact threshold and, if given,
    <= the floor threshold (same column names as parse_agat_introns)."""
    cutoffs = [("Transcripts with introns <= {:g}bp (potential artifacts)".format(artifact_threshold),
                artifact_threshold)]
    if floor_threshold is not None:
        cutoffs.append(("Transcripts with introns <= {:g}bp (floor)".format(floor_threshold), floor_threshold))
    columns = [(LONGEST_PREFIX + label + " (N)", LONGEST_PREFIX + "% " + label, threshold)
               for label, threshold in cutoffs]
    if operation_failed(agat["AGAT introns"]):
        return {column: FAILED for count, pct, _ in columns for column in (count, pct)}
    flagged = {threshold: set() for _, _, threshold in columns}
    with open(agat["AGAT introns"]["outfile"]) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "intron":
                continue
            try:
                length = abs(int(fields[4]) - int(fields[3])) + 1
            except ValueError:
                continue
            parents = None
            for threshold in flagged:
                if length <= threshold:
                    if parents is None:
                        parents = parse_attributes(fields[8]).get("Parent", "").split(",")
                    flagged[threshold].update(parent for parent in parents if parent in longest)
    results = {}
    for count_column, pct_column, threshold in columns:
        results[count_column] = len(flagged[threshold])
        results[pct_column] = percentage(len(flagged[threshold]), len(longest))
    return results


def longest_isoform_stats(longest_isoform, artifact_threshold, coding_longest, agat=None, psauron=None,
                          detenga=None, protein_homology=None, floor_threshold=None):
    """Per-gene (longest isoform) versions of the model-quality metrics of the analyses that were run.

    longest_isoform: report of get_longest_isoform (agat_sp_keep_longest_isoform.pl output).
    coding_longest: IDs of the longest isoforms with a CDS (gffread). As in the per-transcript columns, the
    protein and CDS metrics (DIAMOND, DeTEnGA, PSAURON, AGAT START/STOP) are measured on the coding longest
    isoforms (MEASURED_COLUMN) and the intron metrics on every longest isoform (mRNA).
    The other arguments are the reports of the analyses (None if the analysis was not run)."""
    if "Failed" in longest_isoform["status"] or not longest_isoform["outfile"].is_file():
        return {MEASURED_COLUMN: FAILED}
    longest = read_gff_transcripts(longest_isoform["outfile"])
    coding = longest & set(coding_longest)
    results = {MEASURED_COLUMN: len(coding)}
    if protein_homology is not None:
        results.update(homology_longest(protein_homology, coding))
    if psauron is not None:
        results.update(psauron_longest(psauron, coding))
    if detenga is not None:
        results.update(detenga_longest(detenga, coding))
    if agat is not None:
        results.update(agat_incomplete_longest(agat, coding))
        results.update(introns_longest(agat, longest, artifact_threshold, floor_threshold))
    return results
