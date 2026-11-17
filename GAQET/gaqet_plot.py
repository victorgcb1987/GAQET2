#!/usr/bin/env python

import argparse
import contextlib
import os
import numpy as np
import matplotlib
import matplotlib.colors
import matplotlib.transforms
import matplotlib.pyplot as plt
import pandas as pd
import re
import sys

from argparse import RawTextHelpFormatter
from matplotlib.patches import Wedge, Rectangle
from pathlib import Path


VERSION = "v2.0.1"

# Use a sans-serif font family throughout, suitable for publication figures
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"]


def parse_arguments():
    description = '''\t\t\t#################\n\t\t\t## GAQET  PLOT ##\n\t\t\t#################\n
            Draw an annotation benchmarking plot'''
    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

    help_input = "(Required) GAQET results input"
    parser.add_argument("--input", "-i",
                        help=help_input, required=True)
    
    help_prefix = ("Prefix of the output files (it can include a directory). Two plots and two tables are "
                   "written, with the model-quality metrics measured on all transcripts and on the longest "
                   "isoform of each gene: {prefix}_plot.{format}, {prefix}_longest_isoform_plot.{format}, "
                   "{prefix}_GAS.tsv and {prefix}_longest_isoform_GAS.tsv. Default: GAQET")
    parser.add_argument("--prefix", "-p", type=str,
                        help=help_prefix, default="GAQET")

    help_format = "Image format of the plots: png (default), jpeg or svg"
    parser.add_argument("--format", "-f", type=str, choices=("png", "jpeg", "svg"),
                        help=help_format, default="png")
    
    help_weights = ("Share (percentage) of the Global Annotation Score given to metrics, as key=percentage "
                    "separated by commas (e.g. completeness=40 or busco=40,omark=10). A class key "
                    "(completeness, precision, accuracy) or a role key (homology, busco, omark_completeness, "
                    "omark (consistency), psauron, detenga, agat, introns, omark_fragments) sets the total share of those metrics, split among them as the "
                    "default shares (by analysis). "
                    "Any other key is matched as a substring of one metric name (e.g. swissprot=15, "
                    "embryophyta=20). The remaining percentage is split among the metrics not mentioned in "
                    "proportion to their default shares. Default: the three classes have the same share "
                    "(the GAS is the mean of the completeness, precision and accuracy sub-scores), "
                    "split equally among the analyses of each class (an analysis with several metrics, e.g. "
                    "homology against several databases, counts once)")
    parser.add_argument("--weights", "-w", type=str,
                        help=help_weights, default="")

    help_exclude = ("Comma separated keys (same matching as --weights) of metrics to leave out of "
                    "the Global Annotation Score. They are still drawn in the plot")
    parser.add_argument("--exclude", "-e", type=str,
                        help=help_exclude, default="")

    help_radial = ("Value (%%) at the centre of the radar plot. Default: the lowest plotted value minus 10 "
                   "(rounded down, 0 at least), so differences between annotations are visible and no value is "
                   "clipped. Values below a manual centre are drawn at the centre. Use 0 for the full 0-100%% axis")
    parser.add_argument("--radial_min", "-r", type=float, help=help_radial, default=None)

    help_fill = ("When some annotations lack a value that the other annotations have (e.g. one of ten has no "
                 "PSAURON result), give them 0 for that value and go on (behaviour of earlier versions). "
                 "Without this option GAQET_PLOT lists the missing values and stops. With it, annotations "
                 "without longest-isoform columns make the longest-isoform score be skipped, as before. An "
                 "analysis missing in every annotation is always left out of the score")
    parser.add_argument("--fill_missing_with_zero", action="store_true", help=help_fill)

    help_version = "Print version and exit"
    parser.add_argument("--version", "-v", type=str,
                        help=help_version, default="")
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    return {"input": Path(parser.input),
            "prefix": parser.prefix,
            "format": parser.format,
            "radial_min": parser.radial_min,
            "weights": parser.weights,
            "fill_missing": parser.fill_missing_with_zero,
            "exclude": [key.strip() for key in parser.exclude.split(",") if key.strip()]}


def area_polar_polygon(r_values, theta_values):
    """Calculate area of polygon in polar coordinates (kept for backward compatibility)."""
    r_values = np.array(r_values)
    theta_values = np.array(theta_values)

    # Asegura que el polígono esté cerrado
    if r_values[0] != r_values[-1] or theta_values[0] != theta_values[-1]:
        r_values = np.append(r_values, r_values[0])
        theta_values = np.append(theta_values, theta_values[0])

    area = 0.5 * np.sum(
        r_values[:-1] * r_values[1:] * np.sin(theta_values[1:] - theta_values[:-1])
    )
    return abs(area)

# Gene model descriptors (key, ring label), drawn in the slots of the PSAURON, START/STOP and homology axes
MODEL_DESCRIPTORS = (("gene_models", "Gene Models"), ("transcript_models", "Transcripts models"),
                     ("utr", "Transcripts with\nboth sides UTR"), ("cds_length", "Average CDS\nlength"))

ROLES = ("homology", "busco", "omark_completeness", "omark", "psauron", "detenga", "agat", "introns",
         "omark_fragments")

# The GAS is the mean of three sub-scores, one per class of error, each the mean of its metrics
# (by default each class has 1/3 of the score, split equally among its metrics)
CLASSES = ("completeness", "precision", "accuracy")
CLASS_LABELS = {"completeness": "Completeness", "precision": "Precision", "accuracy": "Accuracy"}
# Okabe-Ito blue, vermilion and bluish green: colour-blind safe and distinct in grayscale (checked with
# the dataviz palette validator, all pairs)
CLASS_COLORS = {"completeness": "#0072B2", "precision": "#D55E00", "accuracy": "#009E73"}


def _mix(color, other, fraction):
    """Blend a hex colour with another one (fraction of the other one, 0-1)."""
    a, b = np.array(matplotlib.colors.to_rgb(color)), np.array(matplotlib.colors.to_rgb(other))
    return matplotlib.colors.to_hex(a * (1 - fraction) + b * fraction)


_C, _F, _A = CLASS_COLORS["completeness"], CLASS_COLORS["precision"], CLASS_COLORS["accuracy"]
# Outer bars: original GAQET colours; bars added later (OMArk consistency and fragments, floor introns)
# take colours not used by the original ones
DESCRIPTOR_COLORS = {
    "busco_S": "#4481cc", "busco_D": "#91b8e8", "busco_F": "#c9d380", "busco_M": "#c0564b",
    "omark_S": "#47CE0D", "omark_D": "#4F9B53D7", "omark_M": "#27C6A9",
    "omark_cons": "#7b52ae", "omark_inco": "#b99ad9", "omark_cont": "#5c3a1e", "omark_unkn": "#c9c9c9",
    "te_related": "#eba434",
    "omark_partial": "#e377c2", "omark_fragmented": "#f4b6dc",
    "floor_introns": "#8c564b",
    "gene_models": "#2d99ba", "transcript_models": "#097c77", "utr": "#D21F22", "cds_length": "#BDD21F",
}
# Original black outlines around bar segments
BAR_EDGE = {"edgecolor": "black"}
# Annotation colours: original tab20 order (extended with tab20b/tab20c for more annotations)
ANNOTATION_COLORS = []
ROLE_CLASS = {"busco": "completeness", "omark_completeness": "completeness",
              "homology": "precision", "omark": "precision", "psauron": "precision",
              "detenga": "precision",
              "agat": "accuracy", "introns": "accuracy", "omark_fragments": "accuracy"}

# Metrics of these roles are left out of the score when their analysis is missing
MISSING_FLAGS = {"omark": "OMARK", "omark_completeness": "OMARK", "omark_fragments": "OMARK",
                 "psauron": "PSAURON", "detenga": "DETENGA", "agat": "STOP", "introns": "INTRONS"}

# OMArk metrics besides consistency (OMArk is always run on the longest isoform of each gene):
# completeness: % of the conserved gene families (HOGs) of the lineage found (single + duplicated);
# fragments: % of the proteins placed in a known family (consistent, inconsistent or contaminant) that
# are neither partial hits nor fragments (< 50% of the family median length). Unknown proteins, which
# OMArk cannot evaluate, are not in the denominator.
OMARK_COMPLETENESS_METRIC = "Conserved gene families found (OMARK completeness)"
OMARK_FRAGMENTS_METRIC = "Placed proteins not fragmented or partial (OMARK fragments)"
AXIS_LABELS = {OMARK_COMPLETENESS_METRIC: "Conserved gene\nfamilies found\n(OMArk)",
               OMARK_FRAGMENTS_METRIC: "Proteins not\nfragmented or\npartial (OMArk)"}


def parse_omark_consistency(text):
    """'Cons:84.92%[P:3.85%;F:2.64%],Inco:1.54%[P:0.59%,F:0.35%],Cont:0.00%,Unkn:13.53%' ->
    {cons, inco, cont, unkn, partial, fragmented} (% of the proteome); NaN if it cannot be parsed."""
    values = {key: float("nan") for key in ("cons", "inco", "cont", "unkn", "partial", "fragmented")}
    if not isinstance(text, str):
        return values
    partial = fragmented = 0.0
    for name, key in (("Cons", "cons"), ("Inco", "inco"), ("Cont", "cont")):
        match = re.search(name + r":([\d.]+)%(?:\[P:([\d.]+)%[;,]\s*F:([\d.]+)%\])?", text)
        if not match:
            return values
        values[key] = float(match.group(1))
        partial += float(match.group(2) or 0)
        fragmented += float(match.group(3) or 0)
    match = re.search(r"Unkn:([\d.]+)%", text)
    values["unkn"] = float(match.group(1)) if match else float("nan")
    values["partial"], values["fragmented"] = partial, fragmented
    return values

INTRON_METRIC = "Transcripts with no artifact (very short) introns (AGAT introns)"

# DeTEnGA classes with TE evidence (TE Pfam domains and/or TE signal in the mRNA): their sum is the
# "Transposon related Transcripts" bar
TE_CLASSES = ["PteM0", "PchM0", "PcpMte", "PteMte", "PchMte", "P0Mte"]

# Columns written by gaqet.py with the model-quality metrics measured on the longest isoform of each
# gene (src/longest_isoform_parsers.py): "Longest isoform: <per-transcript column name>"
LONGEST_PREFIX = "Longest isoform: "
LONGEST_MEASURED = LONGEST_PREFIX + "Measured transcripts (N)"
CODING_COLUMN = "Coding transcripts (N)"
SHORT_INTRON_COLUMN = re.compile(r"^% Transcripts with introns <= [\d.]+bp \(potential artifacts\)$")
FLOOR_INTRON_COLUMN = re.compile(r"^% Transcripts with introns <= ([\d.]+)bp \(floor\)$")

# Input columns of the analyses in the stats table. A column empty in every annotation is treated as an
# analysis that was not run (placeholder, left out of the score); one empty in only some annotations is a
# missing value (check_missing_values)
ANALYSIS_INPUT_COLUMNS = ["DETENGA_FP%", "OMArk Species Composition", "OMArk Completeness Results",
                          "OMArk Consistency Results", "PSAURON SCORE", "Models START missing",
                          "Models STOP missing", "Models START & STOP missing"]
ANALYSIS_INPUT_PREFIXES = ("Annotation_BUSCO_", "ProteinsWith")


def drop_empty_analysis_columns(df):
    """Remove the analysis input columns that have no value in any annotation."""
    empty = [col for col in df.columns
             if (col in ANALYSIS_INPUT_COLUMNS or col.startswith(ANALYSIS_INPUT_PREFIXES)
                 or SHORT_INTRON_COLUMN.match(col) or FLOOR_INTRON_COLUMN.match(col))
             and df[col].isna().all()]
    return df.drop(columns=empty)


def failed_to_missing(df):
    """gaqet.py writes FAILED in the columns of an analysis that failed: make those cells empty, so a failed
    analysis is handled as a missing one (left out of the score if it failed in every annotation, otherwise
    check_missing_values stops), and give back their numeric type to the columns that had them."""
    failed = df.astype(str).apply(lambda column: column.str.strip().str.upper() == "FAILED")
    if not failed.any().any():
        return df
    for column in df.columns[failed.any()]:
        species = df.loc[failed[column], "Species"].astype(str).tolist()
        print("\033[31mWARNING: {} FAILED in {} of {} annotations: {}\033[0m".format(
            column, len(species), len(df), ", ".join(species)))
        df[column] = df[column].mask(failed[column])
        converted = pd.to_numeric(df[column], errors="coerce")
        if converted.notna().sum() == df[column].notna().sum():
            df[column] = converted
    return df


def check_missing_values(df_metrics, fill, measured):
    """Stop if some annotations lack a value that other annotations have. An analysis missing in every
    annotation already has a placeholder and is left out of the score, so any value still empty here is
    missing in only some annotations. With fill=True the empty values are set to 0 (behaviour of earlier
    versions) and a warning is printed."""
    values = df_metrics.drop(columns="Species")
    empty = values.isna()
    if not empty.any().any():
        return
    lines = []
    for column in values.columns[empty.any()]:
        species = df_metrics.loc[empty[column], "Species"].astype(str).tolist()
        lines.append("  {} ({} of {} annotations): {}".format(column.replace("\n", " "), len(species),
                                                             len(df_metrics), ", ".join(species)))
    unit = "the longest isoform of each gene" if measured == "longest" else "all transcripts"
    report = "values missing in some annotations (model-quality metrics on {}):\n{}".format(unit, "\n".join(lines))
    if not fill:
        sys.exit("\033[31mERROR: {}\nCheck the GAQET runs of those annotations (rerun gaqet.py on the same "
                 "Basedir; finished analyses are reused), remove them from the input, or use "
                 "--fill_missing_with_zero to score the missing values as 0\033[0m".format(report))
    print("\033[31mWARNING: {}\nThey are scored as 0 (--fill_missing_with_zero)\033[0m".format(report))
    df_metrics[values.columns] = values.fillna(0.0)


def use_longest_isoform_columns(df):
    """Replace the per-transcript model-quality columns by their longest-isoform versions.

    Returns the column with the number of measured transcripts (one per gene)."""
    longest_columns = [col for col in df.columns if col.startswith(LONGEST_PREFIX)]
    if LONGEST_MEASURED not in df.columns:
        raise ValueError("the '{}' columns are not in the input".format(LONGEST_PREFIX))
    incomplete = df.loc[pd.to_numeric(df[LONGEST_MEASURED], errors="coerce").isna(), "Species"].tolist()
    if incomplete:
        raise ValueError("no longest-isoform values for: {}".format(", ".join(map(str, incomplete))))
    for col in longest_columns:
        if col == LONGEST_MEASURED:
            continue
        base = col[len(LONGEST_PREFIX):]
        # columns of an analysis that was not run stay absent, so the usual placeholder is used
        df[base] = df[col]
    return LONGEST_MEASURED


def metric_role(metric):
    """Return the analysis (role) a metric column belongs to."""
    if metric.startswith(("ProteinsWith", "MissingWithHomology")):
        return "homology"
    if "BUSCO" in metric:
        return "busco"
    for role, tag in (("omark_completeness", "(OMARK completeness)"), ("omark_fragments", "(OMARK fragments)"),
                      ("omark", "(OMARK)"), ("psauron", "(PSAURON)"),
                      ("detenga", "(DeTEnGA)"), ("introns", "(AGAT introns)"), ("agat", "(AGAT)")):
        if tag in metric:
            return role
    return None


def metric_class(metric):
    """Class of error (completeness, precision, accuracy) a metric measures."""
    return ROLE_CLASS.get(metric_role(metric))


def normalize_key(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


CLASS_KEYS = {normalize_key(c): c for c in CLASSES}
CLASS_KEYS[normalize_key("false_positives")] = "precision"  # former name of the precision class
ROLE_KEYS = {normalize_key(r): r for r in ROLES}


def get_available_metrics(metrics, missing_analises):
    """Metrics backed by an analysis that was actually run."""
    available = []
    for metric in metrics:
        if "Missing" in metric:
            continue
        flag = MISSING_FLAGS.get(metric_role(metric))
        if flag is not None and flag in missing_analises:
            continue
        available.append(metric)
    return available


def parse_weights(text):
    """Parse 'busco=40,omark=10' into [("busco", 40.0), ("omark", 10.0)]."""
    spec = []
    for item in [i.strip() for i in text.split(",") if i.strip()]:
        key, sep, value = item.partition("=")
        try:
            share = float(value)
        except ValueError:
            share = None
        if not sep or not key.strip() or share is None or not 0 < share <= 100:
            raise ValueError("Invalid --weights item '{}'. Use key=percentage, with a percentage "
                             "greater than 0 and up to 100 (e.g. busco=40,omark=10). "
                             "To leave a metric out of the score use --exclude".format(item))
        spec.append((key.strip(), share))
    return spec


def match_metrics(key, metrics, available_text):
    """Metrics selected by a class name (completeness, precision, accuracy), a role name
    (busco, omark...) or a substring of a metric name."""
    nkey = normalize_key(key)
    if nkey in CLASS_KEYS:
        matched = [m for m in metrics if metric_class(m) == CLASS_KEYS[nkey]]
        is_role = True
    elif nkey in ROLE_KEYS:
        matched = [m for m in metrics if metric_role(m) == ROLE_KEYS[nkey]]
        is_role = True
    else:
        matched = [m for m in metrics if nkey and nkey in normalize_key(m)]
        is_role = False
    if not matched:
        raise ValueError("'{}' does not match any metric used in this run. Classes: {}. Roles: {}. "
                         "Metrics: {}".format(key, ", ".join(CLASSES), ", ".join(ROLES), available_text))
    return matched, is_role


def default_shares(metrics):
    """Default share (fraction) of each metric: the classes present share the score equally, each class
    share is split equally among its analyses (roles) and each analysis share among its metrics. So an
    analysis with several metrics (e.g. homology against several databases, BUSCO with two lineages)
    counts once in its class, whatever the number of databases or lineages given."""
    by_class = {}
    for metric in metrics:
        by_class.setdefault(metric_class(metric), {}).setdefault(metric_role(metric), []).append(metric)
    return {metric: 1 / len(by_class) / len(roles) / len(members)
            for roles in by_class.values() for members in roles.values() for metric in members}


def resolve_shares(metrics, share_spec=None, exclude=None):
    """Share of the global annotation score (fraction, adding up to 1) of each scored metric.

    GAS = sum(share_i * metric_i). By default the classes (completeness, precision, accuracy)
    have the same share and each class share is split equally among its metrics (default_shares),
    so the GAS is the mean of the three class sub-scores.
    share_spec: [(key, percentage)]. A class key (completeness, precision, accuracy) or a role
    key (busco, omark, psauron, detenga, agat, introns, homology) sets the total share of those
    metrics, split equally among them. Any other key is matched as a substring of one metric name
    and sets the share of that metric alone. The remaining percentage is split among the metrics
    not mentioned in proportion to their default shares.
    exclude: keys (same matching) of metrics left out of the score.
    """
    available_text = "; ".join(metrics)
    excluded = set()
    for key in exclude or []:
        excluded.update(match_metrics(key, metrics, available_text)[0])
    scored = [m for m in metrics if m not in excluded]
    if not scored:
        raise ValueError("All metrics have been excluded from the score")

    shares = {}
    defaults = default_shares(scored)
    for key, percentage in share_spec or []:
        matched, is_role = match_metrics(key, scored, "; ".join(scored))
        if not is_role and len(matched) > 1:
            raise ValueError("'{}' matches several metrics ({}). Use a more specific key "
                             "or a role name".format(key, "; ".join(matched)))
        # a class or role share is split among its metrics as the default shares are (by analysis)
        total_default = sum(defaults[metric] for metric in matched)
        for metric in matched:
            if metric in shares:
                raise ValueError("Metric '{}' receives more than one share".format(metric))
            shares[metric] = percentage * defaults[metric] / total_default

    specified = sum(shares.values())
    remaining = [m for m in scored if m not in shares]
    if specified > 100 + 1e-9:
        raise ValueError("The shares given add up to {:.2f}%, more than 100%".format(specified))
    if remaining:
        if specified >= 100 - 1e-9:
            raise ValueError("The shares given use the whole score, no share left for: {}".format("; ".join(remaining)))
        total_default = sum(defaults[m] for m in remaining)
        for metric in remaining:
            shares[metric] = (100 - specified) * defaults[metric] / total_default
    elif abs(specified - 100) > 1e-6:
        raise ValueError("All metrics have a share, so they must add up to 100% (now {:.2f}%)".format(specified))
    return {metric: shares[metric] / 100 for metric in scored}


def metric_status(metrics, shares, missing_analises):
    """Share (%) and status (scored, excluded, analysis missing) of every metric."""
    available = set(get_available_metrics(metrics, missing_analises))
    status = {}
    for metric in metrics:
        if metric in shares:
            status[metric] = (shares[metric] * 100, "scored")
        elif metric in available:
            status[metric] = (0.0, "excluded")
        else:
            status[metric] = (0.0, "analysis missing")
    return status


def print_share_table(metrics, shares, missing_analises):
    """Print the share of the score of every metric."""
    lines = ["Metric\tClass\tRole\tShare_%\tStatus"]
    for metric, (share, status) in metric_status(metrics, shares, missing_analises).items():
        lines.append("{}\t{}\t{}\t{:.2f}\t{}".format(metric.replace("\n", " "), metric_class(metric),
                                                    metric_role(metric), share, status))
    totals = class_shares(shares)
    lines.append("Classes: " + ", ".join("{} {:.2f}%".format(c, totals.get(c, 0.0) * 100) for c in CLASSES))
    print("Global Annotation Score, share of each metric:\n" + "\n".join(lines))


MEASURED_LABELS = {"all": "all_transcripts", "longest": "longest_isoform"}


def metric_column_name(metric):
    """Short column name of a metric in the score tables (value column; the points column adds ' points')."""
    match = re.match(r"ProteinsWith(.+)Hits \(%\)$", metric)
    if match:
        return "{} hits (%)".format(match.group(1))
    match = re.match(r"BUSCO_C_(.+)$", metric)
    if match:
        return "BUSCO complete {} (%)".format(match.group(1))
    return metric.replace("\n", " ") + " (%)"


def measured_on_label(metric, measured):
    role = metric_role(metric)
    # BUSCO and every OMArk metric (consistency, completeness, fragments) use the longest isoform of each gene
    if role in ("busco", "omark", "omark_completeness", "omark_fragments"):
        return "longest_isoform"
    if role is None:
        return "NA"
    return MEASURED_LABELS[measured]


def write_score_table(df_metrics, df, metrics, shares, missing_analises, measured, measured_column, out_path,
                      weights_path):
    """Score table, one row per annotation (sorted by GAS): GAS, rank, number of measured transcripts
    and, for every metric whose analysis was run, its value (%) and its points in the GAS
    (share x value; the points of an annotation add up to its GAS; 0 for excluded metrics).

    Weights table, one row per metric: role, transcripts it is measured on, share (%) and status
    (scored, excluded, analysis missing). Shares are the same for every annotation."""
    status = metric_status(metrics, shares, missing_analises)
    reported = [metric for metric in metrics if status[metric][1] != "analysis missing"]
    measured_n = dict(zip(df["Species"], df[measured_column])) if measured_column in df.columns else {}
    totals = class_shares(shares)
    header = ["Species", "Measured_transcripts", "Measured_transcripts_N", "GAS", "GAS_rank"]
    # class sub-scores, with the share of the GAS of each class in the column name
    header += ["{} score (share {:.2f}%)".format(CLASS_LABELS[c], totals.get(c, 0.0) * 100) for c in CLASSES]
    for metric in reported:
        name = metric_column_name(metric)
        # the share (weight) of the metric is written in the name of its points column
        header += [name, "{} points (share {:.2f}%)".format(name.replace(" (%)", ""), status[metric][0])]
    ranked = df_metrics.sort_values(by="__Score__", ascending=False).reset_index(drop=True)
    with open(out_path, "w") as fhand:
        fhand.write("\t".join(header) + "\n")
        for rank, row in ranked.iterrows():
            line = [str(row["Species"]), MEASURED_LABELS[measured], str(measured_n.get(row["Species"], "NA")),
                    "{:.2f}".format(row["__Score__"]), str(rank + 1)]
            for error_class in CLASSES:
                value = row["__Class_{}__".format(error_class)]
                line.append("NA" if pd.isna(value) else "{:.2f}".format(value))
            for metric in reported:
                share, state = status[metric]
                points = row[metric] * share / 100 if state == "scored" else 0.0
                line += ["{:.2f}".format(row[metric]), "{:.2f}".format(points)]
            fhand.write("\t".join(line) + "\n")
    with open(weights_path, "w") as fhand:
        fhand.write("\t".join(["Metric", "Class", "Role", "Measured_on", "Share_%", "Status"]) + "\n")
        for metric in metrics:
            share, state = status[metric]
            name = metric_column_name(metric) if state != "analysis missing" else metric.replace("\n", " ")
            fhand.write("\t".join([name, str(metric_class(metric)), str(metric_role(metric)),
                                   measured_on_label(metric, measured), "{:.2f}".format(share), state]) + "\n")


def calculate_weighted_score(row, shares):
    """Global Annotation Score (0-100): sum of share * metric value over the scored metrics."""
    return sum(row[metric] * share for metric, share in shares.items())


def class_shares(shares):
    """Share (fraction) of the GAS of each class: the sum of the shares of its scored metrics."""
    totals = {}
    for metric, share in shares.items():
        totals[metric_class(metric)] = totals.get(metric_class(metric), 0.0) + share
    return totals


def calculate_class_score(row, shares, error_class):
    """Sub-score of one class (0-100): share-weighted mean of its scored metrics (NaN if none).
    GAS = sum over classes of class share * class sub-score."""
    members = {metric: share for metric, share in shares.items() if metric_class(metric) == error_class}
    total = sum(members.values())
    if not total:
        return float("nan")
    return sum(row[metric] * share for metric, share in members.items()) / total


def output_paths(prefix, image_format, measured):
    """Plot, score table, weights table and breakdown plot: {prefix}_plot.{format}, {prefix}_GAS.tsv,
    {prefix}_GAS_weights.tsv, {prefix}_GAS_breakdown.{format} (all transcripts) and the same with
    {prefix}_longest_isoform."""
    base = prefix + ("_longest_isoform" if measured == "longest" else "")
    return (Path(base + "_plot." + image_format), Path(base + "_GAS.tsv"), Path(base + "_GAS_weights.tsv"),
            Path(base + "_GAS_breakdown." + image_format))



def place_class_label(ax, fig, label, start, end, margin=3, centre_clearance=10, outer_limit=98):
    """Move a class label along the bisector of its area to the smallest radius (radar units, 0-100)
    where its box fits between the two edges of the area, measured from the rendered text: narrow
    areas get the label further from the centre. Kept clear of the centre text and of the 100% ring."""
    middle, half_span = (start + end) / 2, (end - start) / 2
    renderer = fig.canvas.get_renderer()
    box = label.get_window_extent(renderer=renderer)
    origin = np.array(ax.transData.transform((middle, 0)))
    pixels_per_unit = np.hypot(*(np.array(ax.transData.transform((middle, 100))) - origin)) / 100
    half_diagonal = np.hypot(box.width, box.height) / 2 / pixels_per_unit + margin
    if half_span <= 0:
        radius = outer_limit - half_diagonal
    else:
        radius = half_diagonal / np.sin(min(half_span, np.pi / 2))
    radius = min(max(radius, half_diagonal + centre_clearance), outer_limit - half_diagonal)
    label.set_position((middle, radius))


def push_axis_labels_out(fig, axis_label_groups, obstacles, step=6, max_steps=30):
    """Move each axis label (and its tool line) outwards along its axis direction, in steps of `step`
    points, until it no longer overlaps any of the obstacles (the outer-bar labels)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    obstacle_boxes = [text.get_window_extent(renderer) for text in obstacles if text.get_text()]
    for angle, tick_label, tool_text in axis_label_groups:
        if not tick_label.get_text():
            continue
        dx, dy = np.cos(angle) * step / 72, np.sin(angle) * step / 72
        base_transform = tick_label.get_transform()
        for n_steps in range(max_steps + 1):
            boxes = [tick_label.get_window_extent(renderer)]
            if tool_text is not None:
                boxes.append(tool_text.get_window_extent(renderer))
            if not any(box.overlaps(other) for box in boxes for other in obstacle_boxes):
                break
            tick_label.set_transform(base_transform + matplotlib.transforms.ScaledTranslation(
                dx * (n_steps + 1), dy * (n_steps + 1), fig.dpi_scale_trans))

def plot_breakdown(df_metrics, shares, out_path, image_format, measured):
    """GAS breakdown: one horizontal bar per annotation (sorted by GAS), split into the points each class
    gives (class share x class sub-score); the sub-score is written inside its segment and the GAS at the end."""
    totals = class_shares(shares)
    classes = [c for c in CLASSES if totals.get(c, 0) > 0]
    ranked = df_metrics.sort_values(by="__Score__", ascending=True).reset_index(drop=True)
    n = len(ranked)
    fig, ax = plt.subplots(figsize=(14, max(3.0, 0.45 * n + 1.8)))
    left = np.zeros(n)
    for error_class in classes:
        scores = ranked["__Class_{}__".format(error_class)].fillna(0.0).to_numpy()
        points = scores * totals[error_class]
        ax.barh(range(n), points, left=left, color=CLASS_COLORS[error_class], edgecolor="white",
                label="{} ({:.1f}% of GAS)".format(CLASS_LABELS[error_class], totals[error_class] * 100))
        for i, (start, width, score) in enumerate(zip(left, points, scores)):
            if width >= 6:
                ax.text(start + width / 2, i, "{:.1f}".format(score), ha="center", va="center",
                        color="white", fontsize=9 if n > 40 else 11, fontweight="bold")
        left += points
    for i, gas in enumerate(ranked["__Score__"]):
        ax.text(gas + 0.8, i, "{:.1f}".format(gas), va="center", fontsize=9 if n > 40 else 11)
    ax.set_yticks(range(n))
    ax.set_yticklabels(ranked["Species"], fontsize=8 if n > 40 else 11)
    ax.set_xlim(0, 108)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Global Annotation Score (points from each class; class sub-score inside each segment)")
    ax.set_title("GAS breakdown, model-quality metrics measured on: {}".format(
        "longest isoform of each gene" if measured == "longest" else "all transcripts"), fontweight="bold",
        pad=30)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(classes), frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, format=image_format, dpi=300, bbox_inches="tight")
    plt.close(fig)


def automatic_radial_min(values):
    """Centre of the radar: the lowest plotted value minus 10 (rounded down), not below 0."""
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return 0.0
    return float(max(0, np.floor(values.min() - 10)))


def main():
    sys.tracebacklimit = 1
    if '--version' in sys.argv or "-v" in sys.argv:
        print(VERSION)
        sys.exit(0)
    arguments = get_arguments()
    if arguments["radial_min"] is not None and not 0 <= arguments["radial_min"] < 100:
        sys.exit("\033[31mERROR: --radial_min must be between 0 and 100 (not included)\033[0m")
    df = pd.read_csv(arguments["input"], sep="\t")
    df = failed_to_missing(df)

    # Longest-isoform columns in only some annotations: stop before writing anything (default)
    if LONGEST_MEASURED in df.columns:
        lacking = df.loc[pd.to_numeric(df[LONGEST_MEASURED], errors="coerce").isna(), "Species"].astype(str)
        if 0 < len(lacking) < len(df) and not arguments["fill_missing"]:
            sys.exit("\033[31mERROR: no longest-isoform values for {} of {} annotations: {}\nRerun gaqet.py on "
                     "the same Basedir of those annotations (finished analyses are reused) or use "
                     "--fill_missing_with_zero to skip the longest-isoform score\033[0m".format(
                         len(lacking), len(df), ", ".join(lacking)))

    # The model-quality metrics (DIAMOND, DeTEnGA, PSAURON, AGAT START/STOP, short introns) are measured
    # on every transcript and on the longest isoform of each gene: one plot and one table for each.
    # BUSCO and OMArk use the longest isoform of each gene in both.
    try:
        df_longest = df.copy()
        measured_column = use_longest_isoform_columns(df_longest)
    except ValueError as error:
        df_longest = None
        longest_error = error
    # Denominator of AGAT START/STOP on all transcripts: coding transcripts (gaqet.py "Coding transcripts (N)";
    # older outputs: AGAT "CDS_Models (N)", the number of mRNAs with a CDS). mRNAs without CDS are not counted
    all_column = next((col for col in (CODING_COLUMN, "CDS_Models (N)") if col in df.columns),
                      "Transcript_Models (N)")
    # missing values of both scores are checked before any output is written
    make_plot(df.copy(), arguments, "all", all_column, check_only=True)
    if df_longest is not None:
        make_plot(df_longest.copy(), arguments, "longest", measured_column, check_only=True)
    make_plot(df.copy(), arguments, "all", all_column)
    if df_longest is None:
        print("\033[31mWARNING: longest-isoform score not made: {}. Rerun gaqet.py on the same Basedir of "
              "those annotations (finished analyses are reused) to add them\033[0m".format(longest_error))
        return
    make_plot(df_longest, arguments, "longest", measured_column)


def make_plot(df, arguments, measured, measured_column="Transcript_Models (N)", check_only=False):
    """GAQET plot and score table of one unit of measure (measured: 'all' or 'longest').
    check_only: only read the metrics and stop if values are missing (nothing written; the messages are
    printed by the real run)."""
    if check_only is True:
        with open(os.devnull, "w") as devnull, contextlib.redirect_stdout(devnull):
            return make_plot(df, arguments, measured, measured_column, check_only="silenced")
    plot_path, table_path, weights_path, breakdown_path = output_paths(arguments["prefix"], arguments["format"],
                                                                       measured)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    print("\n#### Model-quality metrics measured on: {}".format(
        "the longest isoform of each gene" if measured == "longest" else "all transcripts"))

    df = drop_empty_analysis_columns(df)
    missing_analises = []
    ## Data Reading
    # BUSCO
    busco_columns = [col for col in df.columns if col.startswith("Annotation_BUSCO_")]
    if len(busco_columns) > 2:
        print("\033[31mWARNING: more than 2 BUSCO columns found, showing only the first 2 columns\033[0m")
        busco_columns = busco_columns[0:2]
    busco_types = []
    for col in busco_columns:
        match = re.search(r"Annotation_BUSCO_(\w+)", col)
        if match:
            label = match.group(1)
            df[f"BUSCO_C_{label}"] = df[col].str.extract(r"C:([\d.]+)").astype(float)
            df[f"BUSCO_S_{label}"] = df[col].str.extract(r"S:([\d.]+)").astype(float)
            df[f"BUSCO_D_{label}"] = df[col].str.extract(r"D:([\d.]+)").astype(float)
            df[f"BUSCO_F_{label}"] = df[col].str.extract(r"F:([\d.]+)").astype(float)
            df[f"BUSCO_M_{label}"] = df[col].str.extract(r"M:([\d.]+)").astype(float)
            busco_types.append(label)
    #padding if missing busco
   
    missing_busco = 1
    busco_features = ["BUSCO_C_", "BUSCO_S_", "BUSCO_D_",
                      "BUSCO_F_", "BUSCO_M_"]
    while len(busco_columns) < 2:
        print("\033[31mWARNING: Number of BUSCO analysis found is less than 2, placeholders will be generated\033[0m")
        label = "BUSCO_ _{}Missing".format(missing_busco)
        for feature in busco_features:
            df[f"{feature}{label}"] = 0.0
        busco_types.append(label)
        busco_columns.append(label)
        missing_busco += 1
        missing_analises.append(label)

    # DeTEnGA
    try:
        # a class absent from DETENGA_FP% is 0%; an annotation without DETENGA_FP% stays empty (missing value)
        detenga_present = df["DETENGA_FP%"].notna()
        for detenga_class in ["PcpM0"] + TE_CLASSES:
            df[f"DETENGA_{detenga_class}"] = df["DETENGA_FP%"].str.extract(
                r"(?:^|;)\s*{}:\s*([\d.]+)".format(detenga_class))[0].astype(float).fillna(0.0).where(
                    detenga_present)
    except KeyError:
        print("\033[31mWARNING: DeTEnGA analysis not found, placeholder will be generated\033[0m")
        for detenga_class in ["PcpM0"] + TE_CLASSES:
            df[f"DETENGA_{detenga_class}"] = 0.0
        missing_analises.append("DETENGA")
    df["DETENGA_TE_related"] = sum(df[f"DETENGA_{c}"] for c in TE_CLASSES)
    df["Proteins with non-TE Pfam domains and no TE mRNA signal (DeTEnGA)"] = df["DETENGA_PcpM0"]


    # OMArk
    try:
        df["Proteins fitting their Taxonomic rank (OMARK)"] = df["OMArk Species Composition"].str.extract(r":\s*([\d.]+)%").astype(float)
        df["OMARK_S"] = df["OMArk Completeness Results"].str.extract(r"S:([\d.]+)%").astype(float)
        df["OMARK_D"] = df["OMArk Completeness Results"].str.extract(r"D:([\d.]+)%").astype(float)
        df["OMARK_M"] = df["OMArk Completeness Results"].str.extract(r"M:([\d.]+)%").astype(float)
        df[OMARK_COMPLETENESS_METRIC] = df["OMARK_S"] + df["OMARK_D"]
        consistency = pd.DataFrame([parse_omark_consistency(text) for text in df["OMArk Consistency Results"]],
                                   index=df.index)
        placed = consistency["cons"] + consistency["inco"] + consistency["cont"]
        for key in ("cons", "inco", "cont", "unkn"):
            df["OMARK_{}".format(key)] = consistency[key]
        # partial and fragmented proteins as % of the placed proteins
        df["OMARK_partial_placed"] = 100 * consistency["partial"] / placed
        df["OMARK_fragmented_placed"] = 100 * consistency["fragmented"] / placed
        df[OMARK_FRAGMENTS_METRIC] = 100 - df["OMARK_partial_placed"] - df["OMARK_fragmented_placed"]
        # not parsed in any annotation: analysis missing; in only some: missing values (checked below)
        if df[[OMARK_COMPLETENESS_METRIC, OMARK_FRAGMENTS_METRIC]].isna().all(axis=0).any():
            raise KeyError("OMArk results could not be parsed")
    except KeyError:
        print("\033[31mWARNING: Omark analysis not found, placeholder will be generated\033[0m")
        for column in ("Proteins fitting their Taxonomic rank (OMARK)", "OMARK_S", "OMARK_D", "OMARK_M",
                       OMARK_COMPLETENESS_METRIC, OMARK_FRAGMENTS_METRIC, "OMARK_cons", "OMARK_inco",
                       "OMARK_cont", "OMARK_unkn", "OMARK_partial_placed", "OMARK_fragmented_placed"):
            df[column] = 0.0
        missing_analises.append("OMARK")

    # PSAURON
    try:
        df["Transcripts Flagged as coding sequences (PSAURON)"] = df["PSAURON SCORE"]
    except KeyError:
        print("\033[31mWARNING: PSAURON analysis not found, placeholder will be generated\033[0m")
        df["Transcripts Flagged as coding sequences (PSAURON)"] = 0.0
        missing_analises.append("PSAURON")


    # AGAT
    try:
        df["No start/stop codon errors (AGAT)"] = 100 - 100 * (
            df["Models START missing"] + df["Models STOP missing"] + df["Models START & STOP missing"]
        ) / df[measured_column]
    except KeyError:
        print("\033[31mWARNING: START/STOP codons analyses not found, placeholder will be generated\033[0m")
        df["No start/stop codon errors (AGAT)"] = 0.0
        missing_analises.append("STOP")

    # Artifact (very short) introns: transcripts with no intron <= Intron_Artifact_Threshold (longest-
    # isoform score: the longest isoforms, as the "Longest isoform:" columns replace the per-transcript ones)
    intron_columns = [col for col in df.columns if SHORT_INTRON_COLUMN.match(col)]
    if len(intron_columns) > 1:
        # one column per threshold value: rows from runs with different thresholds (e.g. several
        # lineages) are combined, each row has a value in only one of them
        print("\033[33mNOTE: several potential-artifact intron thresholds in the input: {}\033[0m".format(
            ", ".join(intron_columns)))
    if intron_columns:
        short_introns = pd.to_numeric(df[intron_columns].bfill(axis=1).iloc[:, 0], errors="coerce")
        df[INTRON_METRIC] = 100 - short_introns
    if not intron_columns or df[INTRON_METRIC].isna().all():
        print("\033[31mWARNING: short intron metrics not found (older GAQET output?), placeholder will be "
              "generated and the metric is left out of the score\033[0m")
        df[INTRON_METRIC] = 0.0
        missing_analises.append("INTRONS")

    # Outer bar: % transcripts with an intron <= Intron_Floor_Threshold (introns under the floor,
    # which include the potential artifacts); longest-isoform plot: the longest isoforms
    floor_columns = [col for col in df.columns if FLOOR_INTRON_COLUMN.match(col)]
    floor_label = "Transcripts with intron\nlength shorter\nthan floor"
    if floor_columns:
        df["Transcripts with introns under the floor (%)"] = pd.to_numeric(
            df[floor_columns].bfill(axis=1).iloc[:, 0], errors="coerce")
    else:
        print("\033[31mWARNING: floor intron metrics not found (older GAQET output?), the bar is left empty"
              "\033[0m")
        df["Transcripts with introns under the floor (%)"] = 0.0
        floor_label = ""
    try:
        df["Gene_Models"] = df["Gene_Models (N)"]
        df["Transcript_Models"] = df["Transcript_Models (N)"]
        df["Both sides UTR'"] = df["Both sides UTR' (N)"]
    except KeyError:
        print("\033[31mWARNING: AGAT analysis not found, placeholder will be generated\033[0m")
        df["Gene_Models"] = 0
        df["Transcript_Models"] = 0
        df["Both sides UTR' (N)"] = 0
        df["Mean CDS Model Length (bp)"] = 0

        missing_analises.append("MODELS")

    # Homology Hits
    proteins_with_cols = [col for col in df.columns if col.startswith("ProteinsWith")]
    index = 1
    if len(proteins_with_cols) < 2:
        print("\033[31mWARNING: Number of Protein Homology analyses less than 2, placeholder will be generated\033[0m")
    while len(proteins_with_cols) < 2:
        new_col = f"MissingWithHomology_{index}"
        df[new_col] = 0
        proteins_with_cols.append(new_col)
        index += 1
        missing_analises.append(new_col)
    
    if len(proteins_with_cols) > 2:
        print("\033[31mWARNING: more than 2 Homology analyses found, showing only the first 2 columns\033[0m")
        proteins_with_cols = proteins_with_cols[0:2]

    # Axes grouped by class: completeness, precision, accuracy (one shaded sector per class)
    metrics = (
        [f"BUSCO_C_{bt}" for bt in busco_types] + [OMARK_COMPLETENESS_METRIC] +
        [
            "Proteins fitting their Taxonomic rank (OMARK)",
            "Proteins with non-TE Pfam domains and no TE mRNA signal (DeTEnGA)",
        ] +
        # axes without an outer bar of their own (homology, PSAURON, START/STOP) kept together
        proteins_with_cols +
        [
            "Transcripts Flagged as coding sequences (PSAURON)",
            "No start/stop codon errors (AGAT)",
            INTRON_METRIC,
            OMARK_FRAGMENTS_METRIC,
        ]
    )

    ## Generating Readable legends
    metrics_labels = []
    for metric in metrics:
        if metric in AXIS_LABELS:
            metrics_labels.append("" if "OMARK" in missing_analises else AXIS_LABELS[metric])
            continue
        if "BUSCO" in metric:
            if "Missing" in metric:
                taxon = ""
            else:
                taxon = metric.split("_")[-2]
                metric = f"{taxon.capitalize()}\nCompleteness\n(BUSCO)"
        if "Hits" in metric:
            metric = metric.replace("ProteinsWith", "").replace("Hits", "").replace("(%)", "") + "\nHomology Hits" + "\n(DIAMOND)"
        if "(DeTEnGA)" in metric:
            if "DETENGA" in missing_analises:
                metric = ""
            else:
                metric = "Proteins with\nnon-TE Pfam domains,\nno TE mRNA signal\n(DeTEnGA)"
        if "PSAURON" in metric:
            if "PSAURON" in missing_analises:
                metric = ""
            else:
                metric = metric.split()
                metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) + "\n" + metric[-1]
        if "OMARK" in metric:
            if "OMARK" in missing_analises:
                metric = ""
            else:
                metric = metric.split()
                metric = " ".join(metric[0:3]) + "\n" + " ".join(metric[3:-1]) + "\n(OMArk)"
        if metric == INTRON_METRIC:
            metric = "" if "INTRONS" in missing_analises else "Transcripts with no\nartifact (very short)\nintrons (AGAT)"
        elif "(AGAT)" in metric:
            if "STOP" in missing_analises:
                metric = ""
            else:
                metric = metric.split()
                metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) + "\n" + metric[-1]
        if "Missing" in metric:
            metric = ""
        metrics_labels.append(metric)

    # Metrics for columns
    columns_needed = [
        "Species",
        "DETENGA_PcpM0", "DETENGA_TE_related", "Transcripts with introns under the floor (%)",
        "Gene_Models", "Transcript_Models",
        "OMARK_S", "OMARK_D", "OMARK_M", "OMARK_cons", "OMARK_inco", "OMARK_cont", "OMARK_unkn",
        "OMARK_partial_placed", "OMARK_fragmented_placed",
        "Both sides UTR' (N)",
        "Mean CDS Model Length (bp)" 
    ] + metrics
    for bt in busco_types:
        columns_needed += [f"BUSCO_S_{bt}", f"BUSCO_D_{bt}", f"BUSCO_F_{bt}", f"BUSCO_M_{bt}"]
    df_metrics = df[columns_needed].copy()
    check_missing_values(df_metrics, arguments["fill_missing"], measured)
    if check_only:
        return
    

    ## Annotation Colors
    # tab20 first (original colours), then tab20b and tab20c; repeated if there are more annotations
    palette = (ANNOTATION_COLORS + list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) +
               list(plt.get_cmap("tab20c").colors))
    annot_colors = {species: palette[i % len(palette)] for i, species in enumerate(df["Species"].unique())}



    ## Ring Definition

    n_metrics = len(metrics)
    angle_step = 2 * np.pi / n_metrics

    # Bar enter: the bar of each metric sits in the slot before its axis
    bar_angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False)
    # Gene model descriptors (gene and transcript models, UTR, CDS length) use the free slots of the
    # axes without a bar of their own (PSAURON, START/STOP and the two homology databases)
    model_descriptors = "MODELS" not in missing_analises
    descriptor_angles = [bar_angles[metrics.index(metric)] for metric in
                         ["Transcripts Flagged as coding sequences (PSAURON)", "No start/stop codon errors (AGAT)"]
                         + proteins_with_cols[:2]]

    # Radar plot angles
    angles = (bar_angles + angle_step / 2).tolist()
    angles_closed = angles + [angles[0]]

    # 0.7 rad with 8 metrics (89% of a slot); scaled so the gap between bar groups stays the same
    bar_width = 0.7 * angle_step / (2 * np.pi / 8)

    bar_bottom = 105
    cap_height = 6
    outer_radius = 102
    inner_radius = 100


    fig, ax = plt.subplots(figsize=(26, 20), subplot_kw=dict(polar=True))
    fig.subplots_adjust(left=0.20)
    ax.set_ylim(0, 220)



    # Global Annotation Score: equal share per metric unless --weights/--exclude are given
    try:
        shares = resolve_shares(get_available_metrics(metrics, missing_analises),
                                parse_weights(arguments["weights"]), arguments["exclude"])
    except ValueError as error:
        sys.exit("\033[31mERROR: {}\033[0m".format(error))
    print_share_table(metrics, shares, missing_analises)

    df_metrics["__Score__"] = df_metrics.apply(
        lambda row: calculate_weighted_score(row, shares),
        axis=1
    )
    for error_class in CLASSES:
        df_metrics["__Class_{}__".format(error_class)] = df_metrics.apply(
            lambda row: calculate_class_score(row, shares, error_class), axis=1)
    write_score_table(df_metrics, df, metrics, shares, missing_analises, measured, measured_column, table_path,
                      weights_path)

    plot_breakdown(df_metrics, shares, breakdown_path, arguments["format"], measured)

    df_metrics = df_metrics.sort_values(by="__Score__", ascending=True).reset_index(drop=True)

    # Radial axis of the radar from its centre value to 100%: values are drawn at radius_of(value). By
    # default the centre is the lowest plotted value minus 10 (0 at least), so no value is clipped
    radial_min = arguments["radial_min"]
    if radial_min is None:
        radial_min = automatic_radial_min(df_metrics[metrics].to_numpy().ravel())

    def radius_of(value):
        return max(value - radial_min, 0) * 100 / (100 - radial_min)

    # Class sectors: shaded background behind the axes of each class, with its name and share
    totals = class_shares(shares)
    for error_class in CLASSES:
        members = [i for i, metric in enumerate(metrics) if metric_class(metric) == error_class]
        if not members:
            continue
        # from the axis of the first metric of the class to the axis of its last metric
        start, end = angles[members[0]], angles[members[-1]]
        ax.bar((start + end) / 2, 100, width=end - start, bottom=0, color=CLASS_COLORS[error_class],
               alpha=0.10, edgecolor="none", zorder=0)
        for edge in (start, end):
            ax.plot([edge, edge], [0, 100], color=CLASS_COLORS[error_class], linewidth=1.5, alpha=0.6, zorder=1)
        # thin light-gray separator in the middle of the gap before the class
        separator = bar_angles[members[0]]
        ax.plot([separator, separator], [0, 100], color="#bdbdbd", linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
        label = ax.text((start + end) / 2, 50, "{}\n{:.1f}% of GAS".format(CLASS_LABELS[error_class],
                                                                           totals.get(error_class, 0.0) * 100),
                        ha="center", va="center", fontsize=20, fontweight="bold", color=CLASS_COLORS[error_class],
                        alpha=0.85, zorder=5)
        place_class_label(ax, fig, label, start, end)

    for i, row in df_metrics.iterrows():
        values = [radius_of(value) for value in row[metrics].tolist()]
        values_closed = values + [values[0]]
        sample_color = annot_colors[row["Species"]]
        ax.plot(angles_closed, values_closed, linewidth=2.5, color=sample_color, zorder=3)
        if len(df_metrics) == 1:
            # one annotation: light fill; several: lines only, so the fills do not pile up over the classes
            ax.fill(angles_closed, values_closed, alpha=0.08, color=sample_color, zorder=2)
        offset = -((len(df_metrics)-1)/2 - i) * (bar_width / len(df_metrics))
        
        # BUSCO

        for label in busco_types:
            metric_key = f"BUSCO_C_{label}"
            angle = bar_angles[metrics.index(metric_key)]
            s = row[f"BUSCO_S_{label}"] * 0.5
            d = row[f"BUSCO_D_{label}"] * 0.5
            f = row[f"BUSCO_F_{label}"] * 0.5
            m = row[f"BUSCO_M_{label}"] * 0.5
            height = s + d + f + m
            ax.bar(angle + offset, s, width=bar_width/len(df_metrics), bottom=bar_bottom, color=DESCRIPTOR_COLORS["busco_S"], **BAR_EDGE)
            ax.bar(angle + offset, d, width=bar_width/len(df_metrics), bottom=bar_bottom + s, color=DESCRIPTOR_COLORS["busco_D"], **BAR_EDGE)
            ax.bar(angle + offset, f, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d, color=DESCRIPTOR_COLORS["busco_F"], **BAR_EDGE)
            ax.bar(angle + offset, m, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d + f, color=DESCRIPTOR_COLORS["busco_M"], **BAR_EDGE)
            tip_radius = bar_bottom + height + 3 + cap_height / 2
            ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)


        # OMARK consistency: consistent / inconsistent / contaminant / unknown proteins (% of proteome)
        angle = bar_angles[metrics.index("Proteins fitting their Taxonomic rank (OMARK)")]
        bottom = bar_bottom
        for key in ("cons", "inco", "cont", "unkn"):
            color = DESCRIPTOR_COLORS["omark_" + key]
            value = row["OMARK_{}".format(key)] * 0.5
            ax.bar(angle + offset, value, width=bar_width/len(df_metrics), bottom=bottom, color=color, **BAR_EDGE)
            bottom += value
        ax.scatter([angle + offset], [bar_bottom + 50 + 3 + cap_height / 2],
                   color=sample_color, edgecolors='white', s=400, zorder=10)

        # OMARK partial hits and fragments (% of placed proteins)
        angle = bar_angles[metrics.index(OMARK_FRAGMENTS_METRIC)]
        partial, fragmented = row["OMARK_partial_placed"] * 0.5, row["OMARK_fragmented_placed"] * 0.5
        ax.bar(angle + offset, partial, width=bar_width/len(df_metrics), bottom=bar_bottom, color=DESCRIPTOR_COLORS["omark_partial"], **BAR_EDGE)
        ax.bar(angle + offset, fragmented, width=bar_width/len(df_metrics), bottom=bar_bottom + partial,
               color=DESCRIPTOR_COLORS["omark_fragmented"], **BAR_EDGE)
        ax.scatter([angle + offset], [bar_bottom + 50 + 3 + cap_height / 2],
                   color=sample_color, edgecolors='white', s=400, zorder=10)

        # OMARK completeness: single / duplicated / missing conserved families
        angle = bar_angles[metrics.index(OMARK_COMPLETENESS_METRIC)]
        s = row["OMARK_S"] * 0.5
        d = row["OMARK_D"] * 0.5
        m = row["OMARK_M"] * 0.5
        height = s + d + m 
        ax.bar(angle + offset, s, width=bar_width/len(df_metrics), bottom=bar_bottom, color=DESCRIPTOR_COLORS["omark_S"], **BAR_EDGE)
        ax.bar(angle + offset, d, width=bar_width/len(df_metrics), bottom=bar_bottom + s, color=DESCRIPTOR_COLORS["omark_D"], **BAR_EDGE)
        ax.bar(angle + offset, m, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d, color=DESCRIPTOR_COLORS["omark_M"], **BAR_EDGE)
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
           color=sample_color, edgecolors='white',
           s=400, zorder=10)

        # DeTEnGA
        # Transposon related Transcripts: % of transcripts in any DeTEnGA class with TE evidence (TE_CLASSES)
        angle = bar_angles[metrics.index("Proteins with non-TE Pfam domains and no TE mRNA signal (DeTEnGA)")]
        te_related = row["DETENGA_TE_related"] * 0.5   # same scale as the other % bars (100% = full bar)
        height = 50
        ax.bar(angle + offset, te_related, width=bar_width/len(df_metrics), bottom=bar_bottom, color=DESCRIPTOR_COLORS["te_related"], **BAR_EDGE)
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
           color=sample_color, edgecolors='white',
           s=400, zorder=10)

        # Gene model descriptors: counts and lengths relative to the largest
        # annotation of the plot (largest = full bar), UTR as % of transcripts (100% = full bar)
        if model_descriptors:
            descriptor_values = {}
            for key, column in (("gene_models", "Gene_Models"), ("transcript_models", "Transcript_Models"),
                                ("cds_length", "Mean CDS Model Length (bp)")):
                largest = df_metrics[column].max()
                descriptor_values[key] = 0 if largest == 0 else row[column] * 50 / largest
            descriptor_values["utr"] = (0 if row["Transcript_Models"] == 0 else
                                        row["Both sides UTR' (N)"] * 50 / row["Transcript_Models"])
            for angle, (key, _) in zip(descriptor_angles, MODEL_DESCRIPTORS):
                ax.bar(angle + offset, descriptor_values[key], width=bar_width/len(df_metrics), bottom=bar_bottom,
                       color=DESCRIPTOR_COLORS[key], **BAR_EDGE)
                ax.scatter([angle + offset], [bar_bottom + 50 + 3 + cap_height / 2],
                           color=sample_color, edgecolors='white', s=400, zorder=10)

        # Transcripts with introns <= Intron_Floor_Threshold (%), 100% = full bar
        height = 50
        angle = bar_angles[metrics.index(INTRON_METRIC)]
        ax.bar(angle + offset, row["Transcripts with introns under the floor (%)"] * 0.5, width=bar_width/len(df_metrics),
               bottom=bar_bottom, color=DESCRIPTOR_COLORS["floor_introns"], **BAR_EDGE)
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)


    ## spider plot legend
    ring = Wedge((0, 0), outer_radius, 0, 360, width=outer_radius - inner_radius,
                transform=ax.transData._b, color='black', alpha=0.3)
    ax.add_artist(ring)
    span = 100 - radial_min
    tick_step = 20 if span > 60 else (10 if span >= 25 else (5 if span >= 12 else 2))
    # round levels above the centre, skipping one too close to the centre label
    levels = [lvl for lvl in range(100, -1, -tick_step) if lvl - radial_min >= tick_step / 2][::-1]
    if radial_min > 0:
        levels = [radial_min] + levels          # value at the centre, first label of the scale
    ax.set_yticks([radius_of(lvl) for lvl in levels])
    ax.set_yticklabels([f"{lvl:g}%" for lvl in levels], fontsize=18)
    ax.set_rlabel_position(np.degrees(bar_angles[0]))   # tick labels on the 0° line
    ax.yaxis.grid(True, linestyle='dotted', linewidth=0.8, color='gray')
    ax.set_xticks(angles)
    # Axis labels: metric name in bold, the tool that measures it (last "(...)") in regular weight below
    names, tools = [], []
    for text in metrics_labels:
        match = re.match(r"^(.*?)\s*(\([^()]*\))\s*$", text, flags=re.S)
        names.append(match.group(1).strip() if match else text)
        tools.append(match.group(2) if match else "")
    ax.set_xticklabels(names, fontsize=26, fontweight="bold", color="#222222")
    ax.tick_params(axis="x", pad=110)   # axis labels out of the way of the ring labels
    axis_label_groups = []      # (angle, tick label, tool annotation or None)
    for angle, tick_label, tool in zip(angles, ax.get_xticklabels(), tools):
        tool_text = None
        if tool and tick_label.get_text():
            tool_text = ax.annotate(tool, xy=(0.5, 0), xycoords=tick_label, xytext=(0, -2),
                                    textcoords="offset points", ha="center", va="top", fontsize=24,
                                    color="#444444", annotation_clip=False)
        axis_label_groups.append((angle, tick_label, tool_text))
    ax.spines['polar'].set_visible(False)

    radio_texto = bar_bottom + 50 + 3 + cap_height / 2 + 8  # ring labels just beyond the annotation dots

    # Label of the outer bar of each slot (the bars are drawn at the slot of the metric they sit next to)
    slot_labels = {}
    for taxon in busco_types:
        if "Missing" not in taxon:
            slot_labels[metrics.index(f"BUSCO_C_{taxon}")] = f"{taxon.split('_')[0].capitalize()}\nBUSCO"
    if "OMARK" not in missing_analises:
        slot_labels[metrics.index(OMARK_COMPLETENESS_METRIC)] = "OMArk\ncompleteness"
        slot_labels[metrics.index("Proteins fitting their Taxonomic rank (OMARK)")] = "OMArk\nconsistency"
        slot_labels[metrics.index(OMARK_FRAGMENTS_METRIC)] = "OMArk partial /\nfragmented"
    if "MODELS" not in missing_analises:
        if "DETENGA" not in missing_analises:
            slot_labels[metrics.index("Proteins with non-TE Pfam domains and no TE mRNA signal (DeTEnGA)")] = \
                "Transposon related\nTranscripts"
    slot_labels[metrics.index(INTRON_METRIC)] = floor_label
    col_names = [slot_labels.get(i, "") for i in range(len(metrics))]
    ring_angles = list(bar_angles)
    if model_descriptors:
        col_names += [label for _, label in MODEL_DESCRIPTORS]
        ring_angles += list(descriptor_angles)
    # labels moved (in points) to avoid colliding with the metric labels
    label_shift = {}
    ring_texts = []
    for angulo, etiqueta in zip(ring_angles, col_names):
        # anchored just beyond the annotation dots and growing outwards (away from the dots and bars)
        cos, sin = np.cos(angulo), np.sin(angulo)
        ring_texts.append(ax.annotate(
            etiqueta,
            xy=(angulo, radio_texto),
            xytext=label_shift.get(etiqueta, (0, 0)),
            textcoords="offset points",
            fontsize=22,
            ha="left" if cos > 0.35 else ("right" if cos < -0.35 else "center"),
            va="bottom" if sin > 0.35 else ("top" if sin < -0.35 else "center"),
            multialignment="center",
            color='black',
            annotation_clip=False,
            zorder=20
    ))


    # Columns legend
    # Descriptor legend grouped by tool (not by class: the outer bars are not part of any class score)
    legend_groups = [
        ("BUSCO", [("busco_S", "Complete, Single (%)"), ("busco_D", "Complete, Duplicated (%)"),
                   ("busco_F", "Fragmented (%)"), ("busco_M", "Missing (%)")]),
        ("OMArk completeness", [("omark_S", "Single (%)"), ("omark_D", "Duplicated (%)"),
                                ("omark_M", "Missing (%)")]),
        ("OMArk consistency", [("omark_cons", "Consistent (%)"), ("omark_inco", "Inconsistent (%)"),
                               ("omark_cont", "Contaminant (%)"), ("omark_unkn", "Unknown (%)")]),
        ("OMArk fragments", [("omark_partial", "Partial hits (% of placed proteins)"),
                             ("omark_fragmented", "Fragmented (% of placed proteins)")]),
        ("DeTEnGA", [("te_related", "Transposon related Transcripts (%)")]),
        ("AGAT", [("floor_introns", "Transcripts with introns under the floor (%)"),
                  ("gene_models", "Gene models (N, relative to max)"),
                  ("transcript_models", "Transcript models (N, relative to max)"),
                  ("utr", "Both sides UTR' (%)"),
                  ("cds_length", "Mean CDS Model Length (bp, relative to max)")]),
    ]
    metric_handles, header_labels = [], set()
    for header, items in legend_groups:
        metric_handles.append(Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="none", label=header))
        header_labels.add(header)
        metric_handles += [Rectangle((0, 0), 1, 1, facecolor=DESCRIPTOR_COLORS[key], edgecolor="black", label=label)
                           for key, label in items]
    # Reuse already computed scores from df_metrics
    area_map = dict(zip(df_metrics["Species"], df_metrics["__Score__"]))
    areas = df_metrics["__Score__"].tolist()

    def legend_scores(species):
        """GAS and class sub-scores for the legend: GAS: xx (Completeness: xx, Precision: xx, Accuracy: xx)."""
        row = df_metrics.loc[df_metrics["Species"] == species].iloc[0]
        parts = []
        for error_class, name in zip(CLASSES, ("Completeness", "Precision", "Accuracy")):
            value = row["__Class_{}__".format(error_class)]
            parts.append("{}: {}".format(name, "NA" if pd.isna(value) else "{:.1f}".format(value)))
        return "GAS: {:.1f} ({})".format(area_map[species], ", ".join(parts))

    # Anot legends
    reverse_df_metrics = df_metrics.sort_values(by="__Score__", ascending=False).reset_index(drop=True)

    annotation_handles = [
    plt.Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=annot_colors[species],
               markeredgecolor='black',
               markersize=15,
               label=f"{species} - {legend_scores(species)}")
    for species in reverse_df_metrics["Species"]]
    # descriptors legend from the top, annotations legend below it (no overlap whatever their lengths)
    legend1 = ax.legend(handles=metric_handles, title="Descriptors (outer bars)\nNOT used in the GAS",
                        loc='upper left', bbox_to_anchor=(1.20, 1.05), fontsize=26, frameon=False,
                        title_fontproperties={'size': 30, 'weight': 'bold'}, alignment='left')
    for text in legend1.get_texts():
        if text.get_text() in header_labels:
            text.set_fontweight("bold")
            text.set_color("#333333")
    ax.add_artist(legend1)
    fig.canvas.draw()
    box = legend1.get_window_extent().transformed(ax.transAxes.inverted())
    legend2 = ax.legend(handles=annotation_handles, title="Annotations", loc='upper left',
                        bbox_to_anchor=(1.20, box.y0 - 0.03), fontsize=26, frameon=False,
                        title_fontproperties={'size': 30, 'weight': 'bold'}, alignment='left')
    ax.add_artist(legend2)


    plt.tight_layout()
    # after the final layout: axis labels overlapping an outer-bar label are moved outwards
    push_axis_labels_out(fig, axis_label_groups, ring_texts)
    title = fig.suptitle("Model-quality metrics measured on: {}".format(
        "longest isoform of each gene" if measured == "longest" else "all transcripts"),
        fontsize=34, fontweight="bold", y=1.04)
    plt.savefig(plot_path,
        format=arguments["format"],
        dpi=300,
        bbox_inches='tight',
        bbox_extra_artists=[legend1, legend2, title],
        pad_inches=0.5
    )
    plt.close(fig)
    print("Plot: {}\nBreakdown plot: {}\nScore table: {}\nWeights table: {}".format(
        plot_path, breakdown_path, table_path, weights_path))

if __name__ == "__main__":
    main()