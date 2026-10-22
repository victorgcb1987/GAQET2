#!/usr/bin/env python

import argparse
import numpy as np
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
    
    help_outbase = "(Required) GAQET Plot output"
    parser.add_argument("--output", "-o", type=str,
                        help=help_outbase, default="")
    
    help_weights = ("Share (percentage) of the Global Annotation Score given to a metric, as key=percentage "
                    "separated by commas (e.g. busco=40,omark=10). A role key (homology, busco, omark, "
                    "psauron, detenga, agat) sets the total share of that analysis, split equally among its "
                    "metrics. Any other key is matched as a substring of one metric name (e.g. swissprot=15, "
                    "embryophyta=20). The remaining percentage is split equally among the metrics not "
                    "mentioned. Default: all metrics have the same share")
    parser.add_argument("--weights", "-w", type=str,
                        help=help_weights, default="")

    help_exclude = ("Comma separated keys (same matching as --weights) of metrics to leave out of "
                    "the Global Annotation Score. They are still drawn in the plot")
    parser.add_argument("--exclude", "-e", type=str,
                        help=help_exclude, default="")

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
            "output": Path(parser.output),
            "weights": parser.weights,
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

ROLES = ("homology", "busco", "omark", "psauron", "detenga", "agat")

# Metrics of these roles are left out of the score when their analysis is missing
MISSING_FLAGS = {"omark": "OMARK", "psauron": "PSAURON", "detenga": "DETENGA", "agat": "STOP"}


def metric_role(metric):
    """Return the analysis (role) a metric column belongs to."""
    if metric.startswith(("ProteinsWith", "MissingWithHomology")):
        return "homology"
    if "BUSCO" in metric:
        return "busco"
    for role, tag in (("omark", "(OMARK)"), ("psauron", "(PSAURON)"),
                      ("detenga", "(DeTEnGA)"), ("agat", "(AGAT)")):
        if tag in metric:
            return role
    return None


def normalize_key(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


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
    """Metrics selected by a role name (busco, omark...) or by a substring of a metric name."""
    nkey = normalize_key(key)
    if nkey in ROLES:
        matched = [m for m in metrics if metric_role(m) == nkey]
        is_role = True
    else:
        matched = [m for m in metrics if nkey and nkey in normalize_key(m)]
        is_role = False
    if not matched:
        raise ValueError("'{}' does not match any metric used in this run. Roles: {}. "
                         "Metrics: {}".format(key, ", ".join(ROLES), available_text))
    return matched, is_role


def resolve_shares(metrics, share_spec=None, exclude=None):
    """Share of the global annotation score (fraction, adding up to 1) of each scored metric.

    GAS = sum(share_i * metric_i). By default every metric has the same share.
    share_spec: [(key, percentage)]. A role key (busco, omark, psauron, detenga, agat, homology)
    sets the total share of that analysis, split equally among its metrics. Any other key is
    matched as a substring of one metric name and sets the share of that metric alone.
    The remaining percentage is split equally among the metrics not mentioned.
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
    for key, percentage in share_spec or []:
        matched, is_role = match_metrics(key, scored, "; ".join(scored))
        if not is_role and len(matched) > 1:
            raise ValueError("'{}' matches several metrics ({}). Use a more specific key "
                             "or a role name".format(key, "; ".join(matched)))
        for metric in matched:
            if metric in shares:
                raise ValueError("Metric '{}' receives more than one share".format(metric))
            shares[metric] = percentage / len(matched)

    specified = sum(shares.values())
    remaining = [m for m in scored if m not in shares]
    if specified > 100 + 1e-9:
        raise ValueError("The shares given add up to {:.2f}%, more than 100%".format(specified))
    if remaining:
        if specified >= 100 - 1e-9:
            raise ValueError("The shares given use the whole score, no share left for: {}".format("; ".join(remaining)))
        for metric in remaining:
            shares[metric] = (100 - specified) / len(remaining)
    elif abs(specified - 100) > 1e-6:
        raise ValueError("All metrics have a share, so they must add up to 100% (now {:.2f}%)".format(specified))
    return {metric: shares[metric] / 100 for metric in scored}


def write_share_table(metrics, shares, missing_analises, out_path):
    """Print and save the share of the score of every metric."""
    available = set(get_available_metrics(metrics, missing_analises))
    lines = ["Metric\tRole\tShare_%\tStatus"]
    for metric in metrics:
        if metric in shares:
            status, share = "scored", "{:.2f}".format(shares[metric] * 100)
        elif metric in available:
            status, share = "excluded", "0.00"
        else:
            status, share = "analysis missing", "0.00"
        lines.append("{}\t{}\t{}\t{}".format(metric.replace("\n", " "), metric_role(metric), share, status))
    print("Global Annotation Score, share of each metric:\n" + "\n".join(lines))
    if out_path:
        with open(out_path, "w") as fhand:
            fhand.write("\n".join(lines) + "\n")


def calculate_weighted_score(row, shares):
    """Global Annotation Score (0-100): sum of share * metric value over the scored metrics."""
    return sum(row[metric] * share for metric, share in shares.items())


def main():
    sys.tracebacklimit = 1
    if '--version' in sys.argv or "-v" in sys.argv:
        print(VERSION)
        sys.exit(0)
    arguments = get_arguments()
    df = pd.read_csv(arguments["input"], sep="\t")

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
        df["DETENGA_PcpM0"] = df["DETENGA_FP%"].str.extract(r"PcpM0:\s*([\d.]+)").astype(float)
        df["DETENGA_PteMte"] = df["DETENGA_FP%"].str.extract(r"PteMte:\s*([\d.]+)").astype(float)
    except KeyError:
        print("\033[31mWARNING: DeTEnGA analysis not found, placeholder will be generated\033[0m")
        df["DETENGA_PcpM0"] = 0.0
        df["DETENGA_PteMte"] = 0.0
        missing_analises.append("DETENGA")
    df["Proteins with functional domains (DeTEnGA)"] = df["DETENGA_PcpM0"]


    # OMArk
    try:
        df["Proteins fitting their Taxonomic rank (OMARK)"] = df["OMArk Species Composition"].str.extract(r":\s*([\d.]+)%").astype(float)
        df["OMARK_S"] = df["OMArk Completeness Results"].str.extract(r"S:([\d.]+)%").astype(float)
        df["OMARK_D"] = df["OMArk Completeness Results"].str.extract(r"D:([\d.]+)%").astype(float)
        df["OMARK_M"] = df["OMArk Completeness Results"].str.extract(r"M:([\d.]+)%").astype(float)
    except KeyError:
        print("\033[31mWARNING: Omark analysis not found, placeholder will be generated\033[0m")
        df["Proteins fitting their Taxonomic rank (OMARK)"] = 0.0
        df["OMARK_S"] = 0.0
        df["OMARK_D"] = 0.0
        df["OMARK_M"] = 0.0
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
        ) / df["Transcript_Models (N)"]
    except KeyError:
        print("\033[31mWARNING: START/STOP codons analyses not found, placeholder will be generated\033[0m")
        df["No start/stop codon errors (AGAT)"] = 0.0
        missing_analises.append("STOP")
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

    metrics = (
        proteins_with_cols +
        ["Proteins fitting their Taxonomic rank (OMARK)"] +
        [f"BUSCO_C_{bt}" for bt in busco_types] + 
        [
            "Proteins with functional domains (DeTEnGA)",
            "Transcripts Flagged as coding sequences (PSAURON)",
            "No start/stop codon errors (AGAT)",
            
        ]
    )

    ## Generating Readable legends
    metrics_labels = []
    for metric in metrics:
        if "BUSCO" in metric:
            if "Missing" in metric:
                taxon = ""
            else:
                taxon = metric.split("_")[-2]
                metric = f"{taxon.capitalize()} \nCompletness \n(BUSCO)"
        if "Hits" in metric:
            metric = metric.replace("ProteinsWith", "").replace("Hits", "").replace("(%)", "") + "\nHomology Hits" + "\n(DIAMOND)"
        if "(DeTEnGA)" in metric:
            if "DETENGA" in missing_analises:
                metric = ""
            else:
                metric = metric.split()
                metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) +"\n" + metric[-1]
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
                metric = " ".join(metric[0:3]) + "\n" + " ".join(metric[3:-1]) + "\n" + metric[-1]
        if "AGAT" in metric:
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
        "DETENGA_PcpM0", "DETENGA_PteMte",
        "Gene_Models", "Transcript_Models",
        "OMARK_S", "OMARK_D", "OMARK_M",
        "Both sides UTR' (N)",
        "Mean CDS Model Length (bp)" 
    ] + metrics
    for bt in busco_types:
        columns_needed += [f"BUSCO_S_{bt}", f"BUSCO_D_{bt}", f"BUSCO_F_{bt}", f"BUSCO_M_{bt}"]
    df_metrics = df[columns_needed].copy()
    

    ## Annotation Colors
    annot_colors = {
        species: color for species, color in zip(df["Species"].unique(),
            plt.get_cmap("tab20").colors[:len(df["Species"].unique())])
    }



    ## Ring Definition

    n_metrics = len(metrics)
    angle_step = 2 * np.pi / n_metrics

    # Bar enter
    bar_angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False)

    # Radar plot angles
    angles = (bar_angles + angle_step / 2).tolist()
    angles_closed = angles + [angles[0]]

    bar_width = 0.7

    bar_bottom = 105
    cap_height = 6
    outer_radius = 102
    inner_radius = 100


    fig, ax = plt.subplots(figsize=(26, 20), subplot_kw=dict(polar=True))
    fig.subplots_adjust(left=0.20)
    ax.set_ylim(0, 220)


    # Calcula el desplazamiento angular entre columnas
    angle_step = 2 * np.pi / len(bar_angles)

    # Global Annotation Score: equal share per metric unless --weights/--exclude are given
    try:
        shares = resolve_shares(get_available_metrics(metrics, missing_analises),
                                parse_weights(arguments["weights"]), arguments["exclude"])
    except ValueError as error:
        sys.exit("\033[31mERROR: {}\033[0m".format(error))
    share_table = ""
    if str(arguments["output"]):
        share_table = arguments["output"].with_name(arguments["output"].stem + "_GAS_weights.tsv")
    write_share_table(metrics, shares, missing_analises, share_table)

    df_metrics["__Score__"] = df_metrics.apply(
        lambda row: calculate_weighted_score(row, shares),
        axis=1
    )

    df_metrics = df_metrics.sort_values(by="__Score__", ascending=True).reset_index(drop=True)

    for i, row in df_metrics.iterrows():
        values = row[metrics].tolist()
        values_closed = values + [values[0]]
        area = area_polar_polygon(values_closed, angles_closed)
        sample_color = annot_colors[row["Species"]]
        ax.plot(angles_closed, values_closed, linewidth=2, color=sample_color)
        ax.fill(angles_closed, values_closed, alpha=0.05, color=sample_color)
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
            ax.bar(angle + offset, s, width=bar_width/len(df_metrics), bottom=bar_bottom, color="#4481cc", edgecolor='black')
            ax.bar(angle + offset, d, width=bar_width/len(df_metrics), bottom=bar_bottom + s, color="#91b8e8", edgecolor='black')
            ax.bar(angle + offset, f, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d, color="#c9d380", edgecolor='black')
            ax.bar(angle + offset, m, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d + f, color="#c0564b", edgecolor='black')
            tip_radius = bar_bottom + height + 3 + cap_height / 2
            ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)


        # OMARK
        angle = bar_angles[metrics.index("Proteins fitting their Taxonomic rank (OMARK)")]
        s = row["OMARK_S"] * 0.5
        d = row["OMARK_D"] * 0.5
        m = row["OMARK_M"] * 0.5
        height = s + d + m 
        ax.bar(angle + offset, s, width=bar_width/len(df_metrics), bottom=bar_bottom, color="#47CE0D", edgecolor='black')
        ax.bar(angle + offset, d, width=bar_width/len(df_metrics), bottom=bar_bottom + s, color="#4F9B53D7", edgecolor='black')
        ax.bar(angle + offset, m, width=bar_width/len(df_metrics), bottom=bar_bottom + s + d, color="#27C6A9", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
           color=sample_color, edgecolors='white',
           s=400, zorder=10)

        # DeTEnGA
        angle = bar_angles[metrics.index("Proteins with functional domains (DeTEnGA)")]
        ptemte = row["DETENGA_PteMte"]
        height = 50
        ax.bar(angle + offset, ptemte, width=bar_width/len(df_metrics), bottom=bar_bottom, color="#eba434", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
           color=sample_color, edgecolors='white',
           s=400, zorder=10)

        # AGAT
        angle = bar_angles[metrics.index("Transcripts Flagged as coding sequences (PSAURON)")]
        max_gene = df_metrics["Gene_Models"].max()
        if max_gene == 0:
            gene = 0
        else:
            gene = row["Gene_Models"] * 50 / max_gene
        height = 50
        ax.bar(angle + offset, gene, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#2d99ba", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)


        angle = bar_angles[metrics.index("No start/stop codon errors (AGAT)")]
        max_transcripts = df_metrics["Transcript_Models"].max()
        if max_transcripts == 0:
            transcripts = 0
            
        else:
            transcripts = row["Transcript_Models"] * 50 / max_transcripts
        ax.bar(angle + offset, transcripts, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#097c77", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)

        angle = bar_angles[0]
        if row["Transcript_Models"] == 0:
            utr = 0
        else:
            utr = row["Both sides UTR' (N)"] * 50 / row["Transcript_Models"]
        ax.bar(angle + offset, utr, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#D21F22", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)

        angle = bar_angles[1]
        max_cds_value = df_metrics["Mean CDS Model Length (bp)"].max()
        if max_cds_value == 0:
            cds_length = 0
        else:
            cds_length = row["Mean CDS Model Length (bp)"] * 50 / max_cds_value
        ax.bar(angle + offset, cds_length, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#BDD21F", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)


    ## spider plot legend
    ring = Wedge((0, 0), outer_radius, 0, 360, width=outer_radius - inner_radius,
                transform=ax.transData._b, color='black', alpha=0.3)
    ax.add_artist(ring)
    ax.set_yticks([20, 40, 60, 80, 100])
    ax.set_yticklabels([f"{lvl}%" for lvl in [20, 40, 60, 80, 100]], fontsize=18)
    ax.yaxis.grid(True, linestyle='dotted', linewidth=0.8, color='gray')
    ax.set_xticks(angles)
    ax.set_xticklabels(metrics_labels, fontsize=30, fontweight="bold")
    ax.spines['polar'].set_visible(False)

    radio_texto = bar_bottom + 100  # distancia radial por encima de las columnas

    col_names = []
    if "MODELS" not in missing_analises:
        col_names = ["Transcripts with\nboth sides UTR", "Average CDS\nlength"]
        tail = ["Transposons", "Gene Models", "Transcripts models"]
    else:
        col_names = ["",""]
        tail = ["","",""]
    if "OMARK" not in missing_analises:
        col_names += ["OMARK"]
    else:
        col_names += [""]
    for taxon in busco_types:
        if "Missing" not in taxon:
            col_names += [f"{taxon.split('_')[0].capitalize()}\nBUSCO"]
        else:
            col_names += [""]
    if "DETENGA" in missing_analises:
        tail[0] = ""
    col_names += tail
    #col_names = ["Transcripts with\nboth sides UTR", "Average CDS\nlength", "OMArk"] + [f"{taxon.split('_')[0].capitalize()}\nBUSCO" if "Missing" not in taxon else "" for taxon in busco_types] + ["Transposons", "Gene Models", "Transcripts models"]
    for angulo, etiqueta in zip(bar_angles, col_names):
        ax.text(
            angulo,                           
            radio_texto,                      
            etiqueta,      
            fontsize=26,
            ha='center',
            va='center',
            rotation=0,
            rotation_mode=None,
            color='black',
            clip_on=False,
            zorder=20
    )


    # Columns legend
    metric_handles = [
        Rectangle((0,0),1,1,facecolor="#4481cc", edgecolor='black', label='BUSCO Complete, Single (%)'),
        Rectangle((0,0),1,1,facecolor="#91b8e8", edgecolor='black', label='BUSCO Complete, Duplicated (%)'),
        Rectangle((0,0),1,1,facecolor="#c9d380", edgecolor='black', label='BUSCO Fragmented (%)'),
        Rectangle((0,0),1,1,facecolor="#c0564b", edgecolor='black', label='BUSCO Missing (%)'),
        Rectangle((0,0),1,1,facecolor="#e29c2d", edgecolor='black', label='DeTEnGA PteMte (%)'),
        Rectangle((0,0),1,1,facecolor='#47CE0D', edgecolor='black', label='OMARK Single (%)'),
        Rectangle((0,0),1,1,facecolor="#4F9B53D7", edgecolor='black', label='OMARK Duplicated (%)'),
        Rectangle((0,0),1,1,facecolor='#27C6A9', edgecolor='black', label='OMARK Missing (%)'),
        Rectangle((0,0),1,1,facecolor="#2d99ba", edgecolor='black', label='AGAT Gene models (N)'),
        Rectangle((0,0),1,1,facecolor="#097c77", edgecolor='black', label='AGAT Transcript models (N)'),
        Rectangle((0,0),1,1,facecolor="#D21F22", edgecolor='black', label="AGAT Both sides UTR' (%)"),
        Rectangle((0,0),1,1,facecolor="#BDD21F", edgecolor='black', label="AGAT Mean CDS Model Length (bp)")
    ]
    # Reuse already computed scores from df_metrics
    area_map = dict(zip(df_metrics["Species"], df_metrics["__Score__"]))
    areas = df_metrics["__Score__"].tolist()

    # Anot legends
    reverse_df_metrics = df_metrics.sort_values(by="__Score__", ascending=False).reset_index(drop=True)

    annotation_handles = [
    plt.Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=annot_colors[species],
               markeredgecolor='black',
               markersize=15,
               label=f"{species} (Score: {area_map[species]:.1f})")
    for species in reverse_df_metrics["Species"]]
    legend1 = ax.legend(handles=metric_handles, title="Metrics", loc='center left',
                        bbox_to_anchor=(1.20, 0.75), fontsize=28, frameon=False,
                        title_fontproperties={'size': 32, 'weight': 'bold'}, alignment='center')
    legend2 = ax.legend(handles=annotation_handles, title="Annotations", loc='center left',
                        bbox_to_anchor=(1.20, 0.25), fontsize=28, frameon=False,
                        title_fontproperties={'size': 32, 'weight': 'bold'}, alignment='center')
    ax.add_artist(legend1)
    ax.add_artist(legend2)


    # Circle proportions
    area_total = sum(areas)
    area_fractions = [a / area_total for a in areas]

    # Outer ring radios
    ring_inner = 95  # puedes ajustar
    ring_outer = 100
    start_angle = 0

    for i, (frac, row, color) in enumerate(zip(area_fractions, df_metrics.itertuples(), df_metrics["Species"].map(annot_colors))):
        angle_span = frac * 360
        theta_start = start_angle
        theta_end = start_angle + angle_span

        # Actual drawing
        wedge = Wedge(
            center=(0, 0),
            r=ring_outer,
            theta1=theta_start,
            theta2=theta_end,
            width=ring_outer - ring_inner,
            transform=ax.transData._b,
            color=color,
            alpha=0.8,
            linewidth=0
        )
        ax.add_artist(wedge)
        start_angle += angle_span


    plt.tight_layout()
    plt.savefig(arguments["output"],
        dpi=300,
        bbox_inches='tight',
        bbox_extra_artists=[legend1, legend2],
        pad_inches=0.5
    )

if __name__ == "__main__":
    main()