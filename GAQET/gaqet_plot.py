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


VERSION = "v1.3.2"

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
            "output": Path(parser.output)}


def area_polar_polygon(r_values, theta_values):
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

def main():
    sys.tracebacklimit = 1
    if '--version' in sys.argv or "-v" in sys.argv:
        print(VERSION)
        sys.exit(0)
    arguments = get_arguments()
    df = pd.read_csv(arguments["input"], sep="\t")

    ## Data Reading
    # BUSCO
    busco_columns = [col for col in df.columns if col.startswith("Annotation_BUSCO_")]
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

    # DeTEnGA
    df["DETENGA_PcpM0"] = df["DETENGA_FP%"].str.extract(r"PcpM0:\s*([\d.]+)").astype(float)
    df["DETENGA_PteMte"] = df["DETENGA_FP%"].str.extract(r"PteMte:\s*([\d.]+)").astype(float)
    df["Proteins with functional domains (DeTEnGA)"] = df["DETENGA_PcpM0"]


    # OMArk
    df["Proteins fitting their Taxonomic rank (OMARK)"] = df["OMArk Species Composition"].str.extract(r":\s*([\d.]+)%").astype(float)
    df["OMARK_S"] = df["OMArk Completeness Results"].str.extract(r"S:([\d.]+)%").astype(float)
    df["OMARK_D"] = df["OMArk Completeness Results"].str.extract(r"D:([\d.]+)%").astype(float)
    df["OMARK_M"] = df["OMArk Completeness Results"].str.extract(r"M:([\d.]+)%").astype(float)

    # PSAURON
    df["Transcripts Flagged as coding sequences (PSAURON)"] = df["PSAURON SCORE"]

    # AGAT
    df["No start/stop codon errors (AGAT)"] = 100 - 100 * (
        df["Models START missing"] + df["Models STOP missing"] + df["Models START & STOP missing"]
    ) / df["Transcript_Models (N)"]
    df["Gene_Models"] = df["Gene_Models (N)"]
    df["Transcript_Models"] = df["Transcript_Models (N)"]
    df["Both sides UTR'"] = df["Both sides UTR' (N)"]

    # Homology Hits
    proteins_with_cols = [col for col in df.columns if col.startswith("ProteinsWith")]
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
            taxon = metric.split("_")[-2]
            metric = f"{taxon.capitalize()} \nCompletness \n(BUSCO)"
        if "Hits" in metric:
            metric = metric.replace("ProteinsWith", "").replace("Hits", "").replace("(%)", "") + "\nHomology Hits" + "\n(DIAMOND)"
        if "(DeTEnGA)" in metric:
            metric = metric.split()
            metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) +"\n" + metric[-1]
        if "PSAURON" in metric:
            metric = metric.split()
            metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) + "\n" + metric[-1]
        if "OMARK" in metric:
            metric = metric.split()
            metric = " ".join(metric[0:3]) + "\n" + " ".join(metric[3:-1]) + "\n" + metric[-1]
        if "AGAT" in metric:
            metric = metric.split()
            metric = " ".join(metric[0:2]) + "\n" + " ".join(metric[2:-1]) + "\n" + metric[-1]
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

    # Ángulos de las columnas (centro de cada barra)
    bar_angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False)

    # Ángulos del radarplot y etiquetas (entre columnas)
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

    separator_width = 0.01
    separator_color = 'black'
    n_columns = len(df_metrics)

    # Calcula el desplazamiento angular entre columnas
    angle_step = 2 * np.pi / len(bar_angles)
    half_step = angle_step / 2

    # for angle in bar_angles:
    #     # Desplaza la barra separadora a la mitad entre columnas
    #     sep_angle = angle + half_step
    #     ax.bar(sep_angle, ax.get_ylim()[1], width=separator_width, bottom=bar_bottom,
    #         color=separator_color, alpha=0.4, zorder=0)

    ## Adding columns
    # Sort columns by area
    df_metrics["__Area__"] = df_metrics.apply(
    lambda row: area_polar_polygon(row[metrics].tolist() + [row[metrics].tolist()[0]], angles_closed),
    axis=1
    )


    df_metrics = df_metrics.sort_values(by="__Area__", ascending=True).reset_index(drop=True)

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
        transcripts = row["Transcript_Models"] * 50 / max_transcripts
        ax.bar(angle + offset, transcripts, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#097c77", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)

        angle = bar_angles[0]
        max_utr = df_metrics["Both sides UTR' (N)"].max()
        utr = row["Both sides UTR' (N)"] * 50 / max_utr
        ax.bar(angle + offset, utr, width=bar_width/len(df_metrics), bottom=bar_bottom,
        color="#D21F22", edgecolor='black')
        tip_radius = bar_bottom + height + 3 + cap_height / 2
        ax.scatter([angle + offset], [tip_radius],
                        color=sample_color, edgecolors='white',
                        s=400, zorder=10)

        angle = bar_angles[1]
        max_cds_value = df_metrics["Mean CDS Model Length (bp)"].max()
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
    ax.set_yticklabels([f"{lvl}%" for lvl in [20, 40, 60, 80, 100]], fontsize=10)
    ax.yaxis.grid(True, linestyle='dotted', linewidth=0.8, color='gray')
    ax.set_xticks(angles)
    ax.set_xticklabels(metrics_labels, fontsize=24, fontweight="bold")
    ax.spines['polar'].set_visible(False)

    radio_texto = bar_bottom + 100  # distancia radial por encima de las columnas

    
    col_names = ["UTR", "CDS", "OMArk"] + [f"{taxon.split('_')[0].capitalize()}\nBUSCO" for taxon in busco_types] + ["Transposons", "Gene Models", "Transcripts models"]
    for angulo, etiqueta in zip(bar_angles, col_names):
        ax.text(
            angulo,                           # posición angular (alineado con cada columna)
            radio_texto,                      # radio: más alejado del centro
            etiqueta,       # solo la primera línea del texto
            fontsize=20,
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
        Rectangle((0,0),1,1,facecolor="#D21F22", edgecolor='black', label="AGAT Both sides UTR' (N)"),
        Rectangle((0,0),1,1,facecolor="#BDD21F", edgecolor='black', label="AGAT Mean CDS Model Length (bp)")
    ]
    area_map = {}
    for i, row in df_metrics.iterrows():
        values = row[metrics].tolist()
        values_closed = values + [values[0]]
        area = area_polar_polygon(values_closed, angles_closed)
        area_map[row["Species"]] = area

    # Anot legends
    reverse_df_metrics = df_metrics.sort_values(by="__Area__", ascending=False).reset_index(drop=True)

    annotation_handles = [
    plt.Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=annot_colors[species],
               markeredgecolor='black',
               markersize=15,
               label=f"{species} (Area: {area_map[species]:.0f})")
    for species in reverse_df_metrics["Species"]]
    legend1 = ax.legend(handles=metric_handles, title="Metrics", loc='center left',
                        bbox_to_anchor=(1.20, 0.75), fontsize=24, title_fontsize=26, frameon=False)
    legend2 = ax.legend(handles=annotation_handles, title="Annotations", loc='center left',
                        bbox_to_anchor=(1.20, 0.25), fontsize=24, title_fontsize=26, frameon=False)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ## Area Draw 
    areas = []
    for _, row in df_metrics.iterrows():
        values = row[metrics].tolist()
        values_closed = values + [values[0]]
        area = area_polar_polygon(values_closed, angles_closed)
        areas.append(area)


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