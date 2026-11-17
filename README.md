

## Table of Contents
- [GAQET2 - Genome Annotation Quality Evaluation Tool](#gaqet2---genome-annotation-quality-evaluation-tool)
- [Requirements](#requirements)
  - [Python Dependencies](#python-dependencies)
  - [Software dependencies](#software-dependencies)
- [Installation](#installation)
  - [InterproScan](#interproscan)
  - [Manual installation (Not recommended)](#manual-installation-not-recommended)
  - [Conda installation (Strongly recommended)](#conda-installation-strongly-recommended)
- [Usage](#usage)
  - [GAQET Benchmarking](#gaqet-benchmarking)
    - [YAML file](#yaml-file)
    - [GAQET arguments](#gaqet-arguments)
    - [Command usage](#command-usage)
    - [GAQET benchmarking output](#gaqet-benchmarking-output)
- [GAQET metrics explanation](#gaqet-metrics-explanation)
  - [General metrics](#general-metrics)
  - [DeTEnGA metrics explanation](#detenga-metrics-explanation)
  - [OMArk Consistency metrics explanation](#omark-consistency-metrics-explanation)
  - [OMArk Completness metrics explanation](#omark-completness-metrics-explanation)
- [GAQET PLOT](#gaqet-plot)
  - [Usage](#usage-1)
  - [Output figure and interpretation](#output-figure-and-interpretation)


# WARNING! right now OMARK bundled with GAQET is not working with latest LUCA database (MAY 2026). Please, use the previous version of LUCA found here: https://omabrowser.org/All.Jul2024/LUCA.h5 
# GAQET2 - Genome Annotation Quality Evaluation Tool

GAQET2 is a Python-based tool designed to evaluate the quality of genome annotations. Using GFF and FASTA files, GAQET2 generates statistical reports to help identify common errors and artifacts in structural gene annotations.


## Requirements

### Python Dependencies

- Python == 3.10
- ete3
- PyYAML

### Software dependencies
- AGAT == 1.4.1 (https://github.com/NBISweden/AGAT)
- GFFread == 0.12.7 (https://github.com/gpertea/gffread)
- Omamer == 2.1.0 and OMARK == 0.3.1 (https://github.com/DessimozLab/OMArk)
- TEsorter == 1.4.7 (https://github.com/zhangrengang/TEsorter)
- InterproScan == 5.72 (https://github.com/ebi-pf-team/interproscan)
- BUSCO == 5.8.3 (https://gitlab.com/ezlab/busco)
- Diamond == 2.0.14.152 (https://github.com/bbuchfink/diamond)
- PSAURON == 1.0.4 (https://github.com/salzberg-lab/PSAURON)
- DeTEnGA == 1.0 (https://github.com/victorgcb1987/DeTEnGA). Already bundled with GAQET

[🔝 Back to Table of Contents](#table-of-contents)

## Installation
### InterproScan
**Interproscan**: we recommend using the github version instead of any conda installation (https://github.com/ebi-pf-team/interproscan-docs/blob/v5/docs/HowToDownload.rst)
Then, add interproscan.sh to your PATH variable:

```bash
export PATH=$PATH:/path/to/interproscan.sh
```
Even if you use the conda installation, **interproscan should be installed manually**. That's because conda installion doesn't have included the databases needed to run this program.

### Manual installation (Not recommended)
If you are feelling adventurous enough you can install GAQET manually, but we are not offering support for this kind of installation.

First, clone this repository:

```bash
git clone https://github.com/victorgcb1987/GAQET2.git
```
Then, install all the requierements shown in the requirements section. 

Then, go to your installation directory and run 

```bash
python setup.py install
```

### Conda installation (Strongly recommended)

First, install **Interproscan** as shown using the link provided in the Interproscan section.
Then, create a conda enviroment, activate it and install GAQET:

```bash
conda create -n GAQET
conda activate GAQET
conda install -c victorgcb gaqet
```
you can check if GAQET is installed by running:

```bash
GAQET -v
```
If GAQET is properly installed, this command should show you the GAQET version.  
You can check and install updates by using:  
```bash
conda update -c victorgcb gaqet
```

[🔝 Back to Table of Contents](#table-of-contents)

## Usage
### GAQET Benchmarking 
#### YAML file
GAQET uses as a primary input a **YAML** file as follows:
```yaml
ID: "SpeciesName"
Assembly: "/path/to/assembly.fasta"
Annotation: "/path/to/annotation.gff3"
Basedir: "/path/to/GAQET/results"
Threads: N
Analysis:
  - AGAT
  - BUSCO
  - PSAURON
  - DETENGA
  - OMARK
  - PROTHOMOLOGY
OMARK_db: "/path/to/omark_db.h5"
OMARK_taxid: NCBItaxonID
BUSCO_lineages:
  -  clade1_odb10
  -  clade2_odb10
PROTHOMOLOGY_tags:
  - TREMBL: "/path/to/uniprot_trembl_db.dmnd"
  - SWISSPROT: "/path/to/uniprot_swssprot.dmnd"
  - MYDB: "/path/to/mydb.dmnd"
DETENGA_db: "rexdb-plant"   # TEsorter database: rexdb-plant, rexdb-metazoa or rexdb
Intron_Artifact_Threshold: 30
Intron_Floor_Threshold: 50


```

| Parameter     | Description                                  |
|---------------|----------------------------------------------|
| ID            | Name of the species                     |
| Assembly      | FASTA genome file                            |
| Annotation    | GFF3/GTF annotation file                    |
| Basedir       | GAQET analysis and results directory       |
| Threads       | Number of threads       |
| Analysis      | List of analysis to run. All of them are optional      |
| OMARK_db      | Path to omark db. Only needed if OMARK is in Analysis      |
| OMARK_taxid | NCBI taxid for OMARK. Only needed if OMARK is in Analysis     |
| BUSCO_lineages | List of BUSCO clades to run. Only needed if BUSCO is in Analysis      |
| PROTHOMOLOGY_tags | List of name and path to DIAMOND proteins database. Only needed if  PROTHOMOLOGY is in Analysis     |
| DETENGA_db | REXdb database used by TEsorter in the mRNA step of DeTEnGA: `rexdb-plant` (Viridiplantae), `rexdb-metazoa` (Metazoa) or `rexdb` (both, e.g. for fungi). It does not change the protein step: the TE Pfam list is always [GAQET/docs/TE_pfams_a_priori_GAQET.txt](GAQET/docs/TE_pfams_a_priori_GAQET.txt). Only needed if DETENGA is in Analysis    |
| Intron_Artifact_Threshold | (Optional) Introns with length (bp) equal or below this value are in the range of potential artifacts. Intron limits depend on the lineage. Only needed if AGAT is in Analysis. Defaults to 30 if omitted |
| Intron_Floor_Threshold | (Optional) Introns with length (bp) equal or below this value are in the review band (suspicious introns, includes the potential artifacts). Only needed if AGAT is in Analysis. Defaults to 50 if omitted |


#### GAQET arguments
Some YAML config file values can be override by using **GAQET arguments**:


| Parameter     | Description                                  |
|---------------|----------------------------------------------|
| --species, -s  |  Override YAML species ID  |
| --genome, -g          | Override YAML Assembly                     |
| --annotation, -a          | Override YAML Annotation                            |
| --taxid, -t          | Override NCBI taxid                    |
| --outbase, -o   | Override YAML outbase       |
| --intron_artifact   | Override YAML `Intron_Artifact_Threshold` (bp)       |
| --intron_floor   | Override YAML `Intron_Floor_Threshold` (bp)       |
| --compress_intermediate, -c   | When the run finishes, pack all intermediate results in `GAQET_intermediate_results.tar.gz` and delete them. Kept uncompressed in the output directory: `{species}_GAQET.stats.tsv`, `GAQET.log.txt` and the DeTEnGA summary `{species}_TE_summary.csv`. A compressed run can't be resumed. Also available in the YAML as `Compress_intermediate: true`       |

Every one of this arguments **are optional**.


#### Command usage

With the YAML file you can **run GAQET** as follows:
```bash
GAQET --YAML {yaml_file}
```
You can override YAML parameters using the following optional commands. This is useful, for example when you want to reutilize things like databases in the YAML file but you want to change the name of the species or the NCBI taxid:

```bash
GAQET --YAML {yaml_file} -s {species} -g {assembly.fasta} -a annotation.gff -t {NBCI_taxid} -o {outdir}
```
One way to reuse easily a YAML file is to add it to a env variable, for example YAML_PATH:

```bash
export YAML_PATH=/path/to/YAML/file
GAQET --YAML YAML_PATH -s {species} -g {assembly.fasta} -a annotation.gff -t {NCBI_taxid} -o {outdir}
```

#### GAQET benchmarking output
The ouput directories should be similar to this one:


📂 outpudir  
 ├── 📂 input_sequences  
 ├── 📂 AGAT_run  
 ├── 📂 BUSCOCompleteness_run  
 ├── 📂 DETENGA_run  
 ├── 📂 DIAMOND_run  
 ├── 📂 OMARK_run  
 ├── 📂 PSAURON_run  
 ├── 📄 GAQET.log.txt  
 ├── 📄 {species}_GAQET.stats.tsv  
 ├── 📄 {species}_short_introns_transcripts.tsv  
 ├── 📄 {species}_short_introns_genes.tsv  
 ├── 📄 {species}_excluded_from_protein_analyses.tsv  

 The two ```{species}_short_introns_*.tsv``` files (written when AGAT is run) list the transcripts and the genes with at least one intron equal or shorter than ```Intron_Floor_Threshold``` (it includes the introns in the potential artifacts range, ```Intron_Artifact_Threshold```), with the position and length of those introns.
 Each of the ```*_run``` directories contains the output of each analysis run. ```log.txt``` file contains things like run errors or time consumed running analysis. Al programs outputs are parsed and their results are stored in a tsv file, the ```{species}_GAQET.stats.tsv``` file.  
 [🔝 Back to Table of Contents](#table-of-contents)

 # GAQET metrics explanation
 ## General metrics

 | Parameter     | Description                                  |
|---------------|----------------------------------------------|
| Species                        |  Value assigned by the user when running the command  |
| NCBI_TaxID                     | Value assigned by the user when running the command                     |
| Assembly_Version               | Name of the assembly file given by the user (not of the intermediate files GAQET makes)                            |
| Annotation_Version             | Name of the original annotation file given by the user (not of the mRNA file made after splitting it with AGAT)                    |
| Assembly_md5                   | md5sum of the assembly file given by the user, to identify it exactly                            |
| Annotation_md5                 | md5sum of the original annotation file given by the user, to identify it exactly                    |
| Gene_Models (N)                | (AGAT) Number of coding genes found in the annotation     |
| Transcript_Models (N)          | (AGAT) Number of coding transcripts found in the annotation     |
| CDS_Models (N)                 | (AGAT) Number of CDS found in the annotation     |
| UTR5' (N)                      | (AGAT) Number of UTR5' annotated on coding transcripts     |
| UTR3' (N)                      | (AGAT) Number of UTR3' annotated on coding transcripts     |
| Both sides UTR' (N)            | (AGAT) Number of coding transcripts models with both UTR's annotated    |
| Overlapping_Gene_Models (N)    | (AGAT) Number of overlapping coding genes     |
| Single Exon Gene Models (N)    | (AGAT) Number of monoexonic coding genes     |
| Single Exon Transcripts (N)    | (AGAT) Number of monoexonic coding transcripts     |
| Total Gene Space (Mb)          | (AGAT) Coding sequences' total size     |
| Mean Gene Model Length (bp)    | (AGAT) Average coding gene length     |
| Mean CDS Model Length (bp)     | (AGAT) Average CDS length     |
| Mean Exon Length (bp)          | (AGAT) Average coding exon length    |
| Mean Intron Length (bp)        | (AGAT) Average coding gene's intron length    |
| Longest Gene Model Length (bp) | (AGAT) Longest coding gene length     |
| Longest CDS Model Length (bp)  | (AGAT) Longest CDS length     |
| Longest Intron Length (bp))    | (AGAT) Longest coding gene's intron length   |
| Shortest Gene Model Length (bp)| (AGAT) Shortest coding gene length     |
| Shortest CDS Length (bp)         | (AGAT) Shortest CDS length     |
| Shortest intron Length (bp)      | (AGAT) Shortest intron length     |
| Introns (N)                      | (AGAT) Total number of introns found in the annotation     |
| Introns <= {Intron_Artifact_Threshold}bp (potential artifacts) (N) | (AGAT) Number of introns with length equal or below the configured `Intron_Artifact_Threshold` (bp). |
| % Introns <= {Intron_Artifact_Threshold}bp (potential artifacts) | (AGAT) Same, as a percentage of all introns (`Introns (N)`) |
| Transcripts with introns <= {Intron_Artifact_Threshold}bp (potential artifacts) (N) | (AGAT) Number of transcripts (isoforms) with at least one intron in that range |
| % Transcripts with introns <= {Intron_Artifact_Threshold}bp (potential artifacts) | (AGAT) Same, as a percentage of all transcripts in the annotation (single-exon transcripts included) |
| Genes with introns <= {Intron_Artifact_Threshold}bp (potential artifacts) (N) | (AGAT) Number of genes with at least one transcript (isoform) with an intron in that range |
| % Genes with introns <= {Intron_Artifact_Threshold}bp (potential artifacts) | (AGAT) Same, as a percentage of all genes in the annotation |
| Introns <= {Intron_Floor_Threshold}bp (floor) (N) | (AGAT) Number of introns with length equal or below the configured `Intron_Floor_Threshold` (bp). It includes the potential artifacts. |
| % Introns <= {Intron_Floor_Threshold}bp (floor) | (AGAT) Same, as a percentage of all introns (`Introns (N)`) |
| Transcripts with introns <= {Intron_Floor_Threshold}bp (floor) (N) | (AGAT) Number of transcripts (isoforms) with at least one intron in that range |
| % Transcripts with introns <= {Intron_Floor_Threshold}bp (floor) | (AGAT) Same, as a percentage of all transcripts in the annotation (single-exon transcripts included) |
| Genes with introns <= {Intron_Floor_Threshold}bp (floor) (N) | (AGAT) Number of genes with at least one transcript (isoform) with an intron in that range |
| % Genes with introns <= {Intron_Floor_Threshold}bp (floor) | (AGAT) Same, as a percentage of all genes in the annotation |
| Models with early STOP           | (AGAT) Number of coding transcripts with premature (in-frame) stop codons, read from the `_report.txt` of `agat_sp_flag_premature_stop_codons.pl` (before this version it was always 0) |
| Models START missing             | (AGAT) Number of coding transcripts lacking start codon     |
| Models STOP missing              | (AGAT) Number of coding transcripts lacking stop codon     |
| Models START & STOP missing      | (AGAT) Number of coding transcripts lacking stop and start codon     |
| Coding transcripts excluded from protein analyses (N) | (gffread) Coding transcripts left out of DeTEnGA, DIAMOND, PSAURON and OMArk because their CDS has no START, no STOP or an in-frame STOP (`gffread -J`). BUSCO uses them. Listed in `{species}_excluded_from_protein_analyses.tsv` |
| % Coding transcripts excluded from protein analyses | Same, as a percentage of all coding transcripts |
| Excluded transcripts with {reason} (N) | Excluded transcripts by reason: no START, no STOP, in-frame STOP, CDS length not multiple of 3, other. A transcript can have more than one reason |
| Longest isoforms excluded from OMArk (N) | Genes whose longest isoform is excluded, so they are missing from the OMArk input |
| Annotation_BUSCO_{DB} | (BUSCO) Busco proteome completness for the database {DB}. Refer to https://busco.ezlab.org/busco_userguide.html#interpreting-the-results to get an output's explanation| 
| PSAURON SCORE | (PSAURON) Global Annotation's accuracy at detecting ORFs|
| DETENGA_FPV | (DeTEnGA) Number of classified transcripts. See table below for nomeclature explanation|
| DETENGA_FP% | (DeTEnGA) Classified transcripts in percentages. See table below for nomeclature explanation|
| OMArk Consistency Results| (OMARK) Taxonomic consistency results. Check table below for nomenclature explanation|
| OMArk Completeness Results| (OMARK) Taxonomic Completness results. Check table below for nomenclature explanation|
| OMArk Species Composition | (OMARK) Species Composition in percentage|
| ProteinsWith{db}Hits (%) | (DIAMOND) Percentage of proteins with a significant hit on database {db}|  
| Longest isoform: Measured transcripts (N) | Number of longest isoforms (one per gene, from `agat_sp_keep_longest_isoform.pl`, the same ones given to BUSCO and OMArk). Denominator of the `Longest isoform:` columns |
| Longest isoform: {column} | The model-quality columns `ProteinsWith{db}Hits (%)`, `PSAURON SCORE`, `DETENGA_FPV`, `DETENGA_FP%`, `Models START missing`, `Models STOP missing`, `Models START & STOP missing`, `Transcripts with introns <= {Intron_Artifact_Threshold}bp (potential artifacts) (N)`, `Transcripts with introns <= {Intron_Floor_Threshold}bp (floor) (N)` and their `%`, re-counted on the longest isoform of each gene (T and percentages over genes). Computed from the files the analyses already wrote, nothing is run twice. Used by GAQET_PLOT for the longest-isoform score |

Transcripts with no START, no STOP or an in-frame STOP are not given to the protein analyses, and each analysis
treats them differently: DIAMOND percentages and DeTEnGA classes are computed over all coding transcripts (excluded
transcripts count as without hits, and as `invalid_CDS` in DeTEnGA); the PSAURON score and OMArk are computed only
over the transcripts (longest isoforms for OMArk) with a complete CDS. Check the exclusion columns when comparing annotations.

[🔝 Back to Table of Contents](#table-of-contents)

 ## DeTEnGA metrics explanation
Detection of Transposable Elements as Genes on Annotations (DeTEnGA) is an in-house tool created for coding sequences **classification as a Transposable element (TE) or not**. This classification is at protein sequence level using interpro and at mRNA level using TEsorter. A Pfam domain is TE-related if it is in [GAQET/docs/TE_pfams_a_priori_GAQET.txt](GAQET/docs/TE_pfams_a_priori_GAQET.txt): 236 Pfam-A 38.2 families whose name or description refers to TEs, or with published evidence of TE specificity. The same list is used whatever the `DETENGA_db`. Table below describes the nomenclature used in DeTEnGA:

 | Nomeclature     | Description                                  |
|---------------|----------------------------------------------|
 |T|Total number of coding transcripts (includes all calssifications)|
 |PcpM0|Transcripts with non-TEs interpro PFAMs and non-TE mRNA|  
 |PTeM0|Transcripts TEs only interpro PFAMs and non-TE mRNA|
 |PchM0|Transcripts with mixed TEs and non-TEs pfams and non-TE mRNA|
 |PcpMTe|Transcripts with non-TEs interpro PFAMs and TE mRNA|
 |PteMte|Transcripts with only TEs interpro PFAMs and TE mRNA|
 |PchMte|Transcripts with mixed TEs and non-TEs pfams and non-TE mRNA|  
 |Not_analyzed|Transcripts that could not be classified: with_Ns (mRNA with Ns, not given to TEsorter), longer_than_100kb (mRNA of 100 kb or more, not given to TEsorter) and invalid_CDS (no START, no STOP or an in-frame STOP, protein not extracted)|
 
 [🔝 Back to Table of Contents](#table-of-contents)


 ## OMArk Consistency metrics explanation
Consistency is a quality measurement and describes the proportion of protein sequences placed into known gene families from the same lineage.  
Check https://doi.org/10.1038/s41587-024-02147-w for an in-depth explanation
 | Nomeclature     | Description                                  |
|---------------|----------------------------------------------|
 |Cons|Taxonomic consistent hits (%)|
 |Inco| Taxonomic inconsistent hits (%)| 
 |Cont| Contaminantion hits (%)
 |Unkn| Unkown hits (%)|
 |P| Partial hits (%)|  
 |F| Fragmented hits (%)|  
 
 [🔝 Back to Table of Contents](#table-of-contents)


## OMArk Completness metrics explanation
Completeness describes how our proteome overlaps with a conserved ancestral gene set of the species’ lineage.   
The gene set is classified in hierarchical orthologous groups (**HOGs**). Each HOG represents a single ancestral gene. Check https://doi.org/10.1038/s41587-024-02147-w for a in-depth explanation
 | Nomeclature     | Description                                  |
|---------------|----------------------------------------------|
 |{Taxon}-HOGs|Number of HOGs in our species' nearest {taxon}  |
 |S| HOGs hits by a single query protein (%)| 
 |D| HOGs hits by more than one query protein (%)|
 |U| HOGs hits by more than one query protein, unexpected (no HOG duplication evidence exists)(%)|
 |E| HOGs hits by more than one query protein, expected (HOG duplication evidence exists, known HOG subfamilies)(%)|
 |M| HOGs without hit (%)|  
 
 [🔝 Back to Table of Contents](#table-of-contents)


 # GAQET PLOT
 You can get an intuitive representation of the most important metrics described before using GAQET PLOT.
 ## Usage
GAQET PLOT requires as input a TSV file as shown in [in the example included with this repository](GAQET/docs/Arabidopsis_seed_stats.tsv). Included example shows two columns for BUSCO and Protein Homology, but you can include more columns or even remove all columns of these two kinds. Then you can run GAQET_PLOT as:  
```bash
GAQET_PLOT -i {input_tsv} -p {prefix} -f {png|jpeg|svg}
```

| Argument | Description |
|---|---|
| -i, --input | (Required) TSV with one row per annotation (`{species}_GAQET.stats.tsv` files concatenated) |
| -p, --prefix | Prefix of the output files, it can include a directory (default `GAQET`). Two plots and two score tables are written (see below) |
| -f, --format | Image format of the plots: `png` (default), `jpeg` or `svg` |
| -r, --radial_min | Value (%) at the centre of the radar plot. Default: the lowest plotted value minus 10 (rounded down, 0 at least), so differences between annotations are visible and no value is clipped. With a manual value, values below it are drawn at the centre. `0` gives the full 0-100% axis |
| -w, --weights | Share (%) of the score given to metrics, `key=percentage` separated by commas (e.g. `completeness=40` or `busco=40,omark=10`). A class key (`completeness`, `precision`, `accuracy`; `false_positives` is accepted as the former name of `precision`) or a role key (`homology`, `busco`, `omark_completeness`, `omark` (consistency only), `psauron`, `detenga`, `agat`, `introns`, `omark_fragments`) sets the total share of those metrics, split among them as the default shares (by analysis); any other key is matched as a substring of one metric name (e.g. `swissprot=15`). The rest is split among the metrics not mentioned in proportion to their default shares. Default: each class 33.33%, split equally among its analyses and then among the metrics of each analysis |
| -e, --exclude | Keys (same matching as `--weights`) of metrics left out of the score. They are still drawn |

Outputs, for `-p {prefix} -f {format}` (e.g. `-p Athaliana` gives `Athaliana_plot.png`, `Athaliana_longest_isoform_plot.png`, `Athaliana_GAS.tsv` and `Athaliana_longest_isoform_GAS.tsv`):

| File | Content |
|---|---|
| `{prefix}_plot.{format}` | GAQET plot, model-quality metrics measured on every transcript |
| `{prefix}_longest_isoform_plot.{format}` | GAQET plot, model-quality metrics measured on the longest isoform of each gene |
| `{prefix}_GAS.tsv`, `{prefix}_longest_isoform_GAS.tsv` | Score tables: one row per annotation, sorted by GAS |
| `{prefix}_GAS_weights.tsv`, `{prefix}_longest_isoform_GAS_weights.tsv` | Weights tables: one row per metric |
| `{prefix}_GAS_breakdown.{format}`, `{prefix}_longest_isoform_GAS_breakdown.{format}` | GAS breakdown: one horizontal bar per annotation, sorted by GAS, split into the points given by each class (class share × sub-score), with the sub-score inside each segment and the GAS at the end |

Columns of the score tables:

| Column | Description |
|---|---|
| Species | Annotation (`Species` column of the input) |
| Measured_transcripts | `all_transcripts` or `longest_isoform` |
| Measured_transcripts_N | Transcripts on which the model-quality metrics are measured (all transcripts, or one per gene) |
| GAS, GAS_rank | Global Annotation Score of the annotation and its rank in the table (1 = highest) |
| Completeness / Precision / Accuracy score (share {S}%) | Class sub-scores (0-100, higher is better; a high precision score means few false positives), with the share of the GAS of each class in the column name. GAS = sum of share × sub-score |
| {metric} (%) | Value of each metric (radial axis) whose analysis was run, e.g. `SWISSPROT hits (%)`, `BUSCO complete eudicots_odb10 (%)` |
| {metric} points (share {S}%) | Points of the GAS coming from that metric (share × value), with the metric's share in the column name; the points of an annotation add up to its GAS (±0.01 from rounding); share and points 0 if the metric is excluded |

Columns of the weights tables (the shares are the same for every annotation):

| Column | Description |
|---|---|
| Metric | Metric, named as in the score table |
| Class | Class of the metric: `completeness`, `precision` or `accuracy` |
| Role | Analysis it comes from (the role keys of `--weights` and `--exclude`) |
| Measured_on | Transcripts used by this metric: `longest_isoform` always for BUSCO and OMArk |
| Share_% | Share of the score given to the metric (`--weights`); 0 if excluded or missing |
| Status | `scored`, `excluded` (`--exclude`) or `analysis missing` (not in the score table) |

### Global Annotation Score (GAS)
The GAS (0-100) is the weighted mean of the metrics drawn as radial axes: GAS = sum(share × metric), shares adding
up to 100%. Every metric belongs to one class of error, and each class has its own sub-score (the share-weighted
mean of its metrics):

- **Completeness**: are the conserved gene families of the lineage present? (BUSCO, OMArk completeness). Only conserved genes are checked (e.g. 954-3,354 BUSCO genes, ~11,000 OMArk families for vertebrates), not lineage-specific ones
- **Precision**: are the models real protein-coding genes, not spurious ORFs, TEs or contaminants? (DIAMOND homology, OMArk consistency, PSAURON, DeTEnGA)
- **Accuracy**: are the models of real genes correctly built? (AGAT START/STOP, artifact introns, OMArk partial hits and fragments)

By default the three classes have the same share (33.33% each), split equally among the analyses of each class and
then among the metrics of each analysis, so the GAS is the mean of the three sub-scores and no class or analysis weighs
more because it has more metrics: homology counts once in the precision class whatever the number of databases
given (e.g. Swiss-Prot and TrEMBL 4.17% each, or one database 8.33%), and BUSCO counts once whatever the number of
lineages. `--weights` and
`--exclude` also accept the class names. A metric whose analysis was not run (or a column missing from older GAQET
outputs) is drawn empty and left out of the score; the shares of the other metrics are rescaled (a class with no
metric left is dropped and the other classes share the score). The accuracy metrics detect broken CDS ends, very
short introns and, for proteins of known gene families only, truncated/split models (fragments) and partial matches
(e.g. fusions); wrong exon boundaries and errors in genes of unknown families are not measured.

| Metric (axis) | Class | Role | Value | Measured on |
|---|---|---|---|---|
| {db} Homology Hits (DIAMOND), one per database | precision | homology | `ProteinsWith{db}Hits (%)` | model quality: all transcripts / longest isoform |
| Proteins fitting their Taxonomic rank (OMARK) | precision | omark | consistent proteins (%) | longest isoform of each gene (always) |
| {lineage} Completeness (BUSCO), one per lineage | completeness | busco | Complete (%) | longest isoform of each gene (always) |
| Proteins with non-TE Pfam domains and no TE mRNA signal (DeTEnGA) | precision | detenga | PcpM0 (%): only non-TE Pfam domains (e-value ≤ 0.005) and no TE signal in the mRNA (TEsorter); transcripts with no Pfam hit, TE evidence or not analysed (invalid_CDS, with_Ns, longer_than_100kb) count against it | model quality: all transcripts / longest isoform |
| Transcripts flagged as coding sequences (PSAURON) | precision | psauron | PSAURON score | model quality: all transcripts / longest isoform |
| No start/stop codon errors (AGAT) | accuracy | agat | 100 − % models with START and/or STOP missing | model quality: all transcripts / longest isoform |
| Transcripts with no artifact (very short) introns (AGAT) | accuracy | introns | 100 − % transcripts with an intron ≤ `Intron_Artifact_Threshold` (potential artifacts) | model quality: all transcripts / longest isoform |
| Conserved gene families found (OMArk) | completeness | omark_completeness | Single + Duplicated (%) of the conserved gene families (HOGs) of the lineage (`OMArk Completeness Results`) | longest isoform of each gene (always) |
| Proteins not fragmented or partial (OMArk) | accuracy | omark_fragments | 100 − (partial hits + fragments) as % of the proteins placed in a known family (consistent + inconsistent + contaminant; `OMArk Consistency Results`). Fragments: < 50% of the family median length; partial hits: only part of the sequence matches the family. Unknown proteins are not evaluated and are not in the denominator | longest isoform of each gene (always) |

**Unit of measure.** BUSCO and OMArk always receive one protein per gene (the longest isoform): given every isoform
they would report multi-isoform genes as duplicated. The model-quality metrics evaluate each gene model, and
GAQET_PLOT computes the score twice: measured on every transcript (a gene with many isoforms weighs more) and on
the longest isoform of each gene (the `Longest isoform:` columns; every gene counts once, so annotations with and
without alternative isoforms are compared on the same footing). Both come from the same GAQET run. For runs made
before these columns existed only the all-transcripts outputs are written, with a warning: rerun `GAQET` on the same
`Basedir` (finished analyses and compressed intermediate results are reused) to add them.

## Output figure and interpretation
Running this script should generate a figure similar to this one:  
![GAQET_plot](GAQET/docs/Arabidopsis_plot.png)

This figure is called GAQET plot and helps comparing multiple annotations. This plot has 2 levels:   

- **Inner radar plot (the score)**:  
  Represents the value of each metric of the Global Annotation Score, one radial axis per metric, for example, `No start/stop codon errors`. Higher percentages mean better. The axes are grouped in three shaded sectors, one per class (Completeness, Precision, Accuracy), labelled with the share of the GAS of the class (see [Global Annotation Score](#global-annotation-score-gas)). The radial axis starts at the lowest plotted value minus 10 (or `--radial_min`; the first label of the scale), so the radar area is not proportional to the score. The GAS and the three class sub-scores of each annotation are given in the legend (`GAS: xx (Completeness: xx, Precision: xx, Accuracy: xx)`) and in the breakdown plot.

- **Outer Column sets (descriptors, not in the score)**:  
  Represent **in-depth annotation metrics**. Each of these bars sits next to the axis it describes. Percentages use an absolute scale (100% = full bar):  
  - **OMARK completeness / consistency / partial-fragmented**: Single, Duplicated and Missing conserved families; Consistent, Inconsistent, Contaminant and Unknown proteins (% of the proteome); partial hits and fragments (% of the proteins placed in a family).  
  - **BUSCO sections**: Higher % Complete (Single) means better.  
    - If there are no valid biological explanations (e.g., polyploidy), a high % Complete (Duplicated) could indicate assembly problems, such as redundant contigs generated by pooling DNA from multiple individuals.  
    - High % Fragmented could indicate genome assembly continuity problems, but this could also be explained by annotation errors if the average CDS length is low.
  - **Transposon related Transcripts**: Here is represented the % of transcript models with TE evidence in DeTEnGA, the sum of the classes PteM0, PchM0, PcpMte, PteMte, PchMte and P0Mte (TE Pfam domains and/or TE signal in the mRNA). Usually, higher values means worse, but sometimes these sequences are flagged as TEs because they are genes originated from TEs domestication. Check DeTEnGA output file ```{species}_TE_summary.csv```list of sequences IDs for further analysis on these sequences.
  -  **Transcripts with intron length shorter than floor**: Represents the % of transcripts with at least one intron ≤ `Intron_Floor_Threshold` (the review band of suspicious introns, which includes the potential artifacts scored by the "Transcripts with no artifact (very short) introns" axis; in the longest-isoform plot, % of longest isoforms). 100% = full bar. Lower means better, but the floor is close to the shortest real introns of the lineage, so part of these introns can be genuine. The genes and transcripts are listed in `{species}_short_introns_genes.tsv` and `{species}_short_introns_transcripts.tsv`.

- **Gene model descriptors (AGAT, not in the score)**:  
  Four outer bars in the free slots of the axes without a bar of their own (PSAURON, START/STOP and the two homology databases); they do not describe those axes. Gene models, transcript models and mean CDS length are relative to the largest annotation of the plot (largest = full bar), so they change with the set of annotations; UTR is a % (100% = full bar):  
  -  **Gene models and transcript models**: how many models are in each annotation. Extreme differences could be explained by annotation problems, e.g. Helixer doesn't group isoforms into genes, so its number of gene models can be higher than the actual number of genes.
  -  **Transcripts with both sides UTR**: % of transcripts with both 5' and 3' UTR annotated. Higher means better.
  -  **Average CDS length**: average CDS length. Extreme differences could be explained by annotation errors, e.g. genes merged into one model, or fragmented models.

 
