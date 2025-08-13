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
DETENGA_db: "rexdb-plant"


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
| DETENGA_db | DeTEnGA database for interpro checks. Only needed if DETENGA is in Analysis    |


#### GAQET arguments
Some YAML config file values can be override by using **GAQET arguments**:


| Parameter     | Description                                  |
|---------------|----------------------------------------------|
| --species, -s  |  Override YAML species ID  |
| --genome, -g          | Override YAML Assembly                     |
| --annotation, -a          | Override YAML Annotation                            |
| --taxid, -t          | Override NCBI taxid                    |
| --outbase, -o   | Override YAML outbase       |

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

 Each of the ```*_run``` directories contains the output of each analysis run. ```log.txt``` file contains things like run errors or time consumed running analysis. Al programs outputs are parsed and their results are stored in a tsv file, the ```{species}_GAQET.stats.tsv``` file.  
 [🔝 Back to Table of Contents](#table-of-contents)

 # GAQET metrics explanation
 ## General metrics

 | Parameter     | Description                                  |
|---------------|----------------------------------------------|
| Species                        |  Value assigned by the user when running the command  |
| NCBI_TaxID                     | Value assigned by the user when running the command                     |
| Assembly_Version               | Actually, the name of the assembly file. GAQET_REVIEWER provides md5sums of this file for better identification                            |
| Annotation_Version             | Actually, the name of the annotation file. GAQET_REVIEWER provides md5sums of this file for better identification                    |
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
| Models with early STOP (N)       | (AGAT) Number of coding transcripts with premature stop codons     |
| Models START missing             | (AGAT) Number of coding transcripts lacking start codon     |
| Models START & STOP missing      | (AGAT) Number of coding transcripts lacking stop and start codon     |
| Annotation_BUSCO_{DB} | (BUSCO) Busco proteome completness for the database {DB}. Refer to https://busco.ezlab.org/busco_userguide.html#interpreting-the-results to get an output's explanation| 
| PSAURON SCORE | (PSAURON) Global Annotation's accuracy at detecting ORFs|
| DETENGA_FPV | (DeTEnGA) Number of classified transcripts. See table below for nomeclature explanation|
| DETENGA_FP% | (DeTEnGA) Classified transcripts in percentages. See table below for nomeclature explanation|
| OMArk Consistency Results| (OMARK) Taxonomic consistency results. Check table below for nomenclature explanation|
| OMArk Completeness Results| (OMARK) Taxonomic Completness results. Check table below for nomenclature explanation|
| OMArk Species Composition | (OMARK) Species Composition in percentage|
| ProteinsWith{db}Hits (%) | (DIAMOND) Percentage of proteins with a significant hit on database {db}|  

[🔝 Back to Table of Contents](#table-of-contents)

 ## DeTEnGA metrics explanation
Detection of Transposable Elements as Genes on Annotations (DeTEnGA) is an in-house tool created for coding sequences **classification as a Transposable element (TE) or not**. This classification is at protein sequence level using interpro and at mRNA level using TEsorter. Table below describes the nomenclature used in DeTEnGA:

 | Nomeclature     | Description                                  |
|---------------|----------------------------------------------|
 |T|Total number of coding transcripts (includes all calssifications)|
 |PcpM0|Transcripts with non-TEs interpro PFAMs and non-TE mRNA|  
 |PTeM0|Transcripts TEs only interpro PFAMs and non-TE mRNA|
 |PchM0|Transcripts with mixed TEs and non-TEs pfams and non-TE mRNA|
 |PcpMTe|Transcripts with non-TEs interpro PFAMs and TE mRNA|
 |PteMte|Transcripts with only TEs interpro PFAMs and TE mRNA|
 |PchMte|Transcripts with mixed TEs and non-TEs pfams and non-TE mRNA|  
 
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
GAQET_PLOT -i {input_tsv} -o output_figure.{extension}
```
Figure format is defined by the extension. By default is ```jpeg```

## Output figure and interpretation
Running this script should generate a figure similar to this one:  
![GAQET_plot](GAQET/docs/Arabidopsis_plot.png)

This figure is called GAQET plot and helps comparing multiple annotations. This plot has 3 levels:   

- **Inner radar plot**:  
  Represents percentages for each of the metrics marked by the radial axis, for example, `No start/stop codon errors`. Higher percentages mean better.

- **Inner ring**:  
  Represents the area calculated from the radar plot for each annotation represented.  
  - The area value is also shown as a numerical value in the legend.  
  - This is the annotation quality's **global score** — higher means better.

- **Outer Column sets**:  
  Represent **in-depth annotation metrics**:  
  - **BUSCO sections**: Higher % Complete (Single) means better.  
    - If there are no valid biological explanations (e.g., polyploidy), a high % Complete (Duplicated) could indicate assembly problems, such as redundant contigs generated by pooling DNA from multiple individuals.  
    - High % Fragmented could indicate genome assembly continuity problems, but this could also be explained by annotation errors if the average CDS length is low.
  - **Transposons**: Here is represented the % of transcript models classified as PTeMte (potential TEs). Usually, higher values means worse, but sometimes these sequences are flagged as TEs because they are genes originated from TEs domestication. Check DeTEnGA output file ```{species}_TE_summary.csv```list of sequences IDs for further analysis on these sequences.
  -  **Gene models and transcript models**: Represents how many models are in each annotation. This columns are **normalized by the highest number of models present in all annotations**. For example, if you are representing annotation A with 200 gene models and annotation B with 50 gene models, Annotation B's bar will be 4 times lower than A. Higher means better, but **extreme differences could be explained by annotation problems**, e.g. Helixer doesn't group isoforms into genes, so total number of gene models could be potentially higher than the actual number of organism's genes.
  -  **Transcripts with both sides UTR**: Represents the % of each annotation's transcripts that have both UTR'5 and UTR'3 annotated. Higher means better.
  -  **Average CDS length**: Represents the average CDS length. These values are **normalized by the highest Average CDS length**, so if Annotation A is 1000bp and annotation B is 250bp annotation A bar should be 4 times higher. Higher means better, but **extreme differences could be explained by annotation errors**, e.g. incorrect merging of genes as a single gene model.

 
