# GAQET2 - Genome Annotation Quality Evaluation Tool (v2)

GAQET2 is a Python-based tool designed to evaluate the quality of genome annotations. Using GFF and FASTA files, GAQET2 generates statistical reports to help identify common errors and artifacts in gene structural annotations.


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
- BUSCO == 5.8.3 (https://github.com/metashot/busco)
- Diamond == 2.0.14.152 (https://github.com/bbuchfink/diamond)
- PSAURON == 1.0.4 (https://github.com/salzberg-lab/PSAURON)

## Installation

Clone this repository:

```bash
git clone https://github.com/victorgcb1987/GAQET2.git
```

Create a conda enviroment:

```bash
conda create -c bioconda -n GAQET agat
conda activate GAQET
conda install -c bioconda psauron
conda install -c bioconda busco
conda install -c bioconda gffread
conda install -c bioconda tesorter
conda install -c bioconda diamond
conda install -c bioconda busco
conda install python==3.10
pip install psauron
pip install ete3
pip install PyYAML
```

**Interproscan**: we recommend using the github version instead of any conda installation (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
Then, add interproscan.sh to your PATH variable:

`export PATH=$PATH:/path/to/interproscan.sh`

---

## Usage

```yaml
ID: "Arabidopsis thaliana"
Assembly: "Athaliana_447_TAIR10.fa"
Annotation: "Athaliana_447_Araport11.gene_exons.gff3"
Basedir: "athaliana_QC"
Analysis:
  - AGAT
  - BUSCO
  - PSAURON
  - DETENGA
  - OMARK
  - PROTHOMOLOGY
OMARK_db: "/data/shared_dbs/omark/LUCA.h5"
BUSCO_lineages:
  -  viridiplantae_odb10
  -  embryophyta_odb10
OMARK_taxid: 3702
PROTHOMOLOGY_tags:
  - TREMBL: "/data/shared_dbs/swissprot/uniprot_trembl_r2025_01.dmnd"
  - SWISSPROT: "/data/shared_dbs/swissprot/uniprot_sprot_r2025_01.dmnd"
Threads: 40
DETENGA_db: "rexdb-plant"
```

Run the main script from the terminal:

```bash
python gaqet2.py -g annotation.gff -f genome.fasta -o output_folder
```

### Arguments:

| Parameter     | Description                                  |
|---------------|----------------------------------------------|
| `-g`          | GFF file with annotation                     |
| `-f`          | FASTA genome file                            |
| `-o`          | Output folder for reports                    |
| `--summary`   | (optional) Generates summary statistics       |
| `--plot`      | (optional) Generates distribution plots       |

### Example:

```bash
gaqet -g data/annotation.gff -f data/genome.fasta -o results/ --summary --plot
```
