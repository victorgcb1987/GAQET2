# GAQET2 - Genome Annotation Quality Evaluation Tool

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
- BUSCO == 5.8.3 (https://gitlab.com/ezlab/busco)
- Diamond == 2.0.14.152 (https://github.com/bbuchfink/diamond)
- PSAURON == 1.0.4 (https://github.com/salzberg-lab/PSAURON)
- DeTEnGA == 1.0 (https://github.com/victorgcb1987/DeTEnGA). Already bundled with GAQET

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
conda install -c bioconda omark
conda install -c bioconda busco
conda install python==3.10
pip install psauron
pip install ete3
pip install PyYAML
```

**Interproscan**: we recommend using the github version instead of any conda installation (https://github.com/ebi-pf-team/interproscan-docs/blob/v5/docs/HowToDownload.rst)
Then, add interproscan.sh to your PATH variable:

`export PATH=$PATH:/path/to/interproscan.sh`

---

## Usage

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


With the YAML file you can **run GAQET** as follows:

```bash
GAQET --YAML {yaml_file}
```
Some YAML config file values can be override by using **GAQET arguments**:


| Parameter     | Description                                  |
|---------------|----------------------------------------------|
| --species, -s  |  Override YAML species ID  |
| --genome, -g          | Override YAML Assembly                     |
| --annotation, -a          | Override YAML Annotation                            |
| --taxid, -t          | Override NCBI taxid                    |
| --outbase, -o   | Override YAML outbase       |



```bash
GAQET --YAML {yaml_file} -s {species} -g {assembly.fasta} -a annotation.gff -t 3702 -o {outdir}
```
One way to reuse easily a YAML file is to add it to a env variable, for example YAML_PATH:

```bash
export YAML_PATH=/path/to/YAML/file
GAQET --YAML YAML_PATH -s {species} -g {assembly.fasta} -a annotation.gff -t 3702 -o {outdir}
```
