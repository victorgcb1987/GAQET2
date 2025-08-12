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


---

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

 #### GAQET metrics explanation

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
| Shortest CDS Length (bp)| (AGAT) | (AGAT) Shortest CDS length     |
| Shortest intron Length (bp)      | (AGAT) Shortest intron length     |
| Models with early STOP (N)       | (AGAT) Number of coding transcripts with premature stop codons     |
| Models START missing             | (AGAT) Number of coding transcripts lacking start codon     |
| Models START & STOP missing      | (AGAT) Number of coding transcripts lacking stop and start codon     |
| Annotation_BUSCO_{DB} | (BUSCO) Busco proteome completness for the database {DB}. Refer to https://busco.ezlab.org/busco_userguide.html#interpreting-the-results to get an output's explanation.| 

