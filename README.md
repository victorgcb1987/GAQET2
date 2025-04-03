# GAQET2 - Genome Annotation Quality Evaluation Tool (v2)

GAQET2 is a Python-based tool designed to evaluate the quality of genome annotations. Using GFF and FASTA files, GAQET2 generates statistical reports to help identify common errors and artifacts in gene structural annotations.

---

## ✨ Features

- Support for GFF3 annotations and FASTA sequences
- Automated quality evaluation
- Tabular and graphical reports
- Compatible with both prokaryotic and eukaryotic annotations
- Command-line interface

---

## 📆 Requirements

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
- Diamond == 2.1.11 (https://github.com/bbuchfink/diamond)
- PSAURON == 

## ⚙️ Installation

Clone the repository:

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
conda install -c bioconda busco
conda install python==3.10
pip install psauron
pip install ete3
pip install PyYAML
```

---

## 🚀 Usage

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
python gaqet2.py -g data/annotation.gff -f data/genome.fasta -o results/ --summary --plot
```

---

## 🗂️ Repository Structure

```
GAQET2/
├── gaqet2.py               # Main script
├── utils/                  # Helper functions
│   ├── fasta_utils.py
│   └── gff_utils.py
├── data/                   # Example data
├── results/                # Output folder (generated at runtime)
├── requirements.txt        # Python dependencies
├── README.md               # Documentation
└── LICENSE
```

---

## 👥 Credits

Developed by **Vic
