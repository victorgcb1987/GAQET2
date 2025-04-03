# GAQET
Genome Annotation Quality Evaluation Tools

GAQET is a tools configura

## Requirements and installation



**GffRead**: you can download gffread from https://github.com/gpertea/gffread

**TEsorter**: we recommend installing TEsorter using conda as follows: `conda create -n TEsorter -c bioconda tesorter`. Then, update python from this conda installation (python v3.6 uses some python deprecated functions): `conda install python=3.12`.


**Interproscan**: we recommend using the github version instead of any conda installation (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
Then, add interproscan.sh to your PATH variable:

`export PATH=$PATH:/path/to/interproscan.sh`

Recommended installation instructions
`conda create -c bioconda -n GAQET agat  
 conda install python==3.10`

**DeTEnGA**: just clone this respository `git clone https://github.com/victorgcb1987/DeTEnGA.git`


## How to use
In order to use DeTEnga you will need at least an uncompressed FASTA file with the genome assembly and an uncompressed genome annotation file (gff or gtf). You can run it with multiple annotations and assemblies. There is an exemple for running this program:  

``DeTEnGA.py -i fof.txt -o output_dir -t num_threads -s rexdb-plant``  
**--input, -i**:  (Required) file of files used as an input for DeTEnGA    
**--output, -o**: (Required) output directory  
**--threads, -t**: (Optional, default = 1) number of threads  
**--tesorter_database, -s**: (Optional, default = "rexdb-plant") database used with TEsorter)

The file of files is a plain text in tabular format with three columns, being the first one a label for your analyzed annotation, a path for your assembly and the path for the annotation, for example:  
Nicotiana_benthamiana&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/path/to/nicoben/assembly.fasta&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/path/to/nicoben/annotation.gff  
Arabidopsis_thaliana&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/path/to/arathal/assembly.fasta&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/path/to/arathal/annotation.gff  


