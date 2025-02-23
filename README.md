# pipeline_project
Python pipeline project for COMP 383

# Program Description
Will add a program description later

# Programs
1. Python version 3.10.12 (https://www.python.org/)

# Package Dependencies 
1. SPAdes genome assembler v4.0.0 (https://ablab.github.io/spades/index.html)
2. Bowtie 2 v2.4.4 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
3. SRA Toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
3. NCBI Datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
4. BLAST+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

# Retrieving Transcriptomes (rough)
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360

Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363

For each respective transcriptome:
1. Following the link, go to the SRA database and click on the the Data access tab
2. Copy the SRA normalized link
3. In the command line, use 'wget {link}' to download the SRA data
4. Use fasterq-dump to convert the file to FASTQ format using the SRA ID number: 'fasterq-dump SRRXXXXX' 