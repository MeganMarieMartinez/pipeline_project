# pipeline_project
Python pipeline project for COMP 383

# Program Description
Will add a program description later

# Dependencies (Unspecified for now), adding links later
1. SPAdes
2. Bowtie 2
3. SRA Toolkit 
3. NCBI Datasets
4. BLAST+ 

# Retrieving Transcriptomes (rough)
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363

For each respective transcriptome:
1. Following the link, go to the SRA database onto the Data access tab
2. Copy the SRA normalized link
3. In the command line, use 'wget {link}' to download the SRA data
4. Use fasterq-dump to convert the file to FASTQ format