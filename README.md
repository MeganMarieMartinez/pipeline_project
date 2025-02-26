# Pipeline Project
Python pipeline project for COMP 383

## Program Description
This pipeline maps and assembles the viral transcriptome of Human herpesvirus 5 (HCMV) from a patient-donor 2 and 6 days post infection. It then analyzes the assembly and performs BLAST to determine alignment with other virus strains.

## Programs
1. Python version 3.10.12 (https://www.python.org/)

## Package Dependencies 
1. SPAdes genome assembler v4.0.0 (https://ablab.github.io/spades/index.html)
2. Bowtie 2 v2.4.4 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
3. SRA Toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
3. NCBI Datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
4. BLAST+ (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)
5. Biopython (https://biopython.org/wiki/Documentation)

## Retrieving Transcriptomes
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360

Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363

### For each respective transcriptome:
1. Following the link, go to the SRA database, and click on the the Data Access tab
2. Copy the SRA normalized link
3. In the command line, use `wget {link}` to download the SRA data
4. Use fasterq-dump to convert the SRA file to FASTQ format using the SRA ID number:
    * Donor 1 (2dpi): `fasterq-dump SRX2896360`
    * Donor 1 (6dpi): `fasterq-dump SRX2896363`

## Downloading Pipeline

### Github
1. Install dependencies
2. Download this project using git clone or via zip file

### Sample Data
Accessing sample data to quickly test the whole pipeline
1. After downloading the pipeline files, unzip the sample_data.zip file: `unzip sample_data.zip`
2. After unzipping, there will be 2 sets of paired-end FASTQ files totaling in 4 files

## Running the Pipeline
This pipeline can be run using the command: 
`python final_python_wrapper_script.py -i <forward_fastq_1> <reverse_fastq_1> <forward_fastq_2> <reverse_fastq_2>`

The order of the fastq files does not matter, but both pairs of paired-end reads must be included for a total of 4 input files in order for this script to function correctly.

### Example Command Using the Sample Data: 
`python final_python_wrapper_script.py -i sampledata_1.fastq sampledata_2.fastq sampledata2_1.fastq sampledata2_2.fastq`