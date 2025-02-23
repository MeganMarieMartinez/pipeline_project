import os
import sys
import argparse
from Bio import SeqIO
import subprocess 
import time 


#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="Returns the number of reads that meet or surpass the quality score threshold at the specified percentage of bases")
    parser.add_argument("-i",'--input', help = 'input file', nargs = '+',required=True) # nargs  allows you to pass multiple input files
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
inputs = arguments.input # this creates a list of all the input files 

'''
Per the assignment, the only data that the parser is expected to take is the full data restrieved from step 1
(both donors) or the sample data. This next part determines what type of data has been passed through
'''
two_dpi = []
six_dpi = []
if len(inputs) == 4: #four fastq files total for the paired-end reads of the two donors
    for file in inputs:
        if file.startswith('SRR5660030'):
            two_dpi.append(file)
        else:
            six_dpi.append(file)
else:
    sample = inputs

## PRELIMINARY SET UP##
os.system('mkdir PipelineProject_MeganMarie_Martinez') # makes a directory that everything shuold be going into

# moving everything into the pipeline project file
for file in os.listdir(os.getcwd()):
    if file in inputs:
        os.rename(os.path.join(os.getcwd(), file),os.path.join('PipelineProject_MeganMarie_Martinez',file))

os.chdir('PipelineProject_MeganMarie_Martinez')
log_file = open("PipelineProject.log", "x")

## STEP 2: Using Bowtie2 to create a HCMV index and mapping samples to the index ##
datasets_command = 'datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report --filename HCMV_dataset.zip'
os.system(datasets_command) # runs the dataset command to pull down the reference genome FASTA file for HCMV
os.system('unzip HCMV_dataset.zip')
HCMV_fasta = 'ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna' # okay, this should always be the same? hard coding this path in should be fine
# building the index
bowtie_command = f'bowtie2-build {HCMV_fasta} HCMV'
os.system(bowtie_command)

# determining read pairs before filtering
def count_reads(fastq):
    count = 0
    for record in SeqIO.parse(fastq,'fastq'):
        count += 1
    return count

t_forward = ''
t_reverse = ''
if len(inputs) == 4:
    pre_filter_two_dpi = count_reads(two_dpi[0])
    pre_filter_six_dpi = count_reads(six_dpi[0])
    for i in two_dpi:
        if i.endswith('1.fastq'):
            t_forward = i
        else:
            t_reverse = i
    mapping_command_two = f'nohup bowtie2 -x HCMV -1 {t_forward} -2 {t_reverse} --al-conc-gz two_dpi_mapped_%.fq.gz -p 2 &'
    os.system(mapping_command_two)
    for i in six_dpi:
        if i.endswith('1.fastq'):
            t_forward = i
        else:
            t_reverse = i
    mapping_command_six = f'nohup bowtie2 -x HCMV -1 {t_forward} -2 {t_reverse} --al-conc-gz six_dpi_mapped_%.fq.gz =p 2 &'
    os.system(mapping_command_six)
    time.sleep(60 * 5)
    os.system('gunzip *.gz')
    post_filter_two_dpi = count_reads('two_dpi_mapped_1.fq')
    post_filter_six_dpi = count_reads('six_dpi_mapped_1.fq')
    with open("PipelineProject.log", "a") as f:
        f.write(f'Donor 1 (2dpi) had {pre_filter_two_dpi} read pairs before Bowtie2 filtering and {post_filter_two_dpi} read pairs after.' + '\n')
        f.write(f'Donor 1 (6dpi) had {pre_filter_six_dpi} read pairs before Bowtie2 filtering and {post_filter_six_dpi} read pairs after.' + '\n')
else: # the sample was used
    pre_filter_sample = count_reads(sample[0])
    for i in sample:
        if i.endswith('1.fastq'):
            t_forward = i
        else:
            t_reverse = i
    mapping_command_sample = f'nohup bowtie2 -x HCMV -1 {t_forward} -2 {t_reverse} --al-conc-gz sample_mapped_%.fq &'
    time.sleep(60 * 3)
    os.system('gunzip *gz')
    post_filter_sample = count_reads('sample_mapped_1.fq')
    with open("PipelineProject.log",'a') as f:
        f.write(f'The sample input had {pre_filter_sample} read pairs before Bowtie2 filtering and {post_filter_sample} read pairs after.' + '\n')

