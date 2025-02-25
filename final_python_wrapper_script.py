import os
import sys
import argparse

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="Wrapper script that runs a pipeline to analyze a virus genome")
    parser.add_argument("-i",'--input', help = 'input file', nargs = '+',required=True) # nargs  allows you to pass multiple input files
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
inputs = arguments.input # this creates a list of all the input files 

## PRELIMINARY SET UP##
os.system('mkdir PipelineProject_MeganMarie_Martinez') # makes a directory that everything should be going into

# moving everything into the pipeline project file
for file in os.listdir(os.getcwd()):
    if file in inputs:
        os.rename(os.path.join(os.getcwd(), file),os.path.join('PipelineProject_MeganMarie_Martinez',file))

os.chdir('PipelineProject_MeganMarie_Martinez')
log_file = open("PipelineProject.log", "x")
formatted_inputs = " ".join(inputs)

## STEP 2

from Bio import SeqIO

'''
This script takes fastq files and output a fastq file for each donor that was mapped to the ref genome. these mapped reads will be used for assembly in step 3

'''

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
if len(inputs) == 4 and inputs[0].startswith('S'): #four fastq files total for the paired-end reads of the two donors
    for file in inputs:
        if file.startswith('SRR5660030'):
            two_dpi.append(file)
        else:
            six_dpi.append(file)
else:
    sample = inputs

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
if inputs[0].startswith('S'): # if this is the actual data
    pre_filter_two_dpi = count_reads(two_dpi[0])
    pre_filter_six_dpi = count_reads(six_dpi[0])
    for i in two_dpi:
        if i.endswith('1.fastq'):
            t_forward = i
        else:
            t_reverse = i
    mapping_command_two = f' bowtie2 -x HCMV -1 {t_forward} -2 {t_reverse} --al-conc-gz two_dpi_mapped_%.fq.gz -p 2'
    os.system(mapping_command_two)
    for i in six_dpi:
        if i.endswith('1.fastq'):
            t_forward = i
        else:
            t_reverse = i
    mapping_command_six = f'bowtie2 -x HCMV -1 {t_forward} -2 {t_reverse} --al-conc-gz six_dpi_mapped_%.fq.gz -p 2'
    os.system(mapping_command_six)
    os.system('gunzip *.gz')
    post_filter_two_dpi = count_reads('two_dpi_mapped_1.fq')
    post_filter_six_dpi = count_reads('six_dpi_mapped_1.fq')
    with open("PipelineProject.log", "a") as f:
        f.write(f'Donor 1 (2dpi) had {pre_filter_two_dpi} read pairs before Bowtie2 filtering and {post_filter_two_dpi} read pairs after.' + '\n')
        f.write(f'Donor 1 (6dpi) had {pre_filter_six_dpi} read pairs before Bowtie2 filtering and {post_filter_six_dpi} read pairs after.' + '\n')
else: # the sample was used
    t_forward_1 = ''
    t_forward_1 = ''
    t_forward_2 = ''
    t_reverse_2 = ''
    for s in sample:
        if s.startswith('sampledata_'):
            if s.endswith('1.fastq'):
                t_forward_1 = s
            else:
                t_reverse_1 = s
        else:
            if s.endswith('1.fastq'):
                t_forward_2 = s
            else:
                t_reverse_2 = s
    pre_filter_sample_1 = count_reads(t_forward_1)
    pre_filter_sample_2 = count_reads(t_forward_2)
    mapping_command_sample_1 = f'bowtie2 -x HCMV -1 {t_forward_1} -2 {t_reverse_1} --al-conc-gz sample_1_mapped_%.fq.gz -p 2'
    os.system(mapping_command_sample_1)
    mapping_command_sample_2 =  f'bowtie2 -x HCMV -1 {t_forward_2} -2 {t_reverse_2} --al-conc-gz sample2_mapped_%.fq.gz -p 2'
    os.system(mapping_command_sample_2)
    os.system('gunzip *gz')
    post_filter_sample_1 = count_reads('sample_1_mapped_1.fq')
    post_filter_sample_2 = count_reads('sample2_mapped_1.fq')
    with open("PipelineProject.log",'a') as f:
        f.write(f'The first sample input had {pre_filter_sample_1} read pairs before Bowtie2 filtering and {post_filter_sample_1} read pairs after.' + '\n')
        f.write(f'The second sample input had {pre_filter_sample_2} read pairs before Bowtie2 filtering and {post_filter_sample_2} read pairs after.' + '\n')

# since datasets will be used again later, removing the datasets directory and related files
os.system('rm -r ncbi_dataset')
os.system('rm README.md')
os.system('rm md5sum.txt')

## STEP 3

def spades():
    input = []
    for file in os.listdir(os.getcwd()):
        if file.endswith('.fq'): # the other fastq files in the directory right now have a different file extension
            input.append(file)
    paired_ends = pairing(input)
    forwards = [] # these will be the forward reads, this should always just have two 
    for key in paired_ends.keys():
        forwards.append(key)
    # making sure that the correct forward and reverse reads are paired together
    forward_1 = forwards[0]
    reverse_1 = paired_ends[forward_1]
    forward_2 = forwards[1]
    reverse_2 = paired_ends[forward_2]
    spades_command = f'spades.py -k 99 -t 2 --only-assembler --pe1-1 {forward_1} --pe1-2 {reverse_1} --pe2-1 {forward_2} --pe2-2 {reverse_2} -o SPADES_assembly/'
    os.system(spades_command)
    with open("PipelineProject.log",'a') as f:
        f.write(f'SPAdes command used for assembly: {spades_command}' + '\n')

# pairing up the paired end reads via dictionary
def pairing(input):
    pairs = {}
    for i in input:
        for j in input:
            s = i.split('_')
            if j.startswith(s[0]) and j != i:
                if j.endswith('1.fq'):
                    pairs[j] = i
                else:
                    pairs[i] = j
    return pairs

spades() # actually running the script

## STEP 4

def contig_calculations():
    count = 0
    assembly_length = 0
    contigs = './SPADES_assembly/contigs.fasta'  # path to the contigs fasta (current directory when this run should be the PipelineProject_MeganMarie_Martinez)
    for record in SeqIO.parse(open(contigs),'fasta'):
        if len(record.seq) > 1000:
            count += 1
            assembly_length += len(record.seq)
    with open("PipelineProject.log",'a') as f:
        f.write(f'There are {count} contigs > 1000 bp in the assembly' + '\n')
        f.write(f'There are {assembly_length} bp in the assembly' + '\n')

if __name__ == "__main__":
    contig_calculations()

contig_calculations()

## STEP 5 
def blast():
    long_holders = []
    contigs = './SPADES_assembly/contigs.fasta'
    for record in SeqIO.parse(open(contigs),'fasta'):
        if len(record.seq) > 1000:
            long_holders.append(record.seq)
    long_holders.sort(key = len, reverse = True)
    longest_contig = long_holders[0]
    with open('longest_contig.fasta','w') as f:
        f.write(str(longest_contig))
    # setting up the blastdb
    betaherpesvirinae_datasets_command = 'datasets download virus genome taxon betaherpesvirinae --refseq --include genome' # command to pull down the subfamily genomes
    os.system(betaherpesvirinae_datasets_command)
    os.system('unzip ncbi_dataset.zip')
    betaherpesvirinae = './ncbi_dataset/data/*genomic.fna' # path to the fasta w genomes from the Betaherpesvirinae subfamily
    makeblast_db = f'makeblastdb -in {betaherpesvirinae} -out betaherpesvirinae -title beatherpesvirinae -dbtype nucl' # command to make the blastdb
    os.system(makeblast_db)
    # blasting against database
    blastn_command = f'blastn -query longest_contig.fasta -db betaherpesvirinae -max_target_seqs 11 -max_hsps 1  -out blast_output.tsv -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    os.system(blastn_command)
    # writing to .log
    with open("PipelineProject.log",'a') as f:
        f.write('sacc   pident  length  qstart  qend    sstart  send    bitscore    evalue  stitle' + '\n')
        blast_output = open('blast_output.tsv','r')
        blast_output = blast_output.read()
        output_lines = blast_output.splitlines()
        output_lines = output_lines[1:]
        for line in output_lines:
            f.write(line + '\n')
blast()