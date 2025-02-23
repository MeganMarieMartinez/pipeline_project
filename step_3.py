import os


# function to run spades regardless of input type (real data or sample)
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
    spades_command = f'nohup spades.py -k 99 -t 2 --only-assembler --pe1-1 {forward_1} --pe1-2 {reverse_1} --pe2-1 {forward_2} --pe2-2 {reverse_2} -o SPADES_assembly/ &'
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

if __name__ == "__main__":
    spades()