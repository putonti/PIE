import os
import csv
from Bio import SeqIO

# INPUTS

#### THE SCRIPT SHOULD DO THIS & IT SHOULDN'T BE PASSED IN ####
# input: spades assembly of reads to bacterial contigs
bact_contigs_path='/Volumes/TaylorWorkD/Bowtieoutputs/AssembledBowties/C1355contB1354/'
bact_contigs_file='contigs.fasta'
#### THE SCRIPT SHOULD DO THIS & IT SHOULDN'T BE PASSED IN ####

#### THIS HAS TO BE INPUTTED BY THE USER ####
# input: phage sequences for the strain
phage_sequences_path='/Users/taylormiller-ensminger/Desktop/TaylorThesisAll/PhagePredictionsAll/1354PhagePredicts/SecondTry/'
phage_sequences_file='1354PhageAll2.fasta'
#### THIS HAS TO BE INPUTTED BY THE USER ####

# input: bbmap results file to the bacterial contigs
# if you don't already have this, you could include this in the script to actually perform the bbmap
bbmap_path_bact='/Volumes/TaylorWorkD/Bowtieoutputs/CovStats/'
bbmap_file_bact='C1355contB1354covstats.txt'

# input: bbmap results file to the phage predicted sequences
# if you don't already have this, you could include this in the script to actually perform the bbmap
bbmap_path_phage='/Volumes/TaylorWorkD/Bowtieoutputs/CovStats/FixedPhage/'
bbmap_file_phage='C1355Pcont1354covstats.txt'


# OUTPUT FOLDER
output_file='/Users/taylormiller-ensminger/Desktop/1354Cont_output.csv'


# THIS IS MY PATH FOR MY TOOLS -- feel free to remove this here and from the below function
path_of_tools = "/Users/genevieve/Downloads/tools_and_programs/"

# I didn't include paths for the input/output, can easily add those in as an argument
## GJ -- assume that the fastqs are passed in with the paths (that's what I've done elsewhere) -- see comment on line 175
def process_raw_reads(*fastqs, platform='Illumina'):
    if len(fastqs) == 1:
        bbduk_trim = path_of_tools+"bbmap/bbduk.sh -Xmx1G overwrite=t in="+fastqs[0]+" out=trimmed_"+fastqs[0]+" qtrim=rl ftl=15 ftr=135 maq=20 maxns=0 stats="+fastqs[0]+"_read_qualTrimming.stats statscolumns=5 trimq=20"
        spades_assembly = path_of_tools+"SPAdes-3.15.3-Darwin/bin/metaspades.py --only-assembler -12 trimmed_"+fastqs[0]+" -o "+fastqs[0].strip(".fastq.gz")+"_assembly"
        os.system(bbduk_trim)
        os.system(spades_assembly)
        return fastqs[0].strip(".fastq.gz")+"_assembly"
    elif len(fastqs) == 2:
        bbduk_trim = path_of_tools+"bbmap/bbduk.sh -Xmx1G overwrite=t in1="+fastqs[0]+" in2="+fastqs[1]+" out1=trimmed_"+fastqs[0]+" out2=trimmed_"+fastqs[1]+" qtrim=rl ftl=15 ftr=135 maq=20 maxns=0 stats="+fastqs[0]+"_"+fastqs[1]+"_read_qualTrimming.stats statscolumns=5 trimq=20"
        spades_assembly = path_of_tools+"SPAdes-3.15.3-Darwin/bin/metaspades.py --only-assembler -1 trimmed_"+fastqs[0]+" -2 trimmed_"+fastqs[1]+" -o "+fastqs[0].strip(".fastq.gz")+"_assembly"
        os.system(bbduk_trim)
        os.system(spades_assembly)
        return fastqs[0].strip(".fastq.gz")+"_assembly"
    else:
        print("Expected 1 or 2 fastq files. Please check fastq files and try again.")

def trim_assembly(contigs):
    # remove contigs that are less than 1000 bp
    c=list(SeqIO.parse(contigs,'fasta'))
    file_name=contigs[:contigs.rfind(".")]+"_trimmed.fasta"
    with open(file_name,'w') as f:
        for i in c:
            if len(str(i.seq))>=1000:
                f.write('>'+i.id+'\n'+str(i.seq)+'\n')
    return file_name

def categorize_assembled_contigs(contigs,phage_sequences)
    # create blast database of the predicted phage sequences
    phage_name=phage_sequences[:phage_sequences.find('.')]
    command='makeblastdb -in '+phage_sequences+' -out '+phage_name+' -title '+phage_name+' -dbtype nucl'
    os.system(command)  # uncomment to run
    
    # blast assembly against phage database to figure out which ones are phage
    command='blastn -query '+contigs+' -db '+ phage_name+' -max_target_seqs 1 -outfmt="10 qseqid sseqid qcovs pident length evalue bitscore" -out blastn_'+phage_name+'.csv'
    os.system(command)
    
    # read in blast results
    with open('blastn_'+phage_name+'.csv','r') as f:
        fieldnames = ['qseqid','sseqid','qcovs','pident','length','evalue','bitscore']
        reader = csv.DictReader(f,fieldnames=fieldnames)
        blast_results= list(reader)
    
    # grab top hit for each blast result
    top_hits=dict()
    for i in blast_results:
        pair=(i['qseqid'],i['sseqid'])
        if pair in top_hits.keys():
            q,p,b=top_hits[pair]
            if float(i['bitscore'])>float(b):
                top_hits[pair]=(i['qcovs'],i['pident'],i['bitscore'])
        else:
            top_hits[pair]=(i['qcovs'],i['pident'],i['bitscore'])
    
    # toss garbage -- qcov < 90%
    phage_hits=dict()
    for i in top_hits:
        q,p,b=top_hits[i]
        if float(q)>=90:
            phage_hits[i]=top_hits[i]
    
    # make lists of the bacterial contigs that are really bacterial
    bacterial_contigs=list()
    # make list of the phage sequences detected -- don't actually need
    phage_contigs=list()

    x=list(SeqIO.parse(contigs,'fasta'))
    q=[]    # list of bacterial contigs with hits to phage
    for i in phage_hits:
        b,_=i
        q.append(b)
      
    ## CP - needs to be rewritten to return the file names not the list of ids
    for i in x:
        if i.id not in q:
            bacterial_contigs.append(i.id)
        else:
            phage_contigs.append(i.id)
    return bacterial_contigs,phage_contigs

def compute_coverage(bacterial_contigs,phage_sequences,*fastqs):
    # inputs for this are the assemblies for the bacteria, phage sequences and raw data
    
    # compute coverage for bacterial contigs
    command=path_bbmap+'bbmap.sh ref='+bacterial_contigs+' in1='+fastqs[0]+' in2='+fastqs[1]+' out='+path_data+'test.sam basecov='+path_data+'basecov.try'
    print(command)
    #os.system(command)

    b_coverages=dict()
    with open(path_data+'basecov.try','r') as f:
        lines=f.readlines()
    for i in lines[1:]:
        x=i.strip().split('\t')
        if x[0] not in b_coverages.keys():
            b_coverages[x[0]]=[]
            print(x[0])
        b_coverages[x[0]].append(int(x[2]))

    # calculate average coverages for bacteria    
    avg_b_coverages=dict()
    for i in b_coverages:
        avg_b_coverages[i]=sum(b_coverages[i])/len(b_coverages[i])
    
    # compute coverage for phage contigs
    command=path_bbmap+'bbmap.sh ref='+phage_sequences+' in1='+fastqs[0]+' in2='+fastqs[1]+' out='+path_data+'test.sam basecov='+path_data+'basecov.try'
    print(command)
    #os.system(command)

    p_coverages=dict()
    with open(path_data+'basecov.try','r') as f:
        lines=f.readlines()
    for i in lines[1:]:
        x=i.strip().split('\t')
        if x[0] not in p_coverages.keys():
            p_coverages[x[0]]=[]
        p_coverages[x[0]].append(int(x[2]))
        
    # remove phage values that are insignificant (not evenness of coverage)
    # threshold at 10%+300
    revised_phage_list=dict()
    for i in p_coverages:
        length=len(p_coverages[i])
        theshold=length//10+300
        if p_coverages[i].count(0)<theshold:
            revised_phage_list[i]=copy.deepcopy(p_coverages[i])
            p_coverages[i].clear()
    p_coverages.clear()
    
    # calculate average coverages for phages
    avg_p_coverages=dict()
    for i in revised_phage_list:
        avg_p_coverages[i]=sum(revised_phage_list[i])/len(revised_phage_list[i])

    # return bacteria and phage coverages
    return avg_b_coverages,avg_p_coverages

# calls
## all you need to declare here is the fastqs
assembly=process_raw_reads(*fastqs, platform='Illumina')
assmebly=trim_assembly(assembly+'/contigs.fasta')
b_sequence_file,p_sequence_file=categorize_assembled_contigs(assembly,phage_sequence_file)
b_cov,p_cov=compute_coverage(b_sequence_file,p_sequence_file,*fastqs)


