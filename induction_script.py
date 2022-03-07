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


def process_raw_reads(*fastqs,platform='Illumina'):
    # expects 1 or 2 fastq files
    # default Illumina
    
    # trim reads
    ## NEED CODE TO RUN BBDUK
    
    # assemble reads
    ## NEED CODE TO RUN METASPADES
    
    ## NEED TO FIGURE OUT WHERE IT'S GOING TO BE SAVED


def categorize_assembled_contigs()
    # create blast database of the predicted phage sequences
    phage_name=phage_sequences_file[:phage_sequences_file.find('.')]
    command='/Users/taylormiller-ensminger/Desktop/ncbi-blast-2.11.0+/bin/makeblastdb -in '+phage_sequences_path+phage_sequences_file+' -out '+phage_name+' -title '+phage_name+' -dbtype nucl'
    os.system(command)  # uncomment to run
    
    # blast assembly against phage database to figure out which ones are phage
    command='/Users/taylormiller-ensminger/Desktop/ncbi-blast-2.11.0+/bin/blastn -query '+bact_contigs_path+bact_contigs_file+' -db '+ phage_name+' -max_target_seqs 1 -outfmt="10 qseqid sseqid qcovs pident length evalue bitscore" -out /Users/taylormiller-ensminger/Desktop/blastn_'+phage_name+'.csv'
    os.system(command)
    
    # read in blast results
    with open('/Users/taylormiller-ensminger/Desktop/blastn_'+phage_name+'.csv','r') as f:
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

    x=list(SeqIO.parse(bact_contigs_path+bact_contigs_file,'fasta'))
    q=[]    # list of bacterial contigs with hits to phage
    for i in phage_hits:
        b,_=i
        q.append(b)
       
    for i in x:
        if i.id not in q:
            bacterial_contigs.append(i.id)
        else:
            phage_contigs.append(i.id)
    return bacterial_contigs,phage_contigs


def compute_coverage(bacterial_contigs,phage_sequences_file,*fastqs):
    # inputs for this are the assemblies for the bacteria, phage sequences and raw data
    
    # read in metaSPAdes assembly, only select the contigs that are identified in bacterial_contigs, map to them
    
    # read in phage_sequences, map fastqs

   
#### HOT MESS ####
    
# now need bbmap scores - bacteria
with open(bbmap_path_bact+bbmap_file_bact,'r') as f:
    reader = csv.DictReader(f,delimiter='\t')
    temp= list(reader)
bbmap_results_bact=dict()
for i in temp:
    bbmap_results_bact[i['#ID']]=float(i['Avg_fold'])
#print(bbmap_results_bact)    


bacterial_cov=dict()
for i in bacterial_contigs:
    #print(i)
    bacterial_cov[i]=bbmap_results_bact[i]
    #print(bacterial_cov[i])

# phage sequences
x=list(SeqIO.parse(phage_sequences_path+phage_sequences_file,'fasta'))
phage_seqs=[i.id for i in x]

# now need bbmap scores - phage
with open(bbmap_path_phage+bbmap_file_phage,'r') as f:
    reader = csv.DictReader(f,delimiter='\t')
    temp= list(reader)
bbmap_results_phage=dict()
for i in temp:
    bbmap_results_phage[i['#ID']]=float(i['Avg_fold'])
print("hi")

phage_cov=dict()
#print(phage_seqs)
for i in phage_seqs:
    print("hi")
    #I had to add this statement
    #if i in bbmap_results_phage.keys():
    phage_cov[i]=bbmap_results_phage[i]


outf=open(output_file,'w')
outf.write('BACTERIA\n')
for i in bacterial_cov:
    outf.write(i+','+str(bacterial_cov[i])+'\n')
outf.write('PHAGE\n')
for i in phage_cov:
    outf.write(i+','+str(phage_cov[i])+'\n')
outf.close()
