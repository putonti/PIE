import argparse
import os
import csv
from Bio import SeqIO
import copy
import subprocess


def msg(name=None):
    return '''taylor_code.py -p phage_reference_file [read file options] -R Rscript_path -r R_code_path -s path_to_software_tools -o output_path_and_sample_name'''


parser=argparse.ArgumentParser(usage=msg())
parser.add_argument('-o', '--output_path', action="store", metavar='<directory>', help='Directory to store resulting files and name of sample as file prefix (required)')
parser.add_argument('-R', '--Rscript_path', action="store", metavar='<directory>', help='Directory location of local R application (required)')
parser.add_argument('-r', '--R_code_path', action="store", metavar='<directory>', help='Directory location of supplemental R code file induction_R_code.R (required)')
parser.add_argument('-s', '--path_to_software_tools', action="store", metavar='<directory>', help='Directory location of software tools: SPAdes, bbmap, BLAST+ (required)')
parser.add_argument('-t', '--num_threads', action="store", metavar='<int>', type=int, help='Number of processors to use.')
parser.add_argument('-n', '--threshold', action="store", metavar='<int>', help='Threshold for phage coverages.')
parser.add_argument('-f', '--phage_reference_file', action="store", metavar='<filename>', help='Reference file of phage sequences (required)')
group=parser.add_mutually_exclusive_group()
group.add_argument('-i', '--single_read', action="store", metavar='<filename>', help='Single read file')
group.add_argument('-p', '--paired_end_reads', action="store", nargs=2, metavar=('<filename>', '<filename>'), help='Paired-end read files. List both read files with a space between')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args=parser.parse_args()
if args.paired_end_reads is None and args.single_read is None:
    parser.error('Reads must be provided for analysis.')
if args.phage_reference_file is None:
    parser.error('Phage reference file must be provided for analysis.')
if args.output_path is None:
    parser.error('Output path must be provided.')
if args.Rscript_path is None:
    parser.error('Path to local R application must be provided.')
if args.R_code_path is None:
    parser.error('Path to supplemental R code file must be provided.')
if args.path_to_software_tools is None:
    parser.error('Path to local software tools (SPAdes, bbmap, BLAST+) must be provided.')

path_of_tools = args.path_to_software_tools


def process_raw_reads(*fastqs, num_threads=4):
    if len(fastqs) == 1:
        bbduk_trim = path_of_tools+"bbmap/bbduk.sh -Xmx1G overwrite=t in="+fastqs[0]+" out="+\
                     fastqs[0].strip(".fastq")+"_trimmed.fastq"+" qtrim=rl ftl=15 ftr=135 maq=20 maxns=0 stats="\
                     +fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+\
                     "_read_qualTrimming.stats statscolumns=5 trimq=20"
        spades_assembly = path_of_tools+"SPAdes-3.15.3-Darwin/bin/metaspades.py --only-assembler -t "+str(num_threads)+\
            " -12 "+fastqs[0].strip(".fastq")+"_trimmed.fastq"+" -o "+\
            fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+"_assembly"
        os.system(bbduk_trim)
        os.system(spades_assembly)
        return fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+"_assembly"
    elif len(fastqs) == 2:
        bbduk_trim = path_of_tools+"bbmap/bbduk.sh -Xmx1G overwrite=t in1="+fastqs[0]+" in2="+fastqs[1]+\
                     " out1="+fastqs[0].strip(".fastq")+"_trimmed.fastq"+" out2="+\
                     fastqs[1].strip(".fastq")+"_trimmed.fastq qtrim=rl ftl=15 ftr=135 maq=20 maxns=0 " \
                     "stats="+fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+\
                     "_read_qualTrimming.stats statscolumns=5 trimq=20"
        spades_assembly = path_of_tools + "SPAdes-3.15.3-Darwin/bin/metaspades.py --only-assembler -t "+\
            str(num_threads)+" -1 "+fastqs[0].strip(".fastq")+"_trimmed.fastq"+" -2 "+\
            fastqs[1].strip(".fastq")+"_trimmed.fastq"+" -o "+\
            fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+"_assembly"
        os.system(bbduk_trim)
        os.system(spades_assembly)
        return fastqs[0].strip(".fastq").replace("_R1", "").replace("_R2", "")+"_assembly"
    else:
        print("Expected 1 or 2 fastq files. Please check fastq files and try again.")


def trim_assembly(contigs):
    # remove contigs that are less than 1000 bp
    c = list(SeqIO.parse(contigs, 'fasta'))
    file_name = contigs[:contigs.rfind(".")] + "_trimmed.fasta"
    with open(file_name, 'w') as f:
        for i in c:
            if len(str(i.seq)) >= 1000:
                f.write('>' + i.id + '\n' + str(i.seq) + '\n')
    return file_name


## GJ -- you may need to change the version # for blast
def categorize_assembled_contigs(contigs, phage_sequences):
    # create blast database of the predicted phage sequences
    phage_name = phage_sequences[:phage_sequences.rfind('.')]
    command = path_of_tools + 'ncbi-blast-2.12.0+/bin/makeblastdb -in ' + phage_sequences + ' -out ' + phage_name + ' -title ' + phage_name + ' -dbtype nucl'
    os.system(command)  # uncomment to run

    # blast assembly against phage database to figure out which ones are phage
    command = path_of_tools + 'ncbi-blast-2.12.0+/bin/blastn -query ' + contigs + ' -db ' + phage_name + ' -max_target_seqs 1 -outfmt="10 qseqid sseqid qcovs pident length evalue bitscore" -out ' + phage_name + '_blastn.csv'
    os.system(command)

    # read in blast results
    with open(phage_name + '_blastn.csv', 'r') as f:
        fieldnames = ['qseqid', 'sseqid', 'qcovs', 'pident', 'length', 'evalue', 'bitscore']
        reader = csv.DictReader(f, fieldnames=fieldnames)
        blast_results = list(reader)

    # grab top hit for each blast result
    top_hits = dict()
    for i in blast_results:
        pair = (i['qseqid'], i['sseqid'])
        if pair in top_hits.keys():
            q, p, b = top_hits[pair]
            if float(i['bitscore']) > float(b):
                top_hits[pair] = (i['qcovs'], i['pident'], i['bitscore'])
        else:
            top_hits[pair] = (i['qcovs'], i['pident'], i['bitscore'])

    # toss garbage -- qcov < 90%
    phage_hits = dict()
    for i in top_hits:
        q, p, b = top_hits[i]
        if float(q) >= 90:
            phage_hits[i] = top_hits[i]

    # make lists of the bacterial contigs that are really bacterial
    bacterial_contigs = list()
    # make list of the phage sequences detected -- don't actually need
    phage_contigs = list()

    x = list(SeqIO.parse(contigs, 'fasta'))
    q = []  # list of bacterial contigs with hits to phage
    for i in phage_hits:
        b, _ = i
        q.append(b)

    ## return the file of bacterial contigs
    out_b = contigs[:contigs.rfind('/')] + '/filtered_' + contigs[contigs.rfind('/') + 1:]

    with open(out_b, 'w') as f:
        for i in x:
            if i.id not in q:
                bacterial_contigs.append(i.id)
                f.write('>' + i.id + '\n' + str(i.seq) + '\n')
            else:
                phage_contigs.append(i.id)
    return out_b


def compute_coverage(bacterial_contigs, phage_sequences, *fastqs):
    # inputs for this are the assemblies for the bacteria, phage sequences and raw data

    path_data = bacterial_contigs[:bacterial_contigs.rfind('/') + 1]

    # compute coverage for bacterial contigs
    command = path_of_tools + 'bbmap/bbmap.sh -Xmx1G overwrite=t ref=' + bacterial_contigs + ' in1=' + fastqs[
        0] + ' in2=' + fastqs[1] + ' out=' + path_data + 'test.sam basecov=' + path_data + 'basecov.try'
    print(command)
    os.system(command)

    b_coverages = dict()
    with open(path_data + 'basecov.try', 'r') as f:
        lines = f.readlines()
    for i in lines[1:]:
        x = i.strip().split('\t')
        if x[0] not in b_coverages.keys():
            b_coverages[x[0]] = []
            print(x[0])
        b_coverages[x[0]].append(int(x[2]))

    # calculate average coverages for bacteria
    avg_b_names = []
    avg_b_coverages = []
    for i in b_coverages:
        avg_b_names.append(i)
        avg_b_coverages.append(sum(b_coverages[i]) / len(b_coverages[i]))

    # compute coverage for phage contigs
    # IMPORTANT NOTE, PHAGE NAMES MUST BE UNIQUE OTHERWISE BBMAP EXPLODES
    command = path_of_tools + 'bbmap/bbmap.sh -Xmx1G overwrite=t ref=' + phage_sequences + ' in1=' + fastqs[
        0] + ' in2=' + fastqs[1] + ' out=' + path_data + 'test.sam basecov=' + path_data + 'basecov.try'
    print(command)
    os.system(command)

    p_coverages = dict()
    with open(path_data + 'basecov.try', 'r') as f:
        lines = f.readlines()
    for i in lines[1:]:
        x = i.strip().split('\t')
        if x[0] not in p_coverages.keys():
            p_coverages[x[0]] = []
        p_coverages[x[0]].append(int(x[2]))

    # remove phage values that are insignificant (not evenness of coverage)
    # threshold at 10%+300
    revised_phage_list = dict()
    for i in p_coverages:
        length = len(p_coverages[i])
        threshold = length // 10 + 300
        if p_coverages[i].count(0) < threshold:
            revised_phage_list[i] = copy.deepcopy(p_coverages[i])
            p_coverages[i].clear()
    p_coverages.clear()

    # calculate average coverages for phages
    avg_p_names = []
    avg_p_coverages = []
    for i in revised_phage_list:
        avg_p_names.append(i)
        avg_p_coverages.append(sum(revised_phage_list[i]) / len(revised_phage_list[i]))

    # return bacteria and phage coverages
    return avg_b_names, avg_b_coverages, avg_p_names, avg_p_coverages


def write_out(sample_name, bact_name, bact_cov, phage_name, phage_cov):
    final_out = open(sample_name+"_bact_phage_coverages.txt", "w")
    final_out.write("bactMatches:"+str(bact_name).strip("[").strip("]").strip(" ").replace("'", "")+"\n")
    final_out.write("bactCovs:" + str(bact_cov).strip("[").strip("]").strip(" ").replace("'", "")+"\n")
    final_out.write("phageMatches:" + str(phage_name).strip("[").strip("]").strip(" ").replace("'", "")+"\n")
    final_out.write("phageCovs:" + str(phage_cov).strip("[").strip("]").strip(" ").replace("'", "")+"\n")
    final_out.close()


# calls
if args.num_threads is None:
    if args.paired_end_reads is not None:
        assembly = process_raw_reads(args.paired_end_reads[0], args.paired_end_reads[1])
    elif args.single_read is not None:
        assembly = process_raw_reads(args.single_read)
    else:
        parser.error('Reads must be provided for analysis.')
else:
    if args.paired_end_reads is not None:
        assembly = process_raw_reads(args.paired_end_reads[0], args.paired_end_reads[1], num_threads=args.num_threads)
    elif args.single_read is not None:
        assembly = process_raw_reads(args.single_read, num_threads=args.num_threads)
    else:
        parser.error('Reads must be provided for analysis.')

assembly = trim_assembly(assembly + '/contigs.fasta')
b_sequence_file = categorize_assembled_contigs(assembly, args.phage_reference_file)

if args.paired_end_reads is not None:
    b_name, b_cov, p_name, p_cov = compute_coverage(b_sequence_file, args.phage_reference_file, args.paired_end_reads[0], args.paired_end_reads[1])
elif args.single_read is not None:
    b_name, b_cov, p_name, p_cov = compute_coverage(b_sequence_file, args.phage_reference_file, args.single_read)
else:
    parser.error('Reads must be provided for analysis.')

write_out(args.output_path, b_name, b_cov, p_name, p_cov)

subprocess.call([args.Rscript_path, args.R_code_path, args.output_path+"_bact_phage_coverages.txt", "Density Plot", args.threshold, args.output_path])

