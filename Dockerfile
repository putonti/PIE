FROM ubuntu

#Coverage: bbmap
#Assemblers: spades
#Homology:blast+

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y curl vim csh gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget build-essential ncbi-blast+ cmake r-base default-jre
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt install -y python3.9
RUN apt-get install -y python3-pip

ADD BBMap_38.94.tar.gz bbmap
ADD runPIE.py runPIE.py 
ADD calculatePIE.R calculatePIE.R
ADD SPAdes-3.15.3-Linux.tar.gz spades
ADD /inputFiles/ /inputFiles/
ADD /testFiles/ /testFiles/
RUN pip install --upgrade pip
RUN apt-get install -y python3.9-distutils
RUN python3.9 -m pip install biopython


ENV PATH /bbmap/bbmap:/spades/SPAdes-3.15.3-Linux/bin:/blast:$PATH
RUN echo $PATH

#Example Run Command
#CMD ["python3.9", "runPIE.py", "-f", "/inputFiles/phage_reference_file.fasta", "-p" "/inputFiles/R1.fastq /inputFiles/R2.fastq", "-s" "sample_name", "-o", "/output_path/"]
