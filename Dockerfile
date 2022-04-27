FROM ubuntu

#Coverage: bbmap
#Assemblers: spades
#Homology:blast+

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y curl vim csh python3.9 gcc g++ unzip make git bzip2 zlib1g-dev ncurses-dev wget python3-pip build-essential python-pkg-resources python-setuptools ncbi-blast+ cmake r-base default-jre
RUN pip install ipython 
ADD BBMap_38.94.tar.gz bbmap
ADD runPIE.py runPIE.py 
ADD calculatePIE.R calculatePIE.R
ADD SPAdes-3.15.3-Linux.tar.gz spades
ADD /inputFiles/ /inputFiles/
ADD /testFiles/ /testFiles/
RUN pip install --upgrade pip
RUN python3 -m pip install biopython


ENV PATH /bbmap/bbmap:/spades/SPAdes-3.15.3-Linux/bin:/blast:$PATH
RUN echo $PATH

#Example Run Command
#CMD ["python3.9", "runPIE.py", "-f", "phage_reference_file.fasta", "-p" "R1.fastq R2.fastq", "-s" "sample_name", "-o", "/output_path/"]
