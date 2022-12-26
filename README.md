# PIE - Prophage Induction Estimator

## Getting Started
### Option 1:
Clone Project from GitHub
```python
git clone https://github.com/putonti/PIE.git
```
Move input fastq files, as well as the phages reference file, to the inputFiles folder prior to building the docker image.
```python
sudo docker build --tag pie:latest PIE
```
```python
sudo docker run -v ~/pathToLocalFolder/PIE:/pieOutputFolderName -i -t pie
```
### Option 2:
Docker Hub Link: https://hub.docker.com/repository/docker/genevievej16/pie
```python
docker pull genevievej16/pie:latest
```
Move input fastq files, as well as the phages reference file, to your local designated output folder prior to running the docker image. While in the docker you will access the input files from the pieOutputFolderName folder rather than from the inputFiles folder.
```python
sudo docker run -v ~/pathToLocalFolder/PIE:/pieOutputFolderName -i -t genevievej16/pie
```

### Prerequisites

Docker is the only prerequisite for this program to run, all other dependencies are handled by the Dockerfile. If any section of the program causes an error or is unable to run, check that you have enough memory in your Docker resources.

## Command Options

* -o : Directory to store resulting files and name of sample as file prefix (required)
* -s : Name of sample for output file labels (required)
* -t : Number of processors to use (default=4)
* -n : Threshold for phage coverages (default=0.99)
* -f : Reference file of phage sequences (required)
* --version

Single or Paired-End Read Inputs:
* -i : Single read file
* -p : Paired-end read files. List both read files with a space between

### Example Run with Paired-End Reads:
```python
python3.9 runPIE.py -f /inputFiles/phage_reference_file.fasta -p /inputFiles/R1.fastq /inputFiles/R2.fastq -s sample_name -o pieOutputFolderName/sample_output
```
### Example Run with Single Read:
```python
python3.9 runPIE.py -f /inputFiles/phage_reference_file.fasta -i /inputFiles/single_read_file.fastq -s sample_name -o pieOutputFolderName/sample_output
```

If you would prefer to run your own assembly outside of PIE and skip the assembly step, input assembled contig files along with single or paired-end read files:
* -a : Assembled contigs file

### Example Run with Assembled Contigs:
```python
python3.9 runPIE.py -f /inputFiles/phage_reference_file.fasta -a /inputFiles/contigs.fasta -p /inputFiles/R1.fastq /inputFiles/R2.fastq -s sample_name -o pieOutputFolderName/sample_output
```

## Test Data and Example Output
Paired-end reads of a small bacterial community is used for the test data.
The phage reference file is already included in the testFiles folder and contains 10 sample phage sequences ("test_phage_reference.fasta"). We have also included a larger phage dataset, which includes 34, of which the 10 are a subset ("34_phage.fasta").
### Example Run of Test Data:
```python
python3.9 runPIE.py -f /testFiles/test_phage_reference.fasta -p /testFiles/test_data_R1.fastq /testFiles/test_data_R2.fastq -s test_sample -o pieOutputFolderName/test_sample_output
```

## Further Details:
A manuscript describing further details is currently in preparation.

## Authors

* Taylor Miller-Ensminger
* Genevieve Johnson
* Catherine Putonti
* Swarnali Banerjee
