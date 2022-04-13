# PIE - Prophage Induction Estimator

## Getting Started
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
python3 runPIE.py -f inputFiles/phage_reference_file.fasta -p inputFiles/R1.fastq inputFiles/R2.fastq -s sample_name -o pieOutputFolderName/sample_output
```

## Further Details:
A manuscript describing further details is currently in preparation.

## Authors

* Taylor Miller-Ensminger
* Genevieve Johnson
* Catherine Putonti
* Swarnali Banerjee
