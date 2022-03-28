# codename: induce_me

## Getting Started
Clone Project from GitHub
```python
git clone https://github.com/putonti/induce_me.git
```

### Prerequisites

- Biopython must be installed
- bbmap, SPAdes, and BLAST+ must all be located in the same path
- You will need to change the versions of SPAdes and BLAST+ within the python script

## Command Options

* -o : Directory to store resulting files and name of sample as file prefix
* -R : Directory location of local R application
* -r : Directory location of supplemental R code file taylor_code.R
* -s : Directory location of software tools: SPAdes, bbmap, BLAST+
* -t : Number of processors to use
* -n : Threshold for phage coverages
* -f : Reference file of phage sequences
* --version

Single or Paired-End Read Inputs:
* -i : Single read file
* -p : Paired-end read files. List both read files with a space between


### Example Run with Paired-End Reads:
```python
python3 induction_script_with_args.py -p phage_reference_file [read file options] -R Rscript_path -r R_code_path -s path_to_software_tools -o output_path_and_sample_name
```

## Further Details:
A manuscript describing further details is currently in preparation.

## Authors

* Taylor Miller-Ensminger
* Genevieve Johnson
* Catherine Putonti
