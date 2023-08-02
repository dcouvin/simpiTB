# A Python version of simpiTB has been developed to simplify the run

## Requirements
You will need conda. (if not go on : https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
To check if Conda is installed on your machine, you can run the following command:
```
conda --version
```
Java version 6 (or later) programming language must be available in your system to run SpolLineages.
To check if Java is installed on your machine, you can run the following command:
```bash
java -version
```

python (python3) programming language must be available in your system to run simpiTB.
To check if python is installed on your machine, you can run the following command:
```bash
python --version
```
## Starting
First, launch the installer
```bash
sh installer.sh
```
Verify if the "tb-profiler" env is activate. 
If not :
```bash
conda activate tb-profiler
```
## Starting
You can use the following command to display the help
```bash
python3 simpiTB.py --help
```
Command example
```bash
python3 simpiTB.py examplefile.fasta -o outputfile
```
