#!/bin/bash

conda install -c conda-forge mamba

conda install -c bioconda miru-hero
conda install -c bioconda emboss
conda activate tb-profiler

mamba create -n tb-profiler -c conda-forge -c bioconda docxtpl tb-profiler=4.4.2

git clone https://github.com/phglab/MIRUReader.git
git clone https://github.com/matnguyen/SpoTyping.git
git clone https://github.com/dcouvin/SpolLineages.git

pip install -r requirements.txt