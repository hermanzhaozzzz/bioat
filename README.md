# BioinformaticAnalysisTools

## Introduction
A python command line toolkit for Bioinformatics and data science!

\<under development\>

- Author: Hua-nan ZHAO @Tsinghua University
- E-Mail: hermanzhaozzzz@gmail.com

## Installation
```shell
# supported platform: Linux / MacOS / WSL on Windows
pip install --upgrade bioat
```


## usage
```shell
# list commands
bioat list
# check version
bioat version
# check information about bioat
bioat about

# example usage
bioat bam remove_clip --help
samtools view -h test_sorted_n.bam | bioat bam remove_clip | tail
```

## known trouble
- You should run bioat in a rosetta env while using Apple M Silicon-Arm64 (M1/M2) 
    because that one dependent package [`pysam`](https://github.com/pysam-developers/pysam)
    do not support this platform now. See below to use a rosetta env created by conda.
```shell
# create rosetta env
CONDA_SUBDIR=osx-64 conda create -n rosetta python
conda activate rosetta
# in rosetta env
python -c "import platform;print(platform.machine())"  # should print “x86_64”
conda env config vars set CONDA_SUBDIR=osx-64
conda deactivate      # need reactivate rosetta env to enable this var
conda activate rosetta
# now you can do this in this env
pip install --upgrade bioat
bioat version
```
