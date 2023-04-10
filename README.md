# BioinformaticAnalysisTools

## Introduction
A python command line toolkit for Bioinformatics and data science!

\<under development\>

- Author: Hua-nan ZHAO @Tsinghua University
- E-Mail: hermanzhaozzzz@gmail.com

## Installation
```shell
# supported platform: Linux / MacOS (intel & arm64) / WSL on Windows
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
- sometimes, pysam dependent foo.lib can be absent, just `brew install foo.lib` to fix it.
