# BioAT
the Bioinformatic Analysis Tools
## About bioat
Bref:  <u>**bioat**, a python **package** & **command line toolkit** for Bioinformatics and data science!</u>

Home: https://github.com/hermanzhaozzzz/bioat

License: Apache-2.0

| Name | Downloads | Version | Platforms | Test |
| --- | --- | --- | --- | --- |
|[![Conda Recipe](https://img.shields.io/badge/recipe-bioat-green.svg)](https://anaconda.org/conda-forge/bioat) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/bioat.svg)](https://anaconda.org/conda-forge/bioat) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/bioat.svg)](https://anaconda.org/conda-forge/bioat) | [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/bioat.svg)](https://anaconda.org/conda-forge/bioat) | [![Azure Pipelines](https://dev.azure.com/conda-forge/feedstock-builds/_apis/build/status/bioat-feedstock?branchName=main)](https://dev.azure.com/conda-forge/feedstock-builds/_build/latest?definitionId=23719&branchName=main) |

## About author
Name: [赵华男 | Huanan Herman Zhao](https://scholar.google.com/citations?user=ojSVoWQAAAAJ&hl=en) |  EMail: hermanzhaozzzz AT gmail.com | [Zhihu](https://www.zhihu.com/people/hymanzhaozzzz) | [BLOG](http://zhaohuanan.cc)

## Installation
```shell
# supported platform: Linux / MacOS (intel & arm64) / Windows
pip install --upgrade bioat
# use conda
conda install -c conda-forge bioat
```

## Usage
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
[circos plot](docs/demo_circos-plot.ipynb)

## Doc

See [Doc](https://bioat.readthedocs.io/en/latest/)

