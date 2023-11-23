# bioat 
the Bioinformatic Analysis Tools

![macos](https://github.com/hermanzhaozzzz/bioat/actions/workflows/macos.yml/badge.svg)
![linux](https://github.com/hermanzhaozzzz/bioat/actions/workflows/linux.yml/badge.svg)
![windows](https://github.com/hermanzhaozzzz/bioat/actions/workflows/windows.yml/badge.svg)

## about author

> author: [赵华男 | ZHAO Hua-nan](https://scholar.google.com/citations?user=ojSVoWQAAAAJ&hl=en)
>
> email: hermanzhaozzzz@gmail.com
>
> [Zhihu](https://www.zhihu.com/people/hymanzhaozzzz) | [BLOG](http://zhaohuanan.cc)

## Introduction
A python **package** & **command line toolkit** for Bioinformatics and data science!

**\<under development\>!!**

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
[circos plot](docs/demo_circos-plot.ipynb)