# battenberg_generate_probloci

- [1. Instructions for doing test on SNPs, running "PVL" filters and binomial filters on farm](#1-instructions-for-doing-test-on-snps-running-pvl-filters-and-binomial-filters-on-farm)
  - [Running PVL filters](#running-pvl-filters)
  - [Running binomial filters](#running-binomial-filters)
- [2. Python aggregation implementation](#2-python-aggregation-implementation)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Make a virtualenv](#make-a-virtualenv)
    - [Install parse_prob_loci_stats](#install-parse_prob_loci_stats)
  - [Running](#running)
  - [Usage](#usage)
- [Original contents of DW mail](#original-contents-of-dw-mail)
- [LICENCE](#licence)

## Overview

There are two main steps to generating the Probloci file.

1. Generating a list of tested SNPs using filters
2. Aggregating the information to generate `probloci.txt` itself

## 1. Instructions for doing tests on SNPs, running "PVL" filters and binomial filters on farm

### Running PVL filters

run [`bsub_PVL.sh`](scripts/bsub_PVL.sh) on farm, point to directory containing MutCounts files (example in example/sample_input)
[`R/Farm_PVL.R`](R/Farm_PVL.R) will be run over each MutCount file in input directory

point to output directory (example in example/sample_output)

Generates:

- `*_PVL_results_long.csv`  which contains Depth, BAF, and filter pass/fail for each SNP and each sample
- `*_PVL_results_summary.csv` which contains each SNP, and how many samples passed the Depth and BAF filters

submits to compute farm queue, runs very quickly

`./bsub_PVL.sh #input_DIR #output_DI #Path to R directory`

e.g.:

`./bsub_PVL.sh ./example/sample_input ./example/sample_output ~/battenberg_generate_probloci/R`

### Running binomial filters

run [`bsub_farm_binomial.sh`](scripts/bsub_farm_binomial.sh) on farm,  point to directory containing MutCounts files (example in example/sample_input)
[`Farm_binomal.R`](R/Farm_binomal.R) will be run over each MutCount file in input directory

point to output directory (example in example/sample_output)

Generates:

- `*_binomial_results_long.csv`  which contains Depth, BAF, 98% confidence interval of BAF, and filter pass/fail for each SNP and each sample
- `*_binomial_results_summary.csv` which contains each SNP, and how many samples passed the binomial filter

submits to long queue, on MutCount files of 100k SNPs each, over 98 samples runs ~12 hours

`./bsub_PVL.sh #input_DIR #output_DI #Path to R directory`

e.g.:

`./bsub_PVL.sh ./example/sample_input ./example/sample_output ~/battenberg_generate_probloci/R`

## 2. Python aggregation implementation

[Using DW's email](#original-contents-of-dw-mail) we transcribed the [original R](R/legacy) steps into python.
The script `parse_prob_loci_stats.py` takes the output of TB's code and performs the same functions as the provided R code.

### Requirements

- Python 3.7
  - pysam
  - numpy
  - matplotlib

### Installation

#### Make a virtualenv

```bash
python3.7 -m venv /path/to/create/venv
```

#### Install parse_prob_loci_stats

```bash
# Activate venv
source /path/to/create/venv/bin/activate
#change dir into battenberg_generate_probloci/Python
cd battenberg_generate_probloci/Python
python setup.py install
#exit venv
deactivate
```

### Running

```bash
# Activate venv
source /path/to/create/venv/bin/activate
#run script
python parse_prob_loci_stats.py 
```

```bash
usage: parse_prob_loci_stats.py [-h] [-v] -b binomial_.csv -d depth.csv -r
                                path/to/reference.fa -i
                                path/to/ignore_contigs.txt -p
                                /path/to/previous/bad_loci.txt
                                [-o path/to/write/output] [-c 93] [-e 10]
                                [-a 93] [-n 98]
parse_prob_loci_stats.py: error: the following arguments are required: -b/--binomial-file-pattern, -d/--depth-file-pattern, -r/--reference, -i/--ignore-contigs, -p/--previous-bad-loci
```

```bash
#exit venv
deactivate
```

### Usage

```bash
python parse_prob_loci_stats.py -h
usage: parse_prob_loci_stats.py [-h] [-v] -b binomial_.csv -d depth.csv -r
                                path/to/reference.fa -i
                                path/to/ignore_contigs.txt -p
                                /path/to/previous/bad_loci.txt
                                [-o path/to/write/output] [-c 93] [-e 10]
                                [-a 93] [-n 98]

Parse prob loci stats files and output files required for running prob loci
generation code.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version information
  -b binomial_.csv, --binomial-file-pattern binomial_.csv
                        Path to and name pattern of binomial csv input files.
                        Use %s to mark the contig name and %* to mark any
                        number
  -d depth.csv, --depth-file-pattern depth.csv
                        Path to and name pattern of depth and baf csv input
                        files. Use %s to mark the contig name and %* to mark
                        any number
  -r path/to/reference.fa, --reference path/to/reference.fa
                        Path to reference fasta file with associated fasta
                        index.
  -i path/to/ignore_contigs.txt, --ignore-contigs path/to/ignore_contigs.txt
                        Path to file containing list of contigs to ignore. One
                        contig per line
  -p /path/to/previous/bad_loci.txt, --previous-bad-loci /path/to/previous/bad_loci.txt
                        Path to previous bad_loci file
  -o path/to/write/output, --outdir path/to/write/output
                        Path to write output files
  -c 93, --binomial_cutoff 93
                        Cutoff for binomial test. Only results equal to or
                        above this value will be inculded. [default: 93]
  -e 10, --depth_cutoff 10
                        Cutoff for depth test. Only results equal to or above
                        this value will be inculded. [default: 93]
  -a 93, --depth_baf_cutoff 93
                        Cutoff for BAF test in depth file. Only results equal
                        to or above this value will be inculded. [default: 93]
  -n 98, --num-samples 98
                        Number of samples in the normal panel. [default: 98]
```

## Original contents of DW mail

Missing is the first step, in which we took allele counts of 1000G SNPs from a panel of normals (I think we had ~200 samples) and calculated for each sample, which SNPs were heterozygous and which were ‘bad’. The remaining SNPs are assumed to be homozygous.

The categorisation of SNPs was carried out in 2 different ways. The first (called ‘PVL’ because it was developed by Peter) used BAF thresholds, which I think were 0.3-0.7 for heterozygous SNPs and 0.1-0.3 / 0.7-0.9 for bad SNPs. We also used a minimum coverage threshold of 20. The second method (which I developed) used a binomial test of whether BAFs were significantly different from 0.5, 0.01 and 0.99 (I think these were the values used for homozygous SNPs). This method doesn’t need a coverage threshold, since coverage is incorporated into the binomial test.

The various ‘aggregate’ scripts pull together the badloci information across all samples and calculate the fraction of bad samples at each locus.

The ‘SummariseBadLoci’ scripts summarise the results and write lists of bad loci, using a threshold of 5% bad samples.

We found that there were some regions, such as the HLA locus, where few individual SNPs exceeded the 5% threshold, but that across the region a large number of samples had some dodgy SNPs (different SNPs in each samples). We therefore segmented the fraction of samples with bad SNPs to find these regions using the script ‘SegmentBadFraction.R’. ‘SummariseSegmentBadFraction.R’ summarises the results with different gamma parameters.

‘CompareBadLociAndRegions.R’ combines the loci identified from analysis of individual loci and the bad regions, using a gamma threshold for bad regions of 100 (250 on chr 6, so that the whole of the HLA region is excluded).
