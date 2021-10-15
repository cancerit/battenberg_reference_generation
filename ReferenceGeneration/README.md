# Battenberg GRCh38 reference generation

- [Scope](#scope)
- [Background](#background)
  - [Directories](#directories)
  - [Single files](#single-files)
  - [GRCh37 reference generation sources](#grch37-reference-generation-sources)
  - [File generation dependency](#file-generation-dependency)
- [1. Impute file generation](#1-impute-file-generation)
  - [1a. IMPUTE2 reference files](#1a-impute2-reference-files)
  - [1b. Genetic map formatting](#1b-genetic-map-formatting)
  - [1c. impute_info.txt](#1c-impute_infotxt)
- [2. battenberg_wgs_gc_correction_1000g file generation](#2-battenberg_wgs_gc_correction_1000g-file-generation)
- [3. 1000genomesloci file generation](#3-1000genomesloci-file-generation)
- [4. ignore_contigs.txt](#4-ignore_contigstxt)

## Scope

This directory contains the scripts to generate GRCh38 references required by Battenberg, excluding probLoci which is detailed in [the probloci directory](../ProblociGeneration)

## Background

There are 5 type of reference required by Battenberg

``` bash
ls /<ref_path.GRCh37>/battenberg/
1000genomesloci
battenberg_wgs_gc_correction_1000g
ignore_contigs.txt
impute
probloci.txt
```

### Directories

- 1000genomesloci
  - 1000 genomes allele and loci data

    ``` bash
    ls /<ref_path.GRCh37>/battenberg/1000genomesloci/ | sed 's/chr[[:digit:]]*.txt/chr<chr1:23>.txt/' | sort | uniq -c
     23 1000genomesAlleles2012_chr<chr1:23>.txt
     23 1000genomesloci2012_chr<chr1:23>.txt
    ```

- battenberg_wgs_gc_correction_1000g
  - GC correction files

  ``` bash
  /<ref_path.GRCh37>/battenberg/battenberg_wgs_gc_correction_1000g/ | sed 's/chr_[[:digit:]]*.txt/chr_<chr1:23>.txt/' | sort | uniq -c
   23 1000_genomes_GC_corr_chr_<chr1:23>.txt.gz
  ```

- impute
  - Impute reference files

  ``` bash
  ls /<ref_path.GRCh37>/battenberg/impute/ | sed 's/chr[[:digit:]]*_impute/chr<chr1:22>_impute/' | sed 's/chr[[:digit:]]*_combined/chr<chr1:22>_combined/' |sort | uniq -c
   22 ALL_1000G_phase1integrated_v3_chr<chr1:22>_impute.hap.gz
   22 ALL_1000G_phase1integrated_v3_chr<chr1:22>_impute.legend
    1 ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.hap.gz
    1 ALL_1000G_phase1integrated_v3_chrX_nonPAR_impute.legend
    1 ALL_1000G_phase1integrated_v3_chrX_PAR1_impute.hap.gz
    1 ALL_1000G_phase1integrated_v3_chrX_PAR1_impute.legend
    1 ALL_1000G_phase1integrated_v3_chrX_PAR2_impute.hap.gz
    1 ALL_1000G_phase1integrated_v3_chrX_PAR2_impute.legend
   22 genetic_map_chr<chr1:22>_combined_b37.txt
    1 genetic_map_chrX_nonPAR_combined_b37.txt
    1 genetic_map_chrX_PAR1_combined_b37.txt
    1 genetic_map_chrX_PAR2_combined_b37.txt
    1 impute_info.txt
  ```

### Single files

- ignore_contigs.txt
  - Non autosomal or X/Y contigs
- probloci.txt
  - Low quality loci to be ignored

### GRCh37 reference generation sources

[Prerequesites](https://github.com/cancerit/cgpBattenberg#prerequisites)

|File/File set                      |Source    |Source detail                                                                                    |
|-----------------------------------|---------|-------------------------------------------------------------------------------------------------|
|1000genomesloci                  |Script    |[download_generate_bberg_ref_files.pl](https://github.com/cancerit/cgpBattenberg/blob/master/perl/bin/download_generate_bberg_ref_files.pl)  |
|impute                              |Script    |[download_generate_bberg_ref_files.pl](https://github.com/cancerit/cgpBattenberg/blob/master/perl/bin/download_generate_bberg_ref_files.pl)  |
|battenberg_wgs_gc_correction_1000g  |Download |[required-reference-files](https://github.com/Wedge-lab/battenberg#required-reference-files)                                  |
|ignore_contigs.txt                  |User Generated      |List of contigs to ignore during analysis                                                                    |
|probloci.txt                      |Download |[probloci.txt](https://github.com/cancerit/cgpBattenberg/blob/master/perl/share/battenberg/probloci.txt.gz)         |

### File generation dependency

``` bash
├── 1. impute
│   ├── 2. battenberg_wgs_gc_correction_1000g
│   └── 3. 1000genomesloci
├── 4. ignore_contigs.txt
└── 5. probloci.txt
```

## 1. Impute file generation

### 1a. IMPUTE2 reference files

- i. Output format
  - 1) IMPUTE2 file formats
    - a) legend:  [Impute input options](https://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-l)
    - b) hap: [Impute known haps](https://mathgen.stats.ox.ac.uk/impute/input_file_options.html#-known_haps_g)
  - 2) File split
    - a) legend (n=25): chr1-22, chrX_PAR1, chrX_PAR2 & chrX_nonPAR
    - b) hap (n=25): chr1-22, chrX_PAR1, chrX_PAR2 & chrX_nonPAR
  - 3) File name format
    - a) legend: ALL_1000G_phase1integrated_v3_\<split>_impute.legend
    - b) hap: ALL_1000G_phase1integrated_v3_\<split>_impute.hap.gz
- ii. File generation
  - 1) [01a_Impute_file_generation.sh](01a_Impute_file_generation.sh)

### 1b. Genetic map formatting

- i. Output format
  - 1) 3 space delimited columns with header
    - a) position
    - b) COMBINED_rate(cM/Mb)
    - c) Genetic_Map(cM)
  - 2) File split (n=25)
          - a) chr1-22, chrX_PAR1, chrX_PAR2 & chrX_nonPAR
  - 3) File name format
    - genetic_map_\<split>_combined_b37.txt
- ii. File generation
  - 1) [01b_genetic_map_generation.sh](01b_genetic_map_generation.sh)

### 1c. impute_info.txt
  
- i. Output format
  - 1) 7 tabbed delimited columns with no header
    - a) Chr
    - b) legend path
    - c) map path
    - d) hap.gz path
    - e) IMPUTE start position
    - f) IMPUTE end position
    - g) IMPUTE male (1/0 boolean; regions to impute in males)
- ii. File example
  - 1) [01c_impute_info.txt](01c_impute_info.txt)

## 2. battenberg_wgs_gc_correction_1000g file generation

- a. Output format (GC_corr)
  - i. ASCAT GC correction file format
  - ii. File split (n=23)
    - 1) chr1-22 & chrX
  - iii. File name format
    - 1) 1000_genomes_GC_corr_chr_\<split>.txt.gz
- b. File generation
  - i. Process
    - 1) Input
      - a) SnpPositions.tsv (tabbed delimited 3 column file; ID)
    - 2) Code
      - a) [Ascat convert SNPs to SNP GC Corrections](https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv)
  - ii. Script
    - 02a_gcCorrection_file_generation.sh

## 3. 1000genomesloci file generation

- a. Output format
  1. 1000genomesAlleles (3 col tab with header)
       - a. position
       - b. a0 (numeric alleles A=1, C=2, G=3, T=4)
       - c. a1 (numeric alleles A=1, C=2, G=3, T=4)
  2. 1000genomesloci (2 col tab no header)
       - a) Chr
       - b ) position
  3. File split
      - 1000genomesAlleles (n=23): chr1-22 & chrX
      - 1000genomesloci (n=23): chr1-22 & chrX
  4. File name format
      - 1000genomesAlleles: 1000genomesAlleles2012_chr\<split>.txt
      - 1000genomesloci: 1000genomesloci2012_chr\<split>.txt
- b. File generation
  1. Inputs
      - IMPUTE2 legend files
  2. Script
       - 03a_1000genomesloci_file_generation.sh

## 4. ignore_contigs.txt

- a. Format
  - i. No header
  - ii. Single column
    - 1) Contig name
- b. File generation
  - i. Process
    1. Input - refence genome index file (-.fa.fai)
    2. Pseudocode
        - a) Read in from 24th line
        - b) Save first column (contig name) to ignore_contigs.txt
  - ii. Check against current reference

    ``` bash
    # GRCh37d5 
    diff <(cat /<ref_path.GRCh37>/battenberg/ignore_contigs.txt) <(awk 'FNR> 23 {print $1}' /<ref_path.GRCh37>/genome.fa.fai) | diffstat
    0 files changed    
    # GRCh38  
    diff <(cat /<ref_path.GRCh38>/battenberg/ignore_contigs.txt) <(awk 'FNR> 23 {print $1}'  /<ref_path.GRCh38>/genome.fa.fai) | diffstat
    0 files changed
    ```

  - iii. Make ignore_contigs.txt

    ``` bash
    awk 'FNR> 23 {print $1}'  /<ref_path.GRCh38>/genome.fa.fai > GRCh38/battenberg/ignore_contigs.txt 
    ```

## 5. probloci.txt

### Format

- Header
- 2 column (tab separated)
  1. Chr
  1. Pos

### File generation

[probloci](../ProblociGeneration)
