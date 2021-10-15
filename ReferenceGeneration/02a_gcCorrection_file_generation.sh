##########LICENCE##########
# Copyright (c) 2014-2021 Genome Research Ltd.
# 
# Author: Cancer Genome Project cgpit@sanger.ac.uk
# 
# This file is part of battenberg.
# 
# battenberg is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads 'Copyright (c) 2005, 2007-
# 2009, 2011-2012' should be interpreted as being identical to a statement that
# reads 'Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012' and a copyright
# statement that reads "Copyright (c) 2005-2012' should be interpreted as being
# identical to a statement that reads 'Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012'."
##########LICENCE##########

# Aim
# Make GRCh38 reference inputs for Battenberg
# 2. battenberg_wgs_gc_correction_1000g file generation

# Task overview
# Format IMPUTE2 legend files as input to https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv

# Output
# 23 files <split> over chr1:22 & chrX (note double chr in filename)
# 1000_genomes_GC_corr_chr_<1:22 & X>.txt.gz

# Sources
## Conversion script
# https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv

## Reference genome fasta
<genome.fa> = ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Method

### Make SnpPositions.tsv file for each  IMPUTE2 legend file
cat <split>.phased.legend | awk 'NR>1{split($1,a,":"); print $1"\t"a[1]"\t"a[2]}' > <split>_SnpPositions.tsv

### Run over all *.phased.legend.gz files
ascatSnpPanelGcCorrections.pl  <genome.fa> <split>_SnpPositions.tsv > <ref_path>/battenberg_wgs_gc_correction_1000g/1000_genomes_GC_corr_chr_<split>.txt

### ChrX - combine ChrX SnpPositions sections
cat chrX_PAR1_SnpPositions.tsv chrX_nonPAR_SnpPositions.tsv chrX_PAR2_SnpPositions.tsv > <ref_path>/battenberg_wgs_gc_correction_1000g/chrX_SnpPositions.tsv

### Run over combined ChrX
ascatSnpPanelGcCorrections.pl  <genome.fa> chrX_SnpPositions.tsv > <ref_path>/battenberg_wgs_gc_correction_1000g/1000_genomes_GC_corr_chr_chrX.txt

### Gzip files
gzip <ref_path>/battenberg_wgs_gc_correction_1000g/1000_genomes_GC_corr_chr_*.txt
