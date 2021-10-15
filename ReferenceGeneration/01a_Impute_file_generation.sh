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
# 1. IMPUTE file generation

# Sources
## Conversion script
	# https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#scripts
	# replace hg19 pos for GRCh38 > vcf2impute_legend_haps.GRCh38.pl
	# my $chrX_par1_start_hg19 = 10001;
	# my $chrX_par1_end_hg19 = 2781479;
	# my $chrX_par2_start_hg19 = 155701383;
	# my $chrX_par2_end_hg19 = 156030895;

## Phased VCFs - GRCh38 remaps of 1000 genomes
	# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/
	# Download into <20181203_biallelic_SNV_path>

## Genetic map
	# From Broad - https://data.broadinstitute.org/alkesgroup/Eagle/#x1-30002

## output dir 
	#<ref_path>

# Method

## 1. IMPUTE - Legend/hap file generation
### Phased VCF clean-up
	# Current filter is to remove singleton variants (bcftools -c 2)
	# Impute script prefixes "chr" on to outputs, so this command renames the chroms to GRCh37 style (1,2,3..) from GRCh38 (chr1,ch2,chr3...)
	bgzip -c <(bcftools view --no-version -h -c 2 .<20181203_biallelic_SNV_path>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=chr\([0-9XY]\)/##contig=<ID=\1/'; bcftools view --no-version -H -c 2 <20181203_biallelic_SNV_path>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz | sed 's/^chr//') > <20181203_biallelic_SNV_path.chrFix>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.chrFix.gz

### Autosomes
	# Run impute script
	mkdir -p 0<ref_path>/impute/

	for chr in {1..22}; do
		echo ./vcf2impute_legend_haps.pl -vcf <20181203_biallelic_SNV_path.chrFix>/ALL.chr$chr.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.chrFix.gz -leghap <ref_path>/impute/ALL.chr$chr.shapeit2_integrated_v1a.GRCh38.20181129.phased -chr $chr;
	Done

### ChrX
	# Needs to be split over par and non-par regions
	# X: 10001 - 2781479 (PAR1)
	
	./vcf2impute_legend_haps.GRch38.pl -vcf <20181203_biallelic_SNV_path.chrFix>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.chrFix.gz -leghap 0<ref_path>/impute//ALL.chrX_PAR1.shapeit2_integrated_v1a.GRCh38.20181129.phased -chr X -start 10001 -end 2781479

	# X: 155701383 - 156030895 (PAR2)
	
	./vcf2impute_legend_haps.GRch38.pl -vcf <20181203_biallelic_SNV_path.chrFix>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.chrFix.gz -leghap 0<ref_path>/impute//ALL.chrX_PAR2.shapeit2_integrated_v1a.GRCh38.20181129.phased -chr X -start 155701383 -end 156030895
	
	# X 2781480 - 155701382 (non-PAR) - no data with above data source (20181203_biallelic_SNV)
	
./vcf2impute_legend_haps.GRch38.pl -vcf <20181203_biallelic_SNV_path.chrFix>/ALL.chrX.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.chrFix.gz -leghap 0<ref_path>/impute//ALL.chrX_nonPAR.shapeit2_integrated_v1a.GRCh38.20181129.phased -chr X -start 2781480 -end 155701382 -chrX_nonpar
