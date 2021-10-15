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
#	3. 1000genomesloci file generation

# Task overview
# Format and filter (single base subs only) IMPUTE2 legend files into 1000genomesloci files

# Output
	# 1. 1000genomesAlleles (3 col tab with header)
		# a. position
		# b. a0 (numeric alleles A=1, C=2, G=3, T=4)
		# c. a1 (numeric alleles A=1, C=2, G=3, T=4)
	# 2. 1000genomesloci (2 col tab no header)
		# ) Chr
		# b) position
	# 3. File split
		# 1) 1000genomesAlleles (n=23): chr1-22 & chrX 
		# 2) 1000genomesloci (n=23): chr1-22 & chrX
	# 4. File name format
		# 1) 1000genomesAlleles: 1000genomesAlleles2012_chr<split>.txt
		# 2) 1000genomesloci: 1000genomesloci2012_chr<split>.txt
	
# Inputs
# IMPUTE2 legend files  

# Method
## Make 1000genomesAllele files
### Autosomes
# Per autosome
	for i in $(seq 1 22 | sed 's:^:chr:g'); 
	do 
	cat <ref_path>/impute/ALL.v1a.shapeit2_integrated_"$i".GRCh38.20181129.phased.legend \
	| awk 'function alchk(c) {if(($c == "A") || ($c == "C") || ($c == "G") || ($c == "T")) return 1} {if(NR==1) {print} else {if(alchk(3) && alchk(4)) print}}' \
	| awk 'BEGIN{OFS="\t"} function g2n(c) {gsub("A","1",$c);gsub("C","2",$c);gsub("G","3",$c);gsub("T","4",$c)} {g2n(3) g2n(4)} {print $2,$3,$4}' \
	> <ref_path>/1000genomesloci/1000genomesAlleles2012_chr"$i".txt
	done

### ChrX
# ChrX only
	cat </impute/ALL.v1a.shapeit2_integrated_chrX_*.GRCh38.20181129.phased.legend \
	| sort -nk2 \
	| awk 'function alchk(c) {if(($c == "A") || ($c == "C") || ($c == "G") || ($c == "T")) return 1} {if(NR==1) {print} else {if(alchk(3) && alchk(4)) print}}' \
	| awk 'BEGIN{OFS="\t"} function g2n(c) {gsub("A","1",$c);gsub("C","2",$c);gsub("G","3",$c);gsub("T","4",$c)} {g2n(3) g2n(4)} {print $2,$3,$4}' > <ref_path>/1000genomesloci/1000genomesAlleles2012_chrchrX.txt

## Make 1000genomesloci files
### Autosomes
# Per autosome
	for i in $(seq 1 22 | sed 's:^:chr:g'); 
	do 
	cat </impute/ALL.v1a.shapeit2_integrated_"$i".GRCh38.20181129.phased.legend \
		| awk -v s="$i" 'function alchk(c) {if(($c == "A") || ($c == "C") || ($c == "G") || ($c == "T")) return 1} {if(alchk(3) && alchk(4)) print s"\t"$2}' \
	> <ref_path>/1000genomesloci/1000genomesloci2012_chr"$i".txt
	done

### ChrX
# ChrX only
	i="chrX";

	cat </impute/ALL.v1a.shapeit2_integrated_chrX_*.GRCh38.20181129.phased.legend \
	| sort -nk2 \
	| awk -v s="$i" 'function alchk(c) {if(($c == "A") || ($c == "C") || ($c == "G") || ($c == "T")) return 1} {if(alchk(3) && alchk(4)) print s"\t"$2}' \
	> <ref_path>/1000genomesloci/1000genomesloci2012_chrchrX.txt

