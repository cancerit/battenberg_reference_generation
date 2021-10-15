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
# 1b. IMPUTE - genetic map

# Sources
## Genetic map
	# From Broad - provenance unclear
	# https://data.broadinstitute.org/alkesgroup/Eagle/#x1-30002
    # https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz


# Method

## 1. IMPUTE - genetic map file generation
## copy file
    cd <ref_path>/impute
    wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

## Split by chroms
    zcat ./genetic_map_hg38_withX.txt.gz | awk 'FNR==1{hdr=$0;next} {if (!seen[$1]++) print hdr>"genetic_map_chr"$1"_combined_b38.txt"; print>"genetic_map_chr"$1"_combined_b38.txt"}'

## Trim first columns (chr) from files
    find . -mindepth 1 -maxdepth 1 | xargs -I {} sh -c 'cp -f {} {}.old; cat {}.old | cut -d " " -f2-4 > {}; rm {}.old; '

## Split Chr X 
### PAR1 
    cat genetic_map_chr23_combined_b38.txt | awk 'NR==1; NR>1 {if ($1 <= 2781479) print $0}' > genetic_map_chrX_PAR1_combined_b38.txt
### PAR1 
    cat genetic_map_chr23_combined_b38.txt | awk 'NR==1; NR>1 {if ($1 >= 155701383) print $0}' > genetic_map_chrX_PAR2_combined_b38.txt
### non PAR
    cat genetic_map_chr23_combined_b38.txt | awk 'NR==1; NR>1 {if ($1 > 2781479 && $2 < 155701383) print $0}' > genetic_map_chrX_nonPAR_combined_b38.txt

## Clean-up
### Remove genetic_map_chr23_combined_b38.txt
    rm genetic_map_chr23_combined_b38.txt
