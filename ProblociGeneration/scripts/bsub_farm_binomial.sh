# LICENSE
# Copyright (c) 2018-2021 Genome Research Ltd.
# Author: CASM-IT <cgphelp@sanger.ac.uk>
#
#
# This file is part of battenberg.
#
# problocifileparser is free software: you can redistribute it and/or modify it under
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
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
#

#!/bin/bash

INDIR=$1
OUTDIR=$2
PATH_TO_R=$3

# Change to working dir (where error logs will go)
cd $OUTDIR/
mkdir logfiles

MEM=8000
CORE=1
q=long

# Loop through all samples in BAM_DIR
for allele_count in $INDIR/MutCount*;
do
  echo $allele_count
  analysis=${allele_count##*/}
  echo $analysis

# Run R code
  bsub -J $analysis.binom -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -n $CORE -q $q \
          -e $OUTDIR/logfiles/$analysis.stderr -o $OUTDIR/logfiles/$analysis.stdout \
          "Rscript ${PATH_TO_R}/Farm_binomal.R \
          $allele_count \
          $OUTDIR > $OUTDIR/logfiles/$analysis.logfile"
done

