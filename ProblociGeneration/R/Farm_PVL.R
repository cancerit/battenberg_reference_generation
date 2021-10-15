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

# Load Libraries
library(data.table)
library(rlang)
library(dplyr)

# Input Arguments
# 1. input allele count file
# 2. output directory

# Setup input arguments
args = commandArgs(TRUE)
input.file = args[1]
out.dir = args[2]

input.file.name <- gsub(".csv", "", basename(input.file))

# function to get BAFs, pass/fail on depth and BAF filter
baf.fun <- function(x) {
  mutcount.col <- paste0(x, "_MTR")
  depth.col <- paste0(x, "_Depth")
  cat(paste(Sys.time(), x, "\n", sep = "\t"))
  df <- snp.count.tbl %>%
    dplyr::select(chr_pos, depth.col,  mutcount.col) %>%
    dplyr::mutate(BAF = !!sym(mutcount.col) / !!sym(depth.col))
  df_filt <- df %>%
    dplyr::mutate(depth_filt = ifelse(!!sym(depth.col) >= 20, T, F),
                  baf_filt = case_when(BAF <= 0.1 ~ T,
                                       BAF >= 0.9 ~ T,
                                       (BAF >= 0.3 & BAF <= 0.7) ~ T,
                                       is.na(BAF) ~ NA,
                                       TRUE ~ F))
  colnames(df_filt) <- c("chr_pos", depth.col, mutcount.col, paste0(x, "_BAF"), paste0(x, "_Depth_filt"), paste0(x, "_BAF_filt"))
  df_filt
}

# get allele count file
snp.tbl <- tbl_df(fread(input.file, colClasses = "numeric"))

nSNPs = nrow(snp.tbl)
nSamples = (ncol(snp.tbl)-2)/3;

cat(paste0("\nInitial SNP count: ", nSNPs, "\nSample Count: ", nSamples))

# Get list of samples
cols.mutcount <- seq(4, ncol(snp.tbl), by =3)
sample.list <- gsub("_MTR", "",  colnames(snp.tbl[,cols.mutcount]))

# combine chr_pos to use as identifyer
snp.count.tbl <- snp.tbl  %>%
  dplyr::mutate(chr_pos = paste(CHR, POS, sep = "_"))

cat("\nstarting PVL filters\n")
# run binomal tests
baf.out <- lapply(sample.list, baf.fun)

# combine results into single df
baf.df <- Reduce(function(...) merge(..., by = "chr_pos", all.x=TRUE), baf.out)

# count passing samples per snp
filt_df_calc <- baf.df %>%
  mutate(Depth_pass = rowSums(.[,grepl("*Depth_filt", names(.))], na.rm = T),
         BAF_pass = rowSums(.[,grepl("*BAF_filt", names(.))], na.rm = T))
#make short summary tbl
filt_sum <- filt_df_calc %>%
  dplyr::select(chr_pos, Depth_pass, BAF_pass)

#count #SNPs failng filters
ndepth = sum(filt_sum$Depth_pass < 94, na.rm = T)
nbaf = sum(filt_sum$BAF_pass < 94, na.rm = T)
ndpth_baf = sum((filt_sum$BAF_pass < 94 | filt_sum$Depth_pass < 94), na.rm = T)

cat("\n", nbaf, "SNPs failed BAF Filter\n")
cat(ndepth, "SNPs failed Depth Filter\n")
cat(ndpth_baf, "SNPs failed Either Filter\n")

write.csv(filt_df_calc, file = paste0(out.dir, input.file.name, "_PVL_results_long.csv"), quote = F, row.names = F)
write.csv(filt_sum, file = paste0(out.dir, input.file.name, "_PVL_results_summary.csv"), quote = F, row.names = F)


