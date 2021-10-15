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
library(ggplot2)
library(broom)
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

#function for doing binomal test, and summarising results
binom.fun <- function(x) {
  mutcount_col = paste0(x, ("_MTR"))
  depth_col = paste0(x, "_Depth")
  cat(paste(Sys.time(), x, "\n", sep = "\t"))
  snp_depth_filt = filter(snp.count.tbl, !!sym(depth_col) >=1) # removes all sites w/ less than 2X depth
  df <- snp_depth_filt %>%
    group_by(chr_pos) %>%
    do(tidy(binom.test(x=as.numeric(.[,mutcount_col]), n=as.numeric(.[,depth_col]), alternative = "two.sided", p=0.01, conf.level = 0.98)))
  df_join <- left_join(snp.count.tbl[,ncol(snp.count.tbl)], df, by = "chr_pos")
  df_sel <- df_join[,c("chr_pos", "conf.low", "conf.high", "estimate", "parameter")] %>%
    mutate(snp_0.01 = ifelse(conf.low <= 0.01, T, F),
           snp_0.5 = ifelse(conf.low <=0.5 & conf.high >= 0.5, T, F),
           snp_0.99 = ifelse(conf.high >=.99, T,F)) %>%
    mutate(snp_tot = snp_0.01 + snp_0.5 + snp_0.99) %>%
    mutate(binom_filt = ifelse(snp_tot == 1, T, F))
  colnames(df_sel) <- c("chr_pos", paste0(x, "_conf.low"), paste0(x, "_conf.high"), paste0(x, "_BAF"), paste0(x, "_Depth"),
                        paste0(x, "_snp_0.01"), paste0(x, "_snp_0.5"), paste0(x, "_snp_0.99"), paste0(x, "_snp_tot"), paste0(x, "_binom_filt"))
  df_sel
}

# get allele count file
snp.tbl <- tbl_df(fread(input.file, colClasses = "guess"))

nSNPs = nrow(snp.tbl)
nSamples = (ncol(snp.tbl)-2)/3;

cat(paste0("\nInitial SNP count: ", nSNPs, "\nSample Count: ", nSamples))

# Get list of samples
cols.mutcount <- seq(4, ncol(snp.tbl), by =3)
sample.list <- gsub("_MTR", "",  colnames(snp.tbl[,cols.mutcount]))

# combine chr_pos to use as identifyer
snp.count.tbl <- snp.tbl  %>%
  dplyr::mutate(chr_pos = paste(CHR, POS, sep = "_"))

cat("\nstarting binomial tests\n")
# run binomal tests
binom.out <- lapply(sample.list, binom.fun)

# combine results into single df
binom.df <- Reduce(function(...) merge(..., by = "chr_pos", all.x=TRUE), binom.out)

# count passing samples per snp
filt_df_calc <- binom.df %>%
  mutate(binom_pass = rowSums(.[,grepl("*binom_filt", names(.))], na.rm = T))

filt_sum <- filt_df_calc %>%
  dplyr::select(chr_pos, binom_pass)

write.csv(filt_df_calc, file = paste0(out.dir, input.file.name, "_binomal_results_long.csv"), quote = F, row.names = F)
write.csv(filt_sum, file = paste0(out.dir, input.file.name, "_binomal_results_summary.csv"), quote = F, row.names = F)
