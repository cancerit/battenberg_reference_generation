# LICENSE
# Copyright (c) 2018-2021 Genome Research Ltd.
# Author: Cancer Genome Project cgphelp@sanger.ac.uk
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

args = commandArgs(T)
chr.no = as.numeric(args[1])
source("fastPCF.R")
gammas = c(25,50,100,250)
outdir = "segmented_plots"
chrs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
chr=chrs[chr.no]
print(chr)
bad.loci.info = read.table(paste("bad_loci_info_chr",chr,"_March2015.txt",sep=""),header=T,sep=" ",stringsAsFactors=F)
bad.loci.info$fraction.bad = bad.loci.info$no.bad.samples / (bad.loci.info$no.bad.samples + bad.loci.info$no.het.samples)
bad.loci.info$fraction.bad[is.nan(bad.loci.info$fraction.bad)]=0
sdev <- getMad(bad.loci.info$fraction.bad,k=25)
SNP.pos = bad.loci.info$pos
for(gamma in gammas){
	print(gamma)
	res= selectFastPcf(bad.loci.info$fraction.bad,3,gamma*sdev,T)
	segBad = res$yhat
	png(filename = paste(outdir,"/chr",chr,"_segmented_gamma",gamma,".png",sep=""), width = 2000, height = 1000, res = 200)
	par(mar = c(5,5,5,0.5), cex = 0.6, cex.main=3, cex.axis = 2, cex.lab = 2)
	plot(c(min(SNP.pos)/1000000,max(SNP.pos)/1000000),c(0,1),pch=".",type = "n", 
	main = paste("Chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "fraction bad loci")
	points(SNP.pos/1000000,bad.loci.info$fraction.bad,pch=20,col="red",cex=0.25)
	points(SNP.pos/1000000,segBad,pch=20,cex=1,col="green")
	dev.off()

	cum.pos = c(0,cumsum(res$Lengde))
	#we should also have the number of SNPS
	out = cbind(bad.loci.info$pos[cum.pos[-length(cum.pos)]+1],bad.loci.info$pos[cum.pos[-1]],res$Lengde,res$yhat[cum.pos[-1]])
	write.table(out,paste(outdir,"/chr",chr,"_gamma",gamma,"_segmentedBadFraction.txt",sep=""),col.names=c("startpos","endpos","numberOfSNPs","segmentedBadFraction"),quote=F,row.names=F)
}
q(save="no")
