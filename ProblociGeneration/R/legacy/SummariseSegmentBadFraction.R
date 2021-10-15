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

source("fastPCF.R")
source("ProbLociUtils.R")
contigs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
outdir = "segmented_plots"
chrs = c(1:22,"X")
gammas = c(25,50,100,250)
no.bad.SNPs = array(NA,c(23,4))
all.fraction.bad = list()
exclude.regions = list()
for(g in 1:4){
	all.fraction.bad[[g]] = vector(mode="numeric",length=0)
	#exclude.regions[[g]] = data.frame(chr=vector(mode="character",length=0),startpos=vector(mode="integer",length=0),endpos=vector(mode="integer",length=0),noSNPs=vector(mode="integer",length=0))
	exclude.regions[[g]] = array(NA,c(0,4))
}

for(c in 1:23){
	chr=chrs[c]
	for(g in 1:4){
		gamma = gammas[g]
		#bad.regions = read.table(paste(outdir,"/chr",chr,"_gamma",gamma,"_segmentedBadFraction.txt",sep=""),header=T,stringsAsFactors=F)
		bad.regions = read.table(paste(outdir,"/gamma",gamma,"/chr",chr,"_gamma",gamma,"_segmentedBadFraction.txt",sep=""),header=T,stringsAsFactors=F)
		bad.regions$exclude=F
		bad.regions$exclude[bad.regions$segmentedBadFraction>=0.04]=T
		no.bad.SNPs[c,g] = sum(bad.regions$numberOfSNPs[bad.regions$exclude])
		all.fraction.bad[[g]] = c(all.fraction.bad[[g]], bad.regions$segmentedBadFraction)
		if(sum(bad.regions$exclude)>0){
			start.exclude = which(bad.regions$exclude==T & c(F,bad.regions$exclude[-nrow(bad.regions)])==F)
			end.exclude = which(bad.regions$exclude==T & c(bad.regions$exclude[-1],F)==F)
			for(i in 1:length(start.exclude)){
				exclude.regions[[g]] = rbind(exclude.regions[[g]],c(chr,bad.regions$startpos[start.exclude[i]],bad.regions$endpos[end.exclude[i]],sum(bad.regions$numberOfSNPs[start.exclude[i]:end.exclude[i]])))
			}
		}
	}
}
print(no.bad.SNPs)
print(colSums(no.bad.SNPs))

png("bad_fraction_histogram.png",width=1000,height=1000)
par(mfrow=c(2,2))
for(g in 1:4){
	hist(all.fraction.bad[[g]],breaks=seq(0,1,0.005),col="blue",xlab="fraction bad loci",main=paste("gamma = ",gammas[g],sep=""))
	write.table(exclude.regions[[g]],paste("CopyNumberExcludeRegions_gamma",gammas[g],".txt",sep=""),row.names=F,quote=F,sep="\t",col.names=c("chr","startpos","endpos","numberOfSNPs"))
}
dev.off()

#q(save="no")
