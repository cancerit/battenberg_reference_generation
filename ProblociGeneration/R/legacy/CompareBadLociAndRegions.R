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

contigs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
normal.samples = read.table("extended_fresh_frozen_unmatched_normals.txt",header=F,row.names=NULL,stringsAsFactors=F)[,2]

bad.regions = read.table(paste("CopyNumberExcludeRegions_gamma100.txt",sep=""),header=T,sep="\t",row.names=NULL,stringsAsFactors=F)
bad.regions$chrposstart = 1000000000*match(bad.regions$chr,contigs)+bad.regions$startpos
bad.regions$chrposend = 1000000000*match(bad.regions$chr,contigs)+bad.regions$endpos

#100415 for chr6 filter out larger regions, to exclude the whole of the HLA region
bad.regions.250 = read.table(paste("CopyNumberExcludeRegions_gamma250.txt",sep=""),header=T,sep="\t",row.names=NULL,stringsAsFactors=F)
bad.regions.250 = bad.regions.250[bad.regions.250$chr==6,]
bad.regions.250$chrposstart = 1000000000*match(bad.regions.250$chr,contigs)+bad.regions.250$startpos
bad.regions.250$chrposend = 1000000000*match(bad.regions.250$chr,contigs)+bad.regions.250$endpos

previous.bad.loci = read.table("/path/to/probloci.txt",header=T,stringsAsFactors=F)

new.badloci = read.table("badloci_6March2014.txt",sep="\t",header=T,stringsAsFactors=F)
new.badloci.PVL = read.table("badloci_6March2014_PVL.txt",sep="\t",header=T,stringsAsFactors=F)

new.badloci$chrpos = 1000000000*match(new.badloci$chr,contigs)+new.badloci$pos
new.badloci.PVL$chrpos = 1000000000*match(new.badloci.PVL$chr,contigs)+new.badloci.PVL$pos

#040315 exclude loci on chrY
previous.bad.loci = previous.bad.loci[previous.bad.loci$Chr %in% contigs,]
#quicker method to get chr/pos matches
previous.bad.loci$chrpos = 1000000000*match(previous.bad.loci$Chr,contigs)+previous.bad.loci$Pos

previous.bad.loci$in.new.regions = F
new.badloci$in.new.regions = F
new.badloci.PVL$in.new.regions = F
for(r in 1:nrow(bad.regions)){
	previous.bad.loci$in.new.regions[previous.bad.loci$chrpos >= bad.regions$chrposstart[r] & previous.bad.loci$chrpos <= bad.regions$chrposend[r]] = T
	new.badloci$in.new.regions[new.badloci$chrpos >= bad.regions$chrposstart[r] & new.badloci$chrpos <= bad.regions$chrposend[r]] = T
	new.badloci.PVL$in.new.regions[new.badloci.PVL$chrpos >= bad.regions$chrposstart[r] & new.badloci.PVL$chrpos <= bad.regions$chrposend[r]] = T
}

print(paste(sum(previous.bad.loci$in.new.regions), " of ",nrow(previous.bad.loci)," previous badloci in new regions",sep=""))
print(paste(sum(new.badloci$in.new.regions), " of ",nrow(new.badloci)," new badloci in new regions",sep=""))
print(paste(sum(new.badloci.PVL$in.new.regions), " of ",nrow(new.badloci.PVL)," new PVL badloci in new regions",sep=""))

complete.new.bad.loci = new.badloci[,1:2]
for(chr.no in 1:length(contigs)){
	print(chr.no)
	loci.1000G = read.table(paste("path/to/1000genomesloci2012_chr",chr.no,".txt",sep=""))
	names(loci.1000G) = names(complete.new.bad.loci)
	chr.bad.regions = bad.regions[bad.regions$chr==contigs[chr.no],]
	for(r in 1:nrow(chr.bad.regions)){
		complete.new.bad.loci = rbind(complete.new.bad.loci,loci.1000G[loci.1000G[,2]>=chr.bad.regions$startpos[r] & loci.1000G[,2]<=chr.bad.regions$endpos[r],])
	}
	#100415 - additional filtering of larger regions on chr6
	chr.bad.regions = bad.regions.250[bad.regions.250$chr==contigs[chr.no],]
	for(r in 1:nrow(chr.bad.regions)){
		complete.new.bad.loci = rbind(complete.new.bad.loci,loci.1000G[loci.1000G[,2]>=chr.bad.regions$startpos[r] & loci.1000G[,2]<=chr.bad.regions$endpos[r],])
	}
}
complete.new.bad.loci = complete.new.bad.loci[order(match(complete.new.bad.loci$chr,contigs),complete.new.bad.loci$pos),]
dup = which(duplicated(complete.new.bad.loci))
print(paste("# duplicated = ",length(dup),sep=""))
complete.new.bad.loci = complete.new.bad.loci[-dup,]
write.table(complete.new.bad.loci,"probloci_100415.txt",sep="\t",col.names=c("Chr","Pos"),quote=F,row.names=F)
#q(save="no")
