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

#library(hash)
source("ProbLociUtils.R")


#This counts the number of bad samples and sets a threshold to give a similar number of badloci to the previous attempt.
#It does not calculate the fraction of heterozygous SNPs that are bad
summariseBadLoci1 = function(fasta_file, ignored_contigs) {
    chrs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
    #all.loci = NULL
    #chrs = c(1:22,"X")
    #for(c in 1:23){
    #	print(chrs[c])
    #	loci = read.table(paste("bad_loci_info_chr",chrs[c],".txt",sep=""),header=T,stringsAsFactors=F)
    #	if(is.null(all.loci)){
    #		all.loci = loci
    #	}else{
    #		all.loci = rbind(all.loci,loci)
    #	}
    #}
    #TAB = table(all.loci$no.bad.samples)
    #print(TAB)
    #
    #previous.bad.loci = read.table("",header=T,stringsAsFactors=F)
    #nearest.threshold = which.min(abs(cumsum(rev(TAB))-nrow(previous.bad.loci)))
    #thresh = rev(names(TAB))[nearest.threshold]
    #print(paste("nearest threshold = ",thresh,sep=""))
    #filtered.loci = all.loci[all.loci$no.bad.samples>=thresh,]
    #
    #print(warnings())
    #save.image(file="SummariseBadLoci.RData")

    load(file="SummariseBadLoci.RData")
    #filtered.loci$in.previous = sapply(1:nrow(filtered.loci),function(i){sum(filtered.loci$chr[i]==previous.bad.loci$Chr & filtered.loci$pos[i]==previous.bad.loci$Pos)>0})
    #previous.bad.loci$in.new.filter = sapply(1:nrow(previous.bad.loci),function(i){sum(previous.bad.loci$Chr[i]==filtered.loci$chr & previous.bad.loci$Pos[i]==filtered.loci$pos)>0})
    #previous.bad.loci$in.new.1000G = sapply(1:nrow(previous.bad.loci),function(i){sum(previous.bad.loci$Chr[i]==all.loci$chr & previous.bad.loci$Pos[i]==all.loci$pos)>0})
    #040315 exclude loci on chrY
    previous.bad.loci = previous.bad.loci[previous.bad.loci$Chr %in% chrs,]
    #quicker method to get chr/pos matches
    previous.bad.loci$chrpos = 1000000000*match(previous.bad.loci$Chr,chrs)+previous.bad.loci$Pos
    filtered.loci$chrpos = 1000000000*match(filtered.loci$chr,chrs)+filtered.loci$pos
    all.loci$chrpos = 1000000000*match(all.loci$chr,chrs)+all.loci$pos
    filtered.loci$in.previous = filtered.loci$chrpos %in% previous.bad.loci$chrpos
    previous.bad.loci$in.new.filter = previous.bad.loci$chrpos %in% filtered.loci$chrpos
    previous.bad.loci$in.new.1000G = previous.bad.loci$chrpos %in% all.loci$chrpos

    print(paste(sum(filtered.loci$in.previous), " of ",nrow(filtered.loci)," newFilteredLoci in previous badloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.filter), " of ",nrow(previous.bad.loci)," previousBadLoci in newFilteredloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.1000G), " of ",nrow(previous.bad.loci)," previousBadLoci in 1000G list",sep=""))

    q(save="no")
}


summariseBadLoci2 = function(prev_prob_loci_file, fasta_file, ignored_contigs, output_directory) {
    all.loci = NULL
    all.het.loci = NULL
    chrs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
    for(c in 1:length(chrs)){
        print(chrs[c])
        loci = read.table(paste(output_directory,"/bad_loci2_info_chr",chrs[c],".txt",sep=""),header=T,stringsAsFactors=F)
        if(is.null(all.loci)){
            all.loci = loci
        }else{
            all.loci = rbind(all.loci,loci)
        }

        #het.loci = read.table(paste("het_loci_info_chr",chrs[c],".txt",sep=""),header=T,stringsAsFactors=F)
        #if(is.null(all.hetloci)){
        #	all.het.loci = het.loci
        #}else{
        #	all.het.loci = rbind(all.het.loci,het.loci)
        #}
    }

    save.image(file="SummariseBadLoci2.RData")
    #load(file="SummariseBadLoci2.RData")

    all.loci$fraction.bad = all.loci$no.bad.samples/(all.loci$no.bad.samples + all.loci$no.het.samples)

    hist.data = hist(all.loci$fraction.bad,breaks = seq(0,1,0.05))
    print(hist.data)
    png("fraction_bad_histogram.png")
    plot(hist.data,col="blue",xlab="fraction bad loci")
    dev.off()

    filtered.loci = all.loci[all.loci$fraction.bad>=0.05 & !is.na(all.loci$fraction.bad),]

    previous.bad.loci = read.table(prev_prob_loci_file,header=T,stringsAsFactors=F)

    print(warnings())
    #040315 exclude loci on chrY
    previous.bad.loci = previous.bad.loci[previous.bad.loci$Chr %in% chrs,]
    #quicker method to get chr/pos matches
    previous.bad.loci$chrpos = 1000000000*match(previous.bad.loci$Chr,chrs)+previous.bad.loci$Pos
    filtered.loci$chrpos = 1000000000*match(filtered.loci$chr,chrs)+filtered.loci$pos
    all.loci$chrpos = 1000000000*match(all.loci$chr,chrs)+all.loci$pos
    filtered.loci$in.previous = filtered.loci$chrpos %in% previous.bad.loci$chrpos
    previous.bad.loci$in.new.filter = previous.bad.loci$chrpos %in% filtered.loci$chrpos
    previous.bad.loci$in.new.1000G = previous.bad.loci$chrpos %in% all.loci$chrpos

    print(paste(sum(filtered.loci$in.previous), " of ",nrow(filtered.loci)," newFilteredLoci in previous badloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.filter), " of ",nrow(previous.bad.loci)," previousBadLoci in newFilteredloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.1000G), " of ",nrow(previous.bad.loci)," previousBadLoci in 1000G list",sep=""))

    write.table(filtered.loci,paste(output_directory,"/badloci_5March2014.txt"),sep="\t",row.names=F,quote=F)

    q(save="no")

}


summariseBadLoci2_minDepth = function(prev_prob_loci_file, fasta_file, ignored_contigs, output_directory){
    chrs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
    all.loci = NULL
    all.het.loci = NULL
    for(c in 1:length(chrs)){
        print(chrs[c])
        loci = read.table(paste(output_directory,"/bad_loci2_info_chr",chrs[c],"_March2015_minDepth.txt",sep=""),header=T,stringsAsFactors=F)
        if(is.null(all.loci)){
            all.loci = loci
        }else{
            all.loci = rbind(all.loci,loci)
        }

        #het.loci = read.table(paste("het_loci_info_chr",chrs[c],".txt",sep=""),header=T,stringsAsFactors=F)
        #if(is.null(all.hetloci)){
        #	all.het.loci = het.loci
        #}else{
        #	all.het.loci = rbind(all.het.loci,het.loci)
        #}
    }

    save.image(file="SummariseBadLoci2_minDepth.RData")
    #load(file="SummariseBadLoci2.RData")

    all.loci$fraction.bad = all.loci$no.bad.samples/(all.loci$no.bad.samples + all.loci$no.het.samples)

    hist.data = hist(all.loci$fraction.bad,breaks = seq(0,1,0.05))
    print(hist.data)
    png("fraction_bad_histogram.png")
    plot(hist.data,col="blue",xlab="fraction bad loci")
    dev.off()

    filtered.loci = all.loci[all.loci$fraction.bad>=2.5*0.053876 & !is.na(all.loci$fraction.bad),]

    previous.bad.loci = read.table(prev_prob_loci_file,header=T,stringsAsFactors=F)

    print(warnings())
    #040315 exclude loci on chrY
    previous.bad.loci = previous.bad.loci[previous.bad.loci$Chr %in% chrs,]
    #quicker method to get chr/pos matches
    previous.bad.loci$chrpos = 1000000000*match(previous.bad.loci$Chr,chrs)+previous.bad.loci$Pos
    filtered.loci$chrpos = 1000000000*match(filtered.loci$chr,chrs)+filtered.loci$pos
    all.loci$chrpos = 1000000000*match(all.loci$chr,chrs)+all.loci$pos
    filtered.loci$in.previous = filtered.loci$chrpos %in% previous.bad.loci$chrpos
    previous.bad.loci$in.new.filter = previous.bad.loci$chrpos %in% filtered.loci$chrpos
    previous.bad.loci$in.new.1000G = previous.bad.loci$chrpos %in% all.loci$chrpos

    print(paste(sum(filtered.loci$in.previous), " of ",nrow(filtered.loci)," newFilteredLoci in previous badloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.filter), " of ",nrow(previous.bad.loci)," previousBadLoci in newFilteredloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.1000G), " of ",nrow(previous.bad.loci)," previousBadLoci in 1000G list",sep=""))

    write.table(filtered.loci,paste(output_directory,"/badloci_6March2014_PVL.txt"),sep="\t",row.names=F,quote=F)

    q(save="no")
}

summariseBadLoci3 = function(prev_prob_loci_file, fasta_file, ignored_contigs, output_directory){
    chrs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
    all.loci = NULL
    all.het.loci = NULL
    for(c in 1:length(chrs)){
        print(chrs[c])
        loci = read.table(paste(output_directory,"/bad_loci_info_chr",chrs[c],"_March2015.txt",sep=""),header=T,stringsAsFactors=F)
        if(is.null(all.loci)){
            all.loci = loci
        }else{
            all.loci = rbind(all.loci,loci)
        }

        #het.loci = read.table(paste("het_loci_info_chr",chrs[c],".txt",sep=""),header=T,stringsAsFactors=F)
        #if(is.null(all.hetloci)){
        #	all.het.loci = het.loci
        #}else{
        #	all.het.loci = rbind(all.het.loci,het.loci)
        #}
    }

    save.image(file="SummariseBadLoci3.RData")
    #load(file="SummariseBadLoci3.RData")

    all.loci$fraction.bad = all.loci$no.bad.samples/(all.loci$no.bad.samples + all.loci$no.het.samples)

    hist.data = hist(all.loci$fraction.bad,breaks = seq(0,1,0.05))
    print(hist.data)
    png("fraction_bad_histogram_March2015.png")
    plot(hist.data,col="blue",xlab="fraction bad loci")
    dev.off()

    filtered.loci = all.loci[all.loci$fraction.bad>=0.05 & !is.na(all.loci$fraction.bad),]

    previous.bad.loci = read.table(prev_prob_loci_file,header=T,stringsAsFactors=F)

    print(warnings())
    #040315 exclude loci on chrY
    previous.bad.loci = previous.bad.loci[previous.bad.loci$Chr %in% chrs,]
    #quicker method to get chr/pos matches
    previous.bad.loci$chrpos = 1000000000*match(previous.bad.loci$Chr,chrs)+previous.bad.loci$Pos
    filtered.loci$chrpos = 1000000000*match(filtered.loci$chr,chrs)+filtered.loci$pos
    all.loci$chrpos = 1000000000*match(all.loci$chr,chrs)+all.loci$pos
    filtered.loci$in.previous = filtered.loci$chrpos %in% previous.bad.loci$chrpos
    previous.bad.loci$in.new.filter = previous.bad.loci$chrpos %in% filtered.loci$chrpos
    previous.bad.loci$in.new.1000G = previous.bad.loci$chrpos %in% all.loci$chrpos

    print(paste(sum(filtered.loci$in.previous), " of ",nrow(filtered.loci)," newFilteredLoci in previous badloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.filter), " of ",nrow(previous.bad.loci)," previousBadLoci in newFilteredloci list",sep=""))
    print(paste(sum(previous.bad.loci$in.new.1000G), " of ",nrow(previous.bad.loci)," previousBadLoci in 1000G list",sep=""))

    write.table(filtered.loci,paste(output_directory,"/badloci_6March2014.txt"),sep="\t",row.names=F,quote=F)

    q(save="no")
}
