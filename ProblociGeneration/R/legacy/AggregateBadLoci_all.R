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
library(hash)
source("ProbLociUtils.R")
contigs = get_contigs_from_fasta_exclude_ignored(fasta_file, ignored_contigs)
normal.samples = read.table("extended_fresh_frozen_unmatched_normals.txt",header=F,row.names=NULL,stringsAsFactors=F)[,2]


aggregate_bad_loci = function(fasta_file, ignored_contigs, 
                                normal_sample_list_file, output_directory, thousand_genomes_loci_dir, normal.samples, contigs) {
    no.samples = length(normal.samples)
    bad.indices=list()
    for(c in 1:length(contigs)){
        bad.indices[[c]] = vector(length=0,mode="integer")
    }
    het.indices=list()
    for(c in 1:length(contigs)){
        het.indices[[c]] = vector(length=0,mode="integer")
    }
    for(s in 1:no.samples){
        samplename = normal.samples[s]
        print(samplename)
        for(c in 1:length(contigs)){
            bi = read.table(paste(output_directory,"/",samplename,"_badloci_chr",contigs[c],".txt",sep=""),header=T)[,1]
            bad.indices[[c]] = c(bad.indices[[c]],bi)
            hi = read.table(paste(output_directory,"/",samplename,"_hetloci_chr",contigs[c],".txt",sep=""),header=T)[,1]
            het.indices[[c]] = c(het.indices[[c]],hi)
        }
    }

    for(c in 1:length(contigs)){
        print(contigs[c])
        loci = read.table(paste(thousand_genomes_loci_dir,"/","1000genomesloci2012_chr",contigs[c],".txt",sep=""),sep=" ",header=F,row.names=NULL)
        loci$no.bad.samples=0
        loci$no.het.samples=0
        
        #save.image(file="debugAggregateBadLoci.RData")
        
        TAB = table(bad.indices[[c]])
        loci$no.bad.samples[as.numeric(names(TAB))] = TAB
        
        TAB = table(het.indices[[c]])
        loci$no.het.samples[as.numeric(names(TAB))] = TAB

        write.table(loci,paste(output_directory,"/","bad_loci_info_chr",contigs[c],"_March2015.txt",sep=""),quote=F,row.names=F,col.names=c("chr","pos","no.bad.samples","no.het.samples"))

        png(paste(output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_March2015.png",sep=""))
        barplot(table(loci$no.bad.samples),col="blue")
        dev.off()
        
        png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",chrs[c],"_March2015.png",sep=""))
        hist(loci$no.bad.samples/(loci$no.bad.samples+loci$no.het.samples),breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples",main = paste("chr",contigs[c]," bad loci frequency",sep=""))
        dev.off()
    }

    q(save="no")
}

aggregate_bad_loci2_minDepth = function(fasta_file, ignored_contigs, 
                                normal_sample_list_file, output_directory, thousand_genomes_loci_dir, normal.samples, contigs){
    no.samples = length(normal.samples)
    bad.indices=list()
    for(c in 1:length(contigs)){
        bad.indices[[c]] = vector(length=0,mode="integer")
    }
    het.indices=list()
    for(c in 1:length(contigs)){
        het.indices[[c]] = vector(length=0,mode="integer")
    }
    min.depth=20
    for(s in 1:no.samples){
        samplename = normal.samples[s]
        print(samplename)
        for(c in 1:length(contigs)){
            bi = read.table(paste(output_directory,"/",samplename,"_badloci2_chr",contigs[c],".txt",sep=""),header=T)[,1]
            bd = read.table(paste(output_directory,"/",samplename,"_bad_depths2_chr",contigs[c],".txt",sep=""),header=T)[,1]
            bad.indices[[c]] = c(bad.indices[[c]],bi[bd>=min.depth])
            hi = read.table(paste(output_directory,"/",samplename,"_hetloci_chr",contigs[c],".txt",sep=""),header=T)[,1]
            hd = read.table(paste(output_directory,"/",samplename,"_het_depths_chr",contigs[c],".txt",sep=""),header=T)[,1]
            het.indices[[c]] = c(het.indices[[c]],hi[hd>=min.depth])
        }
    }

    for(c in 1:length(contigs)){
        print(contigs[c])
        loci = read.table(paste(thousand_genomes_loci_dir,"/","1000genomesloci2012_chr",contigs[c],".txt",sep=""),sep=" ",header=F,row.names=NULL)
        loci$no.bad.samples=0
        loci$no.het.samples=0
        
        #save.image(file="debugAggregateBadLoci.RData")
        
        TAB = table(bad.indices[[c]])
        loci$no.bad.samples[as.numeric(names(TAB))] = TAB
        
        TAB = table(het.indices[[c]])
        loci$no.het.samples[as.numeric(names(TAB))] = TAB

        write.table(loci,paste(output_directory,"/","bad_loci2_info_chr",contigs[c],"_March2015_minDepth.txt",sep=""),quote=F,row.names=F,col.names=c("chr","pos","no.bad.samples","no.het.samples"))

        png(paste(output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_March2015_2_minDepth.png",sep=""))
        barplot(table(loci$no.bad.samples),col="blue")
        dev.off()
        
        png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",contigs[c],"_March2015_2_minDepth.png",sep=""))
        hist(loci$no.bad.samples/(loci$no.bad.samples+loci$no.het.samples),breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples",main = paste("chr",contigs[c]," bad loci frequency",sep=""))
        dev.off()
    }

    q(save="no")
}

aggregate_bad_loci2_usingHash = function(fasta_file, ignored_contigs, 
                                normal_sample_list_file, output_directory, thousand_genomes_loci_dir, normal.samples, contigs){
    no.samples = length(normal.samples)
    bad.indices=list()
    for(c in 1:length(contigs)){
        bad.indices[c] = hash()
    }
    het.indices=list()
    for(c in 1:length(contigs)){
        het.indices[c] = hash()
    }
    for(s in 1:no.samples){
        samplename = normal.samples[s]
        print(samplename)
        for(c in 1:length(contigs)){
            bi = read.table(paste(output_directory,"/",samplename,"_badloci2_chr",c,".txt",sep=""),header=T)[,1]
            for(key in bi){
                if(has.key(toString(key),bad.indices[[c]])){
                    bad.indices[[c]][[toString(key)]] <- bad.indices[[c]][[toString(key)]] + 1
                }else{
                    bad.indices[[c]][[toString(key)]] <- 1
                }
            }
            hi = read.table(paste(output_directory,"/",samplename,"_hetloci_chr",c,".txt",sep=""),header=T)[,1]
            for(key in hi){
                if(has.key(toString(key),het.indices[[c]])){
                    het.indices[[c]][[toString(key)]] <- het.indices[[c]][[toString(key)]] + 1
                }else{
                    het.indices[[c]][[toString(key)]] <- 1
                }
            }
        }
    }

    for(c in 1:length(contigs)){
        print(contigs[c])
        loci = read.table(paste(thousand_genomes_loci_dir,"/","1000genomesloci2012_chr",contigs[c],".txt",sep=""),sep=" ",header=F,row.names=NULL)
        loci$no.bad.samples=0
        
        vals = values(bad.indices[[c]])
        loci$no.bad.samples[as.numeric(names(vals))] = vals
        
        vals2 = values(het.indices[[c]])
        loci$no.het.samples[as.numeric(names(vals2))] = vals2

        write.table(loci,paste(output_directory,"/","bad_loci2_info_chr",contigs[c],".txt",sep=""),quote=F,row.names=F,col.names=c("chr","pos","no.bad.samples","no.het.samples"))

        png(paste(output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_2.png",sep=""))
        barplot(table(vals),col="blue")
        dev.off()
        
        png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",contigs[c],"_2.png",sep=""))
        hist(vals/vals2,breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples")
        dev.off()
    }

    q(save="no")
}

aggregate_bad_loci2 = function(fasta_file, ignored_contigs, 
                                normal_sample_list_file, output_directory, thousand_genomes_loci_dir, normal.samples, contigs){
    no.samples = length(normal.samples)
    bad.indices=list()
    for(c in 1:length(contigs)){
        bad.indices[[c]] = vector(length=0,mode="integer")
    }
    het.indices=list()
    for(c in 1:length(contigs)){
        het.indices[[c]] = vector(length=0,mode="integer")
    }
    for(s in 1:no.samples){
        samplename = normal.samples[s]
        print(samplename)
        for(c in 1:length(contigs)){
            #we should also check that the depths are >=20 to match PVL method
            bi = read.table(paste(output_directory,"/",samplename,"_badloci2_chr",contigs[c],".txt",sep=""),header=T)[,1]
            bad.indices[[c]] = c(bad.indices[[c]],bi)
            hi = read.table(paste(output_directory,"/",samplename,"_hetloci_chr",contigs[c],".txt",sep=""),header=T)[,1]
            het.indices[[c]] = c(het.indices[[c]],hi)
        }
    }

    for(c in 1:length(contigs)){
        print(contigs[c])
        loci = read.table(paste(thousand_genomes_loci_dir,"/","1000genomesloci2012_chr",contigs[c],".txt",sep=""),sep=" ",header=F,row.names=NULL)
        loci$no.bad.samples=0
        loci$no.het.samples=0
        
        save.image(file="debugAggregateBadLoci.RData")
        
        TAB = table(bad.indices[[c]])
        loci$no.bad.samples[as.numeric(names(TAB))] = TAB
        
        TAB = table(het.indices[[c]])
        loci$no.het.samples[as.numeric(names(TAB))] = TAB

        write.table(loci,paste(output_directory,"/","bad_loci2_info_chr",contigs[c],".txt",sep=""),quote=F,row.names=F,col.names=c("chr","pos","no.bad.samples","no.het.samples"))

        png(paste(output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_2.png",sep=""))
        barplot(table(loci$no.bad.samples),col="blue")
        dev.off()
        
        png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",contigs[c],"_2.png",sep=""))
        hist(loci$no.bad.samples/(loci$no.bad.samples+loci$no.het.samples),breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples",main = paste("chr",contigs[c]," bad loci frequency",sep=""))
        dev.off()
    }

    q(save="no")

}

aggregate_depths = function(fasta_file, ignored_contigs, 
                        normal_sample_list_file, output_directory, thousand_genomes_loci_dir, normal.samples, contigs){
    no.samples = length(normal.samples)
    bad.indices=list()
    for(c in 1:length(contigs)){
        bad.indices[[c]] = vector(length=0,mode="integer")
    }
    het.indices=list()
    for(c in 1:23){
        het.indices[[c]] = vector(length=0,mode="integer")
    }
    min.depth = 20
    for(s in 1:no.samples){
        samplename = normal.samples[s]
        print(samplename)
        for(c in 1:23){
            bi = read.table(paste(output_directory,"/",samplename,"_badloci2_chr",c,".txt",sep=""),header=T)[,1]
            bd = read.table(paste(output_directory,"/",samplename,"_bad_depths_chr",c,".txt",sep=""),header=T)[,1]				
            bad.indices[[c]] = c(bad.indices[[c]],bi[bd>=min.depth])

            hi = read.table(paste(output_directory,"/",samplename,"_hetloci2_chr",c,".txt",sep=""),header=T)[,1]
            hd = read.table(paste(output_directory,"/",samplename,"_het_depths_chr",c,".txt",sep=""),header=T)[,1]
            het.indices[[c]] = c(het.indices[[c]],hi[hd>=min.depth])
        }
    }

    chrs = c(1:22,"X")
    for(c in 1:length(contigs)){
        print(contigs[c])
        loci = read.table(paste(thousand_genomes_loci_dir,"/1000genomesloci2012_chr",contigs[c],".txt",sep=""),sep=" ",header=F,row.names=NULL)
        loci$no.bad.samples=0
        loci$no.het.samples=0
        
        #save.image(file="debugAggregateBadLoci.RData")
        
        TAB = table(bad.indices[[c]])
        loci$no.bad.samples[as.numeric(names(TAB))] = TAB
        
        TAB = table(het.indices[[c]])
        loci$no.het.samples[as.numeric(names(TAB))] = TAB

        write.table(loci,paste(output_directory,"/","bad_loci2_info_chr",contigs[c],"_March2015_minDepth.txt",sep=""),quote=F,row.names=F,col.names=c("chr","pos","no.bad.samples","no.het.samples"))

        png(paste(output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_March2015_minDepth_2.png",sep=""))
        barplot(table(loci$no.bad.samples),col="blue")
        dev.off()
        
        png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",contigs[c],"_March2015_minDepth_2.png",sep=""))
        hist(loci$no.bad.samples/(loci$no.bad.samples+loci$no.het.samples),breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples",main = paste("chr",contigs[c]," bad loci frequency",sep=""))
        dev.off()
    }

}
