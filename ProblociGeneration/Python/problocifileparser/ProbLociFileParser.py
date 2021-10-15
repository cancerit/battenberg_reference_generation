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
import glob
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



class ProbLociFileParser:

    "Class to hold prob loci file parser logic"

    def __init__(self, bin_cutoff, depth_cutoff, baf_cutoff):
        self.binomial_cutoff = bin_cutoff
        self.depth_cutoff = depth_cutoff
        self.baf_cutoff = baf_cutoff

    def get_matching_files(self, pattern, contig):
        """
        Given file pattern and contig, searches for all matching files 
        once the contig is inserted where %s exists.
        Returns a name sorted list of files found to match.
        """
        # Find and replace %* with * for file pattern match.

        star_pattern = re.sub(r"%\*", r"*", pattern)
        new_pattern = re.sub(r"%s", contig, star_pattern)
        #find matching files.
        filelist = glob.glob(new_pattern)
        # Sprintf contig into %s
        return filelist

    def get_file_lists(self, contigs, binomial_file_pattern, depthbaf_file_pattern):
        """
        Takes a list of contigs and the file patterns 
        for binomial output and PVL output files. 
        Returns two dicts of contig -> filelists that were found matching
        the patterns. (binomial_file_dict, depth_baf_file_dict)
        """
        binomial_file_dict = dict()
        depth_baf_file_dict = dict()
        for contig in contigs:
            # binomial
            binomial_file_list = self.get_matching_files(binomial_file_pattern, contig)
            # depth and BAF
            depth_baf_file_list = self.get_matching_files(depthbaf_file_pattern, contig)
            binomial_file_dict[contig] = binomial_file_list
            depth_baf_file_dict[contig] = depth_baf_file_list
        return (binomial_file_dict, depth_baf_file_dict)

    def check_binomial_cutoff(self, val):
        """
        Return True if the binomial value in the line is
        equal to or above the cutoff
        """
        if int(val) >= self.binomial_cutoff:
            return True
        return False

    def check_baf_cutoff(self, val):
        """
        Return True if the binomial value in the line is
        equal to or above the cutoff
        """
        if int(val) >= self.baf_cutoff:
            return True
        return False

    def check_depth_cutoff(self, val):
        """
        Return True if the depth and BAF values in the line are
        equal to or above the cutoff
        """
        if int(val) >= self.depth_cutoff:
            return True
        return False

    def generate_bad_sample_bar_plot(self, output_location, bad_sample_counts_list):
        print (f"plotting barchart {output_location}")
        plt.bar(np.arange(len(bad_sample_counts_list)),height=bad_sample_counts_list)

        print (f"saving barchart")
        plt.savefig(output_location)
        plt.close()
        # output_directory,"/","no_bad_samples_barchart_chr",contigs[c],"_March2015.png",sep=""))
        # barplot(table(loci$no.bad.samples),col="blue")
        # dev.off()

    def generate_bad_sample_hist(self, output_location, bad_sample_props_list):
        print (f"plotting hist {output_location}")
        plt.hist(bad_sample_props_list, range=(0,max(bad_sample_props_list)))
        plt.xlabel("fraction of bad samples")
        print (f"saving hist")
        plt.savefig(output_location)
        plt.close()
        # png(paste(output_directory,"/","fraction_bad_samples_barchart_chr",contigs[c],"_March2015.png",sep=""))
        # hist(loci$no.bad.samples/(loci$no.bad.samples+loci$no.het.samples),breaks = seq(0,1,0.01),col="blue",xlab="fraction of bad samples",main = paste("chr",contigs[c]," bad loci frequency",sep=""))
        # dev.off()
        
    def parse_bad_loci_file(self, path_to_bad_loci):
        bad_loci_list = list()
        with open(path_to_bad_loci) as file:
            for line in file:
                line = line.strip()
                if line == "Chr\tPos": continue
                value = line.replace("\t","_")
                bad_loci_list.append(value)
        return bad_loci_list

