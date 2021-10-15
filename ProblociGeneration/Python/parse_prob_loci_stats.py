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

"""
Main script to parse output from the statistical tests of the allele counts. 
Allows adjusting of the depth and percentage score cutoffs prior to generating
the files required as input to David W's code to generate the prob loci file
"""

# Core imports
import pkg_resources  # part of setuptools
import argparse
import re
import os
import pysam

from problocifileparser.ProbLociFileParser import ProbLociFileParser

VERSION = pkg_resources.require("problocifileparser")[0].version
DEFAULT_BINOMIAL_PASS = 93
DEFAULT_DEPTH_PASS = 93
DEFAULT_BAF_PASS = 93
DEFAULT_NUM_SAMPLES = 98


def setup_args():
    """
    Function to setup arguments using argparse
    """
    parser = argparse.ArgumentParser(
        description="""Parse prob loci stats files and 
                        output files required for running 
                        prob loci generation code.""",
        prog="parse_prob_loci_stats.py",
    )

    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        action="version",
        help="Print version information",
        version="%(prog)s " + VERSION,
    )

    parser.add_argument(
        "-b",
        "--binomial-file-pattern",
        dest="binomial",
        metavar="binomial_.csv",
        help="""Path to and name pattern of binomial csv input files.
                Use %%s to mark the contig name and %%* to mark any number""",
        required=True,
    )

    parser.add_argument(
        "-d",
        "--depth-file-pattern",
        dest="depthbaf",
        metavar="depth.csv",
        help="""Path to and name pattern of depth and baf csv input files.
                Use %%s to mark the contig name and %%* to mark any number""",
        required=True,
    )

    parser.add_argument(
        "-r",
        "--reference",
        dest="ref",
        metavar="path/to/reference.fa",
        help="Path to reference fasta file with associated fasta index.",
        required=True,
    )

    parser.add_argument(
        "-i",
        "--ignore-contigs",
        dest="ignore_contigs",
        metavar="path/to/ignore_contigs.txt",
        help="Path to file containing list of contigs to ignore. One contig per line",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--previous-bad-loci",
        dest="prevbadloci",
        metavar="/path/to/previous/bad_loci.txt",
        help="Path to previous bad_loci file",
        required=True
    )

    parser.add_argument(
        "-o",
        "--outdir",
        metavar="path/to/write/output",
        help="Path to write output files",
        required=False,
        default="./output",
    )

    parser.add_argument(
        "-c",
        "--binomial_cutoff",
        dest="bincutoff",
        metavar="93",
        help=f"""Cutoff for binomial test. 
                Only results equal to or above 
                this value will be inculded.
                [default: {DEFAULT_BINOMIAL_PASS}]""",
        required=False,
        default=DEFAULT_BINOMIAL_PASS,
        type=int,
    )

    parser.add_argument(
        "-e",
        "--depth_cutoff",
        dest="depthcutoff",
        metavar="10",
        help=f"""Cutoff for depth test. 
                Only results equal to or above 
                this value will be inculded.
                [default: {DEFAULT_DEPTH_PASS}]""",
        required=False,
        default=DEFAULT_DEPTH_PASS,
        type=int,
    )

    parser.add_argument(
        "-a",
        "--depth_baf_cutoff",
        dest="depthbafcutoff",
        metavar="93",
        help=f"""Cutoff for BAF test in depth file. 
                Only results equal to or above 
                this value will be inculded.
                [default: {DEFAULT_BAF_PASS}]""",
        required=False,
        default=DEFAULT_BAF_PASS,
        type=int,
    )

    parser.add_argument(
        "-n",
        "--num-samples",
        dest="numsamp",
        metavar="98",
        help=f"""Number of samples in the normal panel.
                [default: {DEFAULT_NUM_SAMPLES}]""",
        required=False,
        default=DEFAULT_NUM_SAMPLES,
        type=int,
    )

    args = parser.parse_args()

    return args


def read_ignore_file(ignore_file):
    """
    Reads ignore file, splitting on newline
    where each line is a contig in the fasta file to ignore.
    Returns a list where each entry is a contig to ignore.
    """
    with open(ignore_file) as file:
        return list(line.strip() for line in file)
    return list()


def get_fasta_contigs(fasta_file):
    """
    Obtain the contigs within the fasta index file associated
    with the given fasta file. 
    Returns a list where each entry is a contig in the fasta file.
    """
    try:
        fa_file = pysam.FastaFile(fasta_file)
        contigs = list(fa_file.references)
        fa_file.close()
    except (ValueError, IOError) as ex:
        raise Exception from ex
    finally:
        fa_file.close()
    return contigs


def get_filtered_contigs(fasta_file, ignore_conts_file):
    """
    Give an ignore list and a fasta file, return a set of contigs
    that are not in the ignore list.
    Returns a set where each entry is a contig in the fasta file.
    """
    ignore_list = read_ignore_file(ignore_conts_file)
    fasta_list = get_fasta_contigs(fasta_file)
    return list(set(fasta_list).difference(set(ignore_list)))


def main():
    args = setup_args()

    #Ensure output folder exists
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    elif not os.path.isdir(args.outdir):
        raise Exception(f"output directory {args.outdir} exists but is not a directory.")
    
    # Read fasta index and ignore contigs file to get list of wanted contigs
    contigs = get_filtered_contigs(args.ref, args.ignore_contigs)

    # Setup file parser
    fileparser = ProbLociFileParser(
        args.bincutoff, args.depthcutoff, args.depthbafcutoff
    )

    # Find all files for list of contigs
    (binomial_file_dict, depth_baf_file_dict) = fileparser.get_file_lists(
        contigs, args.binomial, args.depthbaf
    )

    print(f"Found {len(binomial_file_dict.values())} binomial files and {len(depth_baf_file_dict.values())} depth/BAF files.")

    binomial_all_fail_proprtion = {}
    baf_all_fail_proportion = {}
    depth_all_fail_proportion = {}

    # Iterate through each contig and aggregate
    for contig in contigs:
        print(f"Starting contig {contig}")
        binomial_results_dict = {}
        baf_results_dict = {}
        depth_results_dict = {}
        bad_sample_counts = list()
        bad_sample_proportions = list()
        bad_sample_counts_baf = list()
        bad_sample_proportions_baf = list()

        # Open binomial files and filter according to set value
        # with open(bin_out_per_contig_file, "w") as binout:
        # binout.print("chr pos no.bad.samples no.het.samples")
        for binfile in binomial_file_dict[contig]:
            with open(binfile) as bin:
                for line in bin:
                    if line.startswith("chr_pos,binom_pass"): continue
                    #if fileparser.check_binomial_cutoff(line, ","):
                    line = line.strip()
                    splitline = line.split(",")
                    binomial_results_dict[splitline[0]] = int(splitline[1])

        print(f"Counting binomials for {contig} done. Got {len(binomial_results_dict)}")

        for key,val in binomial_results_dict.items():
            if fileparser.check_binomial_cutoff(val): continue
            num_fail = args.numsamp - val
            bad_sample_counts.append(num_fail)
            prop_fail = 0.0
            if num_fail > 0:
                prop_fail = float(num_fail) / float(args.numsamp)
            bad_sample_proportions.append(prop_fail)
            binomial_all_fail_proprtion[key] = prop_fail

        print(f"Now plotting binomial graphs for {contig} done. Got {len(bad_sample_proportions)} failing locations")

        barplot_filename = f"{args.outdir}/no_bad_samples_binomial_barchart_unfiltered_chr_{contig}.png"
        fileparser.generate_bad_sample_bar_plot(barplot_filename, list(map(lambda x: args.numsamp - x, binomial_results_dict.values())))

        barplot_filename = f"{args.outdir}/no_bad_samples_binomial_barchart_filtered_chr_{contig}.png"
        fileparser.generate_bad_sample_bar_plot(barplot_filename, bad_sample_counts)

        hist_filename = f"{args.outdir}/fraction_bad_samples_binomial_histogram_chr_{contig}.png"
        fileparser.generate_bad_sample_hist(hist_filename, bad_sample_proportions)


        bad_sample_counts.clear()
        bad_sample_proportions.clear()

        for pvlfile in depth_baf_file_dict[contig]:
            with open(pvlfile) as pvl:
                for line in pvl:
                    if line.startswith("chr_pos,Depth_pass,BAF_pass"): continue
                    line = line.strip()
                    splitline = line.split(",")
                    depth_results_dict[splitline[0]] = int(splitline[1])
                    baf_results_dict[splitline[0]] = int(splitline[2])

        # Do BAF
        for key, bafval in baf_results_dict.items():
            if fileparser.check_baf_cutoff(bafval): continue
            baf_num_fail = args.numsamp - bafval
            bad_sample_counts_baf.append(baf_num_fail)
            baf_prop_fail = 0.0
            if baf_num_fail > 0:
                baf_prop_fail = float(baf_num_fail) / float(args.numsamp)
            bad_sample_proportions_baf.append(baf_prop_fail)
            baf_all_fail_proportion[key] = baf_prop_fail

        # do depth 
        for key, depval in depth_results_dict.items():
            if fileparser.check_depth_cutoff(depval): continue
            dep_num_fail = args.numsamp - depval
            bad_sample_counts.append(dep_num_fail)
            dep_prop_fail = 0.0
            if dep_num_fail > 0:
                dep_prop_fail = float(dep_num_fail) / float(args.numsamp)
            bad_sample_proportions.append(dep_prop_fail)
            depth_all_fail_proportion[key] = dep_prop_fail

        barplot_filename = f"{args.outdir}/no_bad_samples_PVL_depth_barchart_chr_{contig}.png"
        fileparser.generate_bad_sample_bar_plot(barplot_filename, bad_sample_counts)

        barplot_filename = f"{args.outdir}/no_bad_samples_PVL_BAF_barchart_chr_{contig}.png"
        fileparser.generate_bad_sample_bar_plot(barplot_filename, bad_sample_counts_baf)

        hist_filename = f"{args.outdir}/fraction_bad_samples_PVL_depth_histogram_chr_{contig}.png"
        fileparser.generate_bad_sample_hist(hist_filename, bad_sample_proportions)

        hist_filename = f"{args.outdir}/fraction_bad_samples_PVL_BAF_histogram_chr_{contig}.png"
        fileparser.generate_bad_sample_hist(hist_filename, bad_sample_proportions_baf)

        bad_sample_counts.clear()
        bad_sample_proportions.clear()
        bad_sample_counts_baf.clear()
        bad_sample_proportions_baf.clear()
        depth_results_dict.clear()
        baf_results_dict.clear()
        binomial_results_dict.clear()
        print(f"Finished contig {contig}")
    # End of iteration through each contig
    hist_filename = f"{args.outdir}/fraction_bad_samples_binomial_histogram_all_chrs.png"
    fileparser.generate_bad_sample_hist(hist_filename, binomial_all_fail_proprtion.values())
    hist_filename = f"{args.outdir}/fraction_bad_samples_PVL_BAF_histogram_all_chrs.png"
    fileparser.generate_bad_sample_hist(hist_filename, baf_all_fail_proportion.values())
    hist_filename = f"{args.outdir}/fraction_bad_samples_PVL_depth_histogram_all_chrs.png"
    fileparser.generate_bad_sample_hist(hist_filename, depth_all_fail_proportion.values())

    # Read in previous bad_loci_file
    prev_bad_loci_list = fileparser.parse_bad_loci_file(args.prevbadloci)
    print(f"Previous bad loci has {len(prev_bad_loci_list)}")

    # Filtered loci onlyprevious prob loci
    print(f"{len(prev_bad_loci_list)} locations in previous prob loci file.")
    # Filtered loci only binomial
    print(f"{len(binomial_all_fail_proprtion.values())} locations in binomial fail list.")
    # Filtered loci only BAF
    print(f"{len(baf_all_fail_proportion.values())} locations in BAF fail list.")
    # Filtered loci only depth
    print(f"{len(depth_all_fail_proportion.values())} locations in depth fail list.\n")

    # Filtered loci only versus previous bad loci
    binomial_and_prev_loci_intersect = list(set(binomial_all_fail_proprtion.keys()) &
                                            set(prev_bad_loci_list))
    print(f"{len(binomial_and_prev_loci_intersect)} matching locations between previous prob loci and binomial fail list.")

    baf_and_prev_loci_intersect = list(set(baf_all_fail_proportion.keys()) &
                                        set(prev_bad_loci_list))
    print(f"{len(baf_and_prev_loci_intersect)} matching locations between previous prob loci and BAF fail list.")

    depth_and_prev_loci_intersect = list(set(depth_all_fail_proportion.keys()) &
                                        set(prev_bad_loci_list))
    print(f"{len(depth_and_prev_loci_intersect)} matching locations between previous prob loci and depth fail list.")

    # Filtered loci only binomial and BAF
    binomial_and_baf_loci_intersect = list(set(binomial_all_fail_proprtion.keys()) &
                                        set(baf_all_fail_proportion.keys()))
    print(f"{len(binomial_and_baf_loci_intersect)} matching locations between binomial and BAF fail list.")

    # Filtered loci only binomial and depth
    binomial_and_depth_loci_intersect = list(set(binomial_all_fail_proprtion.keys()) &
                                        set(depth_all_fail_proportion.keys()))
    print(f"{len(binomial_and_depth_loci_intersect)} matching locations between binomial and depth fail list.")

    # Filtered loci only BAF and depth
    baf_and_depth_loci_intersect = list(set(baf_all_fail_proportion.keys()) &
                                    set(depth_all_fail_proportion.keys()))
    print(f"{len(baf_and_depth_loci_intersect)} matching locations between BAF and depth fail list.")

    # Filtered loci only BAF and depth and binomial
    baf_depth_binomial_intersect = list(set(binomial_all_fail_proprtion.keys()) &
                                    set(baf_all_fail_proportion.keys()) &
                                    set(depth_all_fail_proportion.keys()))
    print(f"{len(baf_depth_binomial_intersect)} matching locations between binomial, BAF and depth fail list.")


if __name__ == "__main__":
    main()
