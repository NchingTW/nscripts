#!/bin/env python

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy as scipy
from scipy import stats
import math
import pickle
import pybedtools
from argparse import ArgumentParser
import os

from collections import OrderedDict, Counter
from collections import defaultdict
from clipper.src import CLIP_analysis
from clipper.src import CLIP_analysis_display
from gscripts.rnaseq import helpers
from gscripts import GO
from gscripts.general import parsers
from gscripts.general import region_helpers
from clipper.src import get_genomic_regions




def make_region_dictionary():
    regions = OrderedDict()
    regions["cds"] = "CDS"
    regions["three_prime_utrs"] = "3' UTR"
    regions["five_prime_utrs"] = "5' UTR"
    regions["proxintron500"] = "Proximal\nIntron"
    regions["distintron500"] = "Distal\nIntron"
    return regions

def pd_read_peaks_file(peak_bed_file):
    """
    Args:
        peak_bed_file: the bed file after normalization, which contains 'chr','start','stop' on the first three columns
    Returns: Dataframe of bed file
    """
    peak_df = pd.read_table(peak_bed_file)
    return peak_df

def make_RBP_name_as_string(RBP):
    RBP = str(RBP)
    return RBP


def save_datafeame_to_bed_file(dataframe, output_folder,RBP):
    save_bed_file = pybedtools.BedTool.from_dataframe(dataframe).saveas(output_folder+'/'+str(RBP)+'.peak.bed')
    return save_bed_file

def read_bed_through_bedtools(bed_file):
    bedtools_read_bed_file = pybedtools.BedTool(bed_file)
    return bedtools_read_bed_file

def make_bed_for_region_assigned(peak_bed_file):
    """
    Args:
        peak_bed_file: the bed file after normalization, which contains 'chr','start','stop'
    Returns:

    """
    peak_df = pd_read_peaks_file(peak_bed_file)
    bed_df = peak_df.iloc[:, 0:3]
    return bed_df


def assigned_peaks_region(peak_bed_file, output_folder, number_of_permutation,RBP):

    """

    Args:
        peak_bed_file: The bed file after normalization, which contains 'chr','start','stop' three columns
        output_folder: output folder for saving  bed file
        number_of_permutation: number of background(random) bed file you want to create
        RBP: RBP name
    Returns: Random and real bed file for 5utr,cds,3utr,proximal Intron,Distal intron and combine all region

    """

    make_region_dictionary()
    bed_df = make_bed_for_region_assigned(peak_bed_file)
    save_datafeame_to_bed_file(bed_df, output_folder,RBP)
    bedtools_read_bed_file = read_bed_through_bedtools(output_folder+'/'+str(RBP)+'.peak.bed')
    assinged_region = CLIP_analysis.assign_to_regions(bedtools_read_bed_file, assigned_dir=output_folder, nrand=number_of_permutation, species="hg19")
    return assinged_region


def main():

    usage = "I don't know what I should put"


    parser = ArgumentParser(usage)

    parser.add_argument(
        "--peaks",
        dest="peak_bed_file",
        help="peak bed file to assign the regions on these peaks, and use these peaks to created the random background."
    )


    parser.add_argument(
        "-o",
        "--output",
        dest="output_folder",
        help="The folder to save the assinged region bed files, and background bed file."
    )

    parser.add_argument(
        "--n",
        dest="number_of_permutation",
        help="set the number ofo permutation. (number of random bed files you want to created)"
    )

    parser.add_argument(
        "--RBP",
        dest="RBP",
        help="RBP name",
        type=str ##not sure this is correct or nor
    )

    opts = parser.parse_args()

    assigned_peaks_region(opts.peak_bed_file, opts.output_folder, opts.number_of_permutation,opts.RBP)

if __name__ == "__main__":
    main()
