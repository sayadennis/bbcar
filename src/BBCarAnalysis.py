import os
import sys
import glob
import numpy as np
import pandas as pd

def CytobandMatrix(cytoband, cna_dir):
    """
    This function allows you to build a cytoband matrix.
    Arguments:
    cytoband - a Pandas DataFrame that contains columns "name", "CONTIG" (chromosome name), "START", and "END".
    cna_dir - an absolute or relative path to the directory that contains all patient's files of CNA.
    Usage:
    from BBCarAnalysis import CytobandMatrix
    cytoband = pd.read_csv("../cytoband_hg19.tsv", sep="\t")
    cyto_mx = CytobandMatrix(cytoband=cytoband, cna_dir="../../data")
    """
    # Generate a list of cytoband names (chromosome number + cytoband name e.g. "1p36.33" for chromosome 1 cytoband p36.33) 
    cytoband_names = []
    for i in cytoband.index:
        chrom = cytoband.loc[i]["#chrom"]
        name = cytoband.loc[i]["name"]
        bandname = str(chrom[3:]) + str(name)
        cytoband_names.append(bandname)

    # Generate a new cytoband dataframe where the indices are cytoband_names 
    new_cytoband = pd.DataFrame(cytoband.iloc[:,:3], index=cytoband_names)
    for col in ["#chrom", "chromStart", "chromEnd"]:
        new_cytoband[col] = cytoband[col].tolist()
    cytoband = new_cytoband

    # Generate an empty CNV matrix with cytoband_names as columns 
    cna_matrix = pd.DataFrame(columns=[cytoband_names])

    # Iterate through files in the cna_dir and start filling in cna_matrix 
    for filepath in glob.glob(cna_dir + "/*.csv"):
        data = pd.read_csv(filepath)
        # CHANGE BELOW TO ADAPT TO VARYING LENGTH OF SAMPLE IDs AND CONTROL/TISSUE 
        sample_name = filepath.split("/")[-1].split(".")[0] # either "<sampleid>_Tissue" or "<sampleid>_Control" 
        # Append new row with all zeros to cna_matrix 
        cna_matrix = cna_matrix.append(
            pd.DataFrame(np.zeros((1, len(cytoband_names)), dtype=int), index=[sample_name], columns=[cytoband_names])
        )
        # loop through rows of data and fill in mean log2 copy ratio to cna_matrix 
        for i in data.index:
            # see what chr the CNA from that row is from 
            chrom = data.loc[i]["CONTIG"]
            # iterate through cytobands in that row
            cyto_chrom = cytoband.loc[cytoband["#chrom"] == chrom][:]
            for cyto_i in cyto_chrom.index:
                # if positional information matches
                if (data.loc[i]["START"] < cyto_chrom.loc[cyto_i]["chromEnd"]) | (data.loc[i]["END"] < cyto_chrom.loc[cyto_i]["chromStart"]):
                    # add the MEAN_LOG2_COPY_RATIO value to the corresponding cytoband column in cna_matrix 
                    meanlog2val = data.loc[i]["MEAN_LOG2_COPY_RATIO"]
                    cna_matrix.loc[sample_name][cyto_i] += meanlog2val
                else:
                    continue
    
    return cna_matrix

