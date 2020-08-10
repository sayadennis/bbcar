import os
import sys
import glob
import numpy as np
import pandas as pd

def CytobandMatrix(cytoband, cna_dir):
    """
    This function allows you to build a cytoband matrix.
    Arguments:
    cytoband ––a Pandas DataFrame that contains columns "name", "CONTIG" (chromosome name), "START", and "END".
    cna_dir ––an absolute or relative path to the directory that contains all patient's files of CNA.

    Usage:
    from BBCarAnalysis import CytobandMatrix
    cytoband = pd.read_csv("../cytoband_hg19.tsv", sep="\t")
    cyto_mx = CytobandMatrix(cytoband=cytoband, cna_dir="../../data")
    """
    cytoband_names = []
    for i in cytoband.index:
        chrom = cytoband.loc[i]["#chrom"]
        name = cytoband.loc[i]["name"]
        bandname = chrom[3:] + name
        cytoband_names.append(bandname)

    new_cytoband = pd.DataFrame(cytoband.iloc[:,:3], index=cytoband_names)
    for col in ["#chrom", "chromStart", "chromEnd"]:
        new_cytoband[col] = cytoband[col].to_list()
    cytoband = new_cytoband

    cna_matrix = pd.DataFrame(columns=[cytoband_names])

    for filepath in glob.glob(cna_dir + "/*"):
        data = pd.read_csv(filepath)
        sample_name = data["name"].values[0][:4] # captures patient ID from the first row of the "name" column in the CNA data table
        cna_matrix = cna_matrix.append(
            pd.DataFrame(np.zeros((1, len(cytoband_names)), dtype=int), index=[sample_name], columns=[cytoband_names])
        )

        for i in data.index:
            # see what chr the CNA from that row is from 
            chrom = data.loc[i]["CONTIG"]
            # iterate through cytobands in that row
            cyto_chrom = cytoband.loc[cytoband["#chrom"] == chrom][:]
            for cyto_i in cyto_chrom.index:
                # if potitional information matches
                if (data.loc[i]["START"] < cyto_chrom.loc[cyto_i]["chromEnd"]) | (data.loc[i]["END"] < cyto_chrom.loc[cyto_i]["chromStart"]):
                    # add the MEAN_LOG2_COPY_RATIO value to the corresponding cytoband column in cna_matrix 
                    meanlog2val = data.loc[i]["MEAN_LOG2_COPY_RATIO"]
                    cna_matrix.loc[sample_name][cyto_i] += meanlog2val
                else:
                    continue
    
    return cna_matrix

