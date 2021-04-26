import os
import sys
import numpy as np
import pandas as pd

"""
The purpose of this module is to provide the tools to filter the BBCar CNV matrices using a given region list.
The CNV matrix (cnv_mx) should:
    - be a Pandas DataFrame with index = subject IDs, and columns = [HUGO gene names | region name]
    - entries should be float (copy number) or int (threshold value)
The region file (regions_fn) should: 
    - be a string specifying an absolute or relative path to the file
    - have 5 columns: #chrom, chromStart, chromEnd, name, amp_del. (Doesn't have to be this name?)
    - amp_del should take values "+", "-", or "."
Reference file (ref_fn) can be two different things depending on whether features="genes" or features="regions:
    - either way, ref_fn should be a string specifying an absolute or relative path to the file
    - if features="genes", ref_fn file should:
        - have any number of columns but the indices should be (when we include the first 'bin' as column):
            - chromosome name: 2
            - gene name: 12
            - start pos: 4
            - end pos: 5
    - if features="regions", ref_fn file should:
        - be the file that documents our BBCar sample region boundaries. specifically, the format should be
            - column 1: chromosome name in string
            - column 2: region start position
            - column 3: region end position
            - column 4: region name e.g. "Amplification Peak 1"
            - column 5: whether this is an amplification or deletion peak 
"""

def amp_del_to_numeric(amp_del):
    if amp_del == "+":
        amp_del = 1
    elif amp_del == "-":
        amp_del = -1
    else:
        amp_del = 0
    return amp_del

def filter_mx(cnv_mx, regions_fn, ref_fn, features="genes"):
    if features=="genes":
        print("Filtering gene features...")
        regions = pd.read_csv(regions_fn, header=0, index_col=None, sep="\t")
        gene_coor = pd.read_csv(ref_fn, sep="\t", header=None, index_col=0)
        gene_coor = gene_coor.iloc[:,np.array([1,3,4,11])]
        gene_coor.columns = ["chrom", "start", "end", "name"]
        # Get a list of gene names that are overlapping with any regions in your regions list
        gene_names = []
        amp_del_dict = {}
        for i in range(regions.shape[0]):
            chrn, start, end, name = regions.iloc[i,0], regions.iloc[i,1], regions.iloc[i,2], regions.iloc[i,3]
            amp_del = amp_del_to_numeric(regions.iloc[i,4])
            gene_coor_sub = gene_coor.iloc[(gene_coor["chrom"]==chrn).values] # subset of the gene_coor where chromosome name matches the region of interest
            start_in_range = ((start > gene_coor_sub["start"]) & (start < gene_coor_sub["end"]))
            end_in_range = ((end > gene_coor_sub["start"]) & (end < gene_coor_sub["end"]))
            gene_coor_sub = gene_coor_sub.iloc[(start_in_range | end_in_range).values]
            gene_names_sub = gene_coor_sub["name"].values
            for gene_name in gene_names_sub:
                if gene_name not in gene_names:
                    gene_names.append(gene_name)
                    amp_del_dict[gene_name] = amp_del # 1, 0, or -1
                else:
                    continue
        print("Number of genes that mapped to the regions specified by file {}: {}".format(regions_fn, len(gene_names)))
        # subset the cnv_mx based on the HUGO gene symbols in the remaining gene_regions
        filtered_cnv = cnv_mx.iloc[:,[x in gene_names for x in list(cnv_mx.columns)]]
        for colname in filtered_cnv.columns:
            if amp_del_dict[colname] == 1:
                for sample in (filtered_cnv[filtered_cnv[colname] < 0]).index:
                    filtered_cnv.loc[sample][colname] = 0
            elif amp_del_dict[colname] == -1:
                for sample in (filtered_cnv[filtered_cnv[colname] > 0]).index:
                    filtered_cnv.loc[sample][colname] = 0
            else:
                continue
        print("Number of gene features in the CNV matrix that survived: {}\n".format(filtered_cnv.shape[1]))
    elif features=="regions":
        # Read the BBCar feature region file
        ref = pd.read_csv(ref_fn, index_col=None, sep="\t")
        regions = pd.read_csv(regions_fn, header=0, index_col=None, sep="\t")
        reg_names = [] # I will store the names of the regions to keep in this list
        amp_del_dict = {}
        for i in range(ref.shape[0]):
            ref_chrn, ref_start, ref_end, ref_name = ref.iloc[i,0], ref.iloc[i,1], ref.iloc[i,2], ref.iloc[i,3]
            for j in range(regions.shape[0]):
                reg_chrn, reg_start, reg_end, reg_name = regions.iloc[j,0], regions.iloc[j,1], regions.iloc[j,2], regions.iloc[j,3]
                amp_del = amp_del_to_numeric(regions.iloc[j,4])
                if ref_chrn==reg_chrn:
                    if (((ref_start > reg_start) & (ref_start < reg_end)) | ((ref_end > reg_start) & (ref_end < reg_end))):
                        if ref_name not in reg_names:
                            reg_names.append(ref_name)
                            amp_del_dict[ref_name]
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
        print("Number of peaks that overlapped with regions specified by file {}: {}".format(regions_fn, len(reg_names)))
        # subset the cnv_mx based on the HUGO gene symbols in the remaining gene_regions
        col_bool = []
        for item in reg_names:
            if np.any([x.startswith(item) for x in cnv_mx.columns]):
                col_bool.append(True)
            else:
                col_bool.append(False)
        filtered_cnv = cnv_mx.iloc[:,]
        print("Number of peaks in the CNV matrix that survived: {}\n".format(filtered_cnv.shape[1]))
        for colname in filtered_cnv.columns:
            if amp_del_dict[colname] == 1:
                for sample in (filtered_cnv[filtered_cnv[colname] < 0]).index:
                    filtered_cnv.loc[sample][colname] = 0
            elif amp_del_dict[colname] == -1:
                for sample in (filtered_cnv[filtered_cnv[colname] > 0]).index:
                    filtered_cnv.loc[sample][colname] = 0
            else:
                continue
    return filtered_cnv

def filter_by_genes(cnv_mx, gene_names_fn):
    with open(gene_names_fn, "r") as f:
        lines = f.readlines()
    genelist = []
    for line in lines:
        genelist.append(line.rstrip())
    print("Number of genes present in file {}: {}".format(gene_names_fn, len(genelist)))
    filtered_cnv = cnv_mx.iloc[:,[x in genelist for x in list(cnv_mx.columns)]]
    print("Number of genes that survived from matrix: {}".format(filtered_cnv.shape[1]))
    return filtered_cnv

