# testing out script by using UCSC's refGene for gene annotation

import os
import sys
import glob
import numpy as np
import pandas as pd

din="/projects/b1122/saya/bbcar_project/gatk_analyzed_segments"
dout="/projects/b1122/saya/bbcar_project/cnv_matrix"
dref="/projects/b1122/saya/bbcar_project/geneanno_refs"

## Store gene-annotation references
chrref = {}
chrnlist = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX']

for chrn in chrnlist:
    chrref[chrn] = pd.read_csv(os.path.join(dref, "GRCh38_gencode.v27.refFlat." + chrn + ".txt"), sep="\t", header=None)

## Define useful functions 
def generate_genedict(filepath):
    genedict = {
        # will fill in data in the format:
        # gene_name : [list, of, log2_copy_ratios]
    }
    data = pd.read_csv(filepath)
    subid = filepath.split("/")[-1].split(".")[0] # find patient ID # FIX THIS PART if necessary
    data = data[data["CALL"]!="0"] # drop all rows where calls are not +/- # i think "0" should be string but double check
    for ix in data.index:
        # find chromosome via "CONTIG" column 
        chrom = data.loc[ix]["CONTIG"]
        segstart = data.loc[ix]["START"]
        segend = data.loc[ix]["END"]
        ref = chrref[chrom]
        # find overlapping gene name
        for j in range(ref.shape[0]):
            genestart = ref.iloc[j,4]
            geneend = ref.iloc[j,5]
            # note - "start" and "end" are not necessarily the starting and finishing point of transcription
            # it can be the reverse order depending on the strand
            # in this context, "start" is simply the position with the smaller genomic coordinates and "end" the larger
            if (
                ( # when start of the segment is within gene range
                    (segstart > genestart) & (segstart < geneend)
                ) or ( # when end of the segment is within gene range
                    (segend > genestart) & (segend < geneend)
                ) or ( # when the gene is fully contained in the segment
                    (segstart < genestart) & (segend > geneend)
                )
            ): # if there is overlap
                ol_genename = ref.iloc[j,0] # overlapping gene name 
                copyratio = data.loc[ix]["MEAN_LOG2_COPY_RATIO"]
                if ol_genename in genedict.keys(): # if gene name already in genelist:
                    genedict[ol_genename].append(copyratio) # add gene to list and store copy ratio value 
                else: # if gene name not in genelist:
                    genedict[ol_genename] = [copyratio] # add gene to list and store copy ratio value 
            else:
                continue
    return subid, genedict

def take_maxabs(genedict): ## Take max/min value of amplification/deletion per gene
    newdict = {}
    for gene in genedict.keys():
        newdict[gene] = genedict[gene][np.argmax(np.absolute(genedict[gene]))]
    return newdict

def write_mx(genedict, subid, mx):
    ## Create a new row with the subject ID
    if subid in mx.index:
        sys.exit("Error: subject {} already in dataframe. Cannot overwrite.".format(subid))
    else:
        mx = pd.concat([mx, pd.DataFrame(np.zeros((1, mx.shape[1])), index=[subid], columns=mx.columns)])
    ## Loop through genes in the dictionary 
    for gene in genedict.keys():
        # if gene doesn't already exist in mx, add a column
        if gene in mx.columns:
            mx.loc[subid][gene] = genedict[gene]
        else:
            mx[gene] = np.zeros((mx.shape[0]))
            mx.loc[subid][gene] = genedict[gene]
    return mx


## Now iterate through segment files and fill in the CNV matrix 

# initialize empty matrix
mx = pd.DataFrame(None)

for filepath in glob.glob(os.path.join(din, "*.csv")):
    # Generate gene dict for this subject's segment file
    subid, genedict = generate_genedict(filepath)
    # Aggregate data by taking the maximum amplification or minimum deletion value 
    genedict = take_maxabs(genedict)
    # Write info from genedict to this matrix
    mx = write_mx(genedict, subid, mx)
    print("Done writing for file: {}".format(filepath))

# Overwrite CNV gene-level feature matrix file 
mx.to_csv(os.path.join(dout, "gatk_cnv_gene_matrix.csv"), header=True, index=True)

# why am i seeing same values repeated multiple times in the resulting genedict? 
# why did the results look different between when i looped through segments first versus when i looped through genes first? 
