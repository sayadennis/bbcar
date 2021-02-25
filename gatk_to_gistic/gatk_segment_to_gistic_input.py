import numpy as np
import pandas as pd
import itertools
import glob

combmx = pd.DataFrame(None, columns=["name", "CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO"])

for fn in glob.glob("/projects/b1122/saya/bbcar_project/gatk_analyzed_segments/*.csv"):
    samplen = fn.split("/")[-1].split(".")[0]
    data = pd.read_csv(fn)
    newdata = data.iloc[:,:6]
    newdata["name"] = list(itertools.repeat(samplen, len(newdata)))
    combmx = pd.concat((combmx, newdata))

combmx["MEAN_LOG2_COPY_RATIO"] = combmx["MEAN_LOG2_COPY_RATIO"] - 1

combmx.to_csv("/projects/b1122/saya/bbcar_project/gistic2_input/combined_gistic_input.tsv", sep="\t", header=False, index=False)
