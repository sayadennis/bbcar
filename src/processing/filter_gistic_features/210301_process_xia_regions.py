import numpy as np
import pandas as pd

with open("/Users/sayadennis/Projects/bbcar_project/bbcar_repository_staging/regions_xia.txt", "r") as f:
    lines = f.readlines()


newlist=[]
for line in lines:
    newlist.append(line.rstrip())

lines = newlist

reg = pd.DataFrame(None, index=[], columns=["#chrom", "chromStart", "chromEnd", "name", "amp_del"])

for line in lines:
    if line.startswith("chr"):
        chrname = line.split(":")[0]
        start = int(line.split(":")[1].split(".")[0].split("-")[0])
        end = int(line.split(":")[1].split(".")[0].split("-")[1])
        name = line[(len(chrname)+1+len(line.split(":")[1].split(".")[0]))+1:]
        # decide amp/del value
        if "amp" in name:
            amp_del = "+"
        elif "del" in name:
            amp_del = "-"
        else:
            amp_del = "."
        # write to dataframe
        reg = reg.append(pd.DataFrame([[chrname, start, end, name, amp_del]], index=[1], columns=reg.columns))
    else:
        chrname = line.split(".")[2].split(":")[0]
        start = int(line.split(".")[2].split(":")[1].split("-")[0])
        end = int(line.split(".")[2].split(":")[1].split("-")[1])
        if len(line.split()) > 1:
            name1 = line.split()[0]
            name2 = line.split()[1].split(".")[0]
            name = name1 + "_" + name2
        else:
            name1 = line.split(".")[0]
            name2 = line.split(".")[1]
            name = name1 + "_" + name2
        # decide amp/del 
        if "amp" in name:
            amp_del = "+"
        elif "del" in name:
            amp_del = "-"
        else:
            amp_del = "."
        # write to dataframe
        reg = reg.append(pd.DataFrame([[chrname, start, end, name, amp_del]], index=[1], columns=reg.columns))

reg.to_csv("/Users/sayadennis/Projects/bbcar_project/bbcar_repository_staging/regions_xia.tsv", sep="\t", header=True, index=False)
