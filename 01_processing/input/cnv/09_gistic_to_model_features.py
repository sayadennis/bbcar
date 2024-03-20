# pylint: disable=missing-module-docstring)

import os

import process_functions

din = "/projects/b1131/saya/bbcar/data/02b_cnv/09_gistic2_out_conf90"
dout = "/projects/b1131/saya/bbcar/data/02b_cnv/10_cleaned_cnv"

cnv_dict = process_functions.generate_gistic_features(din)

if not os.path.isdir(dout):
    os.mkdir(dout)

for key, df in cnv_dict.items():
    df.to_csv(f"{dout}/{key}_conf90.csv", header=True, index=True)
