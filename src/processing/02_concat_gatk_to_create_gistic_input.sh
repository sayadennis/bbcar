#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name='catgistic'
#SBATCH --output=bbcar/out/210127_concat_gistic_input.out

combfile='/projects/b1122/saya/bbcar/gistic2_input/combined_gistic_input.tsv'
header='/projects/b1122/saya/bbcar/gistic2_input/gistic_input_header.tsv'

cat $header >> $combfile

for fn in $(ls /projects/b1122/saya/bbcar_project/gistic2_input/*.csv); do
    cat $fn >> $combfile
done
