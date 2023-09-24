## This script only needs to be run once

module purge all
module load singularity

cd ~/bbcar/repo/01_processing/input/cnv/signatures/

singularity pull docker://opengenomics/battenberg:2.2.9
