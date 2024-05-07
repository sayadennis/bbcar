#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="somvcf"
#SBATCH --output=/projects/b1131/saya/bbcar/out/create_somatic_vcfs.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd /projects/b1131/saya/bbcar/

#### Copy the tissue-normal matched calls (100% somatic) ####

din="/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls"

for pon_source in bbcar; do # 1000g
    dout=/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/vcfs/${pon_source}
    mkdir -p ${dout}
    cp ${din}/tumor_normal/*_DPfiltered_${pon_source}pon.vcf ${dout}/
    #### Filter the tissue-only calls and select predicted-somatic ####
    python ~/bbcar/repo/01_processing/input/mutation/08_feature_matrix/create_somatic_vcfs.py $pon_source
done
