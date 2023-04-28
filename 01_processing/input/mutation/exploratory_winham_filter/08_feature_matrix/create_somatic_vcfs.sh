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

pon_source=bbcar
for filter_type in liberal classical strict; do
    din=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/04_ml_features
    dout=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/07_predicted_somatic/vcfs
    mkdir -p ${dout}
    # First copy the tissue-normal matched calls (100% somatic) _DPfiltered_strictAF_bbcarpon.vcf
    cp /projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/02_variant_calls/tumor_normal/*_DPfiltered_${filter_type}AF_${pon_source}pon.vcf ${dout}/
    #### Filter the tissue-only calls and select predicted-somatic ####
    python ~/bbcar/repo/01_processing/input/mutation/exploratory_winham_filter/08_feature_matrix/create_somatic_vcfs.py $filter_type
done
