#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-307
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=avfts_%a
#SBATCH --output=/projects/b1131/saya/bbcar/out/generate_annovar_features_%a.out

module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/06_generate_ml_features/

# #### To generate 1000g PON input_args.txt file ####
# cd /projects/b1131/saya/bbcar/data/02a_mutation/
# touch sample_names_all_ml_feature_generation_1000gpon.txt # create file 
# > sample_names_all_ml_feature_generation_1000gpon.txt # clear content of file (in case of re-running)
# avdir="/projects/b1131/saya/bbcar/data/02a_mutation/03_annotated_variants/annovar"
# ls ${avdir}/germline_only/*.hg38_multianno.vcf | grep [0-9].hg38_multianno.vcf >> sample_names_all_ml_feature_generation_1000gpon.txt # germline doesn't have PON info in filename
# for category in tumor_normal tumor_only; do 
#     ls ${avdir}/${category}/*_1000gpon.hg38_multianno.vcf >> sample_names_all_ml_feature_generation_1000gpon.txt
# done
# #########################################
# #### To generate BBCAR PON input_args.txt file ####
# cd /projects/b1131/saya/bbcar/data/02a_mutation/
# touch sample_names_all_ml_feature_generation_bbcarpon.txt # create file 
# > sample_names_all_ml_feature_generation_bbcarpon.txt # clear content of file (in case of re-running)
# avdir="/projects/b1131/saya/bbcar/data/02a_mutation/03_annotated_variants/annovar"
# ls ${avdir}/germline_only/*.hg38_multianno.vcf | grep [0-9].hg38_multianno.vcf >> sample_names_all_ml_feature_generation_bbcarpon.txt # germline doesn't have PON info in filename
# for category in tumor_normal tumor_only; do 
#     ls ${avdir}/${category}/*_bbcarpon.hg38_multianno.vcf >> sample_names_all_ml_feature_generation_bbcarpon.txt
# done
# #########################################

for pon_source in bbcar 1000g; do
    IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_all_ml_feature_generation_${pon_source}pon.txt
    python 01_generate_annovar_features.py ${input_args[$SLURM_ARRAY_TASK_ID]}
done

# #####################
# #### Generic PON ####
# #####################

# IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_all_ml_feature_generation.txt

# python 01_generate_annovar_features.py ${input_args[$SLURM_ARRAY_TASK_ID]}

# ###################
# #### BBCAR PON ####
# ###################

# IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_all_ml_feature_generation_bbcarpon.txt

# python 01_generate_annovar_features.py ${input_args[$SLURM_ARRAY_TASK_ID]}
