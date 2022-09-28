#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-307
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=avfts_%a
#SBATCH --output=bbcar/out/generate_annovar_features_%a.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

# #### To generate input_args.txt file ####
# touch /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation.txt # create file 
# > /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation.txt # clear content of file (in case of re-running)
# cd /projects/b1131/saya/bbcar/03_annotated_variants/annovar
# for category in germline_only tumor_normal tumor_only; do ls ${category}/*.hg38_multianno.vcf | grep [0-9t].hg38_multianno.vcf >> /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation.txt ; done
# #########################################
# #### To generate BBCAR PON input_args.txt file ####
# touch /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation_bbcarpon.txt # create file 
# > /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation_bbcarpon.txt # clear content of file (in case of re-running)
# cd /projects/b1131/saya/bbcar/03_annotated_variants/annovar
# ls germline_only/*.hg38_multianno.vcf | grep [0-9t].hg38_multianno.vcf >> /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation_bbcarpon.txt # germline doesn't have "bbcarpon" in filename
# for category in tumor_normal tumor_only; do ls ${category}/*_bbcarpon.hg38_multianno.vcf >> /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation_bbcarpon.txt ; done
# #########################################

#####################
#### Generic PON ####
#####################

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation.txt

python bbcar/codes/somatic_mutations/02_processing/05_generate_ml_features/generate_annovar_features.py ${input_args[$SLURM_ARRAY_TASK_ID]}

###################
#### BBCAR PON ####
###################

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_all_ml_feature_generation_bbcarpon.txt

python bbcar/codes/somatic_mutations/02_processing/05_generate_ml_features/generate_annovar_features.py ${input_args[$SLURM_ARRAY_TASK_ID]}
