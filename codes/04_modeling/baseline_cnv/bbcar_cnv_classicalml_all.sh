#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcar_cnv_classicalml_all"
#SBATCH --output=bbcar/out/bbcar_cnv_classicalml_all.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

inputdir='/projects/b1122/saya/04_cleaned_cnv/all'
labeldir='/projects/b1122/saya'
outdir='/home/srd6051/bbcar/model_performance/results_classicalml/all'
ixdir='/projects/b1122/saya/indices'

#############################
#### Gene-level features ####
#############################

# ## Gene-level copy number
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/gene_copy_conf90_all.csv \
#     --label $labeldir/bbcar_label_studyid.csv \
#     --outfn $outdir/bbcar_genecopy_roc_auc.csv \
#     --indexdir $ixdir \
#     --scoring roc_auc \
#     --savemodel "true"
# #

# ## Gene-level threshold values
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/gene_thres_conf90_all.csv \
#     --label $labeldir/bbcar_label_studyid.csv \
#     --outfn $outdir/bbcar_genethres_roc_auc.csv \
#     --indexdir $ixdir \
#     --scoring roc_auc \
#     --savemodel "true"
# #

# ###############################
# #### Region-level features ####
# ###############################

# ## Region-level copy number
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/reg_copy_conf90_all.csv \
#     --label $labeldir/bbcar_label_studyid.csv \
#     --outfn $outdir/bbcar_regcopy_roc_auc.csv \
#     --indexdir $ixdir \
#     --scoring roc_auc \
#     --savemodel "true"

# # --mem=0

# ## Region-level threshold values
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/reg_thres_conf90_all.csv \
#     --label $labeldir/bbcar_label_studyid.csv \
#     --outfn $outdir/bbcar_regthres_roc_auc.csv \
#     --indexdir $ixdir \
#     --scoring roc_auc \
#     --savemodel "true"
# #

#################################
#### Cytoband-level features ####
#################################

python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_all.csv \
    --label $labeldir/bbcar_label_studyid.csv \
    --outfn $outdir/bbcar_cytocopy_roc_auc.csv \
    --indexdir $ixdir \
    --scoring roc_auc \
    --savemodel "true"
#

# ################################
# #### Pathway-level features ####
# ################################

# inputdir='/projects/b1122/saya/08_pathway_aggregated_cnv'

# for fn in $(ls $inputdir/*.csv); do
#     pathlist=$(echo $fn | tr '/' '\n') # split by "/"
#     patharr=($pathlist)
#     shortfn=${patharr[4]} # get just the filename and not the whole path
#     fout=results_${shortfn}
#     python classical_ml/ClassicalML/run_classical_ml.py \
#         --input $fn \
#         --label $labeldir/bbcar_label_studyid.csv \
#         --outfn $outdir/$fout \
#         --indexdir $ixdir \
#         --scoring roc_auc
# done

# #####################################
# #### Heuristic-filtered features ####
# #####################################

# inputdir='/projects/b1122/saya/07_selected_region_cnv'

# for fn in $(ls $inputdir/*.csv); do
#     pathlist=$(echo $fn | tr '/' '\n') # split by "/"
#     patharr=($pathlist)
#     shortfn=${patharr[4]} # get just the filename and not the whole path
#     fout=results_${shortfn}
#     python classical_ml/ClassicalML/run_classical_ml.py \
#         --input $fn \
#         --label $labeldir/bbcar_label_studyid.csv \
#         --outfn $outdir/$fout \
#         --indexdir $ixdir \
#         --scoring roc_auc
# done
