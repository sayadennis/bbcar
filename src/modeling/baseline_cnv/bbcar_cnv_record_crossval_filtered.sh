#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarcnvfilter"
#SBATCH --output=bbcar/out/bbcar_cnv_record_crossval_filtered.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

inputdir='/projects/b1122/saya/07_selected_region_cnv'
labeldir='/projects/b1122/saya'
outdir='/home/srd6051/bbcar/model_performance/results_classicalml'
ixdir='/projects/b1122/saya/indices'

for fn in $(ls $inputdir/*.csv); do
    pathlist=$(echo $fn | tr '/' '\n') # split by "/"
    patharr=($pathlist)
    shortfn=${patharr[4]} # get just the filename and not the whole path
    fout=results_${shortfn}
    python classical_ml/ClassicalML/run_classical_ml.py \
        --input $fn \
        --label $labeldir/bbcar_label_studyindex.csv \
        --outfn $outdir/$fout \
        --indexdir $ixdir
done

# ## Region-level copy number
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/reg_copy_conf90_studyindex.csv \
#     --label $labeldir/bbcar_label_studyindex.csv \
#     --outfn $outdir/bbcar_regcopy.csv \
#     --indexdir $ixdir

# # --mem=0

# ## Region-level threshold values
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/reg_thres_conf90_studyindex.csv \
#     --label $labeldir/bbcar_label_studyindex.csv \
#     --outfn $outdir/bbcar_regthres.csv \
#     --indexdir $ixdir
# #

# ## Gene-level copy number
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/gene_copy_conf90_studyindex.csv \
#     --label $labeldir/bbcar_label_studyindex.csv \
#     --outfn $outdir/bbcar_genecopy.csv \
#     --indexdir $ixdir
# #

# ## Region-level threshold values
# python classical_ml/ClassicalML/run_classical_ml.py \
#     --input $inputdir/gene_thres_conf90_studyindex.csv \
#     --label $labeldir/bbcar_label_studyindex.csv \
#     --outfn $outdir/bbcar_genethres.csv \
#     --indexdir $ixdir
# #
