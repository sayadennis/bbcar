# Modeling

We have two categories of predictor features (mutations and CNVs from benign breast biopsies) and one binary outcome variable (whether the patient later develops breast cancer). Here, we investigate the predictive value of each of the two categories of predictor features, as well as investigate several fusion approaches: early, late, and mid-level fusion.

## No fusion models

No fusion models are ones where we use only mutation data or only CNV data to predict the outcome. There are four early fusion models from the two categories of input:
* Mutational data ([script](https://github.com/sayadennis/bbcar/blob/master/04_modeling/mutation/bbcar_mutsig_classicalml.sh))
  * 96-element feature vector indicating the accumulated mutational context
  * COSMIC Mutational Signatures derived from the above feature matrix
* CNV data ([script](https://github.com/sayadennis/bbcar/blob/master/04_modeling/cnv/cn_signature_classicalml.sh))
  * 26-element feature vector indicating the accumulated CNV context
  * In-house CNV Signatures derived from the above feature matrix

Note: COSMIC CN Signatures could not be derived for this dataset since matched germline was available for only a subset of samples, and one of the processing software tools (ASCAT) requires matched germline for WXS data. (See: [here](https://github.com/VanLoo-lab/ascat/blob/fd697b443d5063225ee52a8739bf02ae53d3d1a6/README.md#new-features-in-v3) - "_Please note that `ascat.prepareHTS` does require tumour/normal pairs and is unable to process unmatched tumours for now_") 

## Early fusion models

Early fusion models are ones where we simply concatenate the input matrices. There are two possible early fusion models:
* Concatenated mutation + CNV data from the context feature vectors (96 + 26 dimensions)
* Concatenated signatures of mutation + CNV (dimension depends on the number of signatures in each modality)

Training script can be found [here](https://github.com/sayadennis/bbcar/blob/master/04_modeling/mut_cnv_combined_classicalml.sh).

## Late fusion models

Late fusion models take the average of the prediction scores from models trained separately on the two modalities of the data: mutation and CNV. This can be a simple or a weighted average. However, in our case, the optimal weight could not be determined from cross-validation since the performance way always 1.0 on the training set. We are currently proceeding with a simple average (our dataset is small so best to keep the number of parameters minimal), but the overfitting behavior is worth looking into. 

Training script can be found [here](https://github.com/sayadennis/bbcar/blob/master/05_evaluation/combined_mut_cnv_late_fusion.py).

## Mid-level fusion models

### Supervised hybrid non-negative matrix factorization (Supervised hNMF)



