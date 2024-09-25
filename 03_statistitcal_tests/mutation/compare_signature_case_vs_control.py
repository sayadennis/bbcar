import pandas as pd
from scipy.stats import ranksums

projdir = "/projects/b1131/saya/new_bbcar"
mutdir = f"{projdir}/data/02a_mutation/08_feature_matrix"
vardir = f"{projdir}/data/02a_mutation/02_variant_calls"
dout_stat = f"{projdir}/out"
dout_plot = f"{projdir}/plots"

######################################
#### Get case/control assignments ####
######################################

label = pd.read_csv(f"{projdir}/label_all.csv", index_col=0)
sample_ids = pd.read_csv(
    f"{mutdir}/denovo_signature_per_sample_SBS96.csv", index_col=0
).index

###############################################
#### Loop through sigs and perform t-tests ####
###############################################

stat_results = pd.DataFrame(
    columns=["Signature Type", "Signature name", "statistic", "p-value"]
)

for sigtype in ["SBS96", "DBS78", "ID83"]:
    #### De Novo Signatures ####
    signature = pd.read_csv(
        f"{mutdir}/denovo_signature_per_sample_{sigtype}.csv", index_col=0
    )
    signature_cases = signature.iloc[
        [x in label.iloc[label.values == 1, :].index for x in signature.index], :
    ]
    signature_controls = signature.iloc[
        [x in label.iloc[label.values == 0, :].index for x in signature.index], :
    ]
    #
    signature_cases_ratio = signature_cases.divide(
        signature_cases.sum(axis=1), axis="rows"
    )
    signature_controls_ratio = signature_controls.divide(
        signature_controls.sum(axis=1), axis="rows"
    )
    #
    # print('\n#### Results for raw values ####')
    results = ranksums(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        stat_results = pd.concat(
            (
                stat_results,
                pd.DataFrame(
                    {
                        "Signature Type": f"{sigtype} De Novo",
                        "Signature name": signame,
                        "statistic": t,
                        "p-value": p,
                    },
                    index=[0],
                ),
            )
        )
    #
    # print('\n#### Results for ratios ####')
    # results = ranksums(signature_cases_ratio, signature_controls_ratio, axis=0)
    # for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
    #     print(f'{signame}: {t:.2f} (p={p:.4f})')
    # print('\n')
    #### Decomposed to COSMIC Signatures ####
    signature = pd.read_csv(
        f"{mutdir}/cosmic_signature_per_sample_{sigtype}.csv", index_col=0
    )
    signature_cases = signature.iloc[
        [x in label.iloc[label.values == 1, :].index for x in signature.index], :
    ]
    signature_controls = signature.iloc[
        [x in label.iloc[label.values == 0, :].index for x in signature.index], :
    ]
    #
    signature_cases_ratio = signature_cases.divide(
        signature_cases.sum(axis=1), axis="rows"
    )
    signature_controls_ratio = signature_controls.divide(
        signature_controls.sum(axis=1), axis="rows"
    )
    #
    # print('\n#### Results for raw values ####')
    results = ranksums(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        stat_results = pd.concat(
            (
                stat_results,
                pd.DataFrame(
                    {
                        "Signature Type": f"{sigtype} COSMIC",
                        "Signature name": signame,
                        "statistic": t,
                        "p-value": p,
                    },
                    index=[0],
                ),
            )
        )
    #
    # print('\n#### Results for ratios ####')
    # results = ranksums(signature_cases_ratio, signature_controls_ratio, axis=0)
    # for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
    #     print(f'{signame}: {t:.2f} (p={p:.4f})')
    # print('\n')

stat_results.iloc[stat_results["p-value"].values < 0.05, :].to_csv(
    f"{dout_stat}/mutational_signature_exposure_ranksum_significant.csv",
    header=True,
    index=False,
)
