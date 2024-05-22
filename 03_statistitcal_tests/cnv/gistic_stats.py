import pandas as pd
from scipy.stats import ttest_ind

###################
#### Load data ####
###################

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv/10_cleaned_cnv"

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
case_ids = meta.iloc[meta.label.values == 1, :].patient_id.unique()

gistic_data = {}

for feature_name in ["cyto_thres_aber", "cyto_copy", "reg_copy", "reg_thres"]:
    gistic_data[feature_name] = pd.read_csv(
        f"{din}/{feature_name}_conf90.csv", index_col=0
    )

# Note: CNGPLD showed difference in "Deletion Peak 232 - CN values" between cases and controls

###########################################################
#### Find regions with difference between case/control ####
###########################################################

feature_name = "reg_copy"
df = gistic_data[feature_name]

print(f"\n\n#### {feature_name} ####")
for feature in df.columns:
    case_levels = df.iloc[[x in case_ids for x in df.index], :][feature].values
    control_levels = df.iloc[[x not in case_ids for x in df.index], :][feature].values
    t, p = ttest_ind(case_levels, control_levels)
    if p < 0.05:
        print(f"{feature}: t={t:.2f}, p={p:.4f}")
