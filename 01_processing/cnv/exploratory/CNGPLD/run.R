library(io)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)

library(devtools)
load_all()
# Configurations #############################################################

genome <- "hg38"
datadir <- "/projects/b1131/saya/bbcar/data"
dout <- "/projects/b1131/saya/bbcar/data/02b_cnv/CNGPLD"
fits.fn <- paste0(dout, "/cngpld_bbcar.rds")

# Run analysis ###############################################################

if (file.exists(fits.fn)) {
  warning("Warning: Using cached model")
  fits <- qread(fits.fn)
} else {
  # Load segments and case/control labels
  all.segs <- read.table(
    paste0(
      datadir,
      "/02b_cnv/07_called_cn_segs/tissue_only/gistic_input_all.tsv"
    ),
    sep = "\t", header = FALSE
  )
  labels <- read.table(
    paste0(
      datadir,
      "/clinical/bbcar_label_studyid_from_gatk_filenames.csv"
    ),
    sep = ",", header = TRUE, row.names = 1
  )

  # Update column names and change format of chromosome column from "chr1" to "1"
  colnames(all.segs) <- c(
    "sample", "chromosome", "start", "end", "nprobes", "logr"
  )
  all.segs$chromosome <- sub("^chr", "", all.segs$chromosome)

  # Separate cases and controls
  seg.cases <- all.segs[
    all.segs$sample %in% as.numeric(rownames(labels)[labels == 1]),
  ]
  seg.controls <- all.segs[
    all.segs$sample %in% as.numeric(rownames(labels)[labels == 0]),
  ]

  # Fit and save model
  options(mc.cores = 12)
  fits <- compare_segs(seg.cases, seg.controls, genome = genome)

  qwrite(fits, fits.fn)
}

# Examine results ############################################################

# problematic deletion profiles:
cat("Failed deletion profiles:\n")
print(which(unlist(lapply(fits$del, function(x) !is(x$model, "gpldiff")))))

cat("Failed amplification profiles:\n")
print(which(unlist(lapply(fits$amp, function(x) !is(x$model, "gpldiff")))))

# Remove problematic amplification and deletion profiles
fits.orig <- fits
idx <- unlist(lapply(fits$del, function(x) is(x$model, "gpldiff")))
fits$del <- fits$del[idx]
idx <- unlist(lapply(fits$amp, function(x) is(x$model, "gpldiff")))
fits$amp <- fits$amp[idx]
# Significant regions for cases
regions.cases <- summary(fits, genome = genome, lodds.cut = 0.1)
regions.cases <- filter(
  regions.cases,
  end - start + 5 > 1e6,
  abs(ldiff) > 0.10,
  fdr < 0.30,
  n_obs >= 5
)

# Plot significant regions for cases
for (region_name in regions.cases$chromosome) {
  qdraw(
    {
      with(
        fits$amp[[region_name]],
        plot(model, data,
          which = c("response", "latent", "odds"),
          xlab = "position (Mbp)"
        )
      )
    },
    width = 5,
    height = 10,
    file = paste0(dout, "/cngpld_cases_", region_name, ".pdf")
  )
}

# Significant regions for controls
regions.controls <- summary(fits, direction = -1, genome = genome, lodds.cut = 0.1)
regions.controls <- filter(
  regions.controls,
  end - start + 5 > 1e6,
  abs(ldiff) > 0.10,
  fdr < 0.30,
  n_obs >= 5
)

# Plot significant regions for controls
for (region_name in regions.controls$chromosome) {
  qdraw(
    {
      with(
        fits$amp[[region_name]],
        plot(model, data,
          which = c("response", "latent", "odds"),
          xlab = "position (Mbp)"
        )
      )
    },
    width = 5,
    height = 10,
    file = paste0(dout, "/cngpld_controls_", region_name, ".pdf")
  )
}

qwrite(regions.cases, paste0(dout, "/cngpld_sig-regions_cases.tsv"))
qwrite(regions.controls, paste0(dout, "/cngpld_sig-regions_controls.tsv"))
