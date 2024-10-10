# Paths used by many analyses

HMP2_root <- "~/analysis/hmp2_multiomics_paper"
HMP2_data <- file.path(HMP2_root, "data")

if (!dir.exists(HMP2_root)) {
	stop("Fill in HMP2_root in env_config.r to use the HMP2 scripts!")
}
