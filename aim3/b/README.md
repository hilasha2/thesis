**Aim 3.b. – Analyzing and capturing disease signals in Sheba’s larger cohort.**

Figure 21: Initial process of IIRN_plus dataset. 

A) Sequencing depth curve for sample filtration
process – number of species as a function of number of reads per sample.

	script: figure_21/num_features_by_reads/features_as_function_of_reads_for_iirn.Rmd
	metadata: IIRN_metadata_23_Feb_2023_added_pilot_and_biobakery_results.tsv
B) Feature filtration process. Relative abundance heatmaps.

	script: figure_21/filtering/iirn_analysis.ipynb
	metadata: IIRN_metadata_23_Feb_2023_added_pilot_and_biobakery_results.tsv
	data: ecs_relab_unstratified_shorter_headers.tsv
	pathabundance_relab_unstratified_shorter_headers.tsv
	metaphlan_species_profiles.tsv
C) PCoA plot for testing batch effects of combining two studies: IIRN and SOURCE-pilot

	script: figure_21/pcoa/iirn_pcoa.Rmd
	metadata: IIRN_plus/metadata/IIRN_metadata_31_May_2023_samples_above_threshold_only.tsv
	data: IIRN_plus/mgx/metaphlan/metaphlan_species_filtered_1abund_in_1perc_outof10k.tsv

Figure 22:

A) Violin plots of pairwise Spearman correlations between the feature-microbial profiles either from the same subject (within) or other subjects (across).

	script: figure_22/spearman_correlations/iirn_spearman_inter_vs_intra.Rmd
	metadata: IIRN_plus/metadata/IIRN_metadata_31_May_2023_samples_above_threshold_only.tsv
	data: IIRN_plus/mgx/metaphlan/metaphlan_species_filtered_1abund_in_1perc_outof10k.tsv
	IIRN_plus/mgx/humann/pathbundance_filtered_1abund_in_1perc_outof10k.tsv
	IIRN_plus/mgx/humann/ecs_filtered_1abund_in_1perc_outof10k_ecs2name.tsv
	
B) PERMANOVA analysis for testing variation in Canberra distance matrices of microbial relative abundance tables between different variables in metadata.

	script: figure_22/permanova/iirn_permanova.Rmd
	metadata: IIRN_plus/metadata/IIRN_metadata_31_May_2023_samples_above_threshold_only.tsv
	data: IIRN_plus/mgx/metaphlan/metaphlan_species_filtered_1abund_in_1perc_outof10k.tsv
	IIRN_plus/mgx/humann/pathbundance_filtered_1abund_in_1perc_outof10k.tsv
	IIRN_plus/mgx/humann/ecs_filtered_1abund_in_1perc_outof10k_ecs2name.tsv
	
Figure 23: Box plots of top significant associations between IIRN+’s microbial feature abundance tables and their metadata (FDR < 0.25) as identified by MaAsLin2.

	metadata: IIRN_metadata_17_May_2023.tsv
	data:
	metaphlan_species_flare_and_never_only_10kTSS.tsv
	pathabundance_flare_and_never_only_10kTSS.tsv
	ecs_flare_and_never_only_10kTSS_ecs2name.tsv
	
	IIRN_plus/mgx/metaphlan/metaphlan_species_filtered_1abund_in_1perc_outof10k.tsv
	IIRN_plus/mgx/humann/pathabundance_filtered_1abund_in_1perc_outof10k.tsv 
	IIRN_plus/mgx/humann/ecs_filtered_1abund_in_1perc_outof10k_ecs2name.tsv
	
	command examples: 
	~/maaslin2/code/R/Maaslin2.R --normalization=NONE --random_effects="patient_No" --fixed_effects="Group2,Age,Gender" --standardize=FALSE ~/biobakery_workflows/output_data/01_Jan_2023_mgx_iirn/metaphlan/merged_with_pilot/metaphlan_species_flare_and_never_only_10kTSS.tsv ~/analysis/datasets/IIRN_metadata_17_May_2023.tsv ~/maaslin2/output/17_May_2023_iirn_source_pilot_species_flared_vs_never_flared

	~/maaslin2/code/R/Maaslin2.R --normalization=NONE --random_effects="patient_No" --fixed_effects="diagnosis,Age,Gender" --standardize=FALSE ~/biobakery_workflows/output_data/01_Jan_2023_mgx_iirn/humann/merged_with_pilot/ecs_filtered_1abund_in_1perc_outof10k_ecs2name.tsv ~/analysis/datasets/IIRN_metadata_17_May_2023.tsv ~/maaslin2/output/03_Apr_2023_iirn_source_pilot_ecs_ecs2name


Figure 24: ROC curve for a random forest trees on a ECs relative abundance profile of IIRN+’s quiescent patients. 

	script: figure_24/iirn_analysis.ipynb
	metadata: IIRN_metadata_17_May_2023.tsv
	data: IIRN_plus/mgx/humann/ecs_flare_and_never_one_sample_per_patient_10kTSS_ecs2name.tsv




