# **Aim 4: Develop and test new analytic methods for MGX: Investigating if the adjustment of feature abundance tables to taxa mitigates the over-representation of features originating from abundant species.**

Notice: Missing all input gene families files as they are too big to upload.

Figure 27: Heatmaps of DA of IIRN_plus functional profiles between control and CD patients, comparing results in relative abundance tables and count tables

	script: aim_4.ipynb
	metadata: IIRN_plus/metadata/IIRN_metadata_31_May_2023_samples_above_threshold_only.tsv
	gene families:
		genefamilies_relab_unstratified_shorter_headers.tsv
		genefamilies_count_species_normalized.tsv
	ko:
		ko_relab_unstratified_shorter_headers.tsv
		ko_count_species_normalized.tsv
	extra data: extra_data/kotable.txt
	chickens: ./input/data/chickens/genefamilies_relab_unstratified.tsv

