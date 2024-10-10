Aim 3.c. - Identifying correlations between paired MGX and metabolomic (MBX) samples.

figure 25: DA Heatmap of IIRN_plus’s aggregated KO relative abundance table between control and CD patients
	script: aim_4.ipynb (first section, & 'iirn + source pilot', & 'aggregate' sections)
	metadata: ~/IIRN_plus/metadata/IIRN_metadata_31_May_2023_samples_above_threshold_only.tsv
	data: ko_relab_unstratified_shorter_headers.tsv
	extra data: ~/extra_data/kotable.txt
	output data: aim4_IIRN_ko_relab_unstratified_first_samples_1e07of10k_in10perc_diff_CD_vs_ctl_in_lvl3.csv
	
Figure 26: HAllA heatmap of correlations of IIRN_plus’s MBX table vs. IIRN_plus’s MGX table
 MGX table’s features included sum-aggregated KO features which were found significant in Fig. 25.

For matching mbx and mgx samples:
	script: aim4_matching_iirn_mgx_mbx.Rmd
	ID_MATCH: MultiOmics_IDs.csv
	MBX: Stool_Data_MZM10_N2_filtered_v2.csv
	METADATA: Stool_metadata.tsv
	KO_RELAB: ko_relab_unstratified_samples_above_thres_filt_1e08of10k_in10perc.tsv
	KO_RELAB_META: ko_relab_unstratified_samples_above_thres_filt_1e08of10k_in10perc.csv_feature.tsv
	SIG.RELAB.AGG: aim4_IIRN_ko_relab_unstratified_first_samples_1e07of10k_in10perc_diff_CD_vs_ctl_in_lvl3.csv
	output mgx: IIRN_mgx_relab_agg_sig_23_Nov_2023.tsv
	output mbx: IIRN_mbx_22_Nov_2023.tsv
	
HAllA correlations calculations (version 0.8.20):
	data mgx: IIRN_mgx_relab_agg_sig_23_Nov_2023.tsv
	data mbx: IIRN_mbx_22_Nov_2023.tsv
	terminal command:
halla -x IIRN_mbx_22_Nov_2023.tsv -y IIRN_mgx_relab_agg_sig_23_Nov_2023.tsv -o halla0.8.20_IIRN_mbx_mgx_relab_agg_sig_alpha010 --fdr_alpha 0.1 -m spearman --hallagram
