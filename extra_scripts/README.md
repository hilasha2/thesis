# Extra scripts for general usage.

## changing_ecs_names.Rmd: 
changes numbered nomenclature of ECs into their defined names in ECs tables, e.g. 1.1.1.1 -> Alcohol dehydrogenase. 


## creating_metaphlan_species_file.ipynb:
Since MetaPhlAn tables contain taxa from all levels of the taxonomy tree, i.e. kingdom, phyla, etc., this scripts leaves only the species level. 
It also removes extra strings after sample names, like the '_taxonomic_profile' part in 'SAMPLE_01_taxonomic_profile'. 
(Notice that sometimes you need to skip a row and sometimes not when reading the metaphlan_taxonomic_profiles.tsv file). 


## changing_iirn_tables_column_names.Rmd:
This script removes the extra strings in the sample names columns in the relative abundance tables. 



