# Load Bugs and Virus tables
#   bugs.pcl: Full MetaPhlAn2 table (all taxonomic levels), in relative abundances
#   viruses.pcl: VirMAP viral abundance table, in counts


if (!exists("bugs.pcl")) {
    source("./env_config.r")
    source("./common/pcl_utils.r")
    source("./common/fix_metadata.r")

    library(dplyr)

    # Read in raw table
    # it was previously: metadata.rows = " ")
    # bugs.pcl <- pcl.read(file.path(HMP2_data, "mgx", "taxonomic_profiles.pcl.tsv"), metadata.rows=" ")
    bugs.pcl <- pcl.read.mine(data.path = file.path(HMP2_data, "mgx", "taxonomic_profiles.tsv"), 
                              metadata.path = file.path(HMP2_data, "metadata", "hmp2_metadata_2018-08-20.csv"),
                              file.type = "taxonomy")
    
    bugs.pcl$meta <- bugs.pcl$meta %>% filter(data_type == 'metagenomics')
    
    # 11_11_21 - I have noticed that the sample IDs are not ordered the same in the metadata and data
    # however it is crucial for future analyses.  
    bugs.pcl$x <- bugs.pcl$x[order(row.names(bugs.pcl$x)),]
    sample_idx <- match(rownames(bugs.pcl$x), bugs.pcl$meta$External.ID)
    bugs.pcl$meta <- bugs.pcl$meta[sample_idx, ]

    
    ## Filter out low-read depth samples and fix metadata
    
    # stop if there is a sample which has no reads. 
    stopifnot(all(!is.na(bugs.pcl$meta$reads_filtered)))
    
    # keeping only samples which had above 1M reads. 
    bugs.pcl <- bugs.pcl %>% fix_metadata %>%
        pcl.filter.s(reads_filtered >= 8e6)

    # Remove bizarre samples
    bizarre <- names(which(!pcl.apply.s(bugs.pcl, any((x>0) & (x<100)))))
    if (length(bizarre) > 0) {
        warning(sprintf("These samples have only one bug: %s", do.call(paste, c(as.list(bizarre), list(sep=" ")))))
        bugs.pcl <- pcl.filter.s(bugs.pcl, keep=!(rownames(bugs.pcl$x) %in% bizarre))
    }

    #27.21.2021 - i have noticed that the type of bugs.pcl$meta$diagnosis has suddenly converted to 
    # factors so I converted it back to character. Though it would help but nope. 
    #bugs.pcl$meta[] <- lapply(bugs.pcl$meta, as.character)

    # Read in the viromics table
    #viruses.pcl <- pcl.read(file.path(HMP2_data, "mvx", "HMP2.Virome.VirMAP.pcl.tsv"), metadata.rows=" ")

    # Filter samples with low reads and fix metadata
    #viruses.pcl <- viruses.pcl %>%
    #    pcl.filter.s(keep=viruses.pcl$meta$reads_viral >= 10) %>%
    #    # Makes the virus names look more like metaphlan
    #    pcl.map.fnames(chartr(";, ", "|__", gsub("(s__.*?);[^;_]*$", "\\1",
    #        gsub("(super)?(\\w)\\w*=([^;]+)", "\\2__\\3", gsub(";taxId=\\d+", "", Name))))) %>%
    #    fix_metadata

    # Cleanup
    rm(bizarre)
    
    ####################### this section deals with mgx data from metaphlan2
    bugs.pcl2 <- pcl.read.mine(data.path = file.path(HMP2_data, "mgx", "hmp2_mgx_metaphlan2_taxonomic_profiles.tsv"), 
                              metadata.path = file.path(HMP2_data, "metadata", "hmp2_metadata_2018-08-20.csv"),
                              file.type = "taxonomy")
    #bugs.pcl2$x <- bugs.pcl2$x[row.names(bugs.pcl2$x) != 'HSM5MD8B_P',] # this makes pcl.filter.s not work
    #bugs.pcl2$x <- bugs.pcl2$x[row.names(bugs.pcl2$x) != 'PSM6XBRK_P',]
    #bugs.pcl2$x <- bugs.pcl2$x %>% subset(select = -c('HSM5MD8B_P','PSM6XBRK_P'))
    bugs.pcl2$x <- bugs.pcl2$x[order(row.names(bugs.pcl2$x)),]
    bugs.pcl2$meta <- bugs.pcl2$meta %>% filter(data_type == 'metagenomics')
    
    sample_idx <- match(rownames(bugs.pcl2$x), bugs.pcl2$meta$External.ID)
    bugs.pcl2$meta <- bugs.pcl2$meta[sample_idx, ]
    #bugs.pcl2$x <- bugs.pcl2$x[sample_idx, ]
    
    stopifnot(all(!is.na(bugs.pcl2$meta$reads_filtered)))
    bugs.pcl2 <- bugs.pcl2 %>% fix_metadata %>% pcl.filter.s(reads_filtered >= 8e6)

    #bizarre <- names(which(!pcl.apply.s(bugs.pcl2, any((x>0) & (x<100)))))
    #if (length(bizarre) > 0) {
    #  warning(sprintf("These samples have only one bug: %s", do.call(paste, c(as.list(bizarre), list(sep=" ")))))
    ##}
    #rm(bizzare)
    ########################
}


