

#---read data---
combined_path <- '/Users/zhusinan/Downloads/S-MiXcan_code_folder/3pi/bcac2020_result/bcac2020_result_pi3_02.csv'
combined3 = read.csv(combined_path)

combined_path2 <- '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/bcac2020_result_pi2.csv'
combined2 <- read.csv(combined_path2)

#---add cytoband----
ensembl_ref = read.csv('/Users/zhusinan/Downloads/S-MiXcan_code_folder/code_RealData/RealData/GTEx_Data/ensembl38.txt')
setDT(combined2)
combined2[, gene_id_clean := sub("\\..*$", "", gene_id)]  # drop any ENSG version

# Your Ensembl reference table (the one with Gene.stable.ID and Karyotype.band)
setDT(ensembl_ref)

# Build ENSG -> cytoband map (deduplicate just in case)
ensembl_cyto <- unique(
  ensembl_ref[, .(
    ENSG = sub("\\..*$", "", Gene.stable.ID),
    CYTOBAND = Karyotype.band,
    gene_name = Gene.name
  )],
  by = "ENSG"
)

# Join cytoband by ENSG
combined2 <- ensembl_cyto[combined2, on = .(ENSG = gene_id_clean)]
combine2 = combined2[, ENSG := NULL]

# Put CYTOBAND right after gene_id
setcolorder(combined2, c("gene_name", "gene_id", "CYTOBAND",
                           setdiff(names(combine2), c("gene_name","gene_id","CYTOBAND"))))
write.csv(combined2, '/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan/Results/SMiXcanK_results/bcac2020_result_pi2_annotated.csv', row.names = FALSE)



