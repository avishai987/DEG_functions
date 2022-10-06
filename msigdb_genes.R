msigdb_genes = msigdbr(species = "Homo sapiens") %>%
  dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame() %>% pull(gene_symbol) %>% unique
fwrite(x = msigdb_genes %>% as.list,file = "./msigdb_genes.txt",sep = "\n")
