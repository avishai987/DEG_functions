msigdb_genes = msigdbr(species = "Homo sapiens") %>%
  dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame() %>% pull(gene_symbol) %>% unique
fwrite(x = msigdb_genes %>% as.list,file = "./msigdb_genes.txt",sep = "\n")

human.msigdb.genes = fread(input = "../Data/human.msigdb.genes",header = F,select = c(3)) 
human.msigdb.genes.index = c()
for (row_num in 1:nrow(human.msigdb.genes)) {
  row = human.msigdb.genes[row_num,] %>% t() %>% as.vector()
  splitted = strsplit(x = row,split = ",") %>% unlist()
  human.msigdb.genes.index = c(human.msigdb.genes.index,splitted)
}
human.msigdb.genes.index = human.msigdb.genes.index %>% unique()


human2gene = fread(input = "../Data/human2gene.tsv",header = F,select = c(2,7),col.names = c("gene_index","gene_symbol"))
human2gene_dic = human2gene$gene_symbol
names(human2gene_dic)  = human2gene$gene_index

human.msigdb.genes_symbols = human2gene_dic[human.msigdb.genes.index]
human.msigdb.genes_symbols = human.msigdb.genes_symbols %>% unname
fwrite(x = human.msigdb.genes_symbols %>% as.list,file = "./msigdb_homer_genes.txt",sep = "\n")
