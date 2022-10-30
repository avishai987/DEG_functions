library(dplyr)
library(data.table)
library(msigdbr)
library(stringr)
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd() #set wd to script path


#write msigdb genes
msigdb_genes = msigdbr(species = "Homo sapiens") %>%
  distinct(gs_name, gene_symbol) %>% as.data.frame() %>% pull(gene_symbol) %>% unique
fwrite(x = msigdb_genes %>% as.list,file = "./msigdb_genes.txt",sep = "\n")


#write msigdb genes from homer
human.msigdb.genes = fread(input = "./human.msigdb.genes",header = F,select = c(3)) 
human.msigdb.genes.index = c()
for (row_num in 1:nrow(human.msigdb.genes)) {
  row = human.msigdb.genes[row_num,] %>% t() %>% as.vector()
  splitted = strsplit(x = row,split = ",") %>% unlist()
  human.msigdb.genes.index = c(human.msigdb.genes.index,splitted)
}
human.msigdb.genes.index = human.msigdb.genes.index %>% unique()


human2gene = fread(input = "./human2gene.tsv",header = F,select = c(2,7),col.names = c("gene_index","gene_symbol"))
human2gene_dic = human2gene$gene_symbol
names(human2gene_dic)  = human2gene$gene_index

human.msigdb.genes_symbols = human2gene_dic[human.msigdb.genes.index]
human.msigdb.genes_symbols = human.msigdb.genes_symbols %>% unname
fwrite(x = human.msigdb.genes_symbols %>% as.list,file = "./msigdb_homer_genes.txt",sep = "\n")


#write hallmark db from homer
human.msigdb.pathways = fread(input = "./human.msigdb.genes",header = F,select = c(2,3),col.names = c("pathway_name","gene_symbol")) 
human.msigdb.pathways =  human.msigdb.pathways %>% filter(str_detect(pathway_name,"HALLMARK"))
hallmark_pathways_vector = human.msigdb.pathways$pathway_name %>% unique()

homer_hallmark = data.frame()
for (name in hallmark_pathways_vector) {
  pathway_genes = human.msigdb.pathways %>% filter(pathway_name == name) %>% select(gene_symbol) %>% t() %>% as.vector()
  splitted = strsplit(x = pathway_genes,split = ",") %>% unlist()
  pathway_symbols = human2gene_dic[splitted]
  pathway_df = data.frame(gs_name = name, gene_symbol = pathway_symbols)
  homer_hallmark = rbind(homer_hallmark, pathway_df)
}
fwrite(homer_hallmark,file = "homer_hallmark.csv")

a = human.msigdb.genes_symbols %>% as.data.frame()
inter = intersect(rownames(de_genes), human.msigdb.genes_symbols)
b = rownames(de_genes)[!rownames(de_genes) %in% inter]
d = c(inter,b[47])
fwrite(x = d %>% as.list,file = "./background2.txt",sep = "\n")

c = c(b,fdr_genes)
fwrite(x = c %>% as.list,file = "./intersected_bg.txt",sep = "\n")
c2 = c[1:1000]
fwrite(x = c2 %>% as.list,file = "./intersected_bg2.txt",sep = "\n")

undebug(enricher)
result = enrichment_analysis(de_genes, background = c,fdr_Cutoff = 0.05,
                             ident1_name = "pre-treatment",ident2_name = "on-treatment" ,return_df = T, add_bg_to_db = F, add_msigdb_to_bg = T,db = "homer_hallmark")
result = enrichment_analysis(de_genes, background = rownames(de_genes),fdr_Cutoff = 0.05,
                             ident1_name = "pre-treatment",ident2_name = "on-treatment" ,return_df = T, add_bg_to_db = F, add_msigdb_to_bg = T,db = "homer_hallmark")
library(gprofiler2)
background = rownames(de_genes)
gcon = gconvert(query = background, organism = "hsapiens",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)
names = gcon$name

test = c("ACS.1","CDA") %>% as.character()
background = gsub(pattern = "\\..*",replacement = "",x = background)
gcon = gconvert(query = background, organism = "hsapiens",
                target="ENSG", mthreshold = 1, filter_na = TRUE)
names = gcon$name%>% unique()

result = enrichment_analysis(de_genes, background = names,fdr_Cutoff = 0.05,
                             ident1_name = "pre-treatment",ident2_name = "on-treatment" ,return_df = T, add_bg_to_db = F, add_msigdb_to_bg = T,db = "homer_hallmark")

ac2gene = fread(input = "./human2gene.tsv",header = F,select = c(1,7),col.names = c("gene_index","gene_symbol"))
ac2gene_dic = ac2gene$gene_symbol
names(ac2gene_dic)  = ac2gene$gene_index

ac.genes_symbols = ac2gene_dic[background]
ac.genes_symbols = ac.genes_symbols %>% unname
result = enrichment_analysis(de_genes, background = ac.genes_symbols,fdr_Cutoff = 0.05,
                             ident1_name = "pre-treatment",ident2_name = "on-treatment" ,return_df = T, add_bg_to_db = F, add_msigdb_to_bg = T,db = "homer_hallmark")
ac2gene_dic[1:10]
ac2gene_dic["AAQ88605"]
save(ac2gene_dic, file="ac2gene_dic.rds")
a = load("./ac2gene_dic.rds") %>% get()

debug(enrichment_analysis)
result = enrichment_analysis(de_genes, background = background,fdr_Cutoff = 0.05,
                             ident1_name = "pre-treatment",ident2_name = "on-treatment" ,return_df = T, add_bg_to_db = F, add_msigdb_to_bg = T,db = "homer_hallmark"
                             ,convert_background = T)
