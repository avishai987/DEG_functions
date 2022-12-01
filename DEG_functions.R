data_dir =  file.path(.libPaths()[1], "DEG_functions","DEG_functions-main")

# Load Functions ---------------------------------------------------------------
require(msigdbr,quietly = T)
require(fdrtool,quietly = T)
require(enrichR,quietly = T)
require(clusterProfiler,quietly = T)
require(RCurl,quietly = T)

# Split Differential_expression to up regulated and down regulated:
de_split <- function(Differential_expression_genes) {
  up_regulated = Differential_expression_genes[Differential_expression_genes$avg_log2FC >= 0,] 
  down_regulated = Differential_expression_genes[Differential_expression_genes$avg_log2FC < 0,]
  all_regulated_vec <- list(up_regulated = up_regulated, down_regulated = down_regulated)
  return (all_regulated_vec)
}


# enrichment analysis function. Argument: output from de_split, group name.
enrichment_analysis <- function(Differential_expression_genes = NULL, all_regulated = NULL,
                                group_name="",background, fdr_Cutoff = 0.01,ident1_name = "", 
                                ident2_name = "",write_csv = F, ident_name = "", ident_num = 1,
                                return_df = F, pval_cutoff = F, add_bg_to_db = F , db = NULL,
                                add_msigdb_to_bg = F, convert_background = F) {
library(msigdbr,quietly = T)
library(fdrtool,quietly = T)
library(enrichR,quietly = T)
library(clusterProfiler,quietly = T)
library(RCurl,quietly = T)
  
  if ( is.null(all_regulated)){
    if ( is.null(Differential_expression_genes)){
      print("please provide Differential_expression_genes or all_regulated")
      return (NULL)
    }
    all_regulated = de_split(Differential_expression_genes)
  }
  
  
  i = 1; #for choosing the right the names
  all_graphs <- list() #create list to save graphs
  all_results = list()#create list to enrichment results
  all_output_df = list()#create list to enrichment results
  
  for (regulated in all_regulated) { #do twice for downregulated and up

    #take genes from differential expression
    regulated <- tibble::rownames_to_column(regulated, "genes")
    genes = regulated[,c("genes","avg_log2FC")] #get genes and FC
    #calculate fdr
    regulated$fdr<-p.adjust(p = as.vector(regulated$p_val) ,method = "fdr" )
    
    if (pval_cutoff == T){
      fdr_genes = regulated$genes[regulated$p_val<0.05] #take genes less than the cutoff 
      
    }
    if (pval_cutoff == F){
      fdr_genes = regulated$genes[regulated$fdr<fdr_Cutoff] #take genes less than the cutoff 
    }
    
    #set hallmark pathways
    if (is.null(db)){
      hallmark_gene_set = msigdbr(species = "Homo sapiens", category = "H")
      msigdbr_t2g = hallmark_gene_set %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
      
    }else if(db == "homer_hallmark"){
      msigdbr_t2g <- fread(file.path(data_dir,"homer_hallmark.csv"),sep = ",")
      
      
    } else{ msigdbr_t2g = db}

    #Add all genes to msigdb (beacause the background is (TERM2GENE AND universe) (not like the documentation))
    if (!is.null(background) & add_bg_to_db == T){
      all_genes = data.frame(gs_name = "background",gene_symbol = background) 
      msigdbr_t2g = rbind(all_genes, msigdbr_t2g)
    }
    
    if ( add_msigdb_to_bg == T){
      msigdb_genes <- scan(file.path(data_dir,"msigdb_homer_genes.txt"), character(), quote = "",quiet = T)
      all_genes = data.frame(gs_name = "msigdb",gene_symbol = msigdb_genes) 
      msigdbr_t2g = rbind(all_genes, msigdbr_t2g)
    }
    if ( convert_background == T){
      ac2gene_dic = load("./ac2gene_dic.RData") %>% get()
      background = ac2gene_dic[background] %>% unname
    }
    #perform enrichment analysis
    enrichment_result = enricher(
      gene = fdr_genes,
      pvalueCutoff = 0.0001,
      pAdjustMethod = "fdr",
      universe = background,
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.0001,
      TERM2GENE = msigdbr_t2g,
      TERM2NAME = NA
    )
    #take relevant result
    if (!is.null(enrichment_result)) { 
      enrichment_result = enrichment_result@result #take results
      enrichment_result = enrichment_result[ , -which(names(enrichment_result) %in% c("ID","Description","geneID"))] #remove unwanted cols
      output_df  = enrichment_result
      
        #Add pathways that not in the output_df
        all_pathways = unique(msigdbr_t2g$gs_name)
        computed_pathways = rownames(output_df)
        not_computed = all_pathways [!all_pathways %in% computed_pathways]
        zeros = replicate(0, n =length(not_computed))
        ones = replicate(1, n =length(not_computed))
        
        empty_pathways = data.frame(GeneRatio = zeros , BgRatio = zeros ,
                                    pvalue =ones , p.adjust = ones, qvalue = zeros, Count = zeros, row.names = not_computed)
        output_df = rbind(output_df, empty_pathways)
        
      
        all_output_df[[i]] = output_df
        names (all_output_df)[i] = names(all_regulated)[i]
      enrichment_result = enrichment_result[startsWith(x = rownames(enrichment_result),prefix = "HALLMARK"),]
      enrichment_result <- tibble::rownames_to_column(enrichment_result, "pathway_name") #make row names as column
      enrichment_result = enrichment_result[order(enrichment_result$p.adjust),] #order by pvalue
      enrichment_result = enrichment_result[1:10,] #take best 10 pathways
      enrichment_result[, 5] <- -log10(enrichment_result['p.adjust']) #perform -log10 on p.adjust values
    }
    
    else{ #if no genes were significant, plot empty enrichment result
      # pathway_names = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_MYC_TARGETS_V2", 
      #                   "HALLMARK_UV_RESPONSE_DN", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_E2F_TARGETS", 
      #                   "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_G2M_CHECKPOINT", 
      #                   "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_MYC_TARGETS_V1", 
      #                   "HALLMARK_P53_PATHWAY", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
      msigdbr_t2g = hallmark_gene_set %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
      
      pathway_names =  unique(msigdbr_t2g$gs_name)
  
      enrichment_result = data.frame(pathway_name = pathway_names[1:10], p.adjust = 0)
      zeros = replicate(0, n = 50)
      ones = replicate(1, n = 50)
      
      output_df = data.frame(GeneRatio = zeros , BgRatio = zeros ,
                                  pvalue =ones , p.adjust = ones, qvalue = zeros, Count = zeros, row.names = pathway_names)

      all_output_df[[i]] = output_df
      names (all_output_df)[i] = names(all_regulated)[i]
    }
    
    if (write_csv == T){
      file_name = paste0(names(all_regulated)[i],"_",group_name,".csv")
      write.csv(x=enrichment_result, file=file_name)
    }


    ident = ident2_name
    
    #make graph title
    if (names(all_regulated)[i]=="up_regulated"){
      markers_title="negative markers"
      bar_color = "indianred2"
    }
    else if(names(all_regulated)[i]=="down_regulated"){
      markers_title="positive markers"
      bar_color = "dodgerblue"
    }
      
    else {
      print ("error")
      return (0)
    }
    

    
    #create graph
    p = ggplot(data = enrichment_result, aes_string(x = reorder(enrichment_result$pathway_name,
                                                                enrichment_result$p.adjust), y = "p.adjust")) +
      geom_bar(stat = "identity", fill = bar_color) + 
      geom_hline(yintercept = 1.3, colour="black", linetype = "longdash",size = 1) +
      coord_flip() + xlab("Pathway") +
      scale_fill_manual( drop = FALSE) +
      ylab("-log10(p.adjust)") +
      geom_text(aes_string(label = "pathway_name", y = 0),
                size = 4,
                color = "black",
                position = position_dodge(1),
                hjust = 0)+
      theme(axis.title.y= element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())+
      ggtitle(paste(ident, markers_title)) 
    
    all_graphs =  append(all_graphs, list(p))
    all_results[[names(all_regulated)[i]]] = enrichment_result
    
    i = i+1
    
  }
  p3<-all_graphs[[1]]+all_graphs[[2]]
  print(p3)
  

  if (return_df){
    return (all_output_df)
  }
  else{
    return(p3)
    
  }
}


#Volcano plot
#' @param de_genes - value of FindMarkers
#' @param top_genes_text - show names of top n genes  (0 = all significant genes)
#' @param show_gene_names - show specific genes (list of genes)
#' @param fdr_cutoff - consider significant from cutoff
#' @param log2fc_cutoff - consider significant from cutoff
#' @param max_names - maximum names to show
#' @param return_de_genes - return DEG df instead of ggplot
#' @param clac_fdr_from_logfc - deprecated. calculate FDR from genes that has at list x logFC (0= all genes). FDR on less genes should raise number of significant genes.
#' @importFrom stats p.adjust

  volcano_plot<- function(de_genes, top_genes_text=0, title = "" ,show_gene_names = NULL, ident1 = "",
                      ident2 = "" , fdr_cutoff = 0.05 , log2fc_cutoff = 0.6, max_names = 10,
                      return_de_genes = F, show_graph = T, show_legend = T) {
  library(ggrepel,quietly = T)
  library(dplyr,quietly = T)
  names_for_label = c(paste(ident2,"down genes"),paste(ident2,"up genes"))
  #color genes if there are over/under/same expressed
  de_genes$diffexpressed <- "Same" 
  de_genes$avg_log2FC[!is.finite(de_genes$avg_log2FC)] <- 1000 #remove infinite values
  
  # #calculate fdr from genes that has logFC > clac_fdr_from_logfc (clac_fdr_from_logfc = 0 -> all genes)
  # filtered_markers = filter(de_genes, abs(avg_log2FC) > clac_fdr_from_logfc) #take genes with logFC > clac_fdr_from_logfc
  # filtered_markers$fdr <-p.adjust(p = filtered_markers$p_val ,method = "fdr") #calc fdr
  # de_genes$fdr = 1 
  # filtered_markers_indexes = which(rownames(de_genes) %in% rownames(filtered_markers))
  # de_genes[filtered_markers_indexes, "fdr"] = filtered_markers$fdr 
  
  de_genes$fdr = p.adjust(p = de_genes$p_val ,method = "fdr")
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de_genes$diffexpressed[de_genes$avg_log2FC > log2fc_cutoff & de_genes$fdr < fdr_cutoff] <- names_for_label[1]
  
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de_genes$diffexpressed[de_genes$avg_log2FC < -log2fc_cutoff & de_genes$fdr < fdr_cutoff] <- names_for_label[2]
  
  #set name labels for highest genes
  de_genes$delabel <- NA
  de_genes$delabel[0:top_genes_text] <- rownames(de_genes)[0:top_genes_text]
  
  #set name labels for significant genes
  down_genes = de_genes$diffexpressed == names_for_label[1]
  up_genes = de_genes$diffexpressed == names_for_label[2]
  
  down_genes[max_names:length(down_genes)] = F
  up_genes[max_names:length(up_genes)] = F
  
  de_genes$delabel[up_genes] <- rownames(de_genes)[up_genes]
  de_genes$delabel[down_genes] <- rownames(de_genes)[down_genes]
  
  
  #Show genes that specify in show_gene_names
  if (!is.null(show_gene_names)){
    de_genes_index = match(show_gene_names,rownames(de_genes)) #indexes of de_genes that in show_gene_names
    de_genes_index <- de_genes_index[!is.na(de_genes_index)] #remove NA
  
    show_gene_names_index = show_gene_names %in% rownames(de_genes) #indexes of show_gene_names that in de_genes
    de_genes$delabel[de_genes_index] = show_gene_names [show_gene_names_index]
  }
  
  #colors for diff exp genes
  cols <- structure(c("green4", "red", "grey"), .Names = c(names_for_label[2],names_for_label[1], "Same"))
  
  title = paste(title,"Volcano plot- ", ident1,"vs", ident2)
  p = ggplot(data=de_genes, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    scale_color_manual(values=cols) +
    geom_text_repel(na.rm = T,box.padding = 1,max.overlaps = Inf,color = "blue") +
    xlab(paste("avg_log2FC (Positive = up in", ident1,")"))+ 
    ylab("Significance")+ 
    scale_y_continuous(labels = function(x) {parse(text = paste0("10^-",x))})+
    guides(col=guide_legend(title=paste0("Significant DEG\n(FDR<",fdr_cutoff," ,abs(log2fc) >", log2fc_cutoff,")")))+
    {if(show_legend == F) theme(legend.position="none")}+
    {if(show_legend) ggtitle(title) }+
    theme(axis.text  = element_text( color="black", size=14),axis.title = element_text( color="black", size=17))
  
  if(show_graph == T){print(p)}
  
  if (return_de_genes == T){ return(de_genes)}
  else {return (p)}
  
  }
  
genes_vec_enrichment<- function (genes, background, gene_sets = "", homer = F, title, silent = F, 
                                custom_pathways = NULL, return_all = F) {
  library(clusterProfiler,quietly = T)
  
  if (homer == T){
    gene_sets <- fread(file.path(data_dir,"homer_hallmark.csv"), sep = ",") #read hallmark from homer as gene sets
    
    #convert genes and background to with homer's name:
    background = gsub(pattern = "\\..*$", replacement = "",x = background) #remove "." from gene names
    genes = gsub(pattern = "\\..*$", replacement = "", x = genes) #remove "." from gene names
    
    ac2gene_dic = fread(file.path(data_dir,"ac2gene_dic.txt"),sep = "\t", header = F) #read dictionary supplied by homer
      #make dic:
    values = ac2gene_dic %>% pull(2)
    names = ac2gene_dic %>% pull(1)
    ac2gene_dic = values
    names(ac2gene_dic) = names
      #convert:
    background = ac2gene_dic[background] %>% unname
    genes = ac2gene_dic[genes] %>% unname
    
    #add msigdb genes to background (otherwise, genes that are not in hallmark will be ommited)
    msigdb_genes <- scan(file.path(data_dir,"msigdb_homer_genes.txt"),character(), quote = "", quiet = T)
    all_genes = data.frame(gs_name = "msigdb", gene_symbol = msigdb_genes)
    gene_sets = rbind(all_genes, gene_sets)
  }


  if(custom_pathways %>% is.null()== F){
    gene_sets = rbind(gene_sets, custom_pathways)
    
  }
  
  enrichment_result = enricher(gene = genes, pvalueCutoff = 0.05, 
                               pAdjustMethod = "fdr", universe = background, minGSSize = 10, 
                               maxGSSize = 500, qvalueCutoff = 1, TERM2GENE = gene_sets, 
                               TERM2NAME = NA)
  tryCatch({
    enrichment_result = enrichment_result@result
  }, error=function(e){ return (list(plt = plot.new(), mat = data.frame())) 
    exit()
  }) #if enrichment_result is null
  
  enrichment_result = enrichment_result[, -which(names(enrichment_result) %in% 
                                                   c("ID", "Description"))]
  enrichment_result <- tibble::rownames_to_column(enrichment_result, 
                                                  "pathway_name")
  for (pathway in unique(gene_sets$gs_name)) {
    if (!pathway %in% enrichment_result$pathway_name) {
      enrichment_result = enrichment_result %>% add_row(pathway_name = pathway, 
                                                        GeneRatio = "0/0", BgRatio = "0/0", pvalue = 1, 
                                                        p.adjust = 1, qvalue = 1, Count = 0)
    }
  }
  final_result = enrichment_result[order(enrichment_result$p.adjust), 
  ]
  final_result = final_result[1:10, ]
  final_result[, 5] <- -log10(final_result["p.adjust"])
  bar_color = "dodgerblue"
  p <- ggplot(data = final_result, aes_string(x = reorder(final_result$pathway_name, 
                                                          final_result$p.adjust), y = "p.adjust")) + geom_bar(stat = "identity", 
                                                                                                              fill = bar_color) + geom_hline(yintercept = 1.3, colour = "black", 
                                                                                                                                             linetype = "longdash", size = 1) + coord_flip() + xlab("Pathway") + 
    scale_fill_manual(drop = FALSE) + ylab("-log10(p.adjust)") + 
    geom_text(aes_string(label = "pathway_name", y = 0), 
              size = 4, color = "black", position = position_dodge(1), 
              hjust = 0) + theme(axis.title.y = element_blank(), 
                                 axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    ggtitle(title)
  if (silent == F) {
    print(p)
  }
  if (return_all == T){ return (list(plt = p, mat = enrichment_result))}
  return(enrichment_result)
}
  

sig_heatmap <- function(all_patients_result, title,clustering_distance =  "euclidean", annotation = NULL, silent = F) {
  my_fun <- function(p) {                     
    asterisks_vec = p
    # p = c(0.001, 0.01, 0.05,3.314507e-15)
    asterisks_vec[p<=0.05 & p>0.01] = "*"
    asterisks_vec[p<=0.01 & p > 0.001] = "**"
    asterisks_vec[p<=0.001 & p >= 0] = "***"
    asterisks_vec[p>0.05] = ""
    paste(asterisks_vec)
  }
  asterisks = all_patients_result
  asterisks[] <- lapply(all_patients_result, my_fun)
  
  
  
  all_patients_result = -log(all_patients_result)
  all_patients_result[all_patients_result>35] = 35
  paletteFunc <- colorRampPalette(c("white","navy"));
  
  palette <- paletteFunc(100)
  # breaks = seq(0,max(all_patients_result), length.out =6)
  # breaks = round(breaks, digits = 0)
  # breaks_labels = as.character(breaks)
  # breaks_labels[length(breaks_labels)] = "FDR"
  
  
  p<- pheatmap(all_patients_result,
               cluster_rows = T,
               cluster_cols = T,
               show_rownames = TRUE, 
               color = palette,
               # breaks = seq(0,20,0.2), 
               number_color = "grey30",
               main = title,
               display_numbers = asterisks,
               fontsize_row = 8,
               clustering_distance_rows = clustering_distance,
               clustering_distance_cols = clustering_distance,
               annotation_col = annotation[["myannotation"]],
               annotation_colors = annotation[["ann_colors"]],
               border_color = "black",silent = silent
               # legend_breaks =breaks,
               # legend_labels = breaks_labels
  )
  if(silent == F){
    print(p) }
  return (p)
}



#filter_features before findMarkers
#' @param object - seurat object
#' @param ident.1 from seurat findMarkers
#' @param ident.2 from seurat findMarkers
#' @param min.pct from seurat findMarkers
#' @param var_genes_only -filter for only genes that in variable genes
filter_features <- function(object,ident.1,ident.2,min.pct = 0.1, var_genes_only = F) {
  #filter by min.pct>0.1
  fc.results <- FoldChange(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2
  )
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  
  if (var_genes_only == T) {     #intersect with variable genes
    object = FindVariableFeatures(object = object, nfeatures = 15000)
    var_genes = object@assays[["RNA"]]@var.features
    features = intersect(x = var_genes, features)
  }

  return (features)
}
