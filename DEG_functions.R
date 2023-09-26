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


# enrichment analysis function.
enrichment_analysis <- function(markers,
                               background, fdr_Cutoff = 0.01,ident.1, 
                               ident.2 = "",show_by_ident = 1, by_pval = F,custom_pathways = NULL) {
library(msigdbr,quietly = T)
library(fdrtool,quietly = T)
library(enrichR,quietly = T)
library(clusterProfiler,quietly = T)
library(RCurl,quietly = T)
  
  all_regulated = de_split(markers)
  
  
  
  i = 1; #iterator
  all_results = list()#create list to enrichment results
  
  idents = c(ident.1,ident.2)
  
  title_1 = paste("Upregulated in",idents[show_by_ident])
  title_2 = paste("Downregulated in",idents[show_by_ident])
  titles = c(title_1, title_2)
  colors = c("dodgerblue", "indianred2")
  
  if(show_by_ident == 2){
    titles = titles %>% rev()
    colors = colors %>% rev()
  }
  
  
  for (regulated in all_regulated) { #do twice for downregulated and up

    #take genes from differential expression
    regulated <- tibble::rownames_to_column(regulated, "genes")
    genes = regulated[,c("genes","avg_log2FC")] #get genes and FC
    #calculate fdr
    regulated$fdr<-p.adjust(p = as.vector(regulated$p_val) ,method = "fdr" )
    
    if (by_pval == T){
      genes_to_test = regulated$genes[regulated$p_val<0.05] #take genes less than the pval cutoff 
      
    }
    if (by_pval == F){
      genes_to_test = regulated$genes[regulated$fdr<fdr_Cutoff] #take genes less than the fdr cutoff 
    }
    

    
    enrich_res = genes_vec_enrichment(genes = genes_to_test,background = background,homer = T,title = titles[i],bar_color = colors[i],silent = T
                         ,return_all = T,custom_pathways = custom_pathways)
    all_results[[i]] = enrich_res
    i = i+1
  }
  
  p<-all_results[[1]]$plt+all_results[[2]]$plt
  return(p)
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
  
  last_gene_to_show_down = which(de_genes$diffexpressed == names_for_label[1])[max_names]
  last_gene_to_show_up = which(de_genes$diffexpressed == names_for_label[2])[max_names]
  
  down_genes[last_gene_to_show_down:length(down_genes)] = F
  up_genes[last_gene_to_show_up:length(up_genes)] = F
  
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
    theme(axis.text  = element_text( color="black", size=12),axis.title = element_text( color="black", size=12))
  
  if(show_graph == T){print(p)}
  
  if (return_de_genes == T){ return(de_genes)}
  else {return (p)}
  
}
  
  #perform enrichment analysis of vector of genes
genes_vec_enrichment<- function (genes, background, gene_sets = "", homer = F, title, silent = F, 
                                custom_pathways = NULL, return_all = F,bar_color = "dodgerblue" ) {
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
  
  enrichment_result = enrichment_result[, -which(names(enrichment_result) %in% c("ID", "Description"))]
  enrichment_result <- tibble::rownames_to_column(enrichment_result, "pathway_name")
  for (pathway in unique(gene_sets$gs_name)) {
    if (!pathway %in% enrichment_result$pathway_name) { #add missing pathways
      enrichment_result = enrichment_result %>% add_row(pathway_name = pathway, 
                                                        GeneRatio = "0/0", BgRatio = "0/0", pvalue = 1, 
                                                        p.adjust = 1, qvalue = 1, Count = 0)
    }
  }
  final_result = enrichment_result[order(enrichment_result$p.adjust), ] #order by fdr
  final_result = final_result[1:10, ] #take top 10
  final_result[, 5] <- -log10(final_result["p.adjust"])
  bar_color = bar_color
  p <- ggplot(data = final_result, aes_string(x = reorder(final_result$pathway_name, 
                                                          final_result$p.adjust), y = "p.adjust")) + geom_bar(stat = "identity", 
                                                                                                              fill = bar_color) + geom_hline(yintercept = 1.3, colour = "black", 
                                                                                                                                             linetype = "longdash", size = 1) + coord_flip() + xlab("Pathway") + 
    scale_fill_manual(drop = FALSE) + ylab("-log10(p.adjust)") + 
    geom_text(aes_string(label = "pathway_name", y = 0), 
              size = 4, color = "black", position = position_dodge(1), 
              hjust = 0) + 
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
    ggtitle(title)
  
  if (silent == F) {
    print(p)
  }
  if (return_all == T){ return (list(plt = p, mat = enrichment_result))}
  return(enrichment_result)
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



.handle.genesets <- function(genesets) {
  if (is(genesets, "list")) {
    gsets.obj <- gsets$new(genesets, quiet=TRUE)
  }
  else if (is(genesets, "gsets") | is(genesets, "rgsets")) {
    gsets.obj <- genesets
  } 
  else {
    stop("Genesets must be gsets/rgsets object or named list of genesets")
  }
  return(gsets.obj)
}

#from: https://montilab.github.io/hypeR-docs/articles/docs/fgsea.html
hypeR_fgsea <- function(signature, genesets, sample.size=101, min.size=1, max.size=Inf, up_only = T,...) {
  # Save original arguments
  args <- as.list(environment())
  
  # Save gsets object
  gsets.obj <- .handle.genesets(genesets)
  args$genesets <- gsets.obj
  
  # Run fgsea
  results <- fgsea::fgseaMultilevel(stats=signature, 
                                    pathways=gsets.obj$genesets, 
                                    sampleSize=sample.size, 
                                    minSize=min.size, 
                                    maxSize=max.size,
                                    ...)
  
  data <- results %>%
    data.frame() %>%
    plyr::rename(c("pathway"="label", "padj"="fdr", "log2err"="lte", "size"="overlap", "leadingEdge"="le")) %>%
    dplyr::rename_with(tolower) %>%
    mutate(pval=signif(pval, 2)) %>%
    mutate(fdr=signif(fdr, 2)) %>%
    mutate(le=sapply(le, function(x) paste(x, collapse=','))) %>%
    mutate(signature=length(signature)) %>%
    mutate(geneset=sapply(label, function(x) length(gsets.obj$genesets[[x]]))) %>%
    dplyr::select(c("label", "pval", "fdr", "lte", "es", "nes", "signature", "geneset", "overlap", "le"))
  
  data.up <- data %>%
    dplyr::filter(es > 0) %>%
    dplyr::arrange(pval, es)
  
  data.dn <- data %>%
    dplyr::filter(es < 0) %>%
    dplyr::arrange(pval, es)    
  
  # Reproducibility information
  info <- list(fgsea=paste("v", packageVersion("fgsea"), sep=""),
               signature=length(signature), 
               genesets=args$genesets$info())
  
  info <- c(info, args[c("sample.size", "min.size", "max.size")])
  info <- lapply(info, as.character)
  
  # Wrap dataframe in hyp object
  hyp.up <- hyp$new(data=data.up, args=args, info=info)
  hyp.dn <- hyp$new(data=data.dn, args=args, info=info)
  mhyp <- multihyp$new(data=list("up"=hyp.up, "dn"=hyp.dn))
  if(up_only){return(hyp.up)}else{return(mhyp)}
}
