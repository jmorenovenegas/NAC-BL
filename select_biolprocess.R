# Extraer los genes de una función enriquecida en una linkcomm.
# Anotar cuales de entre esos genes estaba presente en el panel, y si se ha activado o inhibido
workingDir <- '.'
setwd(workingDir)
dir.create('Results/selected_functions')

options(stringsAsFactors = FALSE)
set.seed(1234)

# Load required libraries
source("RCode/0_loadLibraries.R")
loadpkg("dplyr") 
loadpkg("ggplot2")
loadpkg("clusterProfiler")

get_info = function(enrichfile, functionID, wgcnafile, linkcommfile, deafile){
  
  dataDir <- "Results/linkcomm"
  RDataDir <- "Data/RData"
  deaDataDir <- "Data"
  wgcnaDataDir <- "Results"
  
  # load the enrichment data
  d = read.csv(file.path(dataDir, enrichfile),nrows = 50)
  
  # load wgcna data
  w = read.csv(file.path(wgcnaDataDir, wgcnafile)) %>% dplyr::select(entrezID,moduleColor,GS.NAC, p.GS.NAC,MM.blue,p.MM.blue)
  names(w) = c("ENTREZID"  ,  "moduleColor","GS.NAC","p.GS.NAC","MM.blue", "p.MM.blue")
  
  # load the linkcomm
  linkcommid = strsplit(tail(strsplit(enrichfile, split = "_")[[1]], 1), split = "\\.")[[1]][1]
  load(file.path(dataDir,linkcommfile))
  
  # load the gene expression matrix
  # I need the DEA data to get info about over(under) regulated genes
  dea_data = read.csv(file.path(deaDataDir, deafile)) %>% dplyr::select(gene.name, log2FoldChange, padj)
  dea_geneid = bitr(dea_data$gene.name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  names(dea_geneid) = c("gene.name" , "ENTREZID")
  dea_data = merge(dea_geneid, dea_data)
  
  
  genes_linkcomm = sort(unique(c((dplyr::filter(lc$edges, cluster==linkcommid))$node1, (dplyr::filter(lc$edges, cluster==linkcommid))$node2)))
  genes = sort(strsplit((d[which(d$ID == functionID),])$geneID, split = "/")[[1]])
  descript = (d[which(d$ID == functionID),])$Description
  
  o = data.frame("ENTREZID" = genes)
  p= merge(o, dea_data, all.x=T ) %>% dplyr::select(ENTREZID,log2FoldChange,padj)
  q = bitr(p$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  r = merge(q,p)
  
  out=merge(r, w, all.x=T)
  
  write.csv(out, file = file.path('Results/selected_functions',paste(gsub(" ", "-", descript),gsub(":","",functionID),"linkcomm",linkcommid,"genedata.csv", sep = "_")))
  
}

get_info(enrichfile = "metabol_ct_enrich-linkcomm_1.csv",
         functionID = "R-HSA-400206",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")

get_info(enrichfile = "metabol_ct_enrich-linkcomm_1.csv",
         functionID = "GO:0001666",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")


get_info(enrichfile = "metabol_ct_enrich-linkcomm_1.csv",
         functionID = "GO:0006979",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")


get_info(enrichfile = "metabol_ct_enrich-linkcomm_3.csv",
         functionID = "R-HSA-1989781",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")

get_info(enrichfile = "metabol_ct_enrich-linkcomm_4.csv",
         functionID = "R-HSA-2990846",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")

get_info(enrichfile = "metabol_ct_enrich-linkcomm_4.csv",
         functionID = "GO:0019216",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")

get_info(enrichfile = "metabol_ct_enrich-linkcomm_4.csv",
         functionID = "hsa05202",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")

get_info(enrichfile = "metabol_ct_enrich-linkcomm_5.csv",
         functionID = "R-HSA-400206",
         wgcnafile = "blue-genesInfo.csv",
         linkcommfile = "metabol_ct_lc.RData",
         deafile= "BL-NAC-DEA.csv")
