## This script contains functions that needed to run analysis on multi-organ proteomics data (calculate protein mean, pathway median)
## Nhung Pham 18 May 2021


###################################################
## Calculate log10 protein mean in each organ, exclude proteins that are measured in less than 10% samples in each tissue
## Input: 
#   - data: a data frame with columns as uniprotID, Genename, patient IDs starting with tissue, for example lung_b11_129N, heart_b14_127N 
#   - tissue_df: a data frame which contains tissue label as in data file without patient IDs, for example lung, heart 
## Output: a list with two objects
#   - all_organ_mean : log10 protein mean with columns as Uniport, genename, organ1,organ2, ect
#   - organ_express_total: log10 protein mean in each organ in a format to plot (2 columns: log10 mean, organ)

protein_mean <- function(data, tissue_df){
  # create empty dataframe that contains only gene ids and names
  all_organ_mean <- data.frame(data[,1:2])
  # prepare dataframe to plot
  organ_express_total <- data.frame()
  
  for (i in 1: length(tissue_df$tissue)){
    organ <- tissue_df$tissue[i]
    organ_data <- data %>% dplyr::select(contains(organ))
    # calculate log10 (mean+1) for protein that are measured in more than 10% of samples
    # proteins that are measured in less than 10% are left NA in mean values
    # log10(mean +1) to shift the plot to non-negative data
    organ_data$mean_value <- apply(organ_data,1, function(x){ ifelse(sum(!is.na(x))/ncol(organ_data) > 0.1, log10(mean(x, na.rm = TRUE)+1), NA)})
    # Combine data with other organs into 1 data frame
    all_organ_mean <- cbind(all_organ_mean,organ_data[,ncol(organ_data)]) 
    names(all_organ_mean)[names(all_organ_mean) == "organ_data[, ncol(organ_data)]"] <- organ
    organ_data_plot <-  data.frame(organ_data[,ncol(organ_data)])
    organ_data_plot$organ <- organ
    names(organ_data_plot) <- c('log10mean','organ')
    organ_express_total <- rbind(organ_express_total,organ_data_plot)
  }
  # remove proteins that are not measured in any tissue
  all_organ_mean <- all_organ_mean[which(rowSums(is.na(all_organ_mean))< length(tissue_df$tissue),),] 
  m <- length(setdiff(data[,1], all_organ_mean[,1]))
  print (paste("There are ", m, "proteins that are not measured in any tissue"))
  newlist <- list(all_organ_mean, organ_express_total)
}

###################################################

# plot protein distribution in all organs
plot_Mean <- function(organ_express_total, savename){
ggplot(as.data.frame(organ_express_total), aes(x = log10mean, y = organ, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "log10mean", option = "C") +
  labs(title = paste('Protein expression in', savename))
ggsave(paste("Figs/", date(),savename,'allTissue_proteindensityMap.png'))
}

###################################################
# generate gene id pathway table from wpid 
# input: list of wpid to get info

getInfo_wpid <- function(pwy) {
  inf <- rWikiPathways::getPathwayInfo(pwy)
}

getXrefs_wpid <- function(pwy) {
  xrefs <- as.data.frame(rWikiPathways::getXrefList(pwy,"L"))
  xrefs$pwy <- pwy
  colnames(xrefs) <- c("gene","pwy")
  return(xrefs)
}

generate_Info_table_from_wpid <- function(pwy){
x <- lapply(pwys, getInfo_wpid)
x1 <- do.call(rbind.data.frame, x)

y <- lapply(pwys, getXrefs_wpid)
y1 <- do.call(rbind.data.frame, y)

added_p <- data.frame(matrix(ncol=5,nrow=nrow(y1))) # 5 collumns as "name"    "version" "wpid"    "org"     "gene"   
colnames(added_p) <- c("name", "version", "wpid", "org", "gene" )
for (i in 1: length(x1$id)){
  index <- which(y1$pwy == x1$id[i])
  added_p$name[index] <- x1$name[i]
  added_p$version[index] <- x1$revision[i]
  added_p$wpid[index] <- x1$id[i]
  added_p$org[index] <- x1$species[i]
  added_p$gene[index] <- y1$gene[index]
}
added_p
}

###################################################
## Combine gene list in Wikipathways and the COVID19 disease map
# input:
# wp_file, dm_file: string text indicates the gmt file name for wikipathways (wp_file) and the disease map (dm_file)
# pwys (optional) : wikipathways id of additional pathways if any
combineWP_DM <- function(wp_file,dm_file, pwys){
  wp2gene <- readPathwayGMT(gsub(" ","",wp_file))
  dm2gene <- readPathwayGMT(gsub(" ","",dm_file))
  if(missing(pwys)){
    pwy2gene <- dplyr::bind_rows(wp2gene, dm2gene)
    } else{
    added_p <- generate_Info_table_from_wpid(pwys)
    pwy2gene <- dplyr::bind_rows(wp2gene, dm2gene, added_p)
  }
}
###################################################

## Extract pathway IDs from WikiPathways and Disease Map
## Input:
#   - all_organ_mean: dataframe that contains Uniprot ID and gene names from the data file (headers : Uniprot.ID)
#   - wp_file, dm_file: gmt files for wikipathways (wp) and the disease map (dm)
#   - pwys (optional) : additional wpids if they are not in the wikipathways gmt file
## Output:
#   - pathway: dataframe with columns as ENTREZID, Uniprot.ID, Gene.name,SYMBOL, name (pathway name), 
#              version, wpid (pathway ID in Wikipathways and Disease Map), org (organism)

pathway_ID <- function(all_organ_mean, wp_file, dm_file, pwys){
  # Map dataset IDs (uniprot) to databases (add NCBI Gene / Entrez Gene column) 
  all_organ_mean <- as.data.frame(all_organ_mean)
  unip2entrez <- clusterProfiler::bitr(all_organ_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
  data <- merge(all_organ_mean, unip2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
  # Remove duplicate that are resulted from multiple matches
  data2 <- data[!duplicated(data$Uniprot.ID),]
  # Remove gene that are not mapped
  data.covid <- data2 %>% tidyr::drop_na(ENTREZID) # 12588 genes
  # Generate pathway gene list from wikipathway and diseasemap 
  pwy2gene <- combineWP_DM(wp_file, dm_file,pwys)
    # merge gene in dataset with pathway-gene list in wikipathways and disease map through entrezid
  pathway <- merge(data.covid,pwy2gene, by.x = "ENTREZID", by.y = 'gene', all.y = TRUE)
}

###################################################
# Loop over all unique pathways -> get measure proteins per pathway and calculate median for each tissue  
# P_median_per_tissue <- function(pathway,tissue_df){
# uP <- data.frame(matrix(NA,nrow = length(unique(pathway$wpid)), ncol = 9))
# colnames(uP) <- c("wpid","geneCount",tissue_df$tissue)
# uP$wpid <- unique(pathway$wpid) # 665 unique pathways
# 
# tissue_median <- data.frame()
# for (i in 1 :length(uP$wpid)){
#   organ_p <- pathway %>% filter(wpid == uP$wpid[i])
#   organ_p <- unique(organ_p)
#   uP$geneCount[i] <- length(unique(organ_p[,1]))
#   for (j in 1:length(tissue_df$tissue)){
#     b <- get(tissue_df$tissue[j],organ_p)
#     uP[i,j+2] <- sum(!is.na(b))
#     tissue_median [i,j] <- median(b, na.rm = TRUE)
#   }
#   names(tissue_median)[1:7] <- tissue_df$tissue
# }
# 
# tissue_specific <- data.frame(uP$wpid,tissue_median)
# Pathway_median_per_tissue <- tissue_specific[order(tissue_specific$lung),]
# 
# }

P_median_per_tissue <- function(pathway,tissue_df){
  # uP <- data.frame(matrix(NA,nrow = length(unique(pathway$wpid)), ncol = 9))
  # colnames(uP) <- c("wpid","geneCount",tissue_df$tissue)
  uP <- unique(pathway$wpid) # 665 unique pathways
  
  tissue_median <- data.frame()
  for (i in 1 :length(uP)){
    organ_p <- pathway %>% filter(wpid == uP[i])
    organ_p <- unique(organ_p)
    for (j in 1:length(tissue_df$tissue)){
      b <- get(tissue_df$tissue[j],organ_p)
      tissue_median [i,j] <- median(b, na.rm = TRUE)
    }
    names(tissue_median)[1:7] <- tissue_df$tissue
  }
  
  tissue_specific <- data.frame(uP,tissue_median)
  Pathway_median_per_tissue <- tissue_specific[order(tissue_specific$lung),]
  
}

###################################################
# plot tissue specific activities 
# Input: pathway_median_per_tissue from P_median_per_tissue function, levels: the order of tissues to plot
# Example: plot_pathway(test_median, tissue_df, 'pathway_test_new_mean_function',
# levels = c("testi","kidney","spleen","liver","heart","lung","thyroid"))

plot_pathway <- function(pathway_median_per_tissue, tissue_df, savename, levels){
pathway_total <- data.frame()
for (i in 1: length(tissue_df$tissue)){
  organ <- tissue_df$tissue[i]
  pathway_data <- pathway_median_per_tissue %>% dplyr::select(contains(organ))
  pathway_data$organ <- organ
  names(pathway_data) <- c('median','organ')
  pathway_total <- rbind(pathway_total,pathway_data)
}
# remove rows with all NA resulted from pathway without genes such as WP5035
pathway_total <- pathway_total[as.logical((rowSums(is.na(pathway_total))-ncol(pathway_total))),]
# reorder tissues to plot
pathway_total$organ <- factor(pathway_total$organ, levels)

ggplot(pathway_total, aes(x = median, y = organ, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "median", option = "C") +
  labs(title = paste('Pathway activies per tissue in', savename))
ggsave(paste("Figs/", date(),savename,'All_tissue_pathway_activities.png'))

}

###################################################
## Statistic table

# Count number of proteins measured in the pathways 
# Input: all_organ_mean from protein_mean function , pwy2gene from combineWP_DM function
count_proteins <- function(all_organ_mean,pwy2gene){
#Convert dataset IDs (uniprot) to database Ids (entrez)
uniprot2entrez <- clusterProfiler::bitr(all_organ_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
data <- merge(all_organ_mean, uniprot2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
# remove duplicate resulted from multiple matches (1 uniprot can map to many entrezid)
data2 <- data[!duplicated(data$Uniprot.ID),]
#remove ids that are not mapped
data.covid <- data2 %>% tidyr::drop_na(ENTREZID) 

database_protein <- unique(pwy2gene$gene) #pwy2gene: all entrezid genes from wikipathways and diseasemap
table <- data.frame(matrix(NA,nrow = length(database_protein), ncol = 10))
colnames(table) <- c("gene","measured",tissue_df$tissue,"Uniprot.ID")
table$gene <- database_protein

for (i in 1:length(database_protein)){
  index <- which(data.covid$ENTREZID == database_protein[i]) #data.covid: gene in the dataset
  if (!is_empty(index)){
    table$measured[i] <- 'yes'
    table[i,3:9] <- data.covid[index[1],3:9]
    table[i,10] <- data.covid$Uniprot.ID[index[1]]
  } else {
    table$measured[i] <- 'no'
    table[i,3:9] <- NA
  }
}
table
}

# the number of pathways each of the measured protein involves in 
# Input: pathway from pathway_ID function
pathways_per_protein <- function(pathway,n,m1,m2){
pathway_data <- pathway[which(rowSums(is.na(pathway[,m1:m2]))<n),] 
unique_count <- data.frame(unique(na.omit(pathway_data$Uniprot.ID)))
colnames(unique_count) <- 'protein'
for (i in 1 :length(unique_count$protein)){
  index <- which(pathway$Uniprot.ID == unique_count$protein[i])
  count_pathway <- length(unique(pathway_data$wpid[index]))
  unique_count$pathway [i] <- count_pathway
}
unique_count
}

# count percent measured proteins per pathway per tissue
# input: pathway from pathway_ID function
prot_measured_per_pathway <- function(pathway,n,m1,m2){
# pathway_data <- pathway[which(rowSums(is.na(pathway[,m1:m2]))<n),] 
unique_wpid <- unique(pathway$wpid)
table2 <- matrix(nrow = length(unique_wpid), ncol=9)
colnames(table2) <- c("wpid","total_genes",colnames(pathway)[m1:m2])
table1 <- table2 
for (j in 1:length(unique_wpid)){
  one_wp <- unique(pathway[which(pathway$wpid == unique_wpid[j]),])
  table2[j,1] <- one_wp$wpid[1]
  table2[j,2] <- length(unique(one_wp$ENTREZID))
  table1[j,1] <- one_wp$wpid[1]
  table1[j,2] <- length(unique(one_wp$ENTREZID))
  for (i in m1:m2){
    table2[j,i-1] <- as.numeric(round(length(which(!is.na(one_wp[,i])))/ length(one_wp$Uniprot.ID), digits = 2))
    table1[j,i-1] <- as.numeric(length(which(!is.na(one_wp[,i]))))
  }
}
newlist <- list(table1,table2)
}

## Extract pathways with at least 30% gene measured and 3 gene measured
## input: tol1: 30 (30% gene measured), tol2: 3 (3 genes measured)
tol_measured_pathway <- function(pathway,n,m1,m2, tol1, tol2){
count_measured_table <- prot_measured_per_pathway(pathway,n,m1,m2)
percent_measured_table <- as.data.frame(count_measured_table[[2]])
percent_measured_table[,c(3:9)] <- sapply(percent_measured_table[,c(3:9)], as.numeric)
pathway_with_tol_measure <- percent_measured_table[which(rowSums(percent_measured_table[,3:9] > tol1/100) > 0),1]
number_measured_table <- as.data.frame(count_measured_table[[1]])
number_measured_table[,c(3:9)] <- sapply(number_measured_table[,c(3:9)], as.numeric)
pathway_with_3_gene <- number_measured_table[which(rowSums(number_measured_table[,3:9] > tol2) > 0),1]
measured_pathways <- unique(intersect(pathway_with_tol_measure,pathway_with_3_gene))
pathway_data <- data.frame()
for (i in 1: length(measured_pathways)){
  index <- which(pathway$wpid == measured_pathways[i])
  combine_data <- pathway[index,]
  pathway_data <- rbind(pathway_data,combine_data)
}
pathway_data
}

# generate statistic table and histogram plot 
# Input:
# - all_organ_mean from protein_mean, pathway_data from measured_pathway function
# - tol: the minimum percentage of measured proteins per pathways (i.e. 30%)
protein_pathway_statistic <- function(all_organ_mean, wp_file, dm_file, pwys, pathway, n, m1, m2,tol1, tol2){
# the total number of protein in the dataset
total_protein <- length(unique(na.omit(all_organ_mean$Uniprot.ID ))) # 11394 

# the number of pathways in WikiPathways and COVID19 Disease Map gene sets
pwy2gene <- combineWP_DM(wp_file, dm_file, pwys)
wiki_pathway <- length(unique(pwy2gene$wpid[grep("^WP[0-9]*$", pwy2gene$wpid)]))   # 701 (10-09-2021) + 4 (length(pwys)) others from the portal
dm_pathway <- length(unique(pwy2gene$wpid[grep("^MINERVA", pwy2gene$wpid)])) #21 (june 2021)
total_pathway <- length(unique(pwy2gene$wpid))  #722 + length(pwys)

# the number of pathways with measured proteins 
pathway_data <- pathway[which(rowSums(is.na(pathway[,m1:m2]))<n),] 
pathway_measured <- length(unique(na.omit(pathway_data$wpid))) # 723

# the number of pathways with measured proteins in COVID19 portal of Wikipathways
pwy <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
in_covidportal <- length(unique(intersect(pwy, pathway_data$wpid)))
not_inportal <- setdiff(pwy,pathway_data$wpid) #"WP4965" in mouse ; "WP5035" does not have gene

# the number of measured proteins
protein_in_pathways <- length(unique(na.omit(pathway_data$Uniprot.ID))) #5235 

# the number of proteins in dataset but not in pathways
protein_not_in_pathways <- length(unique(setdiff(all_organ_mean$Uniprot.ID,pathway$Uniprot.ID)))

# the number of proteins in the pathways but not measured in the dataset
table <- count_proteins(all_organ_mean, pwy2gene)
length(which(table$measured == 'no')) # 2863 

# the number of pathways each of the measured protein involves in
unique_count <- pathways_per_protein(pathway, n,m1,m2)
par(mar=c(4,5,1,2))
# par(mar=c(5.1, 4.1, 4.1, 2.1)) #default
histo <- hist(unique_count$pathway, col ='#E66100',main = 'Number of pathways per protein' , 
              xlab ='Number of Pathways that each protein involves in', 
              ylab='# Protein', freq = TRUE, border = '#E66100',
              xlim = c(0,max(unique_count$pathway)), cex.main = 2, cex.axis = 2, cex.lab = 2)

text(x = histo$mids, y = histo$counts, labels = histo$counts, cex = 1.5, pos = 3)
# count number of pathways with at least 30% gene measured and 3 genes measured in at least one tissue 
# count_measured_table <- prot_measured_per_pathway(pathway,n,m1,m2)
# percent_measured_table <- count_measured_table[[2]]
# pathway_with_tol_measure <- percent_measured_table[which(rowSums(percent_measured_table[,3:9] > tol1/100) > 0),1]
# number_measured_table <- count_measured_table[[1]]
# pathway_with_3_gene <- number_measured_table[which(rowSums(number_measured_table[,3:9] > tol2) >0),1]
# measured_pathways <- intersect(pathway_with_tol_measure,pathway_with_3_gene)
measured_pathways <- tol_measured_pathway(pathway,n,m1,m2,tol1,tol2)
total_measured_pathways_with_tol <- length(unique(measured_pathways$wpid))

print(paste(total_protein, 'proteins in the dataset'))
print(paste(length(unique(pwy2gene$gene)),'genes in wikipathways and disease map'))
print(paste(total_pathway,'pathways in wikipathways and disease map'))
print(paste(wiki_pathway, 'wikipathways'))
print(paste(dm_pathway,'disease map pathways'))
print(paste(pathway_measured,'pathways with measured proteins'))
print(paste(in_covidportal, 'pathways with measured proteins in covid19 portal in wikipathways'))
print(paste(protein_in_pathways, 'proteins in dataset found in pathways'))
print(paste(protein_not_in_pathways, 'proteins in dataset but not in pathways'))
print(paste(total_measured_pathways_with_tol ,'pathways with at least', tol1, '% gene measured and ', tol2, 'genes measured in at least one tissue'))
}
###################################################
# Extract top n most active pathways per tissue

# active_pathway <- function(Pathway_median_per_tissue, n, savename){
# active_pathway_per_tissue <- list()
# for (i in 1: length(tissue_df$tissue)){
#     tissue = tissue_df[i,]
#     active_tissue <- pathway_median[order(-pathway_median[i+1]),]
#     active_tissue_topN <- active_tissue [1:n,]
#     # plot top 10 pathways in lung and in others
#     
#     active_tissue_topN_plot <- data.frame(active_tissue_topN[1], stack(active_tissue_topN[2:ncol(active_tissue_topN)]))
#     active_tissue_topN_plot  <- active_tissue_topN_plot  %>% dplyr::rename(pathway_IDs = uP.wpid, median = values, tissue = ind) 
#     ggplot(active_tissue_topN_plot , aes(x= pathway_IDs, y=median, color=tissue)) + 
#       geom_point() + labs(title = tissue)
#     scale_color_brewer(palette = "Set1")
#     ggsave(paste("Figs/", date(),tissue,'_top', n, '_pathways in', savename,'.png'))
#     active_pathway_per_tissue [[i]] <- active_tissue_topN[,1]
# }
# topN_active_pathway_per_tissue <- do.call(cbind, active_pathway_per_tissue)
# }


active_pathway <- function(Pathway_median_per_tissue, n, savename){
  active_pathway_per_tissue <- list()
  for (i in 1: length(tissue_df$tissue)){
    tissue = tissue_df[i,]
    active_tissue <- Pathway_median_per_tissue[order(-Pathway_median_per_tissue[i+1]),]
    active_tissue_topN <- active_tissue [1:n,]
    # plot top 10 pathways in lung and in others
    
    active_tissue_topN_plot <- data.frame(active_tissue_topN[1], stack(active_tissue_topN[2:ncol(active_tissue_topN)]))
    active_tissue_topN_plot  <- active_tissue_topN_plot  %>% dplyr::rename(pathway_IDs = uP.wpid, median = values, tissue = ind) 
    ggplot(active_tissue_topN_plot , aes(x= pathway_IDs, y=median, color=tissue)) + 
      geom_point() + labs(title = tissue)
    scale_color_brewer(palette = "Set1")
    ggsave(paste("Figs/", date(),tissue,'_top', n, '_pathways in', savename,'.png'))
    active_pathway_per_tissue [[i]] <- active_tissue_topN[,1]
  }
  topN_active_pathway_per_tissue <- do.call(cbind, active_pathway_per_tissue)
}

###################################################
# Cytoscape visualization

cytoscape_visualization <- function(pathway_to_check,savename, pathway_data, min_value, max_value, check_columns,heatmap_colors){
RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check)) 
toggleGraphicsDetails()
loadTableData(pathway_data, data.key.column = "ENSEMBL", table.key.column = "Ensembl")

# apply visual style 
data.values = c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeCustomHeatMapChart(check_columns, data.values, node.colors, range = c(min_value,max_value), zeroLine = TRUE,
                          orientation = "HORIZONTAL",style.name = "WikiPathways",  colors = heatmap_colors)


# Saving output
svg.file <- file.path(getwd(), paste("Figs/",savename,".svg"))
exportImage(svg.file,'SVG')
png.file <- file.path(getwd(), paste("Figs/",savename,".png"))
exportImage(png.file,'PNG', zoom = 500)
cys.file <- file.path(getwd(), paste("Figs/",savename,".cys"))
saveSession(cys.file) 
# RCy3::closeSession(save.before.closing = F)

}

###################################################
# Master function: return the top n active pathway per organ and their median values. 
# Input: data, n: number of active pathways to extrat 

Master_func_active_pathway <- function(data, tissue_df, n, savename){
  ProteinMean <- protein_mean(data,tissue_df)
  plot_Mean(ProteinMean[2],savename)
  pathway <- pathway_ID(ProteinMean[1])
  pathway_median <- P_median_per_tissue(pathway, tissue_df)
  plot_pathway(pathway_median, tissue_df, savename)
  active_pathway_per_tissue <- active_pathway(pathway_median, n, savename)
  newlist <- list(pathway, pathway_median,active_pathway_per_tissue)
}
###################################################
#For heat map
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
