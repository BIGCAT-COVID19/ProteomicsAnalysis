# tets function 

library(clusterProfiler)
library(readr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(ggridges)
library(org.Hs.eg.db)
library(rWikiPathways)
library(RCy3)
library(RColorBrewer)
library(ggpubr)
library(VennDiagram)
library(msigdbr)
library(purrr)
library(gplots)
library(matrixStats)
## set working directory or open R from terminal 
setwd("/home/nhumpham/Desktop/COVID19-ZonMw/ANALYSIS_MULTI_ORGAN_PROTEOMIC")

# load custom-made functions for this analysis
source("MultiOrganProteomicFunc.R")
tissue_df <- data.frame(tissue = c("lung","spleen","liver","heart","kidney","testi","thyroid")) #renal medulla and cortex are grouped to kidney

# gmt file from wikipathways and the disease map
wp_file <- "Data/wikipathways-20210910-gmt-Homo_sapiens.gmt"
dm_file <- "Data/COVID19_DiseaseMap_June2021.gmt"

## Load data and calculate protein mean and pathway median
# All COVID patients
dataname <- "Data/mmc2_covid.csv" # all covid patients
savename <- 'allCOVID'
df_label <- read.csv(dataname) # to extract patient labels
covid_label <- df_label[,2]
df_data <- read.csv('Data/mmc3.csv') # protein matrix
covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", covid_label)])

ProteinMean <- protein_mean(covid_data,tissue_df)
# plot_Mean(ProteinMean[2],savename)
pwys <- c("WP4936","WP5020","WP5021","WP5035") # extra pathways in wikipathways covid portal but not in the gmt file
pathway <- pathway_ID(ProteinMean[1], wp_file, dm_file, pwys)
pathway_median <- P_median_per_tissue(pathway, tissue_df)
# plot_pathway(pathway_median,tissue_df,savename,levels = c("testi","kidney","spleen","liver","heart","lung","thyroid"))

# NON COVID patients
dataname <- "Data/mmc2_noncovid.csv"
savename <- 'nonCOVID' # name to save figures
noncovid_df_label <- read.csv(dataname) # to extract patient labels
non_covid_label <- noncovid_df_label[,2]
non_covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", non_covid_label)])

NC_ProteinMean <- protein_mean(non_covid_data,tissue_df)
#plot_Mean(NC_ProteinMean[2],savename)
NC_pathway <- pathway_ID(NC_ProteinMean[1], wp_file, dm_file, pwys)
NC_pathway_median <- P_median_per_tissue(NC_pathway, tissue_df)
# plot_pathway(NC_pathway_median,tissue_df,savename, levels = c("testi","kidney","spleen","liver","heart","lung","thyroid"))

################################################################################
## Statistic
n <-  length(tissue_df$tissue)
m1 <- 4
m2 <- 10
tol <- 30 # at least 30% gene measured per pathway
all_organ_mean <- ProteinMean[[1]]
protein_pathway_statistic(all_organ_mean, wp_file, dm_file, pwys, pathway, n, m1, m2,tol)

protein_in_pathways <- 5235 # these results are from the above function
protein_not_in_pathways <- 6149
total_protein <- 11384
# plot 
percent_in <- round((protein_in_pathways/total_protein)*100, digit = 2)
percent_out <- 100 - percent_in
# pie(c(protein_in_pathways,length(protein_not_in_pathways)),labels = c("44%", "56%"), col=c('palevioletred4','pink'), cex = 2)
# pie(c(protein_in_pathways,length(protein_not_in_pathways)),labels = c(percent_in, percent_out), col=c('#FFC20A','#0C7BDC'), cex = 2)
pie(c(protein_in_pathways,protein_not_in_pathways),labels = c(percent_in, percent_out), col=c('#E66100','#5D3A9B'), cex = 2)

legend("top",c('Proteins found in WikiPathways and COVID19 Disease Map', 'Proteins not in WikiPathways and COVID19 Disease Map'), 
       fill =c('#E66100','#5D3A9B'), cex = 1.5)


# Compare with other databases (Msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens")
length(unique(all_gene_sets$ensembl_gene)) # 43244
hgcn2entrez <- clusterProfiler::bitr(ProteinMean[[1]]$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
dataset_gene <- merge(ProteinMean[[1]], hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)

msig2uni <- clusterProfiler::bitr(all_gene_sets$ensembl_gene, fromType = "ENSEMBL", toType = c("ENTREZID", "UNIPROT"), OrgDb =  org.Hs.eg.db)
# 16.23% are fail to map

gene_in_others <- intersect(all_gene_sets$ensembl_gene, dataset_gene$ENSEMBL )
length(unique(gene_in_others)) #12518 

uni_gene_in_others <- intersect(msig2uni$UNIPROT,ProteinMean[[1]]$Uniprot.ID)
length(unique(uni_gene_in_others)) #11158 out of 11394


# count gene/pathway in covid and noncovid

uP <- percent_measured_per_pathway(pathway,n,m1,m2)
nC_uP <- percent_measured_per_pathway(NC_pathway,n,m1,m2)

## Compare COVID vs NonCOVID
difference <- data.frame()
for (i in 1: length(pathway_median[,1])){
  index <- which(NC_pathway_median[,1] %in% pathway_median[i,1])
  for (j in 1: length(tissue_df$tissue)){
    difference [i,j] <- pathway_median[i,j+1] - NC_pathway_median[index,j+1] 
  }
  colnames(difference) <- tissue_df$tissue
}

difference <- data.frame(pathway_median[,1], difference)

# combine with measured gene 
tol <- 0.3 #30% gene measured
gene_measured <- data.frame()
for (i in 1: length(difference$pathway_median...1.)){
  index <- which(uP[,1] == difference$pathway_median...1.[i])
  for (j in 1:length(tissue_df$tissue)){
  percent_measured <- uP[index,j+2]
  if (percent_measured > tol){
  gene_measured[i,j] <- 'yes'
  }
  else {gene_measured[i,j] <- 'no'}
}
}
colnames(gene_measured) <- tissue_df$tissue

# plot the differences with the define threshold to select most active pathways per tissue

# cut_off <- 0.07
# for (i in 1: length(tissue_df$tissue)){
#   tissue <- tissue_df$tissue[i]
#   png(paste("Figs/Pathway differences in", tissue, '.png'))
# plot(difference[,i+1], col=factor(gene_measured[,i]), ylim = c(-0.7, 0.7) ,main = paste('Differences in', tissue), xlab='pathway index',ylab='Median difference COVID vs NONCOVID')
# legend('topleft', legend = c('Pathways with more than 30% gene measured', 'pathways with less than 30% gene measured'), col = c('red','black'), pch = 1,cex = 1)
# abline(h=cut_off, col="black", lty = 2, lwd = 1.5)
# abline(h= -cut_off, col='black',lty = 2, lwd = 1.5)
# dev.off()
# 
# }


# apply cut-off to generate overview result table of top active and downactive pathways per tissue

cut_off <- 0.07
most_up <- data.frame()
most_down <- data.frame()
active_table <- data.frame() 
for (i in 1: length(tissue_df$tissue)){
 index <- which(difference[,i+1] > cut_off)
 up <- difference[index,1]
 print(tissue_df$tissue[i])
 most_up[1:length(up),i] <-up
 index2 <- which(difference[,i+1] < -cut_off)
 down <- difference[index2,1]
 most_down[1:length(down),i] <- down
 active_table[1,i] <- length(up)
 active_table[2,i] <- length(down)
 }
colnames(most_up) <- tissue_df$tissue
colnames(most_down) <- tissue_df$tissue
colnames(active_table) <-tissue_df$tissue
row.names(active_table) <- c('up','down')

png("Figs/0.07_cut_off_Number_of_active_pathways.png")
plotdata <- data.matrix(active_table)
barplot(plotdata, main = "The number of up- and down- regulated pathways in COVID19", xlab = 'tissue', col = c('palevioletred4','pink'),beside = TRUE, border= NA)
legend("topleft",
       c("up","down"),
       fill = c('palevioletred4','pink')
)
dev.off()

# plot all on the same page

layout(matrix(c(1,1,2,3,4,5,6,7,8), nrow = 3, ncol = 3, byrow = TRUE))
plotdata <- data.matrix(active_table)
barplot(plotdata, main = "The number of up- and down- regulated pathways in COVID19 patients", xlab = 'Tissue',
        col= c('#2588b3','#920000'),beside = TRUE, border= NA, cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5)
legend("topleft",
       c("up","down"), cex = 1.5,bty = 'n',
       fill = c('#2588b3','#920000'))

# for (i in 1: length(tissue_df$tissue)){
#   tissue <- tissue_df$tissue[i]
#   plot(difference[,i+1], col=factor(gene_measured[,i]), ylim = c(-0.7, 0.7) , pch = 20, cex = 0.9, 
#        cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
#        main = tissue, xlab ='Pathway index', ylab='Differences')
#   abline(h=cut_off, col="black", lty = 2, lwd = 1.5)
#   abline(h= -cut_off, col='black',lty = 2, lwd = 1.5)
# }            # save image as sgv and add legend using inkscape


for (i in 1: length(tissue_df$tissue)){
  tissue <- tissue_df$tissue[i]
  plot(difference[,i+1], col=c('#ADD8E6','#5D3A9B')[factor(gene_measured[,i])], ylim = c(-0.7, 0.7) , pch = 20, cex = 0.9, 
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
       main = tissue, xlab ='Pathway index', ylab='Differences')
  abline(h=cut_off, col="black", lty = 2, lwd = 1.5)
  abline(h= -cut_off, col='black',lty = 2, lwd = 1.5)
}            # save image as sgv and add legend using inkscape

# remove pathways with less than 30% measured genes
tol = 30 # 30% measured genes cutoff
most_active <- data.frame(matrix(NA,nrow = nrow(difference), ncol = ncol(difference)))
most_active [,1] <- difference[,1] # extract wpid
for (i in 1:length(tissue_df$tissue)){
  active_tissue <- difference[,i+1] # extract difference value
  for (j in 1:length(active_tissue)){
    index <- which(uP[,1] == difference[j,1]) 
    percent_measured <- uP[index,i+2]
    if (percent_measured > tol/100) {
      most_active[j,i+1] <- active_tissue[j]
    }
  }
}
colnames(most_active) <- c('wpid',tissue_df$tissue) # most_active from top: the highest positive means the top upregulated and bottom the highest negative means the most downregulated

# get gene count per changed pathway

# most_up
# most_down
# 
# lung_count <- matrix(nrow = length(most_up$lung), ncol = 3)
# for (i in 1:length(most_up$lung)){
#   index <- which(most_up$lung[i] == uP$wpid)
#   lung_count[i,1] <- most_up$lung[i]
#   lung_count[i,2] <- uP$geneCount[index]
#   lung_count[i,3] <- uP$lung[index]
# }
# lung_count <- lung_count[order(-(as.numeric(lung_count[,3]))),]
# lung_count



# get pathway names
for (i in 1: length (tissue_df$tissue)){
print(tissue_df$tissue[i])
pathway_check <- na.omit(most_up[,i])
  for (j in 1: length(pathway_check)){
    index <- which(pathway$wpid == pathway_check[j])
    print(pathway$name[index[1]])
  }
}

# Create heatmap to visualize overlap


mtrx <- matrix(unlist(difference[,2:8]),nrow= nrow(difference), ncol = 7)
colnames(mtrx) <- colnames(difference[,2:8])
# # col.order <- c("testi","kidney","spleen","liver","heart","lung","thyroid")
# col.order <- c("thyroid","lung","heart","liver","spleen","kidney","testi")
# mtrx2 <- mtrx[,col.order]


dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

# heatmap.2(mtrx3, na.color = 'yellow', distfun = dist_no_na, col =c('#5D3A9B','#ADD8E6','#E66100'),
#           breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
#           Colv= FALSE, Rowv = FALSE, trace = "none")
# heatmap.2(mtrx2, na.color = 'yellow', distfun = dist_no_na, col =c('#2588b3','#ADD8E6','#E66100'),
#           breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
#           Colv= FALSE, Rowv = FALSE, trace = "none", labRow = FALSE)
# 
# 
# heatmap.2(mtrx2, na.color = 'yellow', distfun = dist_no_na, col =c('#007792','#ADD8E6','#920000'),
#           breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
#           Colv= FALSE, Rowv = FALSE, trace = "none", labRow = FALSE)
heatmap.2(mtrx, na.color = 'yellow', distfun = dist_no_na, col =c('#2588b3','#ADD8E6','#920000'),
          breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
          Colv= FALSE, Rowv = FALSE, trace = "none", labRow = FALSE)

heatmap.2(mtrx, na.color = 'yellow', distfun = dist_no_na, col =c('#2588b3','#ADD8E6','#920000'),
          breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
          trace = "none", labRow = FALSE)


heatmap.2(mtrx2, na.color = 'yellow', distfun = dist_no_na, col =c('#2588b3','#ADD8E6','#920000'),
          breaks = c(min(mtrx,na.rm=TRUE),-0.07,0.07,max(mtrx,na.rm = TRUE)), 
          Colv= FALSE, Rowv = FALSE, trace = "none", labRow = FALSE,
          margins = c(5,10), sepwidth = c(1, 1))


# remove pathways that are not either up nor downregulated in at least one tissue
# change_count <- 0
# no_count <- 0
# no_change <- matrix (nrow = nrow(mtrx2), ncol = 8)
# change_mtrx <- matrix (nrow = nrow(mtrx2), ncol = 8)
# 
# for (i in 1:nrow(mtrx2)){
#   up_index <- which(mtrx2[i,] > 0.07)
#   down_index <- which(mtrx2[i,] < -0.07)
#   if (is_empty(up_index) & is_empty(down_index)){
#     no_count <- no_count +1
#     no_change [no_count,1] <- difference[i,1]
#     no_change [no_count,2:8] <- mtrx2[i,]
#   } else {
#       change_count <- change_count + 1
#       change_mtrx [change_count,1] <- difference[i,1]
#       change_mtrx [change_count,2:8] <- mtrx2[i,]
#     }
# }
# change_mtrx2 <- matrix(unlist(change_mtrx[1:change_count,2:8]),nrow= change_count, ncol = 7)
# class(change_mtrx2) <- "numeric"
# colnames(change_mtrx2) <- c("thyroid","lung","heart","liver","spleen","kidney","testi")

# remove pathways that are not either up nor downregulated in at least one tissue
change_count <- 0
no_count <- 0
no_change <- matrix (nrow = nrow(difference), ncol = 8)
change_mtrx <- matrix (nrow = nrow(difference), ncol = 8)

for (i in 1:nrow(difference)){
  up_index <- which(difference[i,] > 0.07)
  down_index <- which(difference[i,] < -0.07)
  if (is_empty(up_index) & is_empty(down_index)){
    no_count <- no_count +1
    no_change [no_count,1] <- difference[i,1] # pathway_id
    no_change [no_count,2:8] <- difference[i,2:8] # pathway data
  } else {
    change_count <- change_count + 1
    change_mtrx [change_count,1] <- difference[i,1] 
    change_mtrx [change_count,2:8] <- matrix(unlist(difference[i,2:8]), nrow =1, ncol = 7)
  }
}
colnames(change_mtrx) <- colnames(difference)
change_mtrx2 <- matrix(unlist(change_mtrx[1:change_count,2:8]),nrow= change_count, ncol = 7)
class(change_mtrx2) <- "numeric"
colnames(change_mtrx2) <- colnames(difference[,2:8])
heatmap.2(change_mtrx2, trace = 'none',na.color = 'yellow', distfun = dist_no_na, col =c('#002992','#ADD8E6','#920000'),
          breaks = c(min(change_mtrx2,na.rm=TRUE),-0.07,0.07,max(change_mtrx2,na.rm = TRUE)))


# Remove pathways with less than 30% proteins measured from the change mtrx
tol <- 30 
count_table <- percent_measured_per_pathway(pathway,7,4,10)
change_mtrx3 <- change_mtrx2
for (i in 1: length(change_mtrx2[,1])) {
  index <- which(count_table[,1] == change_mtrx[i,1])
  # check gene count for each tissue
  for (j in 1: length(tissue_df$tissue)){
    if (count_table[index,j+2] < tol/100){
    change_mtrx3[i,j] <- 'NA'
  }
  }
}

class(change_mtrx3) <- "numeric"
heatmap.2(change_mtrx3, na.color = 'yellow', distfun = dist_no_na, col =c('#002992','#ADD8E6','#920000'),
          breaks = c(min(change_mtrx3,na.rm=TRUE),-0.07,0.07,max(change_mtrx3,na.rm = TRUE)),
          Colv= FALSE, Rowv = FALSE, trace = "none", labRow = FALSE)

heatmap.2(change_mtrx3, trace = 'none',na.color = 'yellow', distfun = dist_no_na, col =c('#002992','#ADD8E6','#920000'),
          breaks = c(min(change_mtrx3,na.rm=TRUE),-0.07,0.07,max(change_mtrx3,na.rm = TRUE)))



# Create venn diagram to visualize overlap

most_up2 <- lapply(most_up, function(x) x[!is.na(x)])
most_down2 <- lapply(most_down, function(x) x[!is.na(x)])
xtest <- list(lung = most_up2$lung, spleen = most_up2$spleen,  
              kidney= most_up2$kidney, liver = most_up2$liver, thyroid = most_up2$thyroid)

xtest_down <- list(lung=most_down2$lung, liver = most_down2$liver, spleen = most_down2$spleen,
                   testi = most_down2$testi, thyroid = most_down2$thyroid)


xtest_updown <- list(lung=most_up2$lung, liver = most_down2$liver, spleen = most_down2$spleen,
                   testi = most_down2$testi, thyroid = most_down2$thyroid)
xtest_downup <- list(lung = most_down2$lung, spleen = most_up2$spleen,  
                     kidney= most_up2$kidney, liver = most_up2$liver, thyroid = most_up2$thyroid)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...,ext.text = FALSE)
  grid.draw(venn_object)
}

color = rainbow(5)
color = viridis(5)
color <- brewer.pal(n =5, 'Set2')
display_venn(xtest_downup,lty = 'blank',fill = color,   
             lwd = 2, fontface = 'bold', cex = 2, cat.cex = 3, cat.fontface = 'bold',
             cat.pos = c(-6,8,-120,80,18))

 ## rank overlap degree
# up
most_upStack <- stack(most_up2)
count_overlap <- data.frame()
for (i in 1: length(most_up2$lung)){
  index <- which(most_upStack == most_up2$lung[i])
  count_n <- data.frame(length(index))
  count_overlap <- rbind(count_overlap, count_n)
}

count_up <- data.frame(most_up2$lung,count_overlap)
count_up <- count_up[order(-count_overlap),]
count_up

# down
most_downStack <- stack(most_down2)
count_downoverlap <- data.frame()
for (i in 1: length(most_down2$spleen)){
  index <- which(most_downStack == most_down2$spleen[i])
  count_n <- data.frame(length(index))
  count_downoverlap <- rbind(count_downoverlap, count_n)
}

count_down <- data.frame(most_down2$spleen,count_downoverlap)
count_down <- count_down[order(-count_downoverlap),]
count_down


## check if there is any pathway up but down in other organ
most_up2 <- lapply(most_up, function(x) x[!is.na(x)])
most_down2 <- lapply(most_down, function(x) x[!is.na(x)])
result <- data.frame()
for (i in 1: length(tissue_df$tissue)){
  for (j in 1:length(tissue_df$tissue)){
    upbutdown <- paste('up in', tissue_df$tissue[i], "but down in", tissue_df$tissue[j])
    upbutdown_id <- intersect(most_up2[[i]], most_down2[[j]])
    upbutdown_name<- unique(pathway$name[which(pathway$wpid == upbutdown_id)])
    downbutup <- paste('down in', tissue_df$tissue[i], "but up in", tissue_df$tissue[j])
    downbutup_id <- intersect(most_down2[[i]],most_up2[[j]])
    downbutup_name <- unique(pathway$name[which(pathway$wpid == downbutup_id)])
    info <- c(upbutdown,upbutdown_name,downbutup,downbutup_name)
    write.table(info,file = paste(date(),'pathways with opposite activities among tissues.csv'), sep = ",", append = TRUE)
    }
}

## visualize pathways
# Create overview network








#set input for visualization
min_value <- 0.02
max_value <- 0.8
pathway_data <- data.frame(pathway, NC_pathway)
check_columns <- c(colnames(pathway_data)[i+3], colnames(pathway_data)[i+19]) # visualize colnames(pathway_data)[i+3] :covid vs colnames(pathway_data)[i+19]): covid
heatmap_colors <- c('red', 'yellow', 'blue','grey')

# extract pathway to visualize
i = 1 # tissue
i = i+1 # other tissue

j = 1 # pathway
j = j+1

pathway_to_check <- most_active[j,i]
# pathway_to_check <- 'WP3601'
# pathway_to_check <- 'WP716'
pathway_to_check <- "WP4211"
i <- 1
pathway_to_check <- most_down2$lung[i]
i <- i+1

pathway_to_check <- 'WP4142'

print (pathway_to_check)
pathway$name[which(pathway$wpid == pathway_to_check)]

savename = paste(pathway_to_check,tissue_df$tissue[i])
cytoscape_visualization(pathway_to_check, savename, pathway_data,min_value,max_value,check_columns, heatmap_colors)

#comment following line if you want to manipulate the visualization in Cytoscape
# RCy3::closeSession(save.before.closing = F)

## Generate overview network with gene-pathways 
# combine wp covid19 portal and the dm covid19 map

dm2gene <- readPathwayGMT(gsub(" ","","Data/COVID19_DiseaseMap_June2021.gmt"))
pwy2gene <- dplyr::bind_rows(wp2gene, dm2gene)

pwy <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")

x <- lapply(pwy, getInfo)
x1 <- do.call(rbind.data.frame, x)

y <- lapply(pwy, getXrefs)
y1 <- do.call(rbind.data.frame, y)
added_p <- data.frame(matrix(ncol=ncol(dm2gene),nrow=nrow(y1)))
colnames(added_p) <- colnames(dm2gene)
for (i in 1: length(x1$id)){
  index <- which(y1$pwy == x1$id[i])
  added_p$name[index] <- x1$name[i]
  added_p$version[index] <- x1$revision[i]
  added_p$wpid[index] <- x1$id[i]
  added_p$org[index] <- x1$species[i]
  added_p$gene[index] <- y1$gene[index]
}

wp_portal <- added_p
# combine WikiPathways covid19 portal and the disease map 
wp_portal_dm_combine <- dplyr::bind_rows(wp_portal, dm2gene)
network <- data.frame(wp_portal_dm_combine$wpid, wp_portal_dm_combine$gene)
colnames(network) <- c("source", "target")
RCy3::createNetworkFromDataFrames(edges=network) # export the network table in cytoscape
table <- read.csv('WpPortal_dm.csv')

for (i in 1:length(wp_portal_dm_combine$wpid)){
  index <- which(table$id == wp_portal_dm_combine$wpid[i])
  table$label [index] <- c('Pathway')
}


hgcn2entrez <- clusterProfiler::bitr(ProteinMean[[1]]$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
data <- merge(ProteinMean[[1]], hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
data.covid <- data %>% tidyr::drop_na(ENTREZID) # 12588 genes

table$measured2 <- NA
table$lung[i] <- NA
table$spleen[i] <- NA
table$liver[i] <- NA
table$heart[i] <- NA
table$kidney[i] <- NA
table$testi[i] <- NA
table$thyroid[i] <- NA
for (i in 55:length(table$id)){
  index <- which(data.covid$ENTREZID == table$id[i])
  if (!is_empty(index)){
    table$measured2[i] <- 'yes'
    table$lung[i] <- data.covid$lung[index]
    table$spleen[i] <- data.covid$spleen[index]
    table$liver[i] <- data.covid$liver[index]
    table$heart[i] <- data.covid$heart[index]
    table$kidney[i] <- data.covid$kidney[index]
    table$testi[i] <- data.covid$testi[index]
    table$thyroid[i] <- data.covid$thyroid[index]
  } else {
    table$measured2[i] <- 'no'
  }
}

write_csv(table, '23_09_2021_wpportal_dm')


# combine WikiPathways and COVID19 Disease Map gene sets
wpid2gene <- dm2gene %>% dplyr::select(wpid,gene) #TERM2GENE
colnames(wpid2gene) <- c("source", "target")
RCy3::createNetworkFromDataFrames(edges=wpid2gene)

wp_dm_combine <- pwy2gene
colnames(wp_dm_combine) <- c("name","version","source","org","target")
RCy3::createNetworkFromDataFrames(edges=wp_dm_combine)
write_csv(wp_dm_combine,'wp_dm_combine')
write_csv(pathway,'pathway')

# convert id to gene symbol
hgcn2entrez <- clusterProfiler::bitr(wpid2gene$target, fromType = "ENTREZID",toType = c("UNIPROT","SYMBOL","GENENAME"), OrgDb = org.Hs.eg.db)
dm <- merge(wpid2gene, hgcn2entrez, by.x="target", by.y="ENTREZID", all.x = TRUE)
write_csv(dm,'dm_wp_network')

RCy3::createNetworkFromDataFrames(edges=dm)



pwy <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
wp2gene.covid <- wp2gene[wp2gene$wpid %in% pwy,]
RCy3::createNetworkFromDataFrames(edges=wp2gene.covid)

toggleGraphicsDetails()
loadTableData(pathway, data.key.column = "ENTREZID", table.key.column = "ENTREZID")

min_value <- min(pathway$lung, na.rm = TRUE)
max_value <- max(pathway$lung, na.rm = TRUE)

# apply visual style 
data.values = c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))

setNodeCustomHeatMapChart(c("lung", 'liver'), data.values, node.colors, range = c(0.4, 0.8), zeroLine = TRUE,
                          orientation = "HORIZONTAL",style.name = "WikiPathways",  colors = c('yellow', 'red', 'blue','grey'))



write_csv(pathway, 'allCOVID_pathway')


pwy <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
wp2gene.covid <- wp2gene[wp2gene$wpid %in% pwy,]
colnames(wp2gene.covid$wpid) <- 'source' 

# https://www.statmethods.net/graphs/images/boxplot1.jpg 

new <- merge(pathway, NC_ProteinMean[1], by.pathway = 'Uniprot.ID')

covid_mean <- data.frame(ProteinMean[1])
noncovid_mean <- data.frame(NC_ProteinMean[1])
if(covid_mean$Uniprot.ID == noncovid_mean$Uniprot.ID)
{
  print("Column A and B are identical")
}

combine_protein_mean <- data.frame(covid_mean,noncovid_mean)
hgcn2entrez <- clusterProfiler::bitr(combine_protein_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL"), OrgDb = org.Hs.eg.db)
data <- merge(combine_protein_mean, hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
data.covid <- data %>% tidyr::drop_na(ENTREZID) # 12588 genes
write_csv(data.covid,'combine_protein_mean_covid_and_noncovid')


# prepare table to load to cytoscape

table <- read.csv('covidVsNonCovid_node_table.csv')

# for (i in 1: length(table$shared.name))
# { wp <- table$shared.name[i]
#   index <- which(uP$wpid == wp)
#   if (identical(index,integer(0)) == FALSE)
#   {
#   table$measuredLung[i] <- uP$lung[index]
#   table$measuredSpleen[i] <- uP$spleen[index]
#   table$measuredLiver[i] <- uP$liver[index]
#   table$measuredHeart[i] <- uP$heart[index]
#   table$measuredKidney[i] <- uP$kidney[index]
#   table$measuredTesti[i] <- uP$testi[index]
#   table$measuredThyroid[i] <- uP$thyroid[index]
#   
#   table$NotmeasuredLung[i] <- uP$geneCount[index] - uP$lung[index]
#   table$NotmeasuredSpleen[i] <- uP$geneCount[index] -uP$spleen[index]
#   table$NotmeasuredLiver[i] <- uP$geneCount[index] -uP$liver[index]
#   table$NotmeasuredHeart[i] <- uP$geneCount[index] -uP$heart[index]
#   table$NotmeasuredKidney[i] <- uP$geneCount[index] -uP$kidney[index]
#   table$NotmeasuredTesti[i] <- uP$geneCount[index] -uP$testi[index]
#   table$NotmeasuredThyroid[i] <- uP$geneCount[index] -uP$thyroid[index]
#   }
#   else { table$measuredLung[i] <- ''
#   table$measuredSpleen[i] <- ''
#   table$measuredLiver[i] <- ''
#   table$measuredHeart[i] <- ''
#   table$measuredKidney[i] <- ''
#   table$measuredTesti[i] <- ''
#   table$measuredThyroid[i] <- ''
#   table$NotmeasuredLung[i] <- ''
#   table$NotmeasuredSpleen[i] <- ''
#   table$NotmeasuredLiver[i] <- ''
#   table$NotmeasuredHeart[i] <- ''
#   table$NotmeasuredKidney[i] <- ''
#   table$NotmeasuredTesti[i] <- ''
#   table$NotmeasuredThyroid[i] <- ''
#   }
# }
# 
# write_csv(table, 'covid_noncovid_genecount_incovid')


## the most gene count per pathway table

# in covid
for (i in 1: length(table$shared.name)) {
  wp <- table$shared.name[i]
  index <- which(uP$wpid == wp)
  if (identical(index,integer(0)) == FALSE)
  { table$measured [i] <- max(uP[index,3:9])
  table$notMeasured [i] <- uP$geneCount[index] - max(uP[index,3:9])
  } else {
    table$measured [i]  <- ''
    table$notMeasured [i] <- ''
  }
}

# in noncovid

for (i in 1: length(table$shared.name)) {
  wp <- table$shared.name[i]
  index <- which(nC_uP$wpid == wp)
  if (identical(index,integer(0)) == FALSE)
  { table$nonCovid_measured [i] <- max(nC_uP[index,3:9])
  table$nonCovid_notMeasured [i] <- nC_uP$geneCount[index] - max(nC_uP[index,3:9])
  } else {
    table$nonCovid_measured [i]  <- ''
    table$nonCovid_notMeasured [i] <- ''
  }
}


write_csv(table, 'covid_nonCovid_Maxgenecount_covidvsnoncovid2')

## Add log2fc and pvalue 

FC_data <- read.csv('Data/mmc4_covidvsnoncovid.csv') # protein matrix


for (i in 1: length(table$shared.name)){
  gene <- table$Gene.name[i]
  index <- which(FC_data$Gene.name == gene)
  table$lungLog2FC[i] <- FC_data$Lung [index]
  table$lungpadjustedvalue[i] <- FC_data$X.1[index]
}

write_csv(table, 'lungFC_covid_nonCovid_Maxgenecount_covidvsnoncovid2')

