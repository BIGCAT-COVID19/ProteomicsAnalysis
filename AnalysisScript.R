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

# Load functions
source("MultiOrganProteomicFunc.R")
# Define tissues
tissue_df <- data.frame(tissue = c("lung","spleen","liver","heart","kidney","testi","thyroid")) #renal medulla and cortex are grouped to kidney

# All COVID patients: load data, calculate and plot protein mean and pathways activities
dataname <- "Data/mmc2_covid.csv" 
df_label <- read.csv(dataname) # to extract patient labels
covid_label <- df_label[,2]
df_data <- read.csv('Data/mmc3.csv') # protein matrix of all covid and noncovid patients
covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", covid_label)])

# Calculate protein means
ProteinMean <- protein_mean(covid_data,tissue_df)

# Plot protein means
# savename <- 'allCOVID'
# plot_Mean(ProteinMean[2],savename)

# Map to pathway databases
# Load pathway database: gmt file from wikipathways and the disease map
wp_file <- "Data/wikipathways-20210910-gmt-Homo_sapiens.gmt"
dm_file <- "Data/COVID19_DiseaseMap_June2021.gmt"

pwys <- c("WP4936","WP5020","WP5021","WP5035") # extra pathways in wikipathways covid portal but not in the gmt file
all_pathway <- pathway_ID(ProteinMean[1], wp_file, dm_file, pwys)

## Statistic
n <-  length(tissue_df$tissue)
m1 <- 4 # (the first tissue column in all_pathway)
m2 <- 10 # (the last tissue column in all_pathway)
tol1 <- 30 # at least 30% gene measured per pathway
tol2 <- 3 # at least 3 genes measured per pathway
all_organ_mean <- ProteinMean[[1]]
protein_pathway_statistic(all_organ_mean, wp_file, dm_file, pwys, all_pathway, n, m1, m2,tol1,tol2)

# Remove pathways with less than 30% gene measured and less than 3 genes measured
measured_pathways <- tol_measured_pathway(all_pathway,n,m1,m2,tol1,tol2)
pathway_median <- P_median_per_tissue(measured_pathways, tissue_df)
# savename <- "COVID"
# plot_pathway(pathway_median,tissue_df,savename,levels = c("testi","kidney","spleen","liver","heart","lung","thyroid"))

# nonCOVID patients: load data, calculate and plot protein mean and pathways activities

dataname <- "Data/mmc2_noncovid.csv"
savename <- 'nonCOVID' # name to save figures
noncovid_df_label <- read.csv(dataname) # to extract patient labels
non_covid_label <- noncovid_df_label[,2]
non_covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", non_covid_label)])

NC_ProteinMean <- protein_mean(non_covid_data,tissue_df)
# plot_Mean(NC_ProteinMean[2],savename)
NC_all_pathway <- pathway_ID(NC_ProteinMean[1], wp_file, dm_file, pwys)

NC_measured_pathways <- tol_measured_pathway(NC_all_pathway,n,m1,m2,tol1,tol2)
NC_pathway_median <- P_median_per_tissue(NC_measured_pathways, tissue_df)
# plot_pathway(NC_pathway_median,tissue_df,savename, levels = c("testi","kidney","spleen","liver","heart","lung","thyroid"))

## Compare COVID vs NonCOVID

# any pathway in COVID but not in nonCOVID
setdiff(pathway_median$uP, NC_pathway_median$uP) # "WP4300" "WP2363" "WP4483" "WP382"  "WP2817" "WP268" 

difference <- data.frame()
for (i in 1: length(pathway_median[,1])){
  index <- which(NC_pathway_median[,1] %in% pathway_median[i,1])
  if (!is_empty(index)){
  for (j in 1: length(tissue_df$tissue)){
    difference [i,j] <- pathway_median[i,j+1] - NC_pathway_median[index,j+1] 
  }
  } else {
    difference [i,j] <- pathway_median[i,j+1]
  }
  colnames(difference) <- tissue_df$tissue
}

difference <- data.frame(pathway_median[,1], difference)

# count gene measured 
uP <- prot_measured_per_pathway (measured_pathways,n,m1,m2)
nC_uP <- prot_measured_per_pathway (NC_measured_pathways,n,m1,m2)

number_count <- as.data.frame(uP[[1]])
number_count[,c(3:ncol(number_count))] <- sapply(number_count[,c(3:ncol(number_count))], as.numeric)
percent_count <- as.data.frame(uP[[2]])
percent_count[,c(3:ncol(number_count))] <- sapply(percent_count[,c(3:ncol(number_count))], as.numeric)

NC_number_count <- as.data.frame(nC_uP[[1]])
NC_number_count[,c(3:ncol(NC_number_count))] <- sapply(NC_number_count[,c(3:ncol(NC_number_count))], as.numeric)
NC_percent_count <- as.data.frame(nC_uP[[2]])
NC_percent_count[,c(3:ncol(NC_number_count))] <- sapply(NC_percent_count[,c(3:ncol(NC_number_count))], as.numeric)


# replace number to NA for pathway with less than 30% and 3 genes measured in difference
tol = 0.3 # 30% measured genes cutoff
tol2 = 3 # at least 3 genes measured

difference_tol <- data.frame()
for (i in 1: length(difference$pathway_median...1.)){
  index <- which(number_count[,1] == difference$pathway_median...1.[i])
  index_p <- which(percent_count[,1] == difference$pathway_median...1.[i])
  for (j in 1: length(tissue_df$tissue)) {
  if (number_count[index,j+2] > tol2 & percent_count[index_p,j+2] > tol){
    difference_tol[i,j] <- difference[i,j+1]
  } else {
    difference_tol [i,j] <- "NA"
  }
  }
}

difference_tol <- data.frame(difference[,1], difference_tol)
tissue_col <- c(2:ncol(difference_tol))
difference_tol[, tissue_col] <- sapply(difference_tol[,tissue_col], as.numeric)
colnames(difference_tol) <- c("wpid",tissue_df$tissue)

# print to test
for (i in 1: 10){
print(number_count[which(number_count$wpid == difference_tol$wpid[i]),])
}
difference_tol[1:10,]


# apply cut-off to generate overview result table of top active and downactive pathways per tissue

cut_off <- 0.1
active_table <- data.frame() 
for (i in 1: length(tissue_df$tissue)){
  index <- which(as.numeric(difference_tol[,i+1]) > cut_off)
  if (!is_empty(index)){
  up <- difference_tol[index,1]
  print(tissue_df$tissue[i])
  active_table[1,i] <- length(up)
  } else {
    up <- NA
    print(paste(tissue_df$tissue[i], 'has no up pathways'))
    active_table[1,i] <- 0
  }
  index2 <- which(difference_tol[,i+1] < -cut_off)
  if (!is_empty(index2)){
    down <- difference_tol[index2,1]
    active_table[2,i] <- length(down)
  } else {
    down <- NA
    print(paste(tissue_df$tissue[i], 'has no down pathways'))
    active_table[2,i] <- 0
  }
}
colnames(active_table) <-tissue_df$tissue
row.names(active_table) <- c('up','down')


# Plot
svg(paste(date(),"multiplot_up_down.svg"))
layout(matrix(c(1,1,2,3,4,5,6,7,8), nrow = 3, ncol = 3, byrow = TRUE))
plotdata <- data.matrix(active_table)
barplot(plotdata, main = "The number of up- and down- regulated pathways in COVID19 patients", xlab = 'Tissue',
        col= c('#2588b3','#920000'),beside = TRUE, border= NA, cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5)
legend("topleft",
       c("up","down"), cex = 1.5,bty = 'n',
       fill = c('#2588b3','#920000'))
for (i in 1: length(tissue_df$tissue)){
  tissue <- tissue_df$tissue[i]
  plot(difference[,i+1], col = '#5D3A9B', ylim = c(-0.7, 0.7) , pch = 20, cex = 0.9, 
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
       main = tissue, xlab ='Pathway index', ylab='Differences')
  abline(h=cut_off, col="black", lty = 2, lwd = 1.5)
  abline(h= -cut_off, col='black',lty = 2, lwd = 1.5)
}            # save image as sgv and add legend using inkscape
dev.off()

# extract pathways that change in at least 1 tissue, replace value as NA for tissue that not change
most_active <- data.frame(matrix(NA, nrow = nrow(difference_tol), ncol = ncol(difference_tol)))
most_active [,1] <- difference_tol[,1] # extract wpid

for (i in 1:length(tissue_df$tissue)){
  active_tissue <- difference_tol[,i+1] # extract difference value
  for (j in 1:length(active_tissue)){
    diff <- abs(difference_tol[j,i+1])
    if (!is.na(diff) & diff > cut_off) {
      most_active[j,i+1] <- active_tissue[j]
    }
  }
}
colnames(most_active) <- c('wpid',tissue_df$tissue) 

# Remove pathways that are not changed in any tissue
ind <- apply(most_active[,2:8],1,function(x) all(is.na(x)))
active_list <- most_active[!ind,]
# print out to test
for (i in 1:5) {
print(number_count[which(number_count[,1] == active_list[i,1]),])
print (active_list[i,])
print (difference_tol[which(difference_tol$wpid == active_list[i,1]),])
  
}

# extract pathways that change in at least 1 tissue, keep values for all tissue

n <- 1 
all_active <- data.frame(matrix(NA, nrow = nrow(difference_tol), ncol = ncol(difference_tol)))
for (i in 1: length(difference_tol$wpid)){
  index <- which(abs(difference_tol[i,2:ncol(difference_tol)]) > cut_off)
  if (!is_empty(index)){
    all_active[n,] <- difference_tol[i,]
    n <- n+1
  }
}

ind2 <- apply(all_active[,2:8],1,function(x) all(is.na(x)))
all_active_p <- all_active[!ind2,]
colnames(all_active_p) <- c('wpid',tissue_df$tissue) 

# Print to test
setdiff(active_list$wpid, all_active_p$wpid)

# heatmap

# mtrx <- matrix(unlist(active_list[,2:8]),nrow= nrow(active_list), ncol = 7)
# colnames(mtrx) <- colnames(active_list[,2:8])

mtrx <- matrix(unlist(all_active_p[,2:8]),nrow= nrow(all_active_p), ncol = 7)
colnames(mtrx) <- colnames(all_active_p[,2:8])

# png(paste("Figs/", date(), "heat_map_changed_pathways.png"), res = 300)
# 
# my_color <- c("#002992",'#ADD8E6','#920000')
# 
# heatmap.2(mtrx, na.color = 'grey', distfun = dist_no_na, col = my_color,
#           breaks = c(min(mtrx,na.rm=TRUE),-0.1,0.1,max(mtrx,na.rm = TRUE)), 
#           trace = "none", labRow = FALSE, )
# dev.off()


## gradients
stack_data2 <- data.frame(all_active_p$wpid, stack(all_active_p[2:8]))
rank_data2 <- stack_data2[order(-stack_data2[,2]),]
rank_data2 <- na.omit(rank_data2)
breaks = seq(min(rank_data2$values), max(rank_data2$values), length.out= length(rank_data2$values))

# define colors within zones
gradient1 <- colorpanel(sum(breaks[-1] < -0.1), "#002992", "#51a7db" )
gradient2 <- colorpanel(sum(breaks[-1] >= -0.1 & breaks[-1] <= 0), "#51a7db", "#ffffff")
gradient3 <- colorpanel(sum(breaks[-1] > 0 & breaks[-1] <=0.1), "#ffffff" , "#db8151")
gradient4 <- colorpanel(sum(breaks[-1] > 0.1), "#db8151", '#920000')

hm.colors = c(gradient1, gradient2, gradient3, gradient4)

heatmap.2(mtrx, na.color = 'grey', distfun = dist_no_na, col= hm.colors, breaks = breaks, trace = "none")

## Inspect heat map cluster

top_right <- c(42,21,26,19,8,17,1,43,32,18,40,39,14,44,11,38,34,15)
bottom <- c(36,27,13,7,25,20,41,30,28,16,23,29,35,31,37)

top_red <- all_active_p$wpid[top_right]
up_group <- filter(all_pathway,wpid %in% top_red)
unique(up_group$name)
bottom_blue <- all_active_p$wpid[bottom]
down_blue <- filter(all_pathway,wpid %in%bottom_blue)
unique(down_blue$name)

## rank overlap degree
active_list$rank <- rowSums(!is.na(active_list[,2:8]))
active_list[order(-active_list$rank),]

## Generate top active pathways table

stack_data <- data.frame(active_list$wpid, stack(active_list[2:8]))
rank_data <- stack_data[order(-stack_data[,2]),]
rank_data <- na.omit(rank_data)
plot(rank_data[,2])
# rank_table <-data.frame(rank_data$active_list.wpid[1:10], rank_data$active_list.wpid[73:82])

up_table <- matrix(nrow=10, ncol = 3)
for (i in 1:10){
  up_table[i,1] <- rank_data$active_list.wpid[i]
  up_table[i,2] <- all_pathway$name[which(all_pathway$wpid == rank_data$active_list.wpid[i])[1]]
  up_table[i,3] <- tissue_df$tissue[rank_data$ind[i]]
}
colnames(up_table) <- c("up pathway", "name", "tissue")


down_table <- matrix(nrow=10, ncol = 3)
for (i in 49:58){
  down_table[i-48,1] <- rank_data$active_list.wpid[i]
  down_table[i-48,2] <- all_pathway$name[which(all_pathway$wpid == rank_data$active_list.wpid[i])[1]]
  down_table[i-48,3] <- tissue_df$tissue[rank_data$ind[i]]
}

colnames(down_table) <- c("down pathway", "name", "tissue")
# write_csv(as.data.frame(up_table),paste(date(),'most_up_pathways'))
# write_csv(as.data.frame(down_table),paste(date(),'most_down_pathways'))


## generate visualization for changed pathways in each tissue
for (i in 1: length(tissue_df$tissue)){
  print (tissue_df$tissue[i])
  index <- which(!is.na(active_list[,i+1]))
  for (j in 1: length(index)) {
  pathway_to_check <- active_list[index[j],1]
  if (startsWith(pathway_to_check, "WP")){
  print (pathway_to_check)
  all_pathway$name[which(all_pathway$wpid == pathway_to_check)]
  pathway_data <- data.frame(all_pathway[which(all_pathway$wpid == pathway_to_check),],
                             NC_all_pathway[which(NC_all_pathway$wpid == pathway_to_check),])
  tissue_data <- data.frame(pathway_data[,1:3], pathway_data[,(i+3)],pathway_data[,11:19], pathway_data[,(i+19)], pathway_data[,27:31])
  check_columns <- c(tissue_df$tissue[i],gsub(" ","",paste(tissue_df$tissue[i], ".1")))
  colnames(tissue_data)[4] <- check_columns[1]
  colnames(tissue_data)[14] <- check_columns[2]
  RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check)) 
  toggleGraphicsDetails()
  loadTableData(pathway_data, data.key.column = "ENSEMBL", table.key.column = "Ensembl")
 # write_csv(pathway_data,paste(gsub(" ","", paste("Figs/", tissue_df$tissue[i], "/")),date(),pathway_to_check, tissue_df$tissue[i], ".csv"))
  # setNodeLabelMapping("SYMBOL")
  # heatmap_colors <- c('red', 'yellow', 'blue','grey')
  # setNodeCustomHeatMapChart(check_columns, data.values, node.colors, zeroLine = TRUE,
                            # orientation = "HORIZONTAL",style.name = "WikiPathways",  colors = heatmap_colors)
  setNodeCustomBarChart(check_columns, type = "GROUPED", colors = c("red","blue"), orientation = "HORIZONTAL", style.name = "WikiPathways")
  # Saving output
  path <- gsub(" ","",paste("Figs/",tissue_df$tissue[i]))
  dir.create(path)
  filename <- gsub(" ","",paste(path, "/",pathway_to_check))
  exportImage(filename,'SVG')
  exportImage(filename,'PNG', zoom = 500)
  saveSession(filename) 
  # RCy3::closeSession(save.before.closing = F)
    }
  }
}

# write file for each tissue

for (i in 1:length(tissue_df$tissue)){
  print(tissue_df$tissue)
  index <- which(!is.na(active_list[,i+1]))
  table <- data.frame()
    for (j in 1:length(index)){
    name <- all_pathway$name[which(all_pathway$wpid == active_list[index[j],1])]
    table [j,1] <- active_list[index[j],1]
    table [j,2] <- name[1]
    table [j,3] <- active_list[index[j],i+1]
    }
colnames(table) <- c("id", "name", "difference")
write_csv(table,paste(gsub(" ","", paste("Figs/", tissue_df$tissue[i], "/")),date(),"changed pathways in", tissue_df$tissue[i], ".csv"))

}


## Check log2FC, p value to identify significant difference

dif_data <- read.csv('Data/mmc4_covidvsnoncovid.csv')
dif_data <-  dif_data[2:nrow(dif_data),] # remove header line
# data order: log2FC, pvalue, adjusted pvalue  (significant adjusted p_value < 0.05, |log2FC| > log2[1.2])
lung <- c(3,4,5)
spleen <- c(6,7,8)
liver <- c(9,10,11)
heart <- c(12,13,14)
renal_cortex <- c(15,16,17)
renal_medulla <- c(18,19,20)
kidney <- c(15:20)
testi <- c(21,22,23)
thyroid <- c(24,25,26) 


## check number of proteins measured per pathway and number of protein that are up and down per pathway


inspect_pathway <- function(pathway_to_check, tissue, tissue_col){
  index <- which(number_count$wpid == pathway_to_check)
  total_protein <- number_count[index,2]
  measured_protein <- number_count[index,(tissue+2)]
  pathway_data <- data.frame(all_pathway[which(all_pathway$wpid == pathway_to_check),],
                             NC_all_pathway[which(NC_all_pathway$wpid == pathway_to_check),])
  up_protein <- which(pathway_data[,tissue+3] > pathway_data[,tissue+19])
  down_protein <- which(pathway_data[,tissue+3] < pathway_data[,tissue+19])
  
  tissue_dif_data <- cbind(dif_data$Uniprot.ID, dif_data$Gene.name, dif_data[,tissue_col])
  colnames (tissue_dif_data) <- c("Uniprot.ID", "Gene.name", "log2FC", "pvalue", "adjusted_pvalue")
  check_genes <- intersect(tissue_dif_data$Uniprot.ID, pathway_data$Uniprot.ID)
  not_in_log2fc <- setdiff(pathway_data$Uniprot.ID, tissue_dif_data$Uniprot.ID)
  
  significant_data <- tissue_dif_data[which(abs(as.numeric(tissue_dif_data$log2FC))>1.2 & as.numeric(tissue_dif_data$adjusted_pvalue) < 0.05),]
  sig_genes <- intersect(significant_data$Uniprot.ID, pathway_data$Uniprot.ID)
  sig_gene_names <- intersect(significant_data$Gene.name, pathway_data$Gene.name)
  
  print(paste(pathway_to_check, "has a total of", total_protein, "protein"))
  print(paste("of which", measured_protein, " proteins are detected"))
  print(paste(length(up_protein), "proteins are up in COVID"))
  print(paste(length(down_protein), "proteins are down in COIVD"))
  print(paste("there are", length(check_genes), "proteins in log2FC data"))
  print("these proteins are not in the log2FC data")
  print(not_in_log2fc)
  print(paste(length(unique(significant_data$Uniprot.ID)), "proteins are significantly different in", tissue_df$tissue[tissue]))
  print(paste(length(sig_genes), "significant proteins in", pathway_to_check))
  print("they are")
  print(sig_genes)
  print("gene names")
  print(sig_gene_names)
  }

pathway_to_check <- "WP3599"
all_pathway$name[which(all_pathway$wpid == pathway_to_check)][1]
tissue <- 3
tissue_col <- liver
inspect_pathway(pathway_to_check, tissue, tissue_col)

## All significant proteins that are not in wikipathways and DM collection

not_in_wiki <- setdiff(tissue_dif_data$Uniprot.ID, all_pathway$Uniprot.ID)
length(unique(not_in_wiki)) # 2589 proteins 

which(tissue_dif_data$Uniprot.ID == not_in_wiki[2])
which(all_pathway$Uniprot.ID == not_in_wiki[2])

## check significant data in a pathway
gene_list <- unique(pathway_data$Uniprot.ID)
filter(dif_data, Uniprot.ID %in% gene_list[!is.na(gene_list)])

## Network of all changed pathways 

changed_pathways <- filter(all_pathway,wpid %in% active_list$wpid)
dir.create("Results")
index <- !is.na(changed_pathways$Uniprot.ID)
node_list <- data.frame(changed_pathways$Uniprot.ID[index], changed_pathways$wpid[index])
write_csv(node_list, "Results/node_list")
cys_table <- read.csv("Results/node_list_cytoscape_table.csv")

pathway_list <- unique(changed_pathways$wpid) 
database <- data.frame(matrix(,ncol = 1, nrow = length(pathway_list)))
change_dif <- data.frame(matrix(,ncol = 7, nrow = length(pathway_list)))
colnames(change_dif) <- c("plung","pspleen","pliver","pheart","pkidney","ptesti","pthyroid")
for (i in 1:length(pathway_list)){
  index <- which(active_list$wpid == pathway_list[i])
  change_dif$plung[i] <- active_list[index,2]
  change_dif$pspleen[i] <- active_list[index,3]
  change_dif$pliver[i] <- active_list[index,4]
  change_dif$pheart[i] <- active_list[index,5]
  change_dif$pkidney[i] <- active_list[index,6]
  change_dif$ptesti[i] <- active_list[index,7]
  change_dif$pthyroid[i] <- active_list[index,8]
  if (startsWith(as.character(pathway_list[i]), "WP")){
    database$matrix...ncol...1..nrow...length.pathway_list..[i] <- "wikipathways"
  } else if (startsWith(as.character(pathway_list[i]),"MINERVA")){
    database$matrix...ncol...1..nrow...length.pathway_list..[i] <- "DiseaseMap"}
  }

all_table <- data.frame(pathway_list, database,change_dif)
colnames(all_table)[1:2] <- c("wpid","database")
write_csv(all_table,"Results/all_pathways_with_dif.csv")


# network of changed pathways in each tissue
testi_p <- active_list$wpid[!is.na(active_list$testi)]
testi_p_data <- filter(all_pathway,wpid %in% testi_p)
index <- !is.na(testi_p_data$Uniprot.ID)
testi_node_list <- data.frame(testi_p_data$Uniprot.ID[index], testi_p_data$wpid[index])
write_csv(testi_node_list, "Results/testi_node_list")

# print out pathway to group them 
active_p <- filter(all_pathway,wpid %in% active_list$wpid)
active_p <- data.frame(active_p$wpid, active_p$name)
active_p <- unique(active_p)
write_csv(active_p,"Results/all_changed_pathways")

# Check pathways of intestest
gene <- "IL6ST"
dif_data[which(dif_data$Gene.name == gene),] 
index <- which(all_pathway$Gene.name == gene)
all_pathway[index,]
filter(difference,pathway_median...1.%in% all_pathway[index,15])

pathway_to_check <- "WP5063"
all_pathway$name[which(all_pathway$wpid == pathway_to_check)][1]
filter(difference,pathway_median...1.%in% pathway_to_check)
i<- 1
tissue_col <- lung
inspect_pathway(pathway_to_check, i, tissue_col)

pathway_data <- data.frame(all_pathway[which(all_pathway$wpid == pathway_to_check),],
                           NC_all_pathway[which(NC_all_pathway$wpid == pathway_to_check),])
ggplot(data = pathway_data, aes(x = Gene.name)) +
  geom_point(aes(y = liver))  + geom_point(aes(y = liver.1), color = "red")
+ labs(title = "Protein abundant")





tissue_data <- data.frame(pathway_data[,1:3], pathway_data[,(i+3)],pathway_data[,11:19], pathway_data[,(i+19)], pathway_data[,27:31])
data.frame(tissue_data$liver, tissue_data$liver.1)
check_columns <- c(tissue_df$tissue[i],gsub(" ","",paste(tissue_df$tissue[i], ".1")))
colnames(tissue_data)[4] <- check_columns[1]
colnames(tissue_data)[14] <- check_columns[2]
RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check)) 
toggleGraphicsDetails()
loadTableData(pathway_data, data.key.column = "ENSEMBL", table.key.column = "Ensembl")
setNodeCustomBarChart(check_columns, type = "GROUPED", colors = c("red","blue"), orientation = "HORIZONTAL", style.name = "WikiPathways")
# Saving output
path <- gsub(" ","",paste("Figs/",tissue_df$tissue[i]))
dir.create(path)
filename <- gsub(" ","",paste(path, "/",pathway_to_check))
exportImage(filename,'SVG')
exportImage(filename,'PNG', zoom = 500)
saveSession(filename) 

