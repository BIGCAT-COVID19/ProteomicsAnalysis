# Calculate log2FC for all proteins in the raw data
covid_mean <- ProteinMean[[1]]
non_covid_mean <- NC_ProteinMean[[1]]
colnames(log2FC) <- tissue_df$tissue
log2FC <- data.frame()
for (i in 1:length(covid_mean$Uniprot.ID)){
index <- which(non_covid_mean$Uniprot.ID == covid_mean$Uniprot.ID[i])
FC <- covid_mean[i,3:9] / non_covid_mean[index,3:9]
df <- c(log(FC,2))
log2FC <- rbind(df,log2FC)
}

all_log2FC <- cbind(covid_mean$Uniprot.ID, covid_mean$Gene.name, log2FC)


index <- which(covid_mean$Uniprot.ID == "P35354")
nc_index <- which(non_covid_mean$Uniprot.ID == "P35354")
fc <- non_covid_mean[index,3:9] / covid_mean[index,3:9]
log3 <- log(fc,2)
dif_data[which(dif_data$Uniprot.ID == "P35354"),]


## DESeq2  

install.packages("htmltools")
library(htmltools)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("DESeq2")
library(DESeq2)
metaData <- read.csv("Data/airway_metadata.csv", header = TRUE, sep = ",")
countData <- read.csv("Data/airway_scaledcounts.csv", header = TRUE, sep =",")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~dex, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
head(res)


metaData_c <- read.csv("Data/Metadata.csv", header = TRUE, sep = ",")
head(metaData_c)
countData_c <- read.csv("Data/mmc3.csv", header = TRUE, sep = ",")
head(countData_c)
dds_c <- DESeqDataSetFromMatrix(countData = countData_c, colData = metaData_c, design = ~dex, tidy = TRUE)

a <- setdiff(colnames(countData_c), metaData_c$id)
setdiff(metaData_c$id, colnames(countData_c))
to_check <- intersect(colnames(countData_c), metaData_c$id)
new_countData_c <- countData_c %>% dplyr::select(contains(to_check))
final_countData_c <- data.frame(countData_c$Uniprot.ID, new_countData_c)
 # new_metaData_c <- filter(metaData_c,id %in% to_check)
 new_metaData_c <- unique(metaData_c)

dds_c <- DESeqDataSetFromMatrix(countData = final_countData_c, colData = metaData_c, design = ~dex, tidy = TRUE)

metaData_c[duplicated(metaData_c$id),]

index2 <- which(final_countData_c$countData_c.Uniprot.ID== "P35354")
nc_index <- which(non_covid_mean$Uniprot.ID == "P35354")
fc <- non_covid_mean[index,3:9] / covid_mean[index,3:9]
log3 <- log(fc,2)

organ_data <- final_countData_c %>% dplyr::select(contains('liver'))
c_organ_data <- organ_data %>% dplyr::select(contains(covid_label))
nc_organ_data <- organ_data %>% dplyr::select(contains(non_covid_label))

c_mean <- mean(as.numeric(c_organ_data[index2,]), na.rm= TRUE)
nc_mean <- mean(as.numeric(nc_organ_data[index2,]), na.rm = TRUE)
c_mean/nc_mean
log(3.4, 2)
dif_data[which(dif_data$Uniprot.ID == "P35354"),]

c_data <- final_countData_c %>% dplyr::select(contains(covid_label))
c_data <- data.frame(final_countData_c$countData_c.Uniprot.ID, c_data)
nc_data <- final_countData_c %>% dplyr::select(contains(non_covid_label))
nc_data <- data.frame(final_countData_c$countData_c.Uniprot.ID,nc_data)

count_protein_mean <- function(data, tissue_df){
  
# create empty dataframe that contains only gene ids and names
  all_organ_mean <- data.frame(data[,1])

  for (i in 1: length(tissue_df$tissue)){
    organ <- tissue_df$tissue[i]
    organ_data <- data%>% dplyr::select(contains(organ))
    # calculate log10 (mean+1) for protein that are measured in more than 10% of samples
    # proteins that are measured in less than 10% are left NA in mean values
    # log10(mean +1) to shift the plot to non-negative data
    organ_data$mean_value <- apply(organ_data,1, function(x){mean(x, na.rm = TRUE)})
    # Combine data with other organs into 1 data frame
    all_organ_mean <- cbind(all_organ_mean,organ_data[,ncol(organ_data)]) 
    names(all_organ_mean)[names(all_organ_mean) == "organ_data[, ncol(organ_data)]"] <- organ
    }
  # remove proteins that are not measured in any tissue
  all_organ_mean <- all_organ_mean[which(rowSums(is.na(all_organ_mean))< length(tissue_df$tissue),),] 
  m <- length(setdiff(data[,1], all_organ_mean[,1]))
  print (paste("There are ", m, "proteins that are not measured in any tissue"))
  all_organ_mean
}

covid_mean_count <- count_protein_mean(c_data, tissue_df)
colnames(covid_mean_count)[1] <- "Uniprot.ID"
nc_mean_count <- count_protein_mean(nc_data, tissue_df )
colnames(nc_mean_count)[1] <- "Uniprot.ID"

log2FC <- data.frame()
for (i in 1:length(covid_mean_count$Uniprot.ID)){
  index <- which(nc_mean_count$Uniprot.ID == covid_mean_count$Uniprot.ID[i])
  FC <- covid_mean_count[i,2:8] / nc_mean_count[index,2:8]
  df <- c(log(FC,2))
  log2FC <- rbind(log2FC, df)
}
colnames(log2FC) <- c("FClung","FCspleen", "FCliver","FCheart","FCkidney","FCtesti","FCthyroid")
log2FC$UniprotID <- covid_mean_count$Uniprot.ID
all_log2FC <- cbind(covid_mean_count$Uniprot.ID, log2FC)
i =29
dif_data[which(dif_data$Uniprot.ID == covid_mean_count$Uniprot.ID[i]),]
log2FC[i,]
write_csv(log2FC, "Results/Calculated_log2FC")

gene_ensemble <- data.frame(all_pathway$Uniprot.ID, all_pathway$Gene.name, all_pathway$ENSEMBL, all_pathway$ENTREZID)
gene_identifiers <- gene_ensemble[!is.na(gene_ensemble$all_pathway.Uniprot.ID),]
write_csv(gene_identifiers, "Results/gene_identifiers")
