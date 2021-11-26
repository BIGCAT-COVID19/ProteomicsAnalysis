# test why there are differences when extracting measured genes from pathway_data and from count_protein functions
source("MultiOrganProteomicFunc.R")



count_proteins <- function(all_organ_mean,pwy2gene){
  hgcn2entrez <- clusterProfiler::bitr(all_organ_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
  entrez2uni <- clusterProfiler::bitr(pwy2gene$gene, fromType = "ENTREZID", toType = c("UNIPROT","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
  
  data <- merge(all_organ_mean, hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
  data.covid <- data %>% tidyr::drop_na(ENTREZID) 
  
  
  data <- merge(pwy2gene, entrez2uni, by.x="gene", by.y="ENTREZID", all.x = TRUE)
  data.covid <- data %>% tidyr::drop_na(gene) 
  
  database_protein <- unique(data.covid$UNIPROT) #pwy2gene: all genes from wikipathways and diseasemap
  table <- data.frame(matrix(NA,nrow = length(database_protein), ncol = 10))
  colnames(table) <- c("gene","measured",tissue_df$tissue,"Uniprot.ID")
  table$gene <- database_protein
  
  for (i in 1:length(database_protein)){
    index <- which(all_organ_mean$Uniprot.ID == database_protein[i]) #data.covid: gene in the dataset
    if (!is_empty(index)){
      table$measured[i] <- 'yes'
      table[i,3:9] <- all_organ_mean[index[1],3:9]
      table[i,10] <- all_organ_mean$Uniprot.ID[index[1]]
    } else {
      table$measured[i] <- 'no'
      table[i,3:9] <- NA
    }
  }
  table
}




## Load dataset -> 11394 unique uniprot ids

# All COVID patients
dataname <- "Data/mmc2_covid.csv" # all covid patients
savename <- 'allCOVID'
df_label <- read.csv(dataname) # to extract patient labels
covid_label <- df_label[,2]
df_data <- read.csv('Data/mmc3.csv') # protein matrix
covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", covid_label)])

ProteinMean <- protein_mean(covid_data,tissue_df)

length(unique(covid_data$Uniprot.ID)) - length(unique(ProteinMean[[1]]$Uniprot.ID)) # should be 0 -> yes

# database -> 8092 unique entrezid 
# gmt file from wikipathways and the disease map
wp_file <- "Data/wikipathways-20210910-gmt-Homo_sapiens.gmt"
dm_file <- "Data/COVID19_DiseaseMap_June2021.gmt"

pwy2gene <- combineWP_DM(wp_file,dm_file) # 8057 unique genes
pwys <- c("WP4936","WP5020","WP5021","WP5035") 
added_p <- generate_Info_table_from_wpid(pwys) # add 87 unique genes
pwy2gene <- dplyr::bind_rows(pwy2gene, added_p) # total 8092 unique genes


## Method 1- pathway_data -> map dataset to entrezid

all_organ_mean <- ProteinMean[[1]]
hgcn2entrez <- clusterProfiler::bitr(all_organ_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
# 1.97% input gene are fail to map

data <- merge(all_organ_mean, hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
data.covid <- data %>% tidyr::drop_na(ENTREZID) # 12588 genes
length(unique(data.covid$Uniprot.ID)) # -> 11170 unique uniprot
length(unique(data.covid$ENTREZID)) # -> 11268 unique entrezid 


pathway <- merge(data.covid,pwy2gene, by.x = "ENTREZID", by.y = 'gene', all.y = TRUE)
length(unique(pathway$ENTREZID)) # -> 8092 entrezid genes (all the gene in the database)

length(unique(pathway$Uniprot.ID)) # -> 5242 unique uniprot ids

# remove unmeasured genes based on entrezid 

pathway2 <- data.frame(pathway[,1],pathway[,4:10])
pathway3 <- data.frame(pathway[,4:10])
pathway4 <- pathway2[which(rowSums(is.na(pathway3))< 7,),] 
length(unique(pathway4$pathway...1.)) # 5270 unique proteins 

# remove unmeasured genes based on uniprot

pathway5 <- data.frame(pathway$Uniprot.ID, pathway[,4:10])
pathway6 <- pathway5[which(rowSums(is.na(pathway5))< 7,),] 
length(unique(pathway6$pathway.Uniprot.ID)) #-> 5238 unique uniprot ids 

# the differences are due to the multi-matches of 1 uniprot to many entrez id 
pathway4a <- pathway4[!duplicated(pathway4[,2:8]),]
length(unique(pathway4a$pathway...1.)) # -> 5228 unique entrez id 

pathway6a <- pathway6[!duplicated(pathway6[,2:8]),]
length(unique(pathway6a$pathway.Uniprot.ID)) # -> 5234 unique entrez id 

pathway2a <- pathway2[!duplicated(pathway2[,2:8]),] # -> 5235 unique entrezid
length(pathway2a$pathway...1.)

pathway5 <- pathway5[!duplicated(pathway5[,2:8]),]
length(pathway5$pathway.Uniprot.ID)  # -> 5235 unique uniprot ids  # one empty uniprot here -> in the end the correct number is 5234   

setdiff(pathway_data$Uniprot.ID, table$Uniprot.ID) # "Q13948" "Q5JWF2" "O96033" "P0DMN0" "P0DPB6" "P42167" "P04908" "P62807" "Q9P0M2"
setdiff(table$Uniprot.ID, pathway_data$Uniprot.ID) # "Q9Y6N8" "O75715" "Q6PUV4"

table[which(table$Uniprot.ID == "Q9Y6N8"),]
all_organ_mean[which(all_organ_mean$Uniprot.ID == "Q9Y6N8"),]

### Same entrezid for multiple uniprot ids
all_organ_mean <- as.data.frame(all_organ_mean)
hgcn2entrez <- clusterProfiler::bitr(all_organ_mean$Uniprot.ID, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
data <- merge(all_organ_mean, hgcn2entrez, by.x="Uniprot.ID", by.y="UNIPROT", all.x = TRUE)
# Remove duplicate that are resulted from multiple matches
data2 <- data[!duplicated(data$Uniprot.ID),]
# Remove gene that are not mapped
data.covid <- data2 %>% tidyr::drop_na(ENTREZID) # 12588 genes
dup_entrez <- data.covid$ENTREZID[which(duplicated(data.covid$ENTREZID))] # 13 unique entrezid are duplicated
for (i in 1: length(dup_entrez)){
  print (data.covid[which(data.covid$ENTREZID == dup_entrez[i]),])
}

## Same uniprot to multiple ensemble 
dup_uniprot <- data$Uniprot.ID[which(duplicated(data$Uniprot.ID))] #692 unique uniprot ids are duplicated
i=1
i=1+i
print(data[which(data$Uniprot.ID == dup_uniprot[i]),])
