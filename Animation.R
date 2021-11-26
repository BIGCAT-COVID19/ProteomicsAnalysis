RCy3::cytoscapePing()

# combine WikiPathways and COVID19 Disease Map gene sets
wp2gene <- readPathwayGMT("data/wikipathways-20210710-gmt-Homo_sapiens.gmt")
pwy.covid <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
wp2gene.filt <- wp2gene[ wp2gene$wpid %in% pwy.covid,]

dm2gene <- readPathwayGMT("data/COVID19_DiseaseMap_June2021.gmt")
pwy2gene <- dplyr::bind_rows(wp2gene.filt, dm2gene)
# edge list pathway > gene
edges <- pwy2gene[,c(3,5)]
colnames(edges) <- c("source","target")

pwys <- unique(pwy2gene[,c(3,1)])
colnames(pwys) <- c("id", "name")
pwys$type <- "pathway"

genes <- as.data.frame(unique(pwy2gene[,c(5)]))
colnames(genes) <- c("id")
genes$type <- "gene"

nodes <- dplyr::bind_rows(genes,pwys)

RCy3::createNetworkFromDataFrames(nodes = nodes, edges=unique(edges), title="overlap network")

# in the following lines, you will see how you can get the gene names for the Entrez Gene identifiers and add them to the node table as a label
# map Entrez Gene to HGNC Gene symbols
RCy3::mapTableColumn("id","Human", "Entrez Gene","HGNC")

table <- RCy3::getTableColumns(columns = c("id","name","HGNC"))
table$label <- ifelse(is.na(table$name), table$HGNC, table$name)
table$label <- ifelse(is.na(table$label), table$id, table$label)

RCy3::loadTableData(table[,c(1,4)], data.key.column = "id", table.key.column = "id")

RCy3::analyzeNetwork(directed=TRUE)
##

RCy3::createVisualStyle("my_style")
RCy3::setVisualStyle("my_style")

RCy3::installApp("CyAnimator")

data <- read.csv(file="data/proteomics.txt", header=TRUE, sep="\t")
data[,3:14] <- log10(data[,3:14] +1)

RCy3::loadTableData(data, data.key.column = "Gene.name", table.key.column = "HGNC")

RCy3::copyVisualStyle("my_style", "my_style_lung_covid")
RCy3::setNodeColorMapping(table.column = "lung_C1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_lung_covid")
RCy3::setVisualStyle("my_style_lung_covid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 

RCy3::copyVisualStyle("my_style", "my_style_lung_noncovid")
RCy3::setNodeColorMapping(table.column = "lung_NC1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_lung_noncovid")
RCy3::setVisualStyle("my_style_lung_noncovid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 

RCy3::copyVisualStyle("my_style", "my_style_liver_covid")
RCy3::setNodeColorMapping(table.column = "liver_C1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_liver_covid")
RCy3::setVisualStyle("my_style_liver_covid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 

RCy3::copyVisualStyle("my_style", "my_style_liver_noncovid")
RCy3::setNodeColorMapping(table.column = "liver_NC1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_liver_noncovid")
RCy3::setVisualStyle("my_style_liver_noncovid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 

RCy3::copyVisualStyle("my_style", "my_style_kidney_covid")
RCy3::setNodeColorMapping(table.column = "kidney_C1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_kidney_covid")
RCy3::setVisualStyle("my_style_kidney_covid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 

RCy3::copyVisualStyle("my_style", "my_style_kidney_noncovid")
RCy3::setNodeColorMapping(table.column = "kidney_NC1", table.column.values = c(0,1), colors = c("#FFFFFF","#440154"), mapping.type = "c", style.name = "my_style_kidney_noncovid")
RCy3::setVisualStyle("my_style_kidney_noncovid")
RCy3::commandsRun('cyanimator capture frame interpolate=60 network=\"current\"') 
RCy3::commandsRun('cyanimator list frames network=\"current\"') 

# TODO: UPDATE OUTPUT DIR!
RCy3::commandsRun('cyanimator record outputDir=\"C:/Users/martina.kutmon/Downloads/test/"') 
RCy3::commandsRun('cyanimator record outputDir=\"/home/nhumpham/Desktop"') 

