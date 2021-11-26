RCy3::cytoscapePing()

# combine WikiPathways and COVID19 Disease Map gene sets
wp2gene <- readPathwayGMT("Data/wikipathways-20210710-gmt-Homo_sapiens.gmt")
pwy.covid <- rWikiPathways::getPathwayIdsByCurationTag("Curation:COVID19")
wp2gene.filt <- wp2gene[ wp2gene$wpid %in% pwy.covid,]

dm2gene <- readPathwayGMT("Data/COVID19_DiseaseMap_June2021.gmt")
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

hgcn2entrez <- clusterProfiler::bitr(pwy2gene$gene, fromType = "ENTREZID",toType = c("UNIPROT","SYMBOL","GENENAME"), OrgDb = org.Hs.eg.db)
dm <- merge(pwy2gene, hgcn2entrez, by.x="gene", by.y="ENTREZID", all.x = TRUE)
colnames(dm)[1] <- 'ENTREZID'
write_csv(dm,'09_08_2021_dm_wp_network')

# manually load the table in cytoscape and generate chart visualization
