library('ggplot2')


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


covid_pathway <- data.frame(pathway_median$uP,stack(pathway_median[,2:8]))
ncovid_pathway <- data.frame(NC_pathway_median$uP, stack(NC_pathway_median[,2:8]))

new_data <- covid_pathway[,2:3]
new_data$group <- 'covid'    

nc_new_data <- ncovid_pathway[,2:3]
nc_new_data$group <-'noncovid'
violin_data <- rbind(new_data, nc_new_data)

levels = c("thyroid", "lung", "heart", "liver", "spleen","kidney","testi")
violin_data$ind<- factor(violin_data$ind, levels)

ggplot(violin_data, aes(ind,values, fill = group)) + geom_split_violin() + 
  labs(title = "Tissue-specific pathway activities in COVID-19 patients", x = "Tissue", y = "Median")+ 
  theme_bw(base_size = 30) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



plots <- VlnPlot(object = violin_data, features = violin_data$ind, split.by = "group", group.by = "active.ident", pt.size = 0, combine = FALSE, log=T)

