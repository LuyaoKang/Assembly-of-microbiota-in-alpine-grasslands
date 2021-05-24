# distance-based redundancy analysis (dbRDA) and plot with color gradient for AI
beta_dbrda <- function(phylo, env_df, geo_df) {
  p_list=c("ggplot2", "vegan", "ggrepel")
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    suppressWarnings(suppressMessages(library(p,character.only=T)))}
  otu_tab <- as.matrix(t(otu_table(phylo)))
  otu_tab <- decostand(otu_tab, 'hellinger')
  env_df <- as.matrix(env_df)
  env <- decostand(env_df,"standardize")
  geo <- as.matrix(geo_df)
  all.vars <- data.frame(env, geo)
  ord <- dbrda(otu_tab ~., all.vars)
  mode <- anova(ord, by = 'margin')
  eig <- ord$CCA$eig
  vare_spp_sco <- scores(ord, display = "species")
  vare_sam_sco <- scores(ord, display = "sites")
  vare_env_sco <- scores(ord, display = "bp")
  vare_env_sco <- vare_env_sco[which(mode$`Pr(>F)` < 0.05),]
  points <- data.frame(vare_sam_sco,data.frame(sample_data(phylo))$AI)
  colnames(points)=c("x", "y", 'AI')
  fit_coord_cont <- data.frame(vare_env_sco, variables = rownames(vare_env_sco))
  ggplot(data = points, aes(x = x, y = y)) + 
    geom_point(data = points, aes(colour = AI), size = 1, alpha = 0.8) + 
    scale_colour_gradientn(colours = terrain.colors(10))+
    geom_segment(data = fit_coord_cont,
                 aes(x = 0, xend = dbRDA1*1.8, y = 0, yend = dbRDA2*1.8),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "grey", lwd = 0.5) +
    geom_text(data = fit_coord_cont, aes(x = dbRDA1*1.8, y = dbRDA2*1.8, label = variables), colour = "black", 
              size = 1) +
    labs(x=paste("dbRDA1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("dbRDA2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
    theme(axis.title = element_text(size = 9, colour = "black"),
          axis.text = element_text(size = 7, colour = "black"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          panel.grid = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black"))
}

p1 <- beta_dbrda(bac.phylo, bacteria_env_df, bacteria_geo_df)
p2 <- beta_dbrda(fungi.phylo, fungi_env_df, fungi_geo_df)
p3 <- beta_dbrda(protist.phylo, protists_env_df, protists_geo_df)
p4 <- beta_dbrda(animal.phylo, animal_env_df, animal_geo_df)

library(cowplot)
plot_grid(p1, p2, p3, p4,
          labels = c("(a)", "(b)", "(c)", "(d)"), ncol = 2, nrow = 2, label_x = .01,
          label_y = 1.005,hjust = 0, label_size=8,align = "v")
