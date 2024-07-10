
Omega <- function(alpha){
  S <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  m <- matrix(0,S,1)
  a <- matrix(0,S,1)
  b <- matrix(Inf,S,1)  
  #X <- (t(chol(Sigma)))
  #C <- ltMatrices(X[lower.tri(X, diag=TRUE)], diag=TRUE)
  #out <- lpmvnorm(lower = rep(0,S), upper = rep(Inf,S), mean = rep(0,S), chol=C, M=50000)
  d <- pmvnorm(lower = rep(0,S), upper = rep(Inf,S), mean = rep(0,S), sigma = Sigma)
  #out <- log10(d[1])
  out <- (d[1])
  # out <- d[1]^(1 / S) # species level
  return(out)
}

# incorporating pollinator abundance
Omega.abun <- function(alpha, vis.mat){
  S <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  #abund <- c(rep(Inf, 1), colSums(vis.mat))
  abund <- c(rep(Inf, S-1), rep(sum(vis.mat), 1))
  d <- pmvnorm(lower = rep(0,S), upper = abund, mean = rep(0,S), sigma = Sigma)
  #out <- log10(d[1])
  out <- (d[1])
  # out <- d[1]^(1 / S) # species level
  return(out)
}


# perform the model fitting process for each focal species
names <- c("CLIB", "HCOM", "HHAL")

for(i in names){
  focal <- i
  
  # analyses at the plant population level
  source('Code/population_level_code.R')
  
  # Community-level effects
  source('Code/community_level_code.R')

}

pclib + phcom + phhal + plot_layout(axis_titles="collect") + plot_annotation(tag_levels = 'A')

res <- rbind(res.clib, res.hcom, res.hhal)


ggplot(filter(res, species!="pol"), aes(scenario, value_diff, 
                                                    colour=species, shape=species, group=species)) +  
  geom_hline(yintercept = 0, linetype="dashed", colour="grey50", size=0.5) +
  geom_point(alpha=1, size=4, stroke=1.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=min_diff, ymax=max_diff), width=.1, position = position_dodge(width = 0.5)) +
  theme_bw() +
  ylab(bquote("\u0394p"^"E")) + 
  theme(text=element_text(size=20), axis.title.x = element_blank(), 
        legend.text = element_text(face="italic", size=13),
        axis.text.x = element_text(face="italic", size=13),
        legend.title = element_blank(),
        legend.spacing.y = unit(3, 'cm')) +
  scale_x_discrete(labels=c("Cistus libatonis", 
                            "Halimium calycinum", 
                            "Halimium halimifolium")) +
  scale_color_manual(values=c("dodgerblue3","burlywood4", "pink2"), name = "Plant species", 
                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
  scale_shape_manual(values=c(16, 15, 17), name = "Plant species", 
                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) 

