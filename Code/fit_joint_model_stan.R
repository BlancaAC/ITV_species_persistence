
# used to fit the model with stan
library(cmdstanr)

# compile the stan code
mutualistic_benefits <- cmdstanr::cmdstan_model(
	"Code/joint_model.stan"
)

# estimate the parameters using variational Bayes
fit_vb <- mutualistic_benefits$variational(
	data = stan_data,
	init=list(list(
		taus = rep(0, ncol(X_scaled))
	)),
	iter=1E6,
	algorithm='fullrank'
)



# summary of beta coefficients with scaled predictors
beta_summary <-
	fit_vb$draws() %>%
	subset_draws(variable='betas') %>%
	summary()

# posterior of beta coefficients with scaled predictors
p_scaled_betas <-
	fit_vb$draws() %>%
	subset_draws(variable='betas') %>%
	mcmc_intervals() +
	scale_y_discrete(
		limits=rev,
		labels=rev(pollinators)
	) +
	scale_x_log10()
print(p_scaled_betas)

# convert beta coefficients with scaled predictors to their equivalent in unscaled predictor space
seeds_per_visit_to_flower <-
	fit_vb$draws() %>%
	subset_draws(variable='betas') %>%
	sweep(2,attr(X_scaled, "scaled:scale"),"/")

# posterior of seeds per visit to flower
p_seeds_per_visit <-
	seeds_per_visit_to_flower %>%
	mcmc_intervals() +
	scale_y_discrete(
		limits=rev,
		labels=rev(pollinators)
	) +
	scale_x_log10()
print(p_seeds_per_visit)

# posterior distributions (beta)
post.distrib <- ggplot(gather(as.data.frame(seeds_per_visit_to_flower)), aes(value)) + 
  geom_histogram(bins = 20, alpha=0.5) + 
  facet_wrap(~key, scales = 'free_x') +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2), 
                     breaks = equal_breaks(n = 3)) + 
  theme_bw() +
  theme(strip.background = element_rect(
    color="white", fill="white", size=1.5, linetype="solid"), 
    strip.text= element_text(face="bold"), text=element_text(size=10),
    axis.title = element_text(size=15),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  xlab("Seeds per visit and flower") + ylab("Count")

print(post.distrib)

colour.sp <- ifelse(focal=="CLIB", "dodgerblue3", 
                    ifelse(focal=="HCOM", "burlywood4", "pink2"))

log.contr.plot <- ggplot(gather(as.data.frame(seeds_per_visit_to_flower)), aes(x = value, y = key)) +
  stat_pointinterval(interval_alpha=0.5, color=colour.sp) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  xlab('Number of seeds per visit and per flower')

print(log.contr.plot)

## 

colnames(seeds_per_visit_to_flower) <- pollinators

post.df <- data.frame(coef=as.character(), median=as.numeric())
for (i in 1:ncol(seeds_per_visit_to_flower)){
  post.df[i,]$coef <- colnames(seeds_per_visit_to_flower)[i]
  temp2 <- median(seeds_per_visit_to_flower[,i])
  post.df[i,]$median <- temp2
}

int <- dd %>% select(any_of(c(pollinators, "Plant_id", "Plant_sp", "N_flowers", 
                              "Flowering_min", "degree"))) %>% unique() %>% as.data.frame()
rownames(int) <- int$Plant_id
write.csv(int, paste0("Output/plant_attr_", focal, ".csv"))

# visitation matrix
# estimate total number of interactions per flower received across the season by each pollinator species 
# (visitation rate per flower per min * minutes with flowers)
vis.mat <- (int[,pollinators] * int$Flowering_min)
write.csv(vis.mat, paste0("Output/vis_matrix_", focal, ".csv"))

# beta coeff matrix
# only variation in n seeds per visit per flower (create a matrix plant ind x pol sp where elements depict the contribution of pollinator species (beta), only if they interact) 
binary <- int[,pollinators]
binary[binary != 0] <- 1
beta.mat <- as.data.frame(t(t(binary) * post.df$median))

write.csv(beta.mat, paste0("Output/beta_matrix_", focal, ".csv"))


# Figure 2

# estimate the total n seeds per visit per flower only based in pollinator identity for each plant individual
#beta.mat.total <- beta.mat %>%
#  dplyr::mutate(N_seeds_per_visit_flower=rowSums(across(pollinators)))

beta.mat.total <- beta.mat %>%
  mutate(N_seeds_per_visit_flower = rowSums(across(pollinators)) / rowSums(across(pollinators) != 0))

# N seeds per flower matrix
# total number of visitations per flower across the season by each pollinator species * beta (its contribution to n seeds per visit to a flower)
seeds.fl.mat <- as.data.frame(t(t(vis.mat) * post.df$median))

# estimate the total number of seeds per flower produced based on pollinator visitation from all species
seeds.fl.mat.total <- seeds.fl.mat %>%
  dplyr::mutate(N_seeds_per_flower=rowSums(across(pollinators)))

# another approach for panel A (same beta different visitation rate)
seeds.fl.mat.total.avbeta <- as.data.frame(vis.mat * mean(post.df$median)) %>%
  dplyr::mutate(N_seeds_per_flower_avbeta=rowSums(across(pollinators)))

# another approach for panel A (different beta same visitation rate)
seeds.fl.mat.total.avvis <- as.data.frame(t(t(beta.mat) * colMeans(as.matrix(vis.mat)))) %>%
  dplyr::mutate(N_seeds_per_flower_avvis=rowSums(across(pollinators)))

seeds.fl.mat.total.avvis <- as.data.frame(beta.mat * mean(as.matrix(vis.mat))) %>%
  dplyr::mutate(N_seeds_per_flower_avvis=rowSums(across(pollinators)))


# N seeds per plant individual matrix
# Overall contribution of each pollinator species to plant individuals' seed production
# that is, total n visits per flower along the season * beta (n seeds produced per visit per flower) * n flowers along the season
seeds.mat <- as.data.frame(seeds.fl.mat)

# estimate the total number of seeds per plant  based on pollinator visitation from all species and the number of flowers produced
seeds.mat.total <- seeds.mat %>%
  dplyr::mutate(N_seeds=rowSums(across(pollinators))* int$N_flowers)

# incoporating variation in pollinator contribution, number of visits and number of flowers

total.ben <- cbind(select(beta.mat.total, N_seeds_per_visit_flower), 
                   select(seeds.fl.mat.total, N_seeds_per_flower),
                   select(seeds.mat.total, N_seeds),
                   select(seeds.fl.mat.total.avbeta, N_seeds_per_flower_avbeta),
                   select(seeds.fl.mat.total.avvis, N_seeds_per_flower_avvis))

if (focal=="CLIB"){
plot.plant <- ggplot(total.ben, aes(x=N_seeds)) + 
  geom_histogram(binwidth=150, alpha=0.5, color="#66C2A5", fill="#66C2A5") + 
  geom_density(aes(y=150 * ..count..), alpha=0.6, color="#66C2A5", fill="#66C2A5") + 
  theme_bw() + 
  theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per plant individual)")
plot.fl <- ggplot(total.ben, aes(x=N_seeds_per_flower)) + 
  geom_histogram(binwidth=0.5, alpha=0.5, color="#FC8D62", fill="#FC8D62") + 
  geom_density(aes(y=0.5 * ..count..), alpha=0.6, color="#FC8D62", fill="#FC8D62") + 
  theme_bw() + 
  theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per flower)")
plot.vis.fl <- ggplot(total.ben, aes(x=N_seeds_per_visit_flower)) + 
  geom_histogram(binwidth=0.002, alpha=0.5, color="#8DA0CB", fill="#8DA0CB") + 
  geom_density(aes(y=0.002 * ..count..), alpha=0.6, color="#8DA0CB", fill="#8DA0CB") + 
  theme_bw() + 
  theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per visit)")
}

if (focal=="HCOM"){
  plot.plant <- ggplot(total.ben, aes(x=N_seeds)) + 
    geom_histogram(binwidth=15, alpha=0.5, color="#66C2A5", fill="#66C2A5") + 
    geom_density(aes(y=15 * ..count..), alpha=0.6, color="#66C2A5", fill="#66C2A5") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per plant individual)")
  plot.fl <- ggplot(total.ben, aes(x=N_seeds_per_flower)) + 
    geom_histogram(binwidth=0.04, alpha=0.5, color="#FC8D62", fill="#FC8D62") + 
    geom_density(aes(y=0.04 * ..count..), alpha=0.6, color="#FC8D62", fill="#FC8D62") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per flower)")
  plot.vis.fl <- ggplot(total.ben, aes(x=N_seeds_per_visit_flower)) + 
    geom_histogram(binwidth=0.00005, alpha=0.5, color="#8DA0CB", fill="#8DA0CB") + 
    geom_density(aes(y=0.00005 * ..count..), alpha=0.6, color="#8DA0CB", fill="#8DA0CB") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per visit)")
}

if (focal=="HHAL"){
  plot.plant <- ggplot(total.ben, aes(x=N_seeds)) + 
    geom_histogram(binwidth=250, alpha=0.5, color="#66C2A5", fill="#66C2A5") + 
    geom_density(aes(y=250 * ..count..), alpha=0.6, color="#66C2A5", fill="#66C2A5") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per plant individual)") + 
    scale_x_continuous(breaks = c(0, 4000, 8000, 12000, 16000))
  plot.fl <- ggplot(total.ben, aes(x=N_seeds_per_flower)) + 
    geom_histogram(binwidth=0.6, alpha=0.5, color="#FC8D62", fill="#FC8D62") + 
    geom_density(aes(y=0.6 * ..count..), alpha=0.6, color="#FC8D62", fill="#FC8D62") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per flower)")
  plot.vis.fl <- ggplot(total.ben, aes(x=N_seeds_per_visit_flower)) + 
    geom_histogram(binwidth=0.01, alpha=0.5, color="#8DA0CB", fill="#8DA0CB") + 
    geom_density(aes(y=0.01 * ..count..), alpha=0.6, color="#8DA0CB", fill="#8DA0CB") + 
    theme_bw() + 
    theme(text=element_text(size=15)) + ylab("Density") + xlab("Number of seeds (per visit per flower)")
}

plot.vis.fl /plot.fl / plot.plant + plot_annotation(tag_levels = 'A')


# Create dataset with visitation and benefits data
fitness.con <- cbind(vis.mat, total.ben, select(int, "Plant_id", "N_flowers", "Flowering_min", "degree")) %>% 
  dplyr::mutate(sum_pol= rowSums(across(pollinators)))

# visualize results with PCA
x <- select(fitness.con, "N_flowers", "sum_pol", "N_seeds_per_visit_flower", "degree")
res.pca <- prcomp(x, scale = TRUE)
fviz_eig(res.pca)
var <- summary(res.pca)

# Extract PC axes for plotting
PCAvalues <- data.frame(PC1= (res.pca$x[,1]), PC2= (res.pca$x[,2]), size=fitness.con$N_seeds)

# Extract loadings of the variables
PCAloadings <- data.frame(Variables = rownames((res.pca$rotation)), (res.pca$rotation))

# PCA plot
if (focal!="HCOM"){
pca.plot <-ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*5),
                                       yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "grey40") +
  geom_point(aes(size=size), color="grey40", alpha=0.7) +
  annotate("text", x = (PCAloadings$PC1*4.5), y = (PCAloadings$PC2*4.5), size=4,
           label = c("Flower \nproduction", "Visitation \nrate", "Pollinator \ncontribution per visit", "Plant \ndegree")) + 
  theme_bw() + xlab(paste0("PC1 ", "(", round(var$importance[2,1]*100, 2), "%", ")")) + 
  ylab(paste0("PC2 ", "(", round(var$importance[2,2]*100,2), "%", ")")) + 
  guides(size=guide_legend(title="Number of seeds \n(per plant individual)")) +
  scale_x_reverse() +
  #scale_y_reverse() +
  theme(text=element_text(size=15))
}

if (focal=="HCOM"){
  pca.plot <-ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
    geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*5),
                                         yend = (PC2*5)), arrow = arrow(length = unit(1/2, "picas")),
                 color = "grey40") +
    geom_point(aes(size=size), color="grey40", alpha=0.7) +
    annotate("text", x = (PCAloadings$PC1*4.5), y = (PCAloadings$PC2*4.5), size=4,
             label = c("Flower \nproduction", "Visitation \nrate", "Pollinator \ncontribution per visit", "Plant \ndegree")) + 
    theme_bw() + xlab(paste0("PC1 ", "(", round(var$importance[2,1]*100, 2), "%", ")")) + 
    ylab(paste0("PC2 ", "(", round(var$importance[2,2]*100,2), "%", ")")) + 
    guides(size=guide_legend(title="Number of seeds \n(per plant individual)")) +
    #scale_x_reverse() +
    #scale_y_reverse() +
    theme(text=element_text(size=15))
}


# Figure 2
fitness.plot <- ((plot.vis.fl / plot.fl / plot.plant) | pca.plot) + plot_annotation(tag_levels = 'A')
print(fitness.plot)
