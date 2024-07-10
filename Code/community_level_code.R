
library(anisoFun)
require(mvtnorm)
scale_to <- function(x, val) x * val/mean(x)


# Read data

# matrix with beta coefficients
filenames.beta <- list.files("Output", pattern="beta_matrix_*", full.names=TRUE)
beta.mat.list <- lapply(filenames.beta, read.csv)
for (i in seq_along(beta.mat.list)){
  beta.mat.list[[i]] <- beta.mat.list[[i]] %>% column_to_rownames("X") 
}

# plant attributes
filenames.attr <- list.files("Output", pattern="plant_attr_*", full.names=TRUE)
attr.list <- lapply(filenames.attr, read.csv)
for (i in seq_along(attr.list)){
  attr.list[[i]] <- attr.list[[i]] %>% column_to_rownames("X") %>%
    dplyr::select(Plant_id, Plant_sp, N_flowers, Flowering_min) %>%
    dplyr::mutate(Prop_N_flowers=N_flowers/sum(N_flowers))
}

attr <- do.call("rbind", attr.list)

# matrix with total number of visits per flower across the season 
filenames.vis <- list.files("Output", pattern="vis_matrix_*", full.names=TRUE)
vis.mat.list <- lapply(filenames.vis, read.csv)
for (i in seq_along(vis.mat.list)){
  temp <- vis.mat.list[[i]] %>% column_to_rownames("X") 
  vis.mat.list[[i]] <- temp
}


# matrix with proportional visitation data 
prop.vis.mat.list <- list()
for (i in seq_along(vis.mat.list)){
  prop.vis.mat.list[[i]] <- vis.mat.list[[i]] %>% 
    dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0)) %>% 
    dplyr::mutate(across(where(is.numeric), ~ ./sum(.))) 
}
colSums(prop.vis.mat.list[[1]])


# build the matrix of mutualistic benefits (gamma) with proportional visitation data multiplied by beta coefficients and flower number
# (proportion of visits by each pollinator species to each plant individual)
# I need to use this to parameterize the gamma submatrix of A matrix
# for plant species with many-seed fruits
int.mat <- prop.vis.mat.list[[1]] * attr.list[[1]]$N_flowers * beta.mat.list[[1]] 
colSums(int.mat)
sum(int.mat)

int.mat.list <- list()
for (i in seq_along(beta.mat.list)){
  temp <- beta.mat.list[[i]] * prop.vis.mat.list[[i]] * attr.list[[i]]$N_flowers
  int.mat.list[[i]] <- temp %>% merge(select(attr.list[[i]], Plant_id, Plant_sp), by="row.names") %>% 
    column_to_rownames("Row.names") 
}

int.mat.list[[1]]

# I need the plant id and the plant species as columns

int.mat <- do.call("rbind.fill", int.mat.list) %>% replace(is.na(.), 0)

int.mat.clean <- int.mat  %>% 
  dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
  column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)

# for plant species with four-seed fruits assuming linear response
#int.f.mat <- as.data.frame((1-exp(t(beta[colnames(prop.vis.mat)] * t(prop.vis.mat))))*4)


######
## check to calculate total number of visits along the season (including number of flowers)
# Multiply visitation data per flower along the season, included in the vis.mat matrix, by number of flowers
total.vis.mat.list <- list()
for (i in seq_along(vis.mat.list)){
  total.vis.mat.list[[i]] <- vis.mat.list[[i]] * attr.list[[i]]$N_flowers
}

# calculate proportion of visits per plant species
sum.visits <- list()
for (i in seq_along(total.vis.mat.list)){
  sum.visits[[i]] <- sum(total.vis.mat.list[[i]])
}

visits.per.sp <- do.call("rbind", sum.visits)
colnames(visits.per.sp)[1] <- "visits" 
visits.per.sp %<>% as.data.frame %>% mutate(prop.visits= visits/sum(visits))


# include total visits or visits per flower? CHECK
# pollinator abundance
vis.plant.id <- list()
for (j in seq_along(vis.mat.list)){
  temp <- vis.mat.list[[j]] %>% rownames_to_column("Plant_id")
  vis.plant.id[[j]] <- temp %>% merge(select(attr.list[[j]], Plant_id, Plant_sp), by="Plant_id")
}

abpol <- do.call("rbind.fill", vis.plant.id) %>% replace(is.na(.), 0)

# total abundance of pollinator across all plant species
abpol.total <- abpol %>% select(-Plant_id, -Plant_sp) %>% sum()

# proportion of pollinator species across all plant species
abpol.prop <- (colSums(select(abpol, -Plant_id, -Plant_sp))/sum(select(abpol, -Plant_id, -Plant_sp)))
sum(abpol.prop)



###########

# Generalized-specialized ind mixtures vs specialized individuals


sp <- focal
rand.output[[40]] # random sample of 40 plant individuals 
gen.output[[40]] # 40 most specialized plant individuals

glimpse(int.mat)

Z.list <- list()
PI.list <- list()
PI.a.list <- list()
A.list <- list()
pol.total.abun <- list()
vis.per.sp <- list()
flowers.per.sp <- list()

for (k in seq_along(rand.output[[40]])){
  y <- int.mat %>% filter(Plant_sp !=sp | Plant_id %in% rownames(rand.output[[40]][[k]]))
  
  Z.plant.list <- list()
  for (i in unique(y$Plant_sp)){
    temp <- filter(y, Plant_sp== i)
    temp2 <- matrix(1, nrow(temp), 1)
    rownames(temp2) <- rownames(temp)
    colnames(temp2) <- i
    Z.plant.list[[i]] <- as.data.frame(temp2)
  }
  
  Z.plant <- do.call("rbind.fill", Z.plant.list)
  
  # this is to collapse only plant individuals within the population
  #Z.pol <- diag(1, ncol(int.mat), ncol(int.mat))
  #rownames(Z.pol) <- colnames(int.mat)
  #colnames(Z.pol) <- colnames(int.mat)
  
  #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
  #Z <- Z.df %>% as.matrix()
  
  ## this is to collapse both plant individuals within the population and all pollinator species within the community
  pol <- y %>% select(-Plant_id, -Plant_sp)
  Z.pol <- as.data.frame(matrix(1, ncol(pol), 1))
  rownames(Z.pol) <- colnames(pol)
  colnames(Z.pol) <- "pol"
  
  Z <- rbind.fill(Z.plant, Z.pol) %>% replace(is.na(.), 0) %>% as.matrix()
  Z.list[[k]] <- Z
  
  ##
  
  # create PI P matrix (matrix with proportion of flowers per plant individual within its population)
  X.plant.list <- list()
  for (i in unique(y$Plant_sp)){
    temp <- y %>% filter(Plant_sp==i)
    temp2 <- matrix(1/nrow(temp), nrow=nrow(temp))
    rownames(temp2) <- rownames(temp)
    colnames(temp2) <- i
    
    X.plant.list[[i]] <- as.data.frame(temp2)
  }
  
  X.plant <- do.call("rbind.fill", X.plant.list)
  
  temp <- y %>% select(-Plant_id, -Plant_sp)
  X.pol <- as.data.frame(diag(1, ncol(temp), ncol(temp)))
  rownames(X.pol) <- colnames(temp)
  colnames(X.pol) <- colnames(temp)
  
  PI <- rbind.fill(X.plant, X.pol) %>% replace(is.na(.), 0) %>% as.matrix()
  colSums(PI)
  
  PI.list[[k]] <- PI
  
  
  # create PI matrix for pollinators (relative abundance proportional to number of visits)
  
  # relative pollinator abundance
  pol.subset <- abpol %>% filter(Plant_sp != sp | Plant_id %in% rownames(rand.output[[40]][[k]])) %>%
    dplyr::select(-Plant_id, -Plant_sp)
  abpol.subset <- colSums(pol.subset)/sum(pol.subset)
  
  # total abundance of pollinators
  pol.total.abun[[k]] <- sum(pol.subset)
  
  # relative proportion of visitations per plant species
  pol.plant.sp <- abpol %>% filter(Plant_sp != sp | Plant_id %in% rownames(rand.output[[40]][[k]])) %>%
    dplyr::select(-Plant_id) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
  vis.per.sp[[k]] <- rowSums(pol.plant.sp[,-1])/sum(pol.plant.sp[,-1])
  vis.per.sp[[k]] <- scale_to(rowSums(pol.plant.sp[,-1]), 1)
  
  # relative proportion of flowers
  fl.plant.sp <- attr %>% filter(Plant_sp != sp | Plant_id %in% rownames(rand.output[[40]][[k]])) %>%
    dplyr::select(Plant_sp, N_flowers) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
  flowers.per.sp[[k]] <- rowSums(fl.plant.sp[,-1])/sum(fl.plant.sp[,-1])
  
  
  pol.prop.abun <- as.data.frame(c(abpol.subset))
  colnames(pol.prop.abun) <- "pol"
  
  PI.a <- rbind.fill(unique(Z.plant), pol.prop.abun) %>% replace(is.na(.), 0) %>% as.matrix()
  colSums(PI.a)
  
  PI.a.list[[k]] <- PI.a
  
  # Interaction matrix A
  
  # for plants as rows
  test <- list()
  for (i in unique(y$Plant_sp)){
    
    con <- y %>% filter(Plant_sp == i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    het <- y %>% filter(Plant_sp != i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    # intraspecific competition submatrix
    intra.plant <- matrix(1, nrow=nrow(con), ncol= nrow(con))
    rownames(intra.plant) <- rownames(con)
    colnames(intra.plant) <- rownames(con)
    
    # interspecific competition submatrix
    inter.plant <- matrix(0.1, nrow=nrow(con), ncol= nrow(het))
    rownames(inter.plant) <- rownames(con)
    colnames(inter.plant) <- rownames(het)
    
    # mutualistic benefits received by plants submatrix
    gamma.plant <- y %>% filter(Plant_sp == i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    test[[i]] <- cbind(intra.plant, inter.plant, -gamma.plant)
    
  }
  
  names(test) <- NULL
  plant.submat <- do.call("rbind", test)
  
  # for pollinators as rows
  
  # interspecific competition
  int.mat.clean <- y  %>% 
    dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  inter.pol <- matrix(0.1, nrow=ncol(int.mat.clean), ncol=ncol(int.mat.clean))
  rownames(inter.pol) <- colnames(int.mat.clean)
  colnames(inter.pol) <- colnames(int.mat.clean)
  
  # mutualistic benefits received by pollinators submatrix
  gamma.pol <- matrix(-0.2, nrow=ncol(select(y, -Plant_id, -Plant_sp)), ncol=nrow(y))
  rownames(gamma.pol) <- colnames(int.mat.clean)
  colnames(gamma.pol) <- rownames(int.mat.clean)
  
  pol.submat <- cbind(gamma.pol, inter.pol)
  
  # A matrix
  
  A <- rbind(plant.submat, pol.submat) %>% as.matrix()
  
  A.list[[k]] <- A
  
}


omega.list <- list()
for (i in seq_along(Z.list)){
  Z <- Z.list[[i]]
  PI <- PI.list[[i]]
  PI.a <- PI.a.list[[i]]
  A <- A.list[[i]]
  abpol.total <- pol.total.abun[[i]]
  
  A.int.sp <- MASS::ginv(Z.list[[i]]) %*% A.list[[i]] %*% PI.list[[i]] %*% PI.a.list[[i]]
  A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*vis.per.sp[[i]] 
  #A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*c(0.22, 0.25, 0.53)
  #A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*flowers.per.sp[[i]]*vis.per.sp[[i]]*100
  
  omega <- data.frame(omega= Omega.abun(A.int.sp, abpol.total))
  omega.list[[i]] <- omega
}


omega.all.ran <- do.call("rbind", omega.list) 

ggplot(omega.all.ran, aes(omega)) + geom_density(colour="grey50", fill="grey50", alpha=0.5) + theme_minimal()


prob.clib.list<- list() 

number_Omega_replicates <- 200
number_boot_replicates <- number_Omega_replicates

for(i in seq_along(A.list)){
  A.int.sp <- MASS::ginv(Z.list[[i]]) %*% A.list[[i]] %*% PI.list[[i]] %*% PI.a.list[[i]]
  #A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*10 # change magnitude of mutualistic benefit
  A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*vis.per.sp[[i]]
  #A.int.sp[1:3, 4] <- c(0.22, 0.25, 0.53)
  #A.int.sp[1:3, 4] <- A.int.sp[1:3, 4]*flowers.per.sp[[i]]*vis.per.sp[[i]]*100
  
  value <- prob_extinction_4_int_matrix(A.int.sp, number_Omega_replicates, number_boot_replicates)
  value.df <- as.data.frame(value) %>% dplyr::select(species, prob_excl_mean) %>% mutate(iteration=i)
  prob.clib.list[[i]] <- value.df
}

ep.rand <- do.call("rbind", prob.clib.list) %>% group_by(species) %>% 
  dplyr::summarise(mean=mean(prob_excl_mean), sd=sd(prob_excl_mean)) %>% 
  dplyr::rename(value=mean) 




########## most specialized

y <- int.mat %>% filter(Plant_sp != sp | Plant_id %in% rownames(gen.output[[40]]))

Z.plant.list <- list()
for (i in unique(y$Plant_sp)){
  temp <- filter(y, Plant_sp== i)
  temp2 <- matrix(1, nrow(temp), 1)
  rownames(temp2) <- rownames(temp)
  colnames(temp2) <- i
  Z.plant.list[[i]] <- as.data.frame(temp2)
}

Z.plant <- do.call("rbind.fill", Z.plant.list)

# this is to collapse only plant individuals within the population
#Z.pol <- diag(1, ncol(int.mat), ncol(int.mat))
#rownames(Z.pol) <- colnames(int.mat)
#colnames(Z.pol) <- colnames(int.mat)

#Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
#Z <- Z.df %>% as.matrix()

## this is to collapse both plant individuals within the population and all pollinator species within the community
pol <- y %>% select(-Plant_id, -Plant_sp)
Z.pol <- as.data.frame(matrix(1, ncol(pol), 1))
rownames(Z.pol) <- colnames(pol)
colnames(Z.pol) <- "pol"

Z <- rbind.fill(Z.plant, Z.pol) %>% replace(is.na(.), 0) %>% as.matrix()

##

# create PI P matrix (matrix with proportion of flowers per plant individual within its population)
X.plant.list <- list()
for (i in unique(y$Plant_sp)){
  temp <- y %>% filter(Plant_sp==i)
  temp2 <- matrix(1/nrow(temp), nrow=nrow(temp))
  rownames(temp2) <- rownames(temp)
  colnames(temp2) <- i
  
  X.plant.list[[i]] <- as.data.frame(temp2)
}

X.plant <- do.call("rbind.fill", X.plant.list)

temp <- y %>% select(-Plant_id, -Plant_sp)
X.pol <- as.data.frame(diag(1, ncol(temp), ncol(temp)))
rownames(X.pol) <- colnames(temp)
colnames(X.pol) <- colnames(temp)

PI <- rbind.fill(X.plant, X.pol) %>% replace(is.na(.), 0) %>% as.matrix()
colSums(PI)


# create PI matrix for pollinators (relative abundance proportional to number of visits)
# pollinator abundance

abpol.subset <- abpol %>% filter(Plant_sp !=sp | Plant_id %in% rownames(rand.output[[40]][[k]])) %>%
  dplyr::select(-Plant_id, -Plant_sp)
abpol.subset <- colSums(abpol.subset)/sum(abpol.subset)

pol.prop.abun <- as.data.frame(c(abpol.subset))
colnames(pol.prop.abun) <- "pol"

PI.a <- rbind.fill(unique(Z.plant), pol.prop.abun) %>% replace(is.na(.), 0) %>% as.matrix()
colSums(PI.a)


# Interaction matrix A

# for plants as rows
test <- list()
for (i in unique(y$Plant_sp)){
  
  con <- y %>% filter(Plant_sp == i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  het <- y %>% filter(Plant_sp != i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  # intraspecific competition submatrix
  intra.plant <- matrix(1, nrow=nrow(con), ncol= nrow(con))
  rownames(intra.plant) <- rownames(con)
  colnames(intra.plant) <- rownames(con)
  
  # interspecific competition submatrix
  inter.plant <- matrix(0.1, nrow=nrow(con), ncol= nrow(het))
  rownames(inter.plant) <- rownames(con)
  colnames(inter.plant) <- rownames(het)
  
  # mutualistic benefits received by plants submatrix
  gamma.plant <- y %>% filter(Plant_sp == i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  test[[i]] <- cbind(intra.plant, inter.plant, -gamma.plant)
  
}

names(test) <- NULL
plant.submat <- do.call("rbind", test)

# for pollinators as rows

# interspecific competition
int.mat.clean <- y  %>% 
  dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
  column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)

inter.pol <- matrix(0.1, nrow=ncol(int.mat.clean), ncol=ncol(int.mat.clean))
rownames(inter.pol) <- colnames(int.mat.clean)
colnames(inter.pol) <- colnames(int.mat.clean)

# mutualistic benefits received by pollinators submatrix
gamma.pol <- matrix(-0.2, nrow=ncol(select(y, -Plant_id, -Plant_sp)), ncol=nrow(y))
rownames(gamma.pol) <- colnames(int.mat.clean)
colnames(gamma.pol) <- rownames(int.mat.clean)

pol.submat <- cbind(gamma.pol, inter.pol)

# A matrix

A <- rbind(plant.submat, pol.submat) %>% as.matrix()


MASS::ginv(Z) %*% A %*% PI %*% PI.a
Omega.abun((MASS::ginv(Z) %*% A %*% PI %*% PI.a), abpol.total)

## shape-based metrics
A_int <- MASS::ginv(Z) %*% A %*% PI %*% PI.a
#A_int <- (-1) * diag(c(3,5,7))
number_Omega_replicates <- 200
number_boot_replicates <- number_Omega_replicates

# relative proportion of visitations per plant species
pol.subset <- abpol %>% filter(Plant_sp != sp | Plant_id %in% rownames(gen.output[[40]])) %>%
  dplyr::select(-Plant_id) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
abpol.subset <- rowSums(pol.subset[,-1])/sum(pol.subset[,-1])
abpol.subset <- scale_to(rowSums(pol.subset[,-1]), 1)

# relative proportion of flowers
fl.subset <- attr %>% filter(Plant_sp != sp | Plant_id %in% rownames(gen.output[[40]])) %>%
  dplyr::select(Plant_sp, N_flowers) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
flowers.subset <- rowSums(fl.subset[,-1])/sum(fl.subset[,-1])


A_int[1:3, 4] <- A_int[1:3, 4]*abpol.subset
#A_int[1:3, 4] <- A_int[1:3, 4]*c(0.22, 0.25, 0.53)
#A_int[1:3, 4] <- A_int[1:3, 4]*flowers.subset*abpol.subset*100



ep.spec <- prob_extinction_4_int_matrix(A_int, number_Omega_replicates, number_boot_replicates)

ep.spec %<>% select(species, prob_excl_mean) %>% dplyr::rename(value=prob_excl_mean) %>% 
  mutate(sd=0, scenario="specialized")

ep.rand <- ep.rand %>% dplyr::mutate(scenario="mixture")

ep.all <- rbind(ep.rand, ep.spec)

df <- ep.all %>% mutate(max=value+sd, min=value-sd)

if (sp=="CLIB") {
  res.clib <- df %>%
    group_by(species) %>%
    summarise(value_diff = value[scenario == "mixture"] - value[scenario == "specialized"],
              max_diff = max[scenario == "mixture"] - max[scenario == "specialized"],
              min_diff = min[scenario == "mixture"] - min[scenario == "specialized"]) %>% 
    mutate(species= c("CLIB", "HCOM", "HHAL", "pol"), scenario="mix_CLIB")
}

if (sp=="HCOM") {
  res.hcom <- df %>%
    group_by(species) %>%
    summarise(value_diff = value[scenario == "mixture"] - value[scenario == "specialized"],
              max_diff = max[scenario == "mixture"] - max[scenario == "specialized"],
              min_diff = min[scenario == "mixture"] - min[scenario == "specialized"]) %>% 
    mutate(species= c("CLIB", "HCOM", "HHAL", "pol"), scenario="mix_HCOM")
}

if (sp=="HHAL") {
  res.hhal <- df %>%
    group_by(species) %>%
    summarise(value_diff = value[scenario == "mixture"] - value[scenario == "specialized"],
              max_diff = max[scenario == "mixture"] - max[scenario == "specialized"],
              min_diff = min[scenario == "mixture"] - min[scenario == "specialized"]) %>% 
    mutate(species= c("CLIB", "HCOM", "HHAL", "pol"), scenario="mix_HHAL")
}


#ggplot2::ggsave(filename = "comm_effects_shape.pdf", 
#                plot = comm.eff, 
#                device = cairo_pdf, 
#                dpi = 1200, 
#                width = 22,
#                height = 10, 
#                units = "cm")




