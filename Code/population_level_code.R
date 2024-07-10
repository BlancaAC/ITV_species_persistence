

sp <- focal

beta.mat <- read.csv(paste0("Output/beta_matrix_", sp, ".csv")) %>% column_to_rownames("X")
vis.mat <- read.csv(paste0("Output/vis_matrix_", sp, ".csv")) %>% column_to_rownames("X") 
attr <- read.csv(paste0("Output/plant_attr_", sp, ".csv")) %>% column_to_rownames("X") 
prop.vis.mat <- vis.mat %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0)) %>% 
  dplyr::mutate(across(where(is.numeric), ~ ./sum(.))) 
colSums(prop.vis.mat)
int.mat0 <- beta.mat * prop.vis.mat

# relative pollinator abundance
abpol.prop <- (colSums(vis.mat)/sum(vis.mat))
abpol.total <- sum(vis.mat)

# flower abundance 
flower.prop <- attr %>% select(N_flowers) %>% mutate(prop_flowers=N_flowers/sum(N_flowers))
sum(flower.prop$prop_flowers)

int.mat<- int.mat0 * flower.prop$N_flowers


# How does community feasibility (considering only one plant species) change when increasing randomly the number of plant individuals?

rand.output <- list()
for(i in 2:nrow(int.mat)){
  rand.output[[i]] <- lapply(1:100, function(x) 
    int.mat[sample(nrow(int.mat), size = i, replace = F),]) 
  for (j in seq_along(rand.output[[i]])) {
    #rand.output[[i]][[j]] %<>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
    rand.output[[i]][[j]]$Plant_id <- NULL
    rand.output[[i]][[j]] <- rand.output[[i]][[j]][order(row.names(rand.output[[i]][[j]])), ] 
    rand.output[[i]][[j]] <- rand.output[[i]][[j]][, order(colnames(rand.output[[i]][[j]]))]
  }
}
rand.output[sapply(rand.output, is.null)] <- NULL # remove elements that only contain 1-4 plant individuals

# mapping of within-species groups to populations
# within-species group X species matrix
Z.plant.list <- rand.output
Z.df.list <- rand.output
Z.list <- rand.output
for (i in seq_along(rand.output)) {
  for (j in seq_along(rand.output[[i]])){
    Z.plant <- matrix(1, nrow(rand.output[[i]][[j]]), 1)
    rownames(Z.plant) <- rownames(rand.output[[i]][[j]])
    colnames(Z.plant) <- sp
    Z.plant.list[[i]][[j]] <- Z.plant # i will need to use this below
    
    Z.pol <- diag(1, ncol(rand.output[[i]][[j]]), ncol(rand.output[[i]][[j]]))
    rownames(Z.pol) <- colnames(rand.output[[i]][[j]])
    colnames(Z.pol) <- colnames(rand.output[[i]][[j]])
    
    ## to collapse only plant individuals
    #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% 
    #  column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
    #Z.df.list[[i]][[j]] <- Z.df
    #
    #Z <- Z.df %>% as.matrix()
    
    ## to collapse both plant individuals and pollinator species
    Z.df <- cbind(c(rep(1, nrow(rand.output[[i]][[j]])), rep(0, ncol(rand.output[[i]][[j]]))),
                  c(rep(0, nrow(rand.output[[i]][[j]])), rep(1, ncol(rand.output[[i]][[j]]))))
    Z <- Z.df %>% as.matrix()
    
    Z.list[[i]][[j]] <-  Z
  }
}

# partitioning of within-species groups across populations
# within-population group X population matrix

PI.list <- rand.output
for (i in seq_along(Z.df.list)) {
  for (j in seq_along(Z.df.list[[i]])){
    
    # create PI P matrix
    #plant.id <- rownames(Z.df.list[[i]][[j]])
    #plant.if.fl <- filter(attr, rownames(attr) %in% plant.id)
    X.plant <- matrix(1/nrow(Z.df.list[[i]][[j]]), nrow=nrow(Z.df.list[[i]][[j]]))
    #as.matrix(plant.if.fl$N_flowers/sum(plant.if.fl$N_flowers))
    rownames(X.plant) <- rownames(Z.df.list[[i]][[j]])
    colnames(X.plant) <- sp
    X.pol <- diag(1, ncol(int.mat), ncol(int.mat))
    rownames(X.pol) <- colnames(int.mat)
    colnames(X.pol) <- colnames(int.mat)
    
    # PI matrix for plants (relative abundance proportional to flower number)
    PI <- merge(X.plant, X.pol, by = "row.names", all = TRUE) %>% 
      column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()
    
    # this is to give the same weight to each plant individual within the population
    #PI <- Z.df.list[[i]][[j]] %>% 
    #  mutate(across(c(1), ~ if_else(. ==1, 1/length(Z.plant.list[[i]][[j]]), .))) %>% as.matrix()
    PI.list[[i]][[j]] <- PI
  }
}

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a.list <- rand.output
for (i in seq_along(Z.df.list)) {
  for (j in seq_along(Z.df.list[[i]])){
    PI.a <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))
    PI.a.list[[i]][[j]] <- PI.a
  }}


# mutualistic benefits from pollinators to plant individuals

plant.gamma.list <- rand.output
for (i in seq_along(rand.output)) {
  for (j in seq_along(rand.output[[i]])) {
    temp <- as.matrix(rand.output[[i]][[j]])
    # classical approach (assuming a baseline contribution gamma0 and penalizing by plant degree)
    #gamma0 <- coeff.pol # contribution of each pollinator sp (gamma0)
    #d <- rowSums(temp!=0) # degree of plant individuals
    #delta <- 0.2 # mutualistic trade-off
    #plant.gamma <- temp
    #for(row in rownames(temp)) {
    #for(col in colnames(temp)) {
    #    if(temp[row, col]>0){
    #      # mutualistic benefit received by plant individuals from each pollinator sp
    #      plant.gamma[row, col] <- (gamma0[col])/d^delta 
    #    }
    #}
    #}
    #plant.gamma.list[[i]][[j]] <- plant.gamma
    plant.gamma.list[[i]][[j]] <- temp
  }
}

# interaction matrix
# within-species groups X within-species groups matrix

A.list <- rand.output
for (i in seq_along(PI.list)) {
  for (j in seq_along(PI.list[[i]])){
    # setting within-species groups X within-species groups matrix
    A <- diag(1, nrow(PI.list[[i]][[j]]), nrow(PI.list[[i]][[j]])) 
    rownames(A) <- rownames(PI.list[[i]][[j]])
    colnames(A) <- rownames(PI.list[[i]][[j]])
    
    # plant within-species group competition set to 1 above
    nplant <- nrow(Z.plant.list[[i]][[j]])
    ntotal <- ncol(A)
    npol <- ntotal - nplant
    
    #rmat.plant <- abs(matrix(rnorm(nplant*nplant),nrow=nplant))
    #rmat.plant.pol <- abs(matrix(rnorm(nplant*npol),nrow=npol))
    #rmat.pol <- abs(matrix(rnorm(npol*npol),nrow=npol))
    
    A[1:nplant, 1:nplant] <- 1
    A[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
    A[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
    A[1:nplant, (1+nplant):ntotal] <- -(plant.gamma.list[[i]][[j]]) # mutualistic benefits of pollinators given to plants
    for (row in rownames(A)){
      for (col in colnames(A)){
        if(row==col){
          A[row, col] <- 1 # intraspecific competition set to 1
        }
      }
    }
    
    A.list[[i]][[j]] <- A
  }
}


# Feasibility domain estimation

## example
Z <- Z.list[[1]][[1]]
PI <- PI.list[[1]][[1]]
PI.a <- PI.a.list[[1]][[1]]
A <- A.list[[1]][[1]]
alpha <- (MASS::ginv(Z) %*% A %*% PI %*% PI.a)

S <- nrow(alpha)
Sigma <- solve(t(alpha) %*% alpha)
#abund <- c(rep(Inf, 1), colSums(vis.mat))
abund <- c(rep(Inf, S-1), rep(sum(vis.mat), 1))
##

omega.list <- list()
for (i in seq_along(Z.list)){
  temp <- Z.list
  for (j in seq_along(Z.list[[i]])){
    Z <- Z.list[[i]][[j]]
    PI <- PI.list[[i]][[j]]
    PI.a <- PI.a.list[[i]][[j]]
    A <- A.list[[i]][[j]]
    omega <- data.frame(N_ind= nrow(Z.plant.list[[i]][[j]]), 
                        #N_pol= nrow(Z.list[[i]][[j]]) - nrow(Z.plant.list[[i]][[j]]),
                        omega= Omega.abun(MASS::ginv(Z) %*% A %*% PI %*% PI.a, vis.mat))
    temp[[i]][[j]] <- omega
  }
  omega.list[[i]] <- do.call("rbind", temp[[i]])
}

omega.all.ran <- do.call("rbind", omega.list) 
#%>% mutate(omega=(abs(omega)^(1/N_pol))*sign(omega))

# calculate mean and sd
omega.all <- as.data.frame(
  as.matrix(aggregate(omega ~ N_ind, omega.all.ran, function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T)))))

ggplot(omega.all, aes(x=N_ind, y=omega.mean)) +
  geom_line(aes(x=N_ind, y=omega.mean), color="grey70") +
  geom_ribbon(aes(ymin=omega.mean - omega.sd,
                  ymax=omega.mean + omega.sd), alpha=0.4, fill="grey70") +
  theme_minimal() + ylab("Omega") + xlab("Number of individuals")



# How does community feasibility (considering only one plant species) change when incorporating plant individuals with increasing generalisation in pollinator use (degree)?

x <- specieslevel(vis.mat, index=c("species specificity"), level="lower") %>% rownames_to_column("Plant_id") %>% arrange((-species.specificity.index))

# sorting by degree
gen.output <- list()
int.mat.d <- int.mat %>% mutate(degree= rowSums(.!=0)) %>% 
  rownames_to_column("Plant_id") %>% 
  left_join(x, by="Plant_id") %>% 
  arrange(degree, -species.specificity.index) %>% column_to_rownames("Plant_id")

for(i in seq(from = 2, to = nrow(int.mat.d), by = 1)) { 
  gen.output[[i]] <- int.mat.d %>% select(-degree, -species.specificity.index) %>% dplyr::slice(1:i) 
  #%>%  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))    
  gen.output[[i]] <- gen.output[[i]][order(row.names(gen.output[[i]])), ] 
  gen.output[[i]] <- gen.output[[i]][, order(colnames(gen.output[[i]]))] 
}

gen.output[sapply(gen.output, is.null)] <- NULL # remove elements that only contain 1-4 plant individuals


## still working on this
#k <-1
#for (k in seq_along(unique(int.mat.d$degree))){
#test2 <- list(list())
#for (ind in 1:nrow(filter(int.mat.d, degree==k))){
#  test3 <- list()
#  for (i in 1:5){
#  # randomly sample individuals with degree 1
#  temp <- int.mat.d %>% filter(degree==k) %>% sample_n(ind)
#  test3[[i]] <- temp
#  }
#  test2[[ind]] <- test2
#}
#test[[k]] <- test2
#}
##

# mapping of within-species groups to populations
# within-species group X species matrix

Z.plant.list <- gen.output
Z.df.list <- gen.output
Z.list <- gen.output
for (i in seq_along(gen.output)) {
  Z.plant <- matrix(1, nrow(gen.output[[i]]), 1)
  rownames(Z.plant) <- rownames(gen.output[[i]])
  colnames(Z.plant) <- sp
  Z.plant.list[[i]] <- Z.plant # i will need to use this below
  
  Z.pol <- diag(1, ncol(gen.output[[i]]), ncol(gen.output[[i]]))
  rownames(Z.pol) <- colnames(gen.output[[i]])
  colnames(Z.pol) <- colnames(gen.output[[i]])
  
  # to collapse only plant individuals
  #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% 
  #  column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
  #Z.df.list[[i]] <- Z.df
  #
  #Z <- Z.df %>% as.matrix()
  
  ## to collapse both plant individuals and pollinator species
  Z.df <- cbind(c(rep(1, nrow(gen.output[[i]])), rep(0, ncol(gen.output[[i]]))),
                c(rep(0, nrow(gen.output[[i]])), rep(1, ncol(gen.output[[i]]))))
  Z <- Z.df %>% as.matrix()
  
  Z.list[[i]] <-  Z
}

# partitioning of within-species groups across populations
# within-population group X population matrix

PI.list <- gen.output
for (i in seq_along(Z.df.list)) {
  
  # create PI P matrix
  #plant.id <- rownames(Z.df.list[[i]])
  #plant.if.fl <- filter(attr, rownames(attr) %in% plant.id)
  #X.plant <- as.matrix(plant.if.fl$N_flowers/sum(plant.if.fl$N_flowers))
  #rownames(X.plant) <- rownames(plant.if.fl)
  X.plant <- matrix(1/nrow(Z.df.list[[i]]), nrow=nrow(Z.df.list[[i]]))
  rownames(X.plant) <- rownames(Z.df.list[[i]])
  colnames(X.plant) <- sp
  X.pol <- diag(1, ncol(int.mat), ncol(int.mat))
  rownames(X.pol) <- colnames(int.mat)
  colnames(X.pol) <- colnames(int.mat)
  
  # PI matrix for plants (relative abundance proportional to flower number)
  PI <- merge(X.plant, X.pol, by = "row.names", all = TRUE) %>% 
    column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()
  
  # this is to give the same weight to each plant individual within the population
  #PI <- Z.df.list[[i]] %>% 
  #  mutate(across(c(1), ~ if_else(. ==1, 1/length(Z.plant.list[[i]]), .))) %>% as.matrix()
  
  PI.list[[i]] <- PI
}

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a.list <- gen.output
for (i in seq_along(Z.df.list)) {
  PI.a <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))
  PI.a.list[[i]] <- PI.a
}

# mutualistic benefits from pollinators to plant individuals

plant.gamma.list <- gen.output
for (i in seq_along(gen.output)) {
  temp <- as.matrix(gen.output[[i]])
  plant.gamma.list[[i]] <- temp
}


#visweb(gen.output[[length(gen.output)]], type="nested")
#hist(total.gamma$degree)

# interaction matrix
# within-species groups X within-species groups matrix

A.list <- gen.output
for (i in seq_along(PI.list)) {
  # setting within-species groups X within-species groups matrix
  A <- diag(1, nrow(PI.list[[i]]), nrow(PI.list[[i]])) 
  rownames(A) <- rownames(PI.list[[i]])
  colnames(A) <- rownames(PI.list[[i]])
  
  # plant within-species group competition set to 1 above
  nplant <- nrow(Z.plant.list[[i]])
  ntotal <- ncol(A)
  npol <- ntotal - nplant
  
  rmat.plant <- abs(matrix(rnorm(nplant*nplant),nrow=nplant))
  rmat.plant.pol <- abs(matrix(rnorm(nplant*npol),nrow=npol))
  rmat.pol <- abs(matrix(rnorm(npol*npol),nrow=npol))
  
  A[1:nplant, 1:nplant] <- 1
  A[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
  A[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
  A[1:nplant, (1+nplant):ntotal] <- -(plant.gamma.list[[i]]) # mutualistic benefits of pollinators given to plants
  for (row in rownames(A)){
    for (col in colnames(A)){
      if(row==col){
        A[row, col] <- 1 # intraspecific competition set to 1
      }
    }
  }
  A.list[[i]] <- A
}

# Feasibility domain estimation

omega.list <- list()
for (i in seq_along(Z.list)){
  Z <- Z.list[[i]]
  PI <- PI.list[[i]]
  PI.a <- PI.a.list[[i]]
  A <- A.list[[i]]
  
  value.it <-list()
  for (k in 1:10){
    temp.value <- Omega.abun((MASS::ginv(Z) %*% (A) %*% PI %*% PI.a), vis.mat)
    value.it[[k]] <- temp.value
  }
  obs.value.mean <- mean(unlist(value.it))
  obs.value.sd <- sd(unlist(value.it))
  
  omega <- data.frame(N_ind= nrow(Z.plant.list[[i]]),
                      omega.mean= obs.value.mean,
                      omega.sd = obs.value.sd)
  omega.list[[i]] <- omega
}

omega.all.asc.degree <- do.call("rbind", omega.list) 
#omega.all.asc.d <- do.call("rbind", omega.list) %>% rename(omega.d=omega) %>% mutate(omega.d=(abs(omega.d)^(1/N_pol))*sign(omega.d)) 

# include degree information

k.omega <- omega.all.asc.degree %>% mutate(degree= int.mat.d[2:nrow(int.mat.d),]$degree)

if (focal == "CLIB") {
  k.omega %<>%
    mutate(degree = na_if(degree, lag(degree))) %>% na.omit() %>% 
    filter(degree %in% c(2, 4, 6, 8, 10, 15))
  
pclib <- ggplot() +
  geom_line(data=omega.all, aes(x=N_ind, y=omega.mean), color="grey70") +
  geom_ribbon(data=omega.all, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                  ymax=omega.mean + omega.sd), alpha=0.4, fill="grey70") +
  geom_line(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean), color="dodgerblue3") +
  geom_ribbon(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                  ymax=omega.mean + omega.sd), alpha=0.4, fill="dodgerblue3") +
  geom_point(data=k.omega, aes(x=N_ind, y=omega.mean), color="dodgerblue3", size=2) + 
  geom_text(parse=T, data=k.omega, color="dodgerblue3", size=4,
            aes(x=N_ind, y=omega.mean, label = paste0("  italic(k)==", degree)), hjust = 0.7,  vjust = -1) +
  theme_bw() + ylab("Size of the feasibility domain") + xlab("Number of plant individuals") +
  theme(legend.position = "none", text=element_text(size=16)) + 
  scale_y_log10()

 print(pclib)
}

if (focal == "HCOM") {
  k.omega %<>%
    mutate(degree = na_if(degree, lag(degree))) %>% na.omit() %>% 
    filter(degree %in% c(2, 4, 6, 8))
phcom <- ggplot() +
  geom_line(data=omega.all, aes(x=N_ind, y=omega.mean), color="grey70") +
  geom_ribbon(data=omega.all, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                                  ymax=omega.mean + omega.sd), alpha=0.4, fill="grey70") +
  geom_line(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean), color="burlywood4") +
  geom_ribbon(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                                             ymax=omega.mean + omega.sd), alpha=0.4, fill="burlywood4") +
  geom_point(data=k.omega, aes(x=N_ind, y=omega.mean), color="burlywood4", size=2) + 
  geom_text(parse=T, data=k.omega, color="burlywood4", size=4,
            aes(x=N_ind, y=omega.mean, label = paste0("  italic(k)==", degree)), hjust = 0.8,  vjust = -0.9) +
  theme_bw() + ylab("Size of the feasibility domain") + xlab("Number of plant individuals") +
  theme(legend.position = "none", text=element_text(size=16)) + 
  scale_y_log10()
  
  print(phcom)
}

if (focal == "HHAL") {
  k.omega %<>%
    mutate(degree = na_if(degree, lag(degree))) %>% na.omit() %>% 
    filter(degree %in% c(2, 4, 6, 8, 14))
  
phhal <- ggplot() +
  geom_line(data=omega.all, aes(x=N_ind, y=omega.mean), color="grey70") +
  geom_ribbon(data=omega.all, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                  ymax=omega.mean + omega.sd), alpha=0.4, fill="grey70") +
  geom_line(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean), color="pink2") +
  geom_ribbon(data=omega.all.asc.degree, aes(x=N_ind, y=omega.mean, ymin=omega.mean - omega.sd,
                  ymax=omega.mean + omega.sd), alpha=0.4, fill="pink2") +
  geom_point(data=k.omega, aes(x=N_ind, y=omega.mean), color="pink2", size=2) + 
  geom_text(parse=T, data=k.omega, color="pink2", size=4,
            aes(x=N_ind, y=omega.mean, label = paste0("  italic(k)==", degree)), hjust = 0.9,  vjust = -1) +
  theme_bw() + ylab("Size of the feasibility domain") + xlab("Number of plant individuals") +
  theme(legend.position = "none", text=element_text(size=16)) + 
  scale_y_log10()

  print(phhal)
}

