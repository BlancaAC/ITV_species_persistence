
# select the focal data of interest
dd <- dd.list[[focal]]

# figure out the names of the pollinators that visited this plant
cn <- colnames(dd)
pollinators <- cn[which(!cn %in% c('Plant_id','Plant_sp','Response',
                                   'Plot','N_flowers','Type', 'Flowering_min', 'degree') )]

# aggregate bernoullis into binomials
dd_poi <-
	dd %>%
	filter(Type == 'poisson') %>%
	rowwise() %>%
	dplyr::mutate(Seeds = sum(Response), Flowers = n()) %>%
	select(-Response) %>%
	distinct() %>%
	ungroup()

dd_bin <-
	dd %>%
	filter(Type == 'binomial') %>%
	group_by(Plant_id) %>%
	dplyr::mutate(Successes = sum(Response), Trials = n()) %>%
	select(-Response) %>%
	distinct() %>%
	ungroup()

# combine poisson and binomial observations together again
dd <- bind_rows(dd_poi, dd_bin)

# recreate the model matrix
X <- as.matrix(dd[,pollinators] * dd$Flowering_min)

# scale the X matrix (to help fitting and regularization)
X_scaled <- scale(X, center=FALSE)

# to split responses across response types
poi_responses <- which(dd$Type != 'binomial')
Seeds_poi <- dd$Seeds[poi_responses]
bin_responses <- which(dd$Type == 'binomial')
Success_bin <- dd$Successes[bin_responses]
Trials_bin <- dd$Trials[bin_responses]

# format the data required by the model in stan
stan_data <- list(
	N_POI = length(Seeds_poi),
	N_BIN = length(Success_bin),
	N_POL = ncol(X_scaled),
	SEEDS_POI = Seeds_poi,
	SUCCESS_BIN = Success_bin,
	TRIAL_BIN = Trials_bin,
	X_POI = X_scaled[poi_responses,],
	X_BIN = X_scaled[bin_responses,],
	prior_only = FALSE
)

