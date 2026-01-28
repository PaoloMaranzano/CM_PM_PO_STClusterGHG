
library(tidyr)
library(dplyr)
library(sf)
library(ggplot2)
library(reshape2)
library(eurostat)
library(rnaturalearth)
library(TSclust)
library(units)
library(ClustGeo)
library(gstat)
library(sp)
library(gridExtra)
library(patchwork)

setwd("C:/.....")

load("CM_PM_PO_STClusterGHG_Data_gas22.RData")
load("CM_PM_PO_STClusterGHG_Data_sector22.RData")

####### res_gas22.RData : res_gas22 is a named list typically containing:
# res_gas22[[1]] or res_gas22[["diss"]]: list of dissimilarity matrices for gases plus the spatial matrix, all normalized to [0,1]. Gas time-series dissimilarities are usually DTW; the spatial one is geodetic distance between centroids.
# res_gas22[[2]] or res_gas22[["data"]]: long table with columns NUTS, Unit (gas), and yearly columns 1990…2022 with emissions per km².
# res_gas22[[3]] or res_gas22[["nuts"]]: sf object of NUTS-2 polygons aligned to the data (and later, a “tree” column with cluster labels).
# Optionally "best_alpha", "exp_inertia", "best_tree", "centroids", "joint_inertia" from prior runs contains the results from spatiotemporal hierarchical clustering.

####### res_sector22.RData 
# res_sector22 mirrors res_gas22 but by sector instead of gas:
# res_sector22[["diss"]]: list of sector dissimilarity matrices plus spatial.
# res_sector22[["data"]]: emissions per km² by sector and year.
# res_sector22[["nuts"]]: NUTS-2 geometry aligned to the data.
# Optionally "best_alpha", "exp_inertia", "best_tree", "centroids", "joint_inertia" from prior runs contains the results from spatiotemporal hierarchical clustering.


nuts22<-res_gas22[[3]]
gas22<-res_gas22[[2]]
d_gas22<-res_gas22[[1]]
sector22<-res_sector22[[2]]
d_sector22<-res_sector22[[1]]

gas_name<-unique(gas22$Unit)
sector_name<-unique(sector22$Sector)

##################### Spatial Centroids ##################### 

nuts_shapes <- get_eurostat_geospatial(nuts_level = "all", year = 2021, resolution = "60")
nuts2_shapes <- nuts_shapes %>% filter(LEVL_CODE == 2)
nuts2_filtered <- nuts2_shapes %>% filter(NUTS_ID %in% gas22$NUTS)

# Calculate centroids for NUTS2 regions
centroids_nuts2 <- st_centroid(nuts2_filtered)

# Extract coordinates for NUTS2 centroids
coord2 <- data.frame(NUTS = centroids_nuts2$NUTS_ID,
                     latitude = st_coordinates(centroids_nuts2)[, 2],
                     longitude = st_coordinates(centroids_nuts2)[, 1])

world <- ne_countries(scale = "medium", returnclass = "sf")
nuts2<-nuts2_filtered

####################################################
############## DISSIMILARITY MATRICES ##############
####################################################

gas22 <- gas22[order(gas22$NUTS), ]
sector22 <- sector22[order(sector22$NUTS), ]
coord2 <- coord2[order(coord2$NUTS), ]

# DTW distances across emission time series by gas
d_gas22<-list() #list of dissimilarity matrices by gas
for (gas_i in unique(gas22$Unit)){
  Y<-as.matrix(gas22[gas22$Unit==gas_i,3:35])
  d_gas22[[gas_i]]<-matrix(ncol=nrow(Y),nrow=nrow(Y))
  for(i in 1:nrow(Y)){
    for(j in 1:nrow(Y)){
      d_gas22[[gas_i]][i,j]<-diss.DTWARP(Y[i,],Y[j,])
    }
  }
  cat("Matrice dissimilarità calcolata per",gas_i," \n")
}

# DTW distances across emission time series by sector
d_sector22<-list() #list of dissimilarity matrices by sector
for (sector_i in unique(sector22$Sector)){
  Y<-as.matrix(sector22[sector22$Sector==sector_i,3:35])
  d_sector22[[sector_i]]<-matrix(ncol=nrow(Y),nrow=nrow(Y))
  for(i in 1:nrow(Y)){
    for(j in 1:nrow(Y)){
      d_sector22[[sector_i]][i,j]<-diss.DTWARP(Y[i,],Y[j,])
    }
  }
  cat("Matrice dissimilarità calcolata per",sector_i," \n")
}


# Geodetic distances across Nuts2 regions' centroid
d_sf2 <- st_as_sf(coord2, coords = c("latitude", "longitude")) %>% 
  st_set_crs(4326) # to indicate lat long degree system
d_sp2<-sf::st_distance(d_sf2)

d_gas22[["spatial"]]<-drop_units(d_sp2)
d_sector22[["spatial"]]<-drop_units(d_sp2)

# dissimilarity matrices are normalized with respect to their maximum
for (i in seq_along(d_gas22)){
  d_gas22[[i]]<-d_gas22[[i]]/max(d_gas22[[i]])
}
for (i in seq_along(d_sector22)){
  d_sector22[[i]]<-d_sector22[[i]]/max(d_sector22[[i]])
}


#####################################################################
########## descriptive statistcs: average and perc change ########### 
#####################################################################

yrs <- as.character(1990:2022)

# Build region-year totals by summing across gases (Unit) per NUTS, then derive stats
reg_year_totals <- gas22 %>%
  dplyr::select(NUTS, Unit, dplyr::all_of(yrs)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(yrs),
                      names_to = "year", values_to = "value") %>%
  dplyr::group_by(NUTS, year) %>%
  dplyr::summarise(total = sum(value, na.rm = TRUE), .groups = "drop")

# Annual average over 1990-2022
avg_by_region <- reg_year_totals %>%
  group_by(NUTS) %>%
  summarise(avg_1990_2022 = mean(total, na.rm = TRUE), .groups = "drop")

# Percent change (2022 vs 1990)
chg_by_region <- reg_year_totals %>%
  dplyr::filter(year %in% c("1990", "2022")) %>%
  tidyr::pivot_wider(names_from = year, values_from = total, names_prefix = "y_") %>%
  dplyr::mutate(pct_change_1990_2022 = 100 * (y_2022 - y_1990) / y_1990) %>%
  dplyr::select(NUTS, pct_change_1990_2022)


nuts2_map <- nuts2 %>%
  left_join(avg_by_region,  by = c("NUTS_ID" = "NUTS")) %>%
  left_join(chg_by_region,  by = c("NUTS_ID" = "NUTS"))

# Two maps: (A) annual average 1990–2022
#                 (B) % change 1990 → 2022

p_avg <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2_map, aes(fill = avg_1990_2022), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-14, 35)) +
  scale_fill_gradientn(
    colors = c("#ffFF00", "#8811aa", "#550077", "#440077"),
    name   = "Ton per km²"
  ) +
  theme_minimal() +
  labs(
    title = "Average annual GHG emissions",
    x = "Longitude", y = "Latitude"
  ) +
  theme(legend.position = "right")

# Diverging palette for % change (green = decrease, red = increase)
p_change <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2_map, aes(fill = pct_change_1990_2022), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-14, 35)) +
  scale_fill_gradient2(
    low = "forestgreen", mid = "white", high = "firebrick",
    midpoint = 0, name = "%"
  ) +
  theme_minimal() +
  labs(
    title = "Percent change in GHG emission",
    x = "Longitude", y = "Latitude"
  ) +
  theme(legend.position = "right")

p_avg + p_change


#####################################################################
################# Cluster analysis: alpha and K #####################
#####################################################################

# alpha_comb: give in input delta_alpha and P (number of matrices), we obtain all possible vectors of alpha_p (all possible combinations)
alpha_comb <- function(delta_alpha, num_matrices, current_sum = 0, current_comb = c()) {
  # Base case: if we have selected num_matrices - 1 alpha values, the last one is just 1 - sum of current_comb
  if (length(current_comb) == num_matrices - 1) {
    last_alpha <- 1 - current_sum
    if (last_alpha >= 0 && last_alpha <= 1) {
      return(list(c(current_comb, last_alpha)))
    } else {
      return(list())
    }
  }
  # Recursive case: try all possible values for the next alpha, making sure we don't exceed the sum of 1
  valid_combinations <- list()
  for (alpha in seq(0, 1, by = delta_alpha)) {
    if (current_sum + alpha <= 1) {
      new_comb <- alpha_comb(delta_alpha, num_matrices, current_sum + alpha, c(current_comb, alpha))
      valid_combinations <- c(valid_combinations, new_comb)
    }
  }
  return(valid_combinations)
}  

# sptClu: function for spatiotemporal optimization algorithm
# Input: D_total (list of dissimilarity matrices), delta_alpha, K
# output: alpa (all possible vectors), D_combined (all possible combined diss matrices),
# part (all possible partitions), W (within inertia in each diss matrix in each combination),
# Q (proportion of explained inertia in each dissimilarity matrices, and weighted average proportion explained inertia in the last column)
sptClu<- function(D_total, delta_alpha, K) {
  num_matrices <- length(D_total)
  
  # Generate all valid alpha combinations
  alpha_combinations <- alpha_comb(delta_alpha, num_matrices)
  # Initialize a list to store results
  results <- list()
  # Initialize an array to store within-cluster dissimilarities
  # The dimensions of the array will be based on the number of alpha combinations and number of dissimilarity matrices
  W <- array(0, dim = c(length(alpha_combinations), num_matrices ))  # +1 for the final combined D
  Q<-array(0, dim = c(length(alpha_combinations), num_matrices + 1))
  # Iterate over each combination of alphas
  delta_total<-list() 
  # Compute the linear combination of dissimilarity matrices
  n <- nrow(D_total[[1]])
  for (j in 1:num_matrices) {
    delta_total[[j]] <- as.matrix(wardinit(as.dist(D_total[[j]])))
  }
  for (i in 1:length(alpha_combinations)) {
    # Extract the current combination of alphas
    alpha <- alpha_combinations[[i]]
    # Initialize the combined dissimilarity matrix
    D_combined <- matrix(0, nrow = nrow(delta_total[[1]]), ncol = ncol(delta_total[[1]]))
    # Compute the linear combination of dissimilarity matrices
    for (j in 1:num_matrices) {
      D_combined <- D_combined + alpha[j] * as.matrix(delta_total[[j]])
    }
    # Perform hierarchical clustering on the combined dissimilarity matrix
    tree <- stats::hclust(as.dist(D_combined), method = "ward.D")
    part <- stats::cutree(tree, k = K)  # Cut the tree into K clusters
    # Compute within-cluster dissimilarity for each dissimilarity matrix and the combined matrix
    for (j in 1:num_matrices) {
      W[i, j] <- withindiss(as.dist(D_total[[j]]), part)
    }
    # Store the result
    results[[i]] <- list(alpha = alpha, D_combined = D_combined, part = part)
  }
  t<-list()
  for (j in 1:num_matrices) {
    t[[j]]<-inertdiss(as.dist(D_total[[j]]))
    Q[,j]<-1-W[,j]/t[[j]]
  }
  tot_diss<-0
  for (j in 1:num_matrices) {
    Q[, num_matrices + 1]<-Q[, num_matrices + 1]+ Q[,j]*t[[j]]
    tot_diss<-tot_diss+t[[j]]
  }
  Q[, num_matrices + 1]<-Q[, num_matrices + 1]/(tot_diss)
  
  return(list(results = results, W = W, Q=Q))
}


########## Clustering algorithm by GAS

d_tot<-d_gas22
d_alp<-0.05
k_max<-10
# best alpha vector for each K
best_alpha<-matrix(0,ncol=length(d_tot),nrow=(k_max-1))
colnames(best_alpha)<-names(d_tot)
# partition according to best alpha vector for each K
best_tree<-list()
# proportion of explained inertia in the partition according to best alpha vector for each K
exp_inertia<-matrix(0,ncol=(length(d_tot)+1),nrow=(k_max))
exp_inertia[1,]<-0
colnames(exp_inertia)<-c(names(d_tot),"w_average")
nuts2$tree<-0
for(k in 2:k_max){
  cat("In elaborazione",k,"cluster...")
  spt_clu_risultati<-sptClu(d_tot,d_alp,k)
  max_row_index <- which.max(spt_clu_risultati$Q[, ncol(spt_clu_risultati$Q)])
  best_alpha[k-1,]<-spt_clu_risultati$results[[max_row_index]]$alpha
  best_tree[[k-1]]<-spt_clu_risultati$results[[max_row_index]]$part
  exp_inertia[k,]<-spt_clu_risultati$Q[max_row_index,]
  nuts2$tree<-best_tree[[k-1]]
  cat("Finito! \n")
  cat("Best alpha_i ",best_alpha[k-1,], " \n")
  cat("Inertia ",exp_inertia[k,], " \n")
}

exp_inertia<-res_gas22[["exp_inertia"]]
# plot of First difference in weighted average proportion of Explained inertia
plot(diff(exp_inertia[,ncol(exp_inertia)]), type = "l", col = "blue", 
     xlab = "Clusters", ylab = "", ylim = c(0, 0.35), 
     main = "First difference explained inertia (emissions by gas)", xaxt = "n")
points(diff(exp_inertia[,ncol(exp_inertia)]), col = "blue", pch = 19, cex = 0.7)
abline(v=3,lty=2) 
axis(1, at = 1:(length(diff(exp_inertia[,ncol(exp_inertia)]))), labels = 2:(length(diff(exp_inertia[,ncol(exp_inertia)])) + 1))
# K*=4
# save as diff_inertia_gas 700x400


res_gas22<-list()
res_gas22[["diss"]]<-d_gas22
res_gas22[["data"]]<-gas22
res_gas22[["nuts"]]<-nuts2
res_gas22[["best_alpha"]]<-best_alpha
res_gas22[["exp_inertia"]]<-exp_inertia
res_gas22[["best_tree"]]<-best_tree
save(res_gas22,file="revision/res_gas22_rev.RData")  

## Centroids in each K for emissions by gas
data_gas22<-res_gas22[["data"]]
centroids<-list()
for(k in 2:k_max){
  nuts2$tree<-res_gas22[["best_tree"]][[k-1]]
  gas22<-data_gas22[,1:35]
  nut<-cbind(nuts2$NUTS_ID,nuts2$tree)
  colnames(nut)<-c("NUTS_ID","tree")
  gas22<-merge(gas22,nut,by.x="NUTS",by.y="NUTS_ID",all.x=TRUE)
  centroids[[k-1]]<-list()
  for(gas_i in unique(gas22$Unit)){
    data<-gas22%>%filter(Unit==gas_i)
    #data<-merge(data,nuts[,c(1,14)],by.x=)
    centroids[[k-1]][[gas_i]]<-matrix()
    centroids[[k-1]][[gas_i]]<-aggregate(. ~ tree, data = data[,3:36], mean)
  }
}
res_gas22[["centroids"]]<-centroids
save(res_gas22,file="revision/res_gas22_rev.RData") 

############ MAP and Centroids ############
# input: centroids.RData, nuts1 (or nuts2) with column "tree"

cluster_color_gas <- c("#9579De", "#451070", 
                   "#FF83FA","#d9108B")

centroids<-res_gas22[["centroids"]]
gas22<-res_gas22[["data"]]
best_tree<-res_gas22[[6]]
nuts2<-res_gas22[[3]]

nuts2 <- nuts2[order(nuts2$NUTS_ID), ]
k<-4
# Initialize an empty list to store plots
plots <- list()
plot_map<-list()
# Loop over each gas in centroids[[1]]
for (i in 1:length(gas_name)) {
  # Extract the data matrix for the current gas
  gas_data <- centroids[[k-1]][[gas_name[i]]]
  # Melt the data to convert it from wide to long format
  gas_data_long <- melt(gas_data, id.vars = "tree", variable.name = "Year", value.name = "Observation")
  gas_data_long$Year <- as.numeric(gas_data_long$Year)
  years <- 1990:2022
  tick_years <- seq(1990, 2020, by = 5)
  
  p <- ggplot(gas_data_long, aes(x = Year, y = Observation,color = factor(tree), group = tree)) +
    geom_line(linewidth = 1.3) +
    labs(title = paste("Annual Trend for", gas_name[i]),
         x = "Year", y = "Annual Centroids", color = "Cluster") +
    theme_minimal() +
    scale_color_manual(values = cluster_color_gas) +
    scale_x_continuous(
      breaks = which(years %in% tick_years),     # index positions
      labels = tick_years
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10))
   plots[[i]]<-p
}
nuts2$tree<-best_tree[[k-1]]
plt <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2, aes(fill = factor(tree)), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
  scale_fill_manual(values = cluster_color_gas) +  # use scale_fill_manual for fill
  theme_minimal() +
  labs(title = "Map of spatiotemporal clusters of emissions by gas", x = "Longitude", y = "Latitude", fill = "Clusters")  
print(plt)
# After your plot code
ggsave(
  "map_sp_clu_gas_rev.png",
  plot = plt,
  width = 1000, height = 700, units = "px",
  dpi = 96,                 # match screen DPI so text/lines aren't oversized
  bg = "white",             # avoid transparent background
  device = ragg::agg_png    # optional, for crisp text/lines
)
plot_map[[k-1]]<-plt
# Display all the plots in a grid layout
do.call(grid.arrange, c(plots, ncol = 2))
#save as sp_clu_centroids_gas 900x700



########## Clustering algorithm by SECTOR

d_tot<-d_sector22
d_alp<-0.1
k_max<-10
# best alpha vector for each K
best_alpha<-matrix(0,ncol=length(d_tot),nrow=(k_max-1))
colnames(best_alpha)<-names(d_tot)
# partition according to best alpha vector for each K
best_tree<-list()
# proportion of explained inertia in the partition according to best alpha vector for each K
exp_inertia<-matrix(0,ncol=(length(d_tot)+1),nrow=(k_max))
exp_inertia[1,]<-0
colnames(exp_inertia)<-c(names(d_tot),"w_average")
for(k in 2:k_max){
  cat("In elaborazione",k,"cluster...")
  spt_clu_risultati<-sptClu(d_tot,d_alp,k)
  max_row_index <- which.max(spt_clu_risultati$Q[, ncol(spt_clu_risultati$Q)])
  best_alpha[k-1,]<-spt_clu_risultati$results[[max_row_index]]$alpha
  best_tree[[k-1]]<-spt_clu_risultati$results[[max_row_index]]$part
  exp_inertia[k,]<-spt_clu_risultati$Q[max_row_index,]
  cat("Finito! \n")
  cat("Best alpha_i ",best_alpha[k-1,], " \n")
  cat("Inertia ",exp_inertia[k,], " \n")
}

exp_inertia<-res_sector22[["exp_inertia"]]
# plot of First difference in weighted average proportion of Explained inertia
plot(diff(exp_inertia[,ncol(exp_inertia)]), type = "l", col = "blue", 
     xlab = "Clusters", ylab = "", ylim = c(0, 0.3), 
     main = "First difference explained inertia (emissions by sector)", xaxt = "n")
points(diff(exp_inertia[,ncol(exp_inertia)]), col = "blue", pch = 19, cex = 0.7)
abline(v=4,lty=2) 
axis(1, at = 1:(length(diff(exp_inertia[,ncol(exp_inertia)]))), labels = 2:(length(diff(exp_inertia[,ncol(exp_inertia)])) + 1))
# K*=5
# save as diff_inertia_sector 700x400

res_sector22<-list()
res_sector22[["diss"]]<-d_sector22
res_sector22[["data"]]<-sector22
res_sector22[["nuts"]]<-nuts2
res_sector22[["best_alpha"]]<-best_alpha
res_sector22[["exp_inertia"]]<-exp_inertia
res_sector22[["best_tree"]]<-best_tree
save(res_sector22,file="res_sector22.RData")  


## Centroids in each K for emissions by sector
centroids<-list()
data_sector22<-res_sector22[["data"]]
nuts2<-nuts2[order(nuts2$NUTS_ID),]
nuts2$tree<-0
for(k in 2:k_max){
  sector22<-res_sector22[["data"]][,1:35]
  nuts2$tree<-res_sector22[["best_tree"]][[k-1]]
  nut<-cbind(nuts2$NUTS_ID,nuts2$tree)
  colnames(nut)<-c("NUTS_ID","tree")
  sector22<-merge(sector22,nut,by.x="NUTS",by.y="NUTS_ID",all.x=TRUE)
  centroids[[k-1]]<-list()
  for(sector_i in unique(sector22$Sector)){
    data<-sector22%>%filter(Sector==sector_i)
    centroids[[k-1]][[sector_i]]<-matrix()
    centroids[[k-1]][[sector_i]]<-aggregate(. ~ tree, data = data[,3:36], mean)
  }
}
res_sector22[["centroids"]]<-centroids
save(res_sector22,file="res_sector22.RData") 
save(res_sector22,file="revision/res_sector22_rev.RData")  

############ MAP and Centroids ############
# input: centroids.RData, nuts1 (or nuts2) with column "tree"

cluster_color_sec <- c("#98FB98", "#107030", "#44e5ee", "#0592a0", "#00CD66")

centroids<-res_sector22[["centroids"]]
sector22<-res_sector22[["data"]]
best_tree<-res_sector22[[6]]
nuts2<-res_sector22[[3]]

nuts2 <- nuts2[order(nuts2$NUTS_ID), ]
k<-5
# Initialize an empty list to store plots
plot_map<-list()
plots <- list()
# Loop over each gas in centroids[[1]]
for (i in 1:length(sector_name)) {
  # Extract the data matrix for the current gas
  sector_data <- centroids[[k-1]][[sector_name[i]]]
  # Melt the data to convert it from wide to long format
  sector_data_long <- melt(sector_data, id.vars = "tree", variable.name = "Year", value.name = "Observation")
  sector_data_long$Year <- as.numeric(sector_data_long$Year)
  years <- 1990:2022
  tick_years <- seq(1990, 2020, by = 5)
  
  p <- ggplot(sector_data_long, aes(x = Year, y = Observation,
                                    color = factor(tree), group = tree)) +
    geom_line(linewidth = 1.3) +
    labs(title = paste("Annual Trend for", sector_name[i]),
         x = "Year", y = "Annual Centroids", color = "Cluster") +
    theme_minimal() +
    scale_color_manual(values = cluster_color_sec) +
    scale_x_continuous(
      breaks = which(years %in% tick_years),     # index positions
      labels = tick_years
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10))
  
plots[[i]]<-p
}
nuts2$tree<-best_tree[[k-1]]
plt <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2, aes(fill = factor(tree)), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
  scale_fill_manual(values = cluster_color_sec) +  # use scale_fill_manual for fill
  theme_minimal() +
  labs(title = "Map of spatiotemporal clusters of emissions by sector", x = "Longitude", y = "Latitude", fill = "Clusters")  
print(plt)
# save as map_sp_clu_sector 1000x700
ggsave(
  "map_sp_clu_sector_rev.png",
  plot = plt,
  width = 1000, height = 700, units = "px",
  dpi = 96,                 # match screen DPI so text/lines aren't oversized
  bg = "white",             # avoid transparent background
  device = ragg::agg_png    # optional, for crisp text/lines
)
plot_map[[k-1]]<-plt
# Display all the plots in a grid layout
do.call(grid.arrange, c(plots, ncol = 2))
#save as sp_clu_centroids_sector 900x900



###############################################################
###################### Joint Inertia ##########################
###############################################################


joint_inertia<-function(d_tot,i,alpha_val,K){
  # D_ts = combined matrix with temporal dissimilarities 
  # (or combined matrix without D_i)
  D_ts <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
  Delta_ts <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
  
  # D_sp = combined spatiotemporal diss. matrix (combination of all matrices)
  # according to best_alpha
  Delta_st <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
  D_st <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
  
  delta_tot<-d_tot
  n <- nrow(d_tot[[1]])
  for (j in 1:(length(d_tot))) {
    delta_tot[[j]] <- as.matrix(wardinit(as.dist(d_tot[[j]])))
  }
  
  for (j in 1:(length(delta_tot))) {
    Delta_st <- Delta_st + alpha_val[j] * as.matrix(delta_tot[[j]])
  }
  for (j in 1:(length(d_tot))) {
    D_st <- D_st + alpha_val[j] * as.matrix(d_tot[[j]])
  }
  delta_tot_new<-delta_tot[-i]
  d_tot_new<-d_tot[-i]
  if (length(d_tot)>5){
    d_alp<-0.1
  }else{d_alp<-0.05}
  spt_clu_risultati<-sptClu(d_tot_new,d_alp,K) # spt_clu exluding matrix i
  max_row_index <- which.max(spt_clu_risultati$Q[, ncol(spt_clu_risultati$Q)])
  alpha_new<-spt_clu_risultati$results[[max_row_index]]$alpha # best alpha in the new combination
  for (j in 1:(length(delta_tot_new))) {
    Delta_ts <- Delta_ts + alpha_new[j] * as.matrix(delta_tot_new[[j]])
  }
  for (j in 1:(length(d_tot_new))) {
    D_ts <- D_ts + alpha_new[j] * as.matrix(d_tot_new[[j]])
  }
  # D_sp = spatial dissimilarity matrix (or diss matrix i)
  D_sp<-as.matrix(d_tot[[i]])
  Delta_sp<-as.matrix(delta_tot[[i]])
  # Total inertia for spatial and temporal component
  T_sp<-inertdiss(as.dist(D_sp))
  # Total inertia for temporal component (or when excluding matrix i)
  T_ts<-inertdiss(as.dist(D_ts))
  T_ts_sum<-0
  for (j in 1:(length(d_tot_new))) {
    #D_ts <- D_ts + alpha_new[j] * as.matrix(d_tot_new[[j]])
    T_ts_sum<-T_ts_sum + inertdiss(as.dist(d_tot_new[[j]]))
  }
  # time-series clustering (clustering without matrix i)
  tree_ts <- stats::hclust(as.dist(Delta_ts), method = "ward.D")
  part_ts <- stats::cutree(tree_ts, k = K)  # Cut the tree into K clusters
  # spatial clustering (clustering only with matrix i)
  tree_sp <- stats::hclust(as.dist(Delta_sp), method = "ward.D")
  part_sp <- stats::cutree(tree_sp, k = K)  # Cut the tree into K clusters
  # spatiotemporal clustering 
  tree_st <- stats::hclust(as.dist(Delta_st), method = "ward.D")
  part_st <- stats::cutree(tree_st, k = K)  # Cut the tree into K clusters
  # W_x_y = within cluster inertia of dissimilarity matrix x 
  # considering clustering partition obtained using y
  W_sp_sp<-withindiss(as.dist(D_sp), part_sp) # W(D_sp)(P^sp_k)
  W_sp_st<-withindiss(as.dist(D_sp), part_st) # W(D_sp)(P^a_k)
  W_ts_ts<-withindiss(as.dist(D_ts), part_ts) # W(D(-sp))(P^(a-sp)_k)
  W_ts_ts_sum<-0
  for (j in 1:(length(d_tot_new))) {
    #D_ts <- D_ts + alpha_new[j] * as.matrix(d_tot_new[[j]])
    W_ts_ts_sum<-W_ts_ts_sum+withindiss(as.dist(d_tot_new[[j]]), part_ts)
  }
  W_ts_st<-withindiss(as.dist(D_ts), part_st) # W(D(-sp))(P^a_k)
  W_ts_st_sum<-0
  for (j in 1:(length(d_tot_new))) {
    #D_ts <- D_ts + alpha_new[j] * as.matrix(d_tot_new[[j]])
    W_ts_st_sum<-W_ts_st_sum+withindiss(as.dist(d_tot_new[[j]]), part_st)
  }
  # Q_x_y = proportion of explained inertia from matrix x 
  # considering clustering partition obtained using y, 
  # with respect to the total inertia in matrix x
  Q_sp_sp<-1-W_sp_sp/T_sp 
  Q_sp_st<-1-W_sp_st/T_sp
  Q_ts_ts<-1-W_ts_ts/T_ts
  Q_ts_st<-1-W_ts_st/T_ts
  Q_ts_ts_av<-1-W_ts_ts_sum/T_ts_sum
  Q_ts_st_av<-1-W_ts_st_sum/T_ts_sum
  Qnorm_ts<-Q_ts_st/Q_ts_ts
  Qnorm_ts_av<-Q_ts_st_av/Q_ts_ts_av
  Qnorm_sp<-Q_sp_st/Q_sp_sp
  # resulting joint inertia
  J_in_old<-Qnorm_sp+Qnorm_ts-1
  J_in<-Qnorm_sp+Qnorm_ts_av-1
  return(c(J_in, Qnorm_sp,Qnorm_ts_av, Q_sp_st, alpha_new))
}

## Joint inertia in ghg emission by GAS clustering
k<-4
d_alp<-0.05
j_inertia<-list()
j_inertia[["proportion"]]<-matrix(0,nrow=length(res_gas22[[1]]),ncol=4)
colnames(j_inertia[["proportion"]])<-c("joint_inertia","QnormD","Qnorm_D","QD")
j_inertia[["new_alpha"]]<-matrix(0,nrow=length(res_gas22[[1]]),ncol=(length(res_gas22[[1]])-1))
as.numeric(res_gas22[[4]][k-1,])
alpha_val<-as.numeric(res_gas22[[4]][k-1,])
d_tot<-res_gas22[[1]]
K<-4

for (i in 1:5){
  ji<-joint_inertia(d_tot=res_gas22[[1]],i=i,alpha_val=alpha_val,K=k)
  print(i)
  print(ji)
  j_inertia[["proportion"]][i,]<-ji[1:4]
  j_inertia[["new_alpha"]][i,]<-ji[5:length(ji)]
}
res_gas22[["joint_inertia"]]<-j_inertia


## Joint inertia in ghg emission by SECTOR clustering
k<-5
d_alp<-0.1
j_inertia<-list()
j_inertia[["proportion"]]<-matrix(0,nrow=length(res_sector22[[1]]),ncol=4)
colnames(j_inertia[["proportion"]])<-c("joint_inertia","QnormD","Qnorm_D","QD")
j_inertia[["new_alpha"]]<-matrix(0,nrow=length(res_sector22[[1]]),ncol=(length(res_sector22[[1]])-1))
for (i in 1:7){
  ji<-joint_inertia(res_sector22[[1]],i,res_sector22[[4]][k-1,],k)
  print(i)
  print(ji)
  j_inertia[["proportion"]][i,]<-ji[1:4]
  j_inertia[["new_alpha"]][i,]<-ji[5:length(ji)]
}
res_sector22[["joint_inertia"]]<-j_inertia


###############################################################
############## APPENDIX Temporal clustering ###################
###############################################################

## ghg emission by GAS clustering
k<-4
d_tot<-res_gas22[["diss"]]
data_gas22<-res_gas22[["data"]]
alpha<-res_gas22[["joint_inertia"]][["new_alpha"]][length(d_tot),]
d_tot<-d_tot[-length(d_tot)]
D <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
n<-nrow(d_tot[[1]])
delta_j <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
for (j in 1:(length(d_tot))) {
  delta_j <- as.matrix(wardinit(as.dist(d_tot[[j]])))
  D <- D + alpha[j] * as.matrix(delta_j)
}


tree <- stats::hclust(as.dist(delta_j), method = "ward.D")
part <- stats::cutree(tree, k = k)  # Cut the tree into K clusters

nuts2 <- nuts2[order(nuts2$NUTS_ID), ]
nuts2$tree<-part
gas22<-data_gas22[,1:35]
nut<-cbind(nuts2$NUTS_ID,nuts2$tree)
colnames(nut)<-c("NUTS_ID","tree")
gas22<-merge(gas22,nut,by.x="NUTS",by.y="NUTS_ID",all.x=TRUE)
centroids<-list()
for(gas_i in unique(gas22$Unit)){
  data<-gas22%>%filter(Unit==gas_i)
  #data<-merge(data,nuts[,c(1,14)],by.x=)
  centroids[[gas_i]]<-matrix()
  centroids[[gas_i]]<-aggregate(. ~ tree, data = data[,3:36], mean)
}


cluster_color_gas <- c("#9579De", "#451070", 
                       "#FF83FA","#d9108B")
k<-4
# Initialize an empty list to store plots
plots <- list()
#plot_map<-list()
# Loop over each gas in centroids[[1]]
for (i in 1:length(gas_name)) {
  # Extract the data matrix for the current gas
  gas_data <- centroids[[gas_name[i]]]
  # Melt the data to convert it from wide to long format
  gas_data_long <- melt(gas_data, id.vars = "tree", variable.name = "Year", value.name = "Observation")
  gas_data_long$Year <- as.numeric(gas_data_long$Year)
  years <- 1990:2022
  tick_years <- seq(1990, 2020, by = 5)
  p <- ggplot(gas_data_long, aes(x = Year, y = Observation, color = factor(tree), group = tree)) +
    geom_line(linewidth = 1.3) +  # Thicker lines
    labs(title = paste("Annual Trend for", gas_name[i]), x = "Year", y = "Annual Centroids", color = "Cluster") +
    theme_minimal() +
    scale_color_manual(values = cluster_color_gas) +
    scale_x_continuous(
      breaks = which(years %in% tick_years),     # index positions
      labels = tick_years
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10))

  plots[[i]]<-p
}
plt <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2, aes(fill = factor(tree)), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
  scale_fill_manual(values = cluster_color_gas) +  # use scale_fill_manual for fill
  theme_minimal() +
  labs(title = "Map of temporal clusters of emissions by gas", x = "Longitude", y = "Latitude", fill = "Clusters")  
print(plt)

# save as map_ts_clu_gas 1000x700
ggsave(
  "map_ts_clu_gas_rev.png",
  plot = plt,
  width = 1000, height = 700, units = "px",
  dpi = 96,                 # match screen DPI so text/lines aren't oversized
  bg = "white",             # avoid transparent background
  device = ragg::agg_png    # optional, for crisp text/lines
)
# Display all the plots in a grid layout
do.call(grid.arrange, c(plots, ncol = 2))
#save as ts_clu_centroids_gas 900x700



## ghg emission by SECTOR clustering

k<-5
d_tot<-res_sector22[["diss"]]
data_sector22<-res_sector22[["data"]]
alpha<-res_sector22[["joint_inertia"]][["new_alpha"]][length(d_tot),]
d_tot<-d_tot[-length(d_tot)]
D <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
n<-nrow(d_tot[[1]])
delta_j <- matrix(0, nrow = nrow(d_tot[[1]]), ncol = ncol(d_tot[[1]]))
for (j in 1:(length(d_tot))) {
  delta_j <-as.matrix(wardinit(as.dist(d_tot[[j]])))
  D <- D + alpha[j] * as.matrix(delta_j)
}
tree <- stats::hclust(as.dist(delta_j), method = "ward.D")
part <- stats::cutree(tree, k = k)  # Cut the tree into K clusters

nuts2 <- nuts2[order(nuts2$NUTS_ID), ]
nuts2$tree<-part
sector22<-data_sector22[,1:35]
nut<-cbind(nuts2$NUTS_ID,nuts2$tree)
colnames(nut)<-c("NUTS_ID","tree")
sector22<-merge(sector22,nut,by.x="NUTS",by.y="NUTS_ID",all.x=TRUE)
centroids<-list()
for(sector_i in unique(sector22$Sector)){
  data<-sector22%>%filter(Sector==sector_i)
  #data<-merge(data,nuts[,c(1,14)],by.x=)
  centroids[[sector_i]]<-matrix()
  centroids[[sector_i]]<-aggregate(. ~ tree, data = data[,3:36], mean)
}

cluster_color_sec <- c("#98FB98", "#107030", "#44e5ee", "#0592a0", "#00CD66")

k<-5
plots <- list()
for (i in 1:length(sector_name)) {
  # Extract the data matrix for the current gas
  sector_data <- centroids[[sector_name[i]]]
  # Melt the data to convert it from wide to long format
  sector_data_long <- melt(sector_data, id.vars = "tree", variable.name = "Year", value.name = "Observation")
  sector_data_long$Year <- as.numeric(sector_data_long$Year)

  years <- 1990:2022
  tick_years <- seq(1990, 2020, by = 5)
  p <- ggplot(sector_data_long, aes(x = Year, y = Observation, color = factor(tree), group = tree)) +
    geom_line(linewidth = 1.3) +  # Thicker lines
    labs(title = paste("Annual Trend for", sector_name[i]), x = "Year", y = "Annual Centroids", color = "Cluster") +
    theme_minimal() +
    scale_color_manual(values = cluster_color_sec) +
    scale_x_continuous(
      breaks = which(years %in% tick_years),     # index positions
      labels = tick_years
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10))

  plots[[i]]<-p
}
plt <- ggplot() +
  geom_sf(data = world, fill = "white", color = "gray50") +
  geom_sf(data = nuts2, aes(fill = factor(tree)), color = "gray50") +
  coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
  scale_fill_manual(values = cluster_color_sec) +  # use scale_fill_manual for fill
  theme_minimal() +
  labs(title = "Map of temporal clusters of emissions by sector", x = "Longitude", y = "Latitude", fill = "Clusters")  
print(plt)

# save as map_ts_clu_sector 1000x700
ggsave(
  "map_ts_clu_sector_rev.png",
  plot = plt,
  width = 1000, height = 700, units = "px",
  dpi = 96,                 # match screen DPI so text/lines aren't oversized
  bg = "white",             # avoid transparent background
  device = ragg::agg_png    # optional, for crisp text/lines
)
# Display all the plots in a grid layout
do.call(grid.arrange, c(plots, ncol = 2))
#save as ts_clu_centroids_sector 900x900


####################################################################
######### APPENDIX : CLUSTERING EACH GAS and each SECTOR ###########
####################################################################

##### CLUSTERING GAS

## Inputs assumed:
## - d_tot: list of dissimilarity matrices (gases first, SPATIAL LAST)
## - sptClu(): your function from above
## - d_alp: step for alpha grid (e.g., 0.05)

K_star <- 4
d_alp<-0.05
d_gas22<-res_gas22[["diss"]]
P <- length(d_gas22)                     # total number of matrices
spatial_idx <- P                       # last is spatial
gas_idx <- seq_len(P - 1)              # all but last

re_single_gas22 <- vector("list", length(gas_idx))
names(re_single_gas22) <- names(d_tot)[gas_idx]

for (g in gas_idx) {
  cat("Gas:", names(d_gas22)[g], " + spatial → optimizing...\n")
  
  # run the optimizer using only two matrices: current gas + spatial
  res <- sptClu(D_total = list(d_gas22[[g]], d_gas22[[spatial_idx]]),
                delta_alpha = d_alp,
                K = K_star)
  # pick the alpha combination with the highest weighted explained inertia
  best_row <- which.max(res$Q[, ncol(res$Q)])
  best_res <- res$results[[best_row]]
  
  # store results for this gas
  re_single_gas22[[g]] <- list(
    gas              = names(d_gas22)[g],
    matrices_used    = c(names(d_gas22)[g], names(d_gas22)[spatial_idx]),
    best_alpha       = best_res$alpha,                # length-2 vector (gas, spatial)
    partition        = best_res$part,                 # cluster labels (K = 4)
    D_combined       = best_res$D_combined,           # combined dissimilarity at best alpha
    W_best           = res$W[best_row, ],             # within-cluster diss for each matrix
    Q_best           = res$Q[best_row, ],             # [gas, spatial, weighted average]
    Qnorm_best       =  c(res$Q[best_row,1 ]/(res$Q[nrow(res$Q),1 ]), res$Q[best_row, 2]/res$Q[1, 2]),
    Q_all            = res$Q,                         # all Q rows (optional, for inspection)
    alpha_grid       = do.call(rbind, lapply(res$results, `[[`, "alpha"))  # tried alphas
  )
  
re_single_gas22[[g]][["joint_inertia"]]<-sum(re_single_gas22[[g]][["Qnorm_best"]])-1
cat("  Best alpha (gas, spatial):",
    paste(round(re_single_gas22[[g]]$best_alpha, 3), collapse = ", "),
    "  |  Q_wavg:",
    round(re_single_gas22[[g]]$Q_best[ncol(res$Q)], 3), 
    "  |  joint_inertia:",
    round(re_single_gas22[[g]]$joint_inertia, 3), "\n")

  }

# palette for K=4 (same family, well separated)
cluster_color_gas <- c("#9579DE", "#451070", "#FF83FA", "#D9108B")

# gas labels = everything except the spatial matrix (assumes names match gas22$Unit)
gas_labels <- names(d_gas22)[seq_len(length(d_gas22) - 1)]

# make sure spatial data order matches your dissimilarities (as you did before)
nuts2 <- nuts2[order(nuts2$NUTS_ID), ]

# wide data with columns: NUTS, Unit, 1990..2022
gas22_df <- res_gas22[["data"]]

# build all six maps and keep them in a list
maps <- vector("list", length(gas_labels))
names(maps) <- gas_labels

for (g in seq_along(gas_labels)) {
  gas_name_i <- gas_labels[g]
  part_i <- re_single_gas22[[g]]$partition  # best partition for this sector (K = 4)
  
  nuts_plot <- nuts2
  nuts_plot$tree <- as.integer(part_i)
  
  map_g <- ggplot() +
    geom_sf(data = world, fill = "white", color = "gray70") +
    geom_sf(data = nuts_plot, aes(fill = factor(tree)), color = "gray50", linewidth = 0.2) +
    coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
    scale_fill_manual(values = cluster_color_gas, drop = FALSE) +
    labs(title = paste("Map of spatiotemporal clusters of emissions of", gas_name_i),
         x = "Longitude", y = "Latitude", fill = "Cluster") +
    theme_minimal()
  
  maps[[g]] <- map_g
}

# arrange into a single picture: 2 columns × 3 rows
grid_out <- do.call(grid.arrange, c(maps, ncol = 2, nrow = 2))
# save as map_sp_single_gas 1100x900


##### CLUSTERING EACH Sector

## Inputs assumed:
## - d_tot: list of dissimilarity matrices (gases first, SPATIAL LAST)
## - sptClu(): your function from above
## - d_alp: step for alpha grid (e.g., 0.05)

K_star <- 5
d_alp<-0.05
d_sector22<-res_sector22[["diss"]]
P <- length(d_sector22)                     # total number of matrices
spatial_idx <- P                       # last is spatial
sector_idx <- seq_len(P - 1)              # all but last

re_single_sector22 <- vector("list", length(sector_idx))
names(re_single_sector22) <- names(d_sector22)[sector_idx]

for (g in sector_idx) {
  cat("sector:", names(d_sector22)[g], " + spatial → optimizing...\n")
  
  # run the optimizer using only two matrices: current gas + spatial
  res <- sptClu(D_total = list(d_sector22[[g]], d_sector22[[spatial_idx]]),
                delta_alpha = d_alp,
                K = K_star)
  # pick the alpha combination with the highest weighted explained inertia
  best_row <- which.max(res$Q[, ncol(res$Q)])
  best_res <- res$results[[best_row]]
  
  # store results for this gas
  re_single_sector22[[g]] <- list(
    gas              = names(d_sector22)[g],
    matrices_used    = c(names(d_sector22)[g], names(d_sector22)[spatial_idx]),
    best_alpha       = best_res$alpha,                # length-2 vector (gas, spatial)
    partition        = best_res$part,                 # cluster labels (K = 4)
    D_combined       = best_res$D_combined,           # combined dissimilarity at best alpha
    W_best           = res$W[best_row, ],             # within-cluster diss for each matrix
    Q_best           = res$Q[best_row, ],             # [gas, spatial, weighted average]
    Qnorm_best       =  c(res$Q[best_row,1 ]/(res$Q[nrow(res$Q),1 ]), res$Q[best_row, 2]/res$Q[1, 2]),
    Q_all            = res$Q,                         # all Q rows (optional, for inspection)
    alpha_grid       = do.call(rbind, lapply(res$results, `[[`, "alpha"))  # tried alphas
  )
  
re_single_sector22[[g]][["joint_inertia"]]<-sum(re_single_sector22[[g]][["Qnorm_best"]])-1
cat("  Best alpha (sector, spatial):",
    paste(round(re_single_sector22[[g]]$best_alpha, 3), collapse = ", "),
    "  |  Q_wavg:",
    round(re_single_sector22[[g]]$Q_best[ncol(res$Q)], 3), 
    "  |  joint_inertia:",
    round(re_single_sector22[[g]]$joint_inertia, 3), "\n")

  }

# palette for K=5 
cluster_color_sec <- c("#98FB98", "#107030", "#44e5ee", "#0592a0", "#00CD66")

# gas labels = everything except the spatial matrix (assumes names match gas22$Unit)
sector_labels <- names(d_sector22)[seq_len(length(d_sector22) - 1)]

# make sure spatial data order matches your dissimilarities (as you did before)
nuts2 <- nuts2[order(nuts2$NUTS_ID), ]

# wide data with columns: NUTS, Unit, 1990..2022
sector22_df <- res_sector22[["data"]]

# a place to store the plots
plots_single_sector22 <- vector("list", length(sector_labels))
names(plots_single_sector22) <- sector_labels

# build all six maps and keep them in a list
maps <- vector("list", length(sector_labels))
names(maps) <- sector_labels

for (g in seq_along(sector_labels)) {
  sector_name_i <- sector_labels[g]
  part_i <- re_single_sector22[[g]]$partition  # best partition for this sector (K = 4)
  
  nuts_plot <- nuts2
  nuts_plot$tree <- as.integer(part_i)
  
  map_g <- ggplot() +
    geom_sf(data = world, fill = "white", color = "gray70") +
    geom_sf(data = nuts_plot, aes(fill = factor(tree)), color = "gray50", linewidth = 0.2) +
    coord_sf(ylim = c(34, 70), xlim = c(-25, 40)) +
    scale_fill_manual(values = cluster_color_sec, drop = FALSE) +
    labs(title = paste("Map of spatiotemporal clusters of emissions in", sector_name_i),
         x = "Longitude", y = "Latitude", fill = "Cluster") +
    theme_minimal()
  
  maps[[g]] <- map_g
}

# arrange into a single picture: 2 columns × 3 rows
grid_out <- do.call(grid.arrange, c(maps, ncol = 2, nrow = 3))
# save as map_sp_single_sector 1100x1300



####################################################################
##################### APPENDIX : Moran's I #########################
####################################################################



# --------------------------- 1) Geometry & data prep ---------------------------
nuts_pts <- centroids_nuts2 %>%
  sf::st_as_sf() %>%
  sf::st_transform(3035) %>%
  dplyr::select(NUTS_ID, geometry)

gas_long <- data_gas22 %>%
  ungroup() %>%
  pivot_longer(
    cols = matches("^\\d{4}$"),
    names_to = "Year",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  mutate(Year = as.integer(Year))

nuts_gas_all <- nuts_pts %>%
  left_join(gas_long, by = c("NUTS_ID" = "NUTS"))

#Ks <- c(6, 10, 15)   # or Ks <- 6
Ks <- 6

compute_moran_group <- function(g, k) {
  g <- g %>% filter(!is.na(value))
  n <- nrow(g)
  if (n < 3 || sd(g$value) == 0) {
    return(tibble(
      N = n, k_used = ifelse(n > 1, min(k, n - 1), NA_integer_),
      Moran_I = NA_real_, E_I = NA_real_, Var_I = NA_real_,
      p_asympt = NA_real_, I_perm = NA_real_, p_perm = NA_real_
    ))
  }
  k_used <- min(k, n - 1)
  coords <- st_coordinates(g)
  nb <- knn2nb(knearneigh(coords, k = k_used, longlat = FALSE))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  mt <- moran.test(g$value, lw, zero.policy = TRUE)
  
  set.seed(123)
  mc <- moran.mc(g$value, lw, nsim = 999, zero.policy = TRUE)
  
  tibble(
    N = n,
    k_used = k_used,
    Moran_I = unname(mt$estimate["Moran I statistic"]),
    E_I     = unname(mt$estimate["Expectation"]),
    Var_I   = unname(mt$estimate["Variance"]),
    p_asympt = mt$p.value,
    I_perm   = as.numeric(mc$statistic),
    p_perm   = mc$p.value
  )
}

set.seed(123)

res_moran_gas <- map_dfr(Ks, function(k) {
  nuts_gas_all %>%
    group_by(Unit, Year) %>%
    group_modify(~ compute_moran_group(.x, k)) %>%
    ungroup() %>%
    mutate(k = k)
}) %>%
  # --- add 95% null band and significance flags
  mutate(
    ci_low   = E_I - 1.96 * sqrt(Var_I),
    ci_high  = E_I + 1.96 * sqrt(Var_I),
    sig_perm = p_perm   < 0.05,
    sig_asym = p_asympt < 0.05
  )


# Plot: time series
res_moran_gas <- res_moran_gas %>%
  mutate(sig_perm = p_perm < 0.05)

k_sel <- 6

p_gas <- res_moran_gas %>%
  filter(k == k_sel, !is.na(Moran_I)) %>%
  ggplot(aes(x = Year, y = Moran_I, color = Unit, group = Unit)) +
  # 95% null band computed inline: E[I] ± 1.96*SE
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = sig_perm), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     name = "Signif. alpha=0.05",
                     labels = c("No", "Yes")) +
  scale_x_continuous(breaks = seq(1990, 2022, by = 4)) +
  labs(
    title = paste0("Global Moran's I over time by gas "), #(k = ", k_sel, ")
    # subtitle = "Grey ribbon = asymptotic 95% null band [E(I) ± 1.96·SE]; filled points = significant by permutation",
    y = "Moran's I", x = NULL, color = "Gas"
  ) +
  theme_minimal()

p_gas
# save moran_gas 700x400

for (gas_i in gas_name){
  res_moran_gas_i<-res_moran_gas%>%filter(Unit==gas_i)
  mean_moran_gas_i<-mean(res_moran_gas_i$Moran_I)
  print(gas_i)
  print(mean_moran_gas_i)
}


###### MORAN I by Sector

# NUTS2 centroids in a projected CRS (meters)

nuts_pts <- centroids_nuts2 %>%
  sf::st_as_sf() %>%
  sf::st_transform(3035) %>%
  dplyr::select(NUTS_ID, geometry)

# Long table: NUTS,Sector, Year, value
sector_long <- data_sector22 %>%
  ungroup() %>%
  pivot_longer(
    cols = matches("^\\d{4}$"),               # year columns 1990..2022
    names_to = "Year",
    values_to = "value",
    values_drop_na = FALSE
  ) %>%
  mutate(Year = as.integer(Year))

# Join values to centroids (one row per NUTS_ID × Unit × Year)
nuts_sector_all <- nuts_pts %>%
  left_join(sector_long, by = c("NUTS_ID" = "NUTS"))

#--------------------------- 2) Moran per sector × year ---------------------------
# Choose the k values you want to evaluate
#Ks <- c(6, 10, 15)   # or just Ks <- 6
Ks <- 6

compute_moran_group <- function(g, k) {
  # g is an sf data frame for one gas × year (has geometry + 'value')
  g <- g %>% filter(!is.na(value))
  n <- nrow(g)
  if (n < 3 || sd(g$value) == 0) {
    return(tibble(
      N = n, k_used = ifelse(n > 1, min(k, n - 1), NA_integer_),
      Moran_I = NA_real_, E_I = NA_real_, Var_I = NA_real_,
      p_asympt = NA_real_, I_perm = NA_real_, p_perm = NA_real_
    ))
  }
  k_used <- min(k, n - 1)
  coords <- st_coordinates(g)                 # projected coords (meters)
  nb <- knn2nb(knearneigh(coords, k = k_used, longlat = FALSE))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  mt <- moran.test(g$value, lw, zero.policy = TRUE)
  
  # Permutation test
  mc <- moran.mc(g$value, lw, nsim = 999, zero.policy = TRUE)
  
  tibble(
    N = n,
    k_used = k_used,
    Moran_I = unname(mt$estimate["Moran I statistic"]),
    E_I     = unname(mt$estimate["Expectation"]),
    Var_I   = unname(mt$estimate["Variance"]),
    p_asympt = mt$p.value,
    I_perm   = as.numeric(mc$statistic),
    p_perm   = mc$p.value
  )
}

set.seed(123)  # for permutation tests reproducibility

res_moran_sector <- map_dfr(Ks, function(k) {
  nuts_sector_all %>%
    group_by(Sector, Year) %>%                  # each gas × year
    group_modify(~ compute_moran_group(.x, k)) %>%
    ungroup() %>%
    mutate(k = k)
})


# Pick which k to visualize (if multiple were computed)

res_moran_sector <- res_moran_sector %>%
  mutate(sig_perm = p_perm < 0.05)

k_sel <- 6

p_sector <- res_moran_sector %>%
  filter(k == k_sel, !is.na(Moran_I)) %>%
  ggplot(aes(x = Year, y = Moran_I, color = Sector, group = Sector)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  geom_point(aes(shape = sig_perm), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     name = "Signif. alpha=0.05",
                     labels = c("Yes", "Yes")) +
  scale_x_continuous(breaks = seq(1990, 2022, by = 4)) +
  labs(
    title = paste0("Global Moran's I over time by sector"), # (k = ", k_sel, ")
    # subtitle = "Grey ribbon = asymptotic 95% null band [E(I) ± 1.96·SE]; filled points = significant by permutation",
    y = "Moran's I", x = NULL, color = "Sector"
  ) +
  theme_minimal()

p_sector
# save moran_sector 700x400

for (sector_i in sector_name){
  res_moran_sector_i<-res_moran_sector%>%filter(Sector==sector_i)
  mean_moran_sector_i<-mean(res_moran_sector_i$Moran_I)
  print(sector_i)
  print(mean_moran_sector_i)
}

