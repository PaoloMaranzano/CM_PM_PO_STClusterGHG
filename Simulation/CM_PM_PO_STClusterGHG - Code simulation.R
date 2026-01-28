library(MASS)
library(ggplot2)
library(ClustGeo)
library(forecast)
library(gridExtra)
library(dplyr)
library(caret)
library(mclust)
library(clue)
sim_res <- data.frame(matrix(ncol = 23, nrow=0))
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max","adj_rand_ch","adj_rand_max") # 2 correct class rate
# Set the seed for reproducibility
set.seed(123)
for (d in c(0,1/3,2/3,1)){
  # Define the number of clusters and their probabilities
  n_obs<-100
  num_clusters <- 4
  cluster_probs <- c(0.25, 0.25, 0.25, 0.25)
  # Define the means and covariance matrices for each cluster
  means_space <- list(
    c(d, d),
    c(d, (-d)),
    c((-d), d),
    c((-d), (-d))
  )
  v_sp<-0.4
  covs_space <- list(
    matrix(c(1, 0, 0, 1)*v_sp, ncol=2),
    matrix(c(1, 0, 0, 1)*v_sp, ncol=2),
    matrix(c(1, 0, 0, 1)*v_sp, ncol=2),
    matrix(c(1, 0., 0, 1)*v_sp, ncol=2)
  )
  # Define the means and variances for var1 in each cluster
  var1_means <- c(2,4,6,8)
  v<-0.4 #variablity of var1
  var1_sds <- c(1, 1, 1, 1)*v
  
  x <- numeric(0)
  y <- numeric(0)
  var1 <- numeric(0)
  cluster_labels <- numeric(0)
  
  # Initialize the dataframe
  data <- data.frame(x = numeric(0), y = numeric(0), var1 = numeric(0), var2 = numeric(0), cluster = integer(0))
  
  # Simulation loop
  while (nrow(data) < n_obs) {
    # Step 1: Determine cluster assignment
    cluster <- sample(num_clusters, 1, prob = cluster_probs)
    
    # Step 2: Simulate position in space
    position <- mvrnorm(1, means_space[[cluster]], covs_space[[cluster]])
    
    # Step 3: Check if position is within [0, 1] x [0, 1]
    # if (all(position >= 0 & position <= 1)) {
    # Step 4: Simulate var1 and var2
    mean_vars <- var1_means[[cluster]]
    sd_vars <- var1_sds[[cluster]]
    var1 <- rnorm(1, mean = mean_vars, sd = sd_vars)
    
    # Step 5: Append to the dataframe
    data <- rbind(data, data.frame(x = position[1], y = position[2], var1 = var1, cluster = cluster))
    #}
  }
  data$cluster <- as.factor(data$cluster)
  cluster_colors <- c("blue", "green", "orange", "purple")
  ggplot(data, aes(x=x, y=y, color=cluster, size=var1)) +
    geom_point() +
    scale_color_manual(values=cluster_colors) +
    labs(title="Scatter Plot of Clusters",
         x="X Coordinate",
         y="Y Coordinate",
         color="Cluster",
         size="Value of var1") +
    theme_minimal()
  
  d0<-dist(data[,3:4])
  d0<-d0/max(d0)
  d1<-dist(data[,1:2])
  d1<-d1/max(d1)
  T1 <- inertdiss(d0)
  T2 <- inertdiss(d1)
  
  
  cr<-choicealpha(as.dist(d0),as.dist(d1),range.alpha = seq(0,1,0.05),K=4,graph = FALSE)
  a<-as.data.frame(cr$Qnorm)
  b<-as.data.frame(cr$Q)
  a$tot<-a$Q0norm+a$Q1norm
  b$tot<-(b$Q0*T1+b$Q1*T2)/(T1+T2)
  a_ch<-cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  a_max<-cr$range.alpha[which.max(b$tot)]  
  
  
  tree_ch<-hclustgeo(as.dist(d0),as.dist(d1), alpha=a_ch)
  clu_ch<-cutree(tree_ch,k=4)
  data<-cbind(data,clu_ch)
  
  tree_max<-hclustgeo(as.dist(d0),as.dist(d1), alpha=a_max)
  clu_max<-cutree(tree_max,k=4)
  data<-cbind(data,clu_max)
  
  tree_d0<-hclust(as.dist(d0))
  clu_d0<-cutree(tree_d0,k=4)
  data<-cbind(data,clu_d0)
  
  tree_d1<-hclust(as.dist(d1))
  clu_d1<-cutree(tree_d1,k=4)
  data<-cbind(data,clu_d1)
  
  #table(data$cluster,data$clu_ch)
  #table(data$cluster,data$clu_max)
  
  c_ch <- aggregate(. ~ clu_ch, data = data[,c(3,5)], mean) 
  c_max <- aggregate(. ~ clu_max, data = data[,c(3,6)], mean) 
  c_ch <- c_ch %>%
    mutate(clu_ch2 = ntile(var1, 4))
  c_max <- c_max %>%
    mutate(clu_max2 = ntile(var1, 4))
  
  data <- data %>%
    mutate(clu_max = replace(clu_max, clu_max %in% c_max$clu_max, c_max$clu_max2[match(clu_max, c_max$clu_max)]))
  data <- data %>%
    mutate(clu_ch = replace(clu_ch, clu_ch %in% c_ch$clu_ch, c_ch$clu_ch2[match(clu_ch, c_ch$clu_ch)]))
  data2<-data
  
  
  plot_ch<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_ch), fill=as.factor(cluster), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Clustering with benchmark method",
         x="X Coordinate",
         y="Y Coordinate",
         size="Value of variable Z") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  plot_max<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_max), fill=as.factor(cluster), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Clustering with optimal average inertia",
         x="X Coordinate",
         y="Y Coordinate",
         size="Value of variable Z") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  grid.arrange(plot_ch, plot_max, ncol = 2)
  
  acc_ch<-sum(data2$cluster==data2$clu_ch)
  acc_max<-sum(data2$cluster==data2$clu_max)
  
  i_ch<-as.numeric(which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2])))
  i_max<-as.numeric(which.max(b$tot))
  #proportion of inertia and weighted inertia
  q_d0_ch<-b$Q0[i_ch]
  q_d0_max<-b$Q0[i_max]
  
  q_d1_ch<-b$Q1[i_ch]
  q_d1_max<-b$Q1[i_max]
  
  wei_ch<-b$tot[i_ch]
  wei_max<-b$tot[i_max]
  
  #normalized inertia, d0 d1 total
  qn_d0_ch<-a$Q0norm[i_ch]
  qn_d0_max<-a$Q0norm[i_max]
  
  qn_d1_ch<-a$Q1norm[i_ch]
  qn_d1_max<-a$Q1norm[i_max]
  
  qn_ch<-a$tot[i_ch]
  qn_max<-a$tot[i_max]
  
  adj_rand_ch <- adjustedRandIndex(data$clu_ch, data$cluster)
  adj_rand_max <- adjustedRandIndex(data$clu_max, data$cluster)
  
  sim_result<-c(n_obs,num_clusters,d,v_sp,v, #input
                a_ch,a_max, # 2 alpha
                q_d0_ch,q_d0_max,q_d1_ch,q_d1_max, # 4 prop inertia
                qn_d0_ch,qn_d0_max,qn_d1_ch,qn_d1_max,# 4 norm prop inertia
                wei_ch,wei_max,qn_ch,qn_max, # 4 weighted and tot norm inertia
                acc_ch,acc_max,
                adj_rand_ch,adj_rand_max) # 2 correct classification rate
  sim_res<-rbind(sim_res,sim_result)
}
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max","adj_rand_ch","adj_rand_max")

#######################################################################
############################# MONTE CARLO #############################
#######################################################################


sample_clu_prob <- function(min_val, max_val, length_out, target_sum) {
  while(TRUE) {
    vec <- runif(length_out, min_val, max_val)
    vec <- vec / sum(vec) * target_sum
    if(all(vec >= min_val & vec <= max_val)) {
      return(vec)
    }
  }
}
cluster_colors <- c("blue", "green", "orange", "purple")


# Add new columns for Adjusted Rand Index, Precision, Recall, and F1 score
sim_res <- data.frame(matrix(ncol = 39, nrow = 0))
colnames(sim_res) <- c("n_obs", "n_clu", "d", "v_sp", "v", # 5 input
                       "clu_pr1","clu_pr2","clu_pr3","clu_pr4",#clusters probabilities
                       "a_ch", "a_max", # 2 alpha
                       "q_d0_ch", "q_d0_max", "q_d1_ch", "q_d1_max", # 4 prop inertia
                       "qn_d0_ch", "qn_d0_max", "qn_d1_ch", "qn_d1_max", # 4 norm prop inertia
                       "wei_ch", "wei_max", "qn_ch", "qn_max", # 4 weighted and tot norm inertia
                       "mean_a","mean_b",
                       "qn_d0_0","qn_d0_1", "qn_d1_0", "qn_d1_1",
                       "acc_ch", "acc_max", # 2 correct class rate
                       "adj_rand_ch", "adj_rand_max", # Adjusted Rand Index
                       "prec_ch", "recall_ch", "f1_ch", # Precision, Recall, F1 for Chavent
                       "prec_max", "recall_max", "f1_max") # Precision, Recall, F1 for Max Inertia
set.seed(123)
# Same simulation setup as before
n_mc<-500
n_obs <- 100
num_clusters <- 4
sim_list <- list()
v_max<-0.6
v_min<-0.1
v_sp<-0.4
v<-0.4
clu_prob_min<-0.15
clu_prob_max<-0.35
due_variabili<-FALSE
plt <- FALSE
values <- seq(0, 1, by = 0.02)

for (i in 1:n_mc) {
  # Define random parameters for each simulation run
  d <- runif(1, 0, 1)

  cluster_probs <- sample_clu_prob(clu_prob_min, clu_prob_max, 4, 1)
  means_space <- list(c(d, d), c(d, -d), c(-d, d), c(-d, -d))

  covs_space <- list(matrix(c(1, 0, 0, 1) * v_sp, ncol = 2),
                     matrix(c(1, 0, 0, 1) * v_sp, ncol = 2),
                     matrix(c(1, 0, 0, 1) * v_sp, ncol = 2),
                     matrix(c(1, 0, 0, 1) * v_sp, ncol = 2))
  if (due_variabili==TRUE){
    means_vars <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
  }else{
    means_vars <- list(2,4,6,8)
  }
  
  v <- runif(1, min = v_min, max = v_max)
  vars_vars <- list(1 * v, 1 * v, 1 * v, 1 * v)
  
  # Initialize data for this simulation
  data <- data.frame(x = numeric(0), y = numeric(0), var1 = numeric(0), var2 = numeric(0), cluster = integer(0))
  
  while (nrow(data) < n_obs) {
    cluster <- sample(num_clusters, 1, prob = cluster_probs)
    position <- mvrnorm(1, means_space[[cluster]], covs_space[[cluster]])
    mean_vars <- means_vars[[cluster]]
    sd_vars <- sqrt(vars_vars[[cluster]])
    var1 <- rnorm(1, mean = mean_vars, sd = sd_vars)
    if (due_variabili==TRUE){
      var2 <- rnorm(1, mean = mean_vars[2], sd = sd_vars[2])
      data <- rbind(data, data.frame(x = position[1], y = position[2], var1 = var1, var2 = var2, cluster = cluster))
    }else{
      data <- rbind(data, data.frame(x = position[1], y = position[2], var1 = var1, cluster = cluster))
    }
  }
  data$cluster <- as.factor(data$cluster)
  
  # Distances
  if (due_variabili==TRUE){
    d0 <- dist(data[, 3:4]) / max(dist(data[, 3:4]))
  }else{
    d0 <- dist(data[, 3]) / max(dist(data[, 3]))
  }
  d1 <- dist(data[, 1:2]) / max(dist(data[, 1:2]))
  T1 <- inertdiss(d0)
  T2 <- inertdiss(d1)
  
  # Calculate alpha for both methods
  cr <- choicealpha(as.dist(d0), as.dist(d1), range.alpha = seq(0, 1, 0.05), K = 4, graph = FALSE)
  
  a<-as.data.frame(cr$Qnorm)
  b<-as.data.frame(cr$Q)
  a$tot<-a$Q0norm+a$Q1norm
  b$tot<-(b$Q0*T1+b$Q1*T2)/(T1+T2)
  a_ch<-cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  a_max<-cr$range.alpha[which.max(b$tot)]  
  # Cluster using both Chavent and Maximum Inertia methods
  tree_ch <- hclustgeo(as.dist(d0), as.dist(d1), alpha = a_ch)
  clu_ch <- cutree(tree_ch, k = 4)
  
  tree_max <- hclustgeo(as.dist(d0), as.dist(d1), alpha = a_max)
  clu_max <- cutree(tree_max, k = 4)
  
  ###############
  
  # Create confusion matrix
  conf_matrix <- table(data$cluster, clu_ch)
  # Use the Hungarian algorithm to find the optimal assignment
  assignment <- solve_LSAP(conf_matrix, maximum = TRUE)
  # Reassign clusters in clu_ch based on the optimal assignment
  clu_ch_new <- clu_ch
  for (i in 1:length(assignment)) {
    clu_ch_new[clu_ch == assignment[i]] <- i
  }
  data$clu_ch <- clu_ch_new
  
  
  conf_matrix <- table(data$cluster, clu_max)
  # Use the Hungarian algorithm to find the optimal assignment
  assignment <- solve_LSAP(conf_matrix, maximum = TRUE)
  # Reassign clusters in clu_ch based on the optimal assignment
  clu_max_new <- clu_max
  for (i in 1:length(assignment)) {
    clu_max_new[clu_max == assignment[i]] <- i
  }
  data$clu_max <- clu_max_new
  
  # Calculate Adjusted Rand Index for both methods
  adj_rand_ch <- adjustedRandIndex(data$clu_ch, data$cluster)
  adj_rand_max <- adjustedRandIndex(data$clu_max, data$cluster)
  
  # Precision, Recall, and F1 Score for Chavent Method
  conf_matrix_ch <- confusionMatrix(as.factor(data$clu_ch), data$cluster)
  prec_ch <- mean(conf_matrix_ch[["byClass"]][,3]) # precision = Pos Pred Value
  recall_ch <- mean(conf_matrix_ch[["byClass"]][,1]) #recall = sensitivity
  f1_ch <- 2 * (prec_ch * recall_ch) / (prec_ch + recall_ch)
  
  # Precision, Recall, and F1 Score for Maximum Inertia Method
  conf_matrix_max <- confusionMatrix(as.factor(data$clu_max), data$cluster)
  prec_max <- mean(conf_matrix_max[["byClass"]][,3]) # precision = Pos Pred Value
  recall_max <- mean(conf_matrix_max[["byClass"]][,1]) #recall = sensitivity
  f1_max <- 2 * (prec_max * recall_max) / (prec_max + recall_max)
  
  # Compute other metrics as in your original code (e.g., weighted and normalized inertia)
  i_ch<-as.numeric(which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2])))
  i_max<-as.numeric(which.max(b$tot))
  #proportion of inertia and weighted inertia
  q_d0_ch<-b$Q0[i_ch]
  q_d0_max<-b$Q0[i_max]
  
  q_d1_ch<-b$Q1[i_ch]
  q_d1_max<-b$Q1[i_max]
  
  wei_ch<-b$tot[i_ch]
  wei_max<-b$tot[i_max]
  
  #normalized inertia, d0 d1 total
  qn_d0_ch<-a$Q0norm[i_ch]
  qn_d0_max<-a$Q0norm[i_max]
  
  qn_d1_ch<-a$Q1norm[i_ch]
  qn_d1_max<-a$Q1norm[i_max]
  
  qn_ch<-a$tot[i_ch]
  qn_max<-a$tot[i_max]
  
  mean_a<-mean(a$tot)
  mean_b<-mean(b$tot)
  
  qn_d0_0<-a$Q0norm[1]
  qn_d0_1<-a$Q0norm[length(a$Q0norm)]
  qn_d1_0<-a$Q1norm[1]
  qn_d1_1<-a$Q1norm[length(a$Q1norm)]
  
  acc_ch <- sum(data$cluster == data$clu_ch) / n_obs
  acc_max <- sum(data$cluster == data$clu_max) / n_obs
  if(plt==TRUE){
    plot_ch<-ggplot(data, aes(x=x, y=y)) +
      geom_point(aes(color=as.factor(clu_ch), fill=as.factor(cluster), size=var1), shape=21, stroke=1.5) +
      scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
      scale_fill_manual(values=cluster_colors, name="Real Cluster") +
      labs(title="Clustering with Chavent method",
           x="X Coordinate",
           y="Y Coordinate",
           size="Value of var1") +
      theme_minimal() +
      #theme(legend.position = "right")
      theme(legend.position = "none")
    plot_max<-ggplot(data, aes(x=x, y=y)) +
      geom_point(aes(color=as.factor(clu_max), fill=as.factor(cluster), size=var1), shape=21, stroke=1.5) +
      scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
      scale_fill_manual(values=cluster_colors, name="Real Cluster") +
      labs(title="Clustering maximum average inertia",
           x="X Coordinate",
           y="Y Coordinate",
           size="Value of var1") +
      theme_minimal() +
      theme(legend.position = "right")
    #theme(legend.position = "none")
    grid.arrange(plot_ch, plot_max, ncol = 2)
  }
  # Store the results
  sim_result <- c(n_obs, num_clusters, d, v_sp, v, 
                  cluster_probs, a_ch, a_max,
                  q_d0_ch, q_d0_max, q_d1_ch, q_d1_max,
                  qn_d0_ch, qn_d0_max, qn_d1_ch, qn_d1_max,
                  wei_ch, wei_max, qn_ch, qn_max,
                  mean_a,mean_b,
                  qn_d0_0,qn_d0_1, qn_d1_0, qn_d1_1,
                  acc_ch, acc_max,
                  adj_rand_ch, adj_rand_max,
                  prec_ch, recall_ch, f1_ch,
                  prec_max, recall_max, f1_max)
  
  sim_res <- rbind(sim_res, sim_result)
}


colnames(sim_res) <- c("n_obs", "n_clu", "d", "v_sp", "v", # 5 input
                       "clu_pr1","clu_pr2","clu_pr3","clu_pr4",#clusters probabilities
                       "a_ch", "a_max", # 2 alpha
                       "q_d0_ch", "q_d0_max", "q_d1_ch", "q_d1_max", # 4 prop inertia
                       "qn_d0_ch", "qn_d0_max", "qn_d1_ch", "qn_d1_max", # 4 norm prop inertia
                       "wei_ch", "wei_max", "qn_ch", "qn_max", # 4 weighted and tot norm inertia
                       "mean_a","mean_b",
                       "qn_d0_0","qn_d0_1", "qn_d1_0", "qn_d1_1",
                       "acc_ch", "acc_max", # 2 correct class rate
                       "adj_rand_ch", "adj_rand_max", # Adjusted Rand Index
                       "prec_ch", "recall_ch", "f1_ch", # Precision, Recall, F1 for Chavent
                       "prec_max", "recall_max", "f1_max") # Precision, Recall, F1 for Max Inertia
#################################################################
sim_res <- sim_res[order(sim_res$d),]

# 4 plot in one image, representing the accuracy,sensitivity,recall and f1 for both methodology. 
color_method<-c("green","#FF33FF")


par(mfrow = c(2, 2))
# 1. Accuracy plot
plot(sim_res$d, sim_res$acc_ch, col = color_method[1], pch = 16, cex=0.9,
     xlab = "d = overlapping parameter", ylab = "rate", main = "Accuracy")
points(sim_res$d, sim_res$acc_max, col = color_method[2], pch = 16, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$acc_ch, span = 0.75), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$acc_max, span = 0.75), col = color_method[2], lwd = 2)
#legend("bottomright", legend = c("Chavent (2018)", "Morelli (2024)"), col = color_method, pch = 16, lwd = 2)

# 2. Precision plot
plot(sim_res$d, sim_res$prec_ch, col = color_method[1], pch = 16, 
     xlab = "d = overlapping parameter", ylab = "rate", main = "Precision", cex=0.9)
points(sim_res$d, sim_res$prec_max, col = color_method[2], pch = 16, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$prec_ch, span = 0.75), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$prec_max, span = 0.75), col = color_method[2], lwd = 2)
#legend("bottomright", legend = c("Chavent (2018)", "Morelli (2024)"), col = color_method, pch = 16, lwd = 2)

# 3. Recall plot
plot(sim_res$d, sim_res$recall_ch, col = color_method[1], pch = 16, 
     xlab = "d = overlapping parameter", ylab = "rate", main = "Sensitivity", cex=0.9)
points(sim_res$d, sim_res$recall_max, col = color_method[2], pch = 16, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$recall_ch, span = 0.75), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$recall_max, span = 0.75), col = color_method[2], lwd = 2)
#legend("bottomright", legend = c("Chavent (2018)", "Morelli (2024)"), col = color_method, pch = 16, lwd = 2)

# 4. Adjusted Rand Index
plot(sim_res$d, sim_res$adj_rand_ch, col = color_method[1], pch = 16, 
     xlab = "d = overlapping parameter", ylab = "rate", main = "Adjusted Rand Index", cex=0.9)
points(sim_res$d, sim_res$adj_rand_max, col = color_method[2], pch = 16, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$adj_rand_ch, span = 0.75), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$adj_rand_max, span = 0.75), col = color_method[2], lwd = 2)
#legend("bottomright", legend = c("Chavent (2018)", "Morelli (2024)"), col = color_method, pch = 16, lwd = 2)

# Reset layout to default
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))

plot(sim_res$d, sim_res$a_ch, col = color_method[1], pch = 1, 
     xlab = "d = overlapping parameter", ylab = "alpha", main = "Alpha",ylim=c(0,1), cex=0.9)
points(sim_res$d, sim_res$a_max, col = color_method[2], pch = 1, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$a_ch, span = 0.5), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$a_max, span = 0.5), col = color_method[2], lwd = 2)
legend("topright", legend = c("Chavent (2018)", "Morelli (2024)"), 
       col = color_method, pch = 1,pt.lwd = 2)


plot(sim_res$d, (sim_res$qn_ch-1), col = color_method[1], pch = 1, 
     xlab = "d = overlapping parameter", ylab = "Joint Inertia", main = "Jount Inertia",ylim=c(0,1), cex=0.9)
points(sim_res$d, (sim_res$qn_max-1), col = color_method[2], pch = 1, cex=0.9)
lines(loess.smooth(sim_res$d, sim_res$qn_ch-1, span = 1), col = color_method[1], lwd = 2)
lines(loess.smooth(sim_res$d, sim_res$qn_max-1, span = 1), col = color_method[2], lwd = 2)
legend("bottomright", legend = c("Chavent (2018)", "Morelli (2024)"), 
       col = color_method, pch = 16, lwd = 2)

par(mfrow = c(1,1))



### summary
stats_df <- data.frame(
  Variable = c("acc_ch", "acc_max", "prec_ch", "prec_max", 
               "recall_ch", "recall_max", "adj_rand_ch", "adj_rand_max"),
  Mean = c(
    mean(sim_res$acc_ch, na.rm = TRUE),
    mean(sim_res$acc_max, na.rm = TRUE),
    mean(sim_res$prec_ch, na.rm = TRUE),
    mean(sim_res$prec_max, na.rm = TRUE),
    mean(sim_res$recall_ch, na.rm = TRUE),
    mean(sim_res$recall_max, na.rm = TRUE),
    mean(sim_res$adj_rand_ch, na.rm = TRUE),
    mean(sim_res$adj_rand_max, na.rm = TRUE)
  ),
  Median = c(
    median(sim_res$acc_ch, na.rm = TRUE),
    median(sim_res$acc_max, na.rm = TRUE),
    median(sim_res$prec_ch, na.rm = TRUE),
    median(sim_res$prec_max, na.rm = TRUE),
    median(sim_res$recall_ch, na.rm = TRUE),
    median(sim_res$recall_max, na.rm = TRUE),
    median(sim_res$adj_rand_ch, na.rm = TRUE),
    median(sim_res$adj_rand_max, na.rm = TRUE)
  ),
  SD = c(
    sd(sim_res$acc_ch, na.rm = TRUE),
    sd(sim_res$acc_max, na.rm = TRUE),
    sd(sim_res$prec_ch, na.rm = TRUE),
    sd(sim_res$prec_max, na.rm = TRUE),
    sd(sim_res$recall_ch, na.rm = TRUE),
    sd(sim_res$recall_max, na.rm = TRUE),
    sd(sim_res$adj_rand_ch, na.rm = TRUE),
    sd(sim_res$adj_rand_max, na.rm = TRUE)
  )
)

##############################################################
############## The role of the spatial component #############
##############################################################

sim_res <- data.frame(matrix(ncol = 21, nrow=0))
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max") # 2 correct class rate
# Set the seed for reproducibility
set.seed(12)

for (d in c(0,1)){
  # Define the number of clusters and their probabilities
  n_obs<-150
  num_clusters <- 2
  cluster_probs <- c(0.5,0.5)
  
  means_space <- list(
    c(d, d),
    c((-d), (-d))
  )
  v_sp<-0.3
  covs_space <- list(
    matrix(c(1, 0, 0, 1)*v_sp, ncol=2),
    matrix(c(1, 0., 0, 1)*v_sp, ncol=2)
  )
  
  # Define the means and variances for var1 in each cluster
  var1_means <- c(1, 2)
  v<-0.32 #variablity of var1
  var1_sds <- c(1, 1)*v
  
  x <- numeric(0)
  y <- numeric(0)
  var1 <- numeric(0)
  cluster_labels <- numeric(0)
  
  # Initialize the dataframe
  data <- data.frame(x = numeric(0), y = numeric(0), var1 = numeric(0), var2 = numeric(0), cluster = integer(0))
  
  # Simulation loop
  while (nrow(data) < n_obs) {
    # Step 1: Determine cluster assignment
    cluster <- sample(num_clusters, 1, prob = cluster_probs)
    
    # Step 2: Simulate position in space
    position <- mvrnorm(1, means_space[[cluster]], covs_space[[cluster]])
    # position<-runif(n=2,min=0,max=1)
    
    # Step 3: Check if position is within [0, 1] x [0, 1]
    # if (all(position >= 0 & position <= 1)) {
    # Step 4: Simulate var1 and var2
    mean_vars <- var1_means[[cluster]]
    sd_vars <- var1_sds[[cluster]]
    var1 <- rnorm(1, mean = mean_vars, sd = sd_vars)
    
    # Step 5: Append to the dataframe
    data <- rbind(data, data.frame(x = position[1], y = position[2], var1 = var1, cluster = cluster))
    #}
  }
  data$cluster <- as.factor(data$cluster)
  cluster_colors <- c("blue", "orange")
  ggplot(data, aes(x=x, y=y, color=cluster, size=var1)) +
    geom_point() +
    scale_color_manual(values=cluster_colors) +
    labs(title="Scatter Plot of Clusters",
         x="X Coordinate",
         y="Y Coordinate",
         color="Cluster",
         size="Value of var1") +
    theme_minimal()
  
  d0<-dist(data[,3:4])
  d0<-d0/max(d0)
  d1<-dist(data[,1:2])
  d1<-d1/max(d1)
  T1 <- inertdiss(d0)
  T2 <- inertdiss(d1)
  
  
  cr<-choicealpha(as.dist(d0),as.dist(d1),range.alpha = seq(0,1,0.05),K=4,graph = FALSE)
  a<-as.data.frame(cr$Qnorm)
  b<-as.data.frame(cr$Q)
  a$tot<-a$Q0norm+a$Q1norm
  b$tot<-(b$Q0*T1+b$Q1*T2)/(T1+T2)
  a_ch<-cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  a_max<-cr$range.alpha[which.max(b$tot)]  
  
  
  tree_ch<-hclustgeo(as.dist(d0),as.dist(d1), alpha=a_ch)
  clu_ch<-cutree(tree_ch,k=2)
  data<-cbind(data,clu_ch)
  
  tree_max<-hclustgeo(as.dist(d0),as.dist(d1), alpha=a_max)
  clu_max<-cutree(tree_max,k=2)
  data<-cbind(data,clu_max)
  
  tree_d0<-hclust(as.dist(d0))
  clu_d0<-cutree(tree_d0,k=2)
  data<-cbind(data,clu_d0)
  
  tree_d1<-hclust(as.dist(d1))
  clu_d1<-cutree(tree_d1,k=2)
  data<-cbind(data,clu_d1)
  
  
  c_ch <- aggregate(. ~ clu_ch, data = data[,c(3,5)], mean) 
  c_max <- aggregate(. ~ clu_max, data = data[,c(3,6)], mean) 
  c_ch <- c_ch %>%
    mutate(clu_ch2 = ntile(var1, 2))
  c_max <- c_max %>%
    mutate(clu_max2 = ntile(var1, 2))
  
  c_d0 <- aggregate(. ~ clu_d0, data = data[,c(3,7)], mean) 
  c_d1 <- aggregate(. ~ clu_d1, data = data[,c(3,8)], mean) 
  c_d0 <- c_d0 %>%
    mutate(clu_d02 = ntile(var1, 2))
  c_d1 <- c_d1 %>%
    mutate(clu_d12 = ntile(var1, 2))
  
  data <- data %>%
    mutate(clu_max = replace(clu_max, clu_max %in% c_max$clu_max, c_max$clu_max2[match(clu_max, c_max$clu_max)]))
  data <- data %>%
    mutate(clu_ch = replace(clu_ch, clu_ch %in% c_ch$clu_ch, c_ch$clu_ch2[match(clu_ch, c_ch$clu_ch)]))
  data2<-data
  
  data <- data %>%
    mutate(clu_d1 = replace(clu_d1, clu_d1 %in% c_d1$clu_d1, c_d1$clu_d12[match(clu_d1, c_d1$clu_d1)]))
  data <- data %>%
    mutate(clu_d0 = replace(clu_d0, clu_d0 %in% c_d0$clu_d0, c_d0$clu_d02[match(clu_d0, c_d0$clu_d0)]))
  data2<-data
  
  plot_clu<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(cluster), fill=as.factor(cluster), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="True Clustering partition",
         x="X Coordinate",
         y="Y Coordinate",
         size="Mean of var1") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  plot_clu
  plot_ch<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_ch), fill=as.factor(clu_ch), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Clustering with Chavent method",
         x="X Coordinate",
         y="Y Coordinate",
         size="Mean of var1") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  plot_max<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_max), fill=as.factor(clu_max), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Combined Clustering ",
         x="X Coordinate",
         y="Y Coordinate",
         size="Mean of var1") +
    theme_minimal() +
    #  theme(legend.position = "right")
    theme(legend.position = "none")
  
  plot_d0<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_d0), fill=as.factor(clu_d0), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Non spatial clustering",
         x="X Coordinate",
         y="Y Coordinate",
         size="Mean of var1") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  plot_d1<-ggplot(data2, aes(x=x, y=y)) +
    geom_point(aes(color=as.factor(clu_d1), fill=as.factor(clu_d1), size=var1), shape=21, stroke=1.5) +
    scale_color_manual(values=cluster_colors, name="Resulting Cluster") +
    scale_fill_manual(values=cluster_colors, name="Real Cluster") +
    labs(title="Only spatial clustering",
         x="X Coordinate",
         y="Y Coordinate",
         size="Mean of var1") +
    theme_minimal() +
    #theme(legend.position = "right")
    theme(legend.position = "none")
  
  grid.arrange(plot_clu,nrow =1 )
  
  acc_ch<-sum(data2$cluster==data2$clu_ch)
  acc_max<-sum(data2$cluster==data2$clu_max)
  
  i_ch<-as.numeric(which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2])))
  i_max<-as.numeric(which.max(b$tot))
  #proportion of inertia and weighted inertia
  q_d0_ch<-b$Q0[i_ch]
  q_d0_max<-b$Q0[i_max]
  
  q_d1_ch<-b$Q1[i_ch]
  q_d1_max<-b$Q1[i_max]
  
  wei_ch<-b$tot[i_ch]
  wei_max<-b$tot[i_max]
  
  #normalized inertia, d0 d1 total
  qn_d0_ch<-a$Q0norm[i_ch]
  qn_d0_max<-a$Q0norm[i_max]
  
  qn_d1_ch<-a$Q1norm[i_ch]
  qn_d1_max<-a$Q1norm[i_max]
  
  qn_ch<-a$tot[i_ch]
  qn_max<-a$tot[i_max]
  
  sim_result<-c(n_obs,num_clusters,d,v_sp,v, #input
                a_ch,a_max, # 2 alpha
                q_d0_ch,q_d0_max,q_d1_ch,q_d1_max, # 4 prop inertia
                qn_d0_ch,qn_d0_max,qn_d1_ch,qn_d1_max,# 4 norm prop inertia
                wei_ch,wei_max,qn_ch,qn_max, # 4 weighted and tot norm inertia
                acc_ch,acc_max) # 2 correct classification rate
  sim_res<-rbind(sim_res,sim_result)
}
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max")
######################
sim_res <- data.frame(matrix(ncol = 21, nrow=0))
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max") # 2 correct class rate
# Set the seed for reproducibility
set.seed(12)
for (d in c(0,1)){
  # Define the number of clusters and their probabilities
  n_obs<-150
  num_clusters <- 2
  cluster_probs <- c(0.5,0.5)

  # Define the means and covariance matrices for each cluster
  means_space <- list(
    c(d, d),
    c((-d), (-d))
  )
  v_sp<-0.25
  covs_space <- list(
    matrix(c(1, 0, 0, 1)*v_sp, ncol=2),
    matrix(c(1, 0., 0, 1)*v_sp, ncol=2)
  )
  
  # Define the means and variances for var1 in each cluster
  var1_means <- c(1, 2)
  v<-0.32 #variablity of var1
  var1_sds <- c(1, 1)*v
  
  x <- numeric(0)
  y <- numeric(0)
  var1 <- numeric(0)
  cluster_labels <- numeric(0)
  
  # Initialize the dataframe
  data <- data.frame(x = numeric(0), y = numeric(0), var1 = numeric(0), var2 = numeric(0), cluster = integer(0))
  
  # Simulation loop
  while (nrow(data) < n_obs) {
    # Step 1: Determine cluster assignment
    cluster <- sample(num_clusters, 1, prob = cluster_probs)
    cluster_sp <- sample(num_clusters, 1, prob = cluster_probs)
    #cluster_sp<-cluster
    # Step 2: Simulate position in space
    position <- mvrnorm(1, means_space[[cluster_sp]], covs_space[[cluster_sp]])
    
    # Step 3: Check if position is within [0, 1] x [0, 1]
    # if (all(position >= 0 & position <= 1)) {
    # Step 4: Simulate var1 and var2
    mean_vars <- var1_means[[cluster]]
    sd_vars <- var1_sds[[cluster]]
    var1 <- rnorm(1, mean = mean_vars, sd = sd_vars)
    #var2 <- rnorm(1, mean = mean_vars[2], sd = sd_vars[2])
    
    # Step 5: Append to the dataframe
    data <- rbind(data, data.frame(x = position[1], y = position[2], var1 = var1, cluster = cluster))
    #}
  }
  data$cluster <- as.factor(data$cluster)
  cluster_colors <- c("blue", "orange")
  ggplot(data, aes(x = x, y = y, color = cluster, size = var1)) +
    geom_point() +
    scale_color_manual(values = cluster_colors) +
    labs(
      title = "Plot of Clusters",
      x = "X Coordinate",
      y = "Y Coordinate",
      color = "Cluster"
    ) +
    guides(size = "none") +  # Remove the size legend
    theme_minimal()
  
}
colnames(sim_res)<-c("n_obs","n_clu","d","v_sp","v",# 5 input
                     "a_ch","a_max", # 2 alpha
                     "q_d0_ch","q_d0_max","q_d1_ch","q_d1_max", # 4 prop inertia
                     "qn_d0_ch","qn_d0_max","qn_d1_ch","qn_d1_max",# 4 norm prop inertia
                     "wei_ch","wei_max", "qn_ch","qn_max", # 4 weighted and tot norm inertia
                     "acc_ch","acc_max")

