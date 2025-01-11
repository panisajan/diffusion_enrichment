library(diffuStats)
library(igraph)
library(pROC)
library(PRROC)

# Load pre-computed Laplacian matrix and dataset files
load("K_laplacian.Rdata")
GO_data <- read.csv("go_pathway_withgene.csv")[,-1]
ppi_data <- read.csv("ppi_with_genes.csv")[,-c(1, 2)]

# Filter GO data based on level
GO_data <- subset(GO_data, GO_data$Level <= 3)

# Prepare PPI data and create the graph
ppi_data$weight <- ppi_data$combined_score
ppi_data <- ppi_data[, -1]
g <- graph_from_data_frame(d = ppi_data, directed = FALSE, vertices = NULL)
E(g)$weight <- ppi_data$weight

# Set number of iterations for the simulation
n_iterations <- 100

# Loop over each pathway in GO_data
for (k in 1:nrow(GO_data)) {
  print(paste("Processing pathway:", k))
  
  # Extract pathway information
  pathway_id <- GO_data$go_id[k]
  namestrue <- strsplit(GO_data$gene_name[k], ",\\s*")[[1]]
  
  # Create a binary test vector indicating which genes are in the current pathway
  test <- ifelse(V(g)$name %in% namestrue, 1, 0)
  names(test) <- V(g)$name
  num_test_data <- sum(test)
  
  # Pre-allocate a dataframe to store results
  results <- data.frame(
    index = integer(num_test_data),
    mean_auc_diff = numeric(num_test_data),
    var_auc_diff = numeric(num_test_data),
    mean_aupr_diff = numeric(num_test_data),
    var_aupr_diff = numeric(num_test_data),
    mean_auc_neigh = numeric(num_test_data),
    var_auc_neigh = numeric(num_test_data),
    mean_aupr_neigh = numeric(num_test_data),
    var_aupr_neigh = numeric(num_test_data),
    mean_auc_neigh2 = numeric(num_test_data),
    var_auc_neigh2 = numeric(num_test_data),
    mean_aupr_neigh2 = numeric(num_test_data),
    var_aupr_neigh2 = numeric(num_test_data),
    percentage_of_seed = numeric(num_test_data),
    num_test_data = integer(num_test_data)
  )
  
  # Loop over the range of gene set sizes
  for (j in seq(from = 1, to = num_test_data)) {
    # Initialize vectors for AUC and AUPR values
    auc_diff_values <- numeric(n_iterations)
    auc_neigh_values <- numeric(n_iterations)
    auc_neigh2_values <- numeric(n_iterations)
    aupr_diff_values <- numeric(n_iterations)
    aupr_neigh_values <- numeric(n_iterations)
    aupr_neigh2_values <- numeric(n_iterations)
    
    # Perform the simulation for the specified number of iterations
    for (i in 1:n_iterations) {
      seed_gene <- sample(namestrue, j)  # Randomly select seed genes
      y <- ifelse(V(g)$name %in% seed_gene, 1, 0)
      
      # Perform diffusion calculation
      F_diffusion <- K %*% y
      F_dif <- F_diffusion[, 1]
      
      # Calculate Neighborhood Lists (1st and 2nd level neighbors)
      neighbors_list <- unique(unlist(lapply(seed_gene, function(gene) {
        V(g)$name[neighbors(g, gene, mode = "all")]
      })))
      neighbors_list2 <- unique(unlist(lapply(neighbors_list, function(gene) {
        V(g)$name[neighbors(g, gene, mode = "all")]
      })))
      
      # Assign neighbor scores (1st and 2nd level)
      neigh_scores <- ifelse(V(g)$name %in% c(seed_gene, neighbors_list), 1, 0)
      neigh_scores2 <- ifelse(V(g)$name %in% c(seed_gene, neighbors_list, neighbors_list2), 1, 0)
      
      # Calculate AUC and AUPR for diffusion scores
      roc_diff <- roc(test, F_dif)
      auc_diff_values[i] <- auc(roc_diff)
      
      pr_dif <- pr.curve(scores.class0 = F_dif[test == 1], scores.class1 = F_dif[test == 0], curve = TRUE)
      aupr_diff_values[i] <- pr_dif$auc.integral
      
      # Calculate AUC and AUPR for neighbor scores (1st level)
      roc_neigh <- roc(test, neigh_scores)
      auc_neigh_values[i] <- auc(roc_neigh)
      
      pr_neigh <- pr.curve(scores.class0 = neigh_scores[test == 1], scores.class1 = neigh_scores[test == 0], curve = TRUE)
      aupr_neigh_values[i] <- pr_neigh$auc.integral
      
      # Calculate AUC and AUPR for neighbor scores (2nd level)
      roc_neigh2 <- roc(test, neigh_scores2)
      auc_neigh2_values[i] <- auc(roc_neigh2)
      
      pr_neigh2 <- pr.curve(scores.class0 = neigh_scores2[test == 1], scores.class1 = neigh_scores2[test == 0], curve = TRUE)
      aupr_neigh2_values[i] <- pr_neigh2$auc.integral
    }
    
    # Store results for current gene set size (j)
    results[j, ] <- data.frame(
      index = j,
      mean_auc_diff = mean(auc_diff_values),
      var_auc_diff = var(auc_diff_values),
      mean_aupr_diff = mean(aupr_diff_values),
      var_aupr_diff = var(aupr_diff_values),
      mean_auc_neigh = mean(auc_neigh_values),
      var_auc_neigh = var(auc_neigh_values),
      mean_aupr_neigh = mean(aupr_neigh_values),
      var_aupr_neigh = var(aupr_neigh_values),
      mean_auc_neigh2 = mean(auc_neigh2_values),
      var_auc_neigh2 = var(auc_neigh2_values),
      mean_aupr_neigh2 = mean(aupr_neigh2_values),
      var_aupr_neigh2 = var(aupr_neigh2_values),
      percentage_of_seed = (j / num_test_data) * 100,
      num_test_data = num_test_data
    )
  }
  
  # Save the results to a CSV file
  pathway_id_no_colon <- gsub(":", "_", pathway_id)
  full_filename <- paste0(pathway_id_no_colon, '_results_GO.csv')
  write.csv(results, full_filename, row.names = FALSE)
}
