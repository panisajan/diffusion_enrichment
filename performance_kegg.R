# Load required libraries
library(diffuStats)
library(igraph)
library(pROC)
library(PRROC)

# Load data
load("K_laplacian.Rdata")

kegg <- read.csv("kegg_pathway.csv", header = TRUE)
ppi_data <- read.csv("ppi_with_genes.csv")[, -c(1, 2)]  # Remove unnecessary columns
ppi_data$weight <- ppi_data$combined_score
ppi_data <- ppi_data[, -1]  # Drop the original `combined_score` column

# Create graph
g <- graph_from_data_frame(d = ppi_data, directed = FALSE)
E(g)$weight <- ppi_data$weight

# Parameters
n_iterations <- 100

# Main loop over pathways
for (k in 1:length(kegg$X)) {  # Adjust range as needed
  print(k)
  
  # Extract pathway data
  pathway_id <- kegg$pathway_id[k]
  gene_symbols <- strsplit(kegg$gene_symbols[k], ",\\s*")[[1]]
  namestrue <- intersect(gene_symbols, V(g)$name)
  test <- as.numeric(V(g)$name %in% namestrue)
  names(test) <- V(g)$name
  num_test_data <- sum(test)
  
  if (num_test_data == 0) next  # Skip if no matching genes
  
  # Set iteration limits
  max_iterations <- if (num_test_data > 100) 100 else num_test_data
  percentages <- round(seq(1, 100, length.out = max_iterations) * num_test_data / 100)
  
  # Initialize results dataframe
  results <- data.frame(
    index = integer(max_iterations),
    mean_auc_diff = numeric(max_iterations),
    var_auc_diff = numeric(max_iterations),
    mean_aupr_diff = numeric(max_iterations),
    var_aupr_diff = numeric(max_iterations),
    mean_auc_neigh = numeric(max_iterations),
    var_auc_neigh = numeric(max_iterations),
    mean_aupr_neigh = numeric(max_iterations),
    var_aupr_neigh = numeric(max_iterations),
    mean_auc_neigh2 = numeric(max_iterations),
    var_auc_neigh2 = numeric(max_iterations),
    mean_aupr_neigh2 = numeric(max_iterations),
    var_aupr_neigh2 = numeric(max_iterations),
    percentage_of_seed = numeric(max_iterations),
    num_test_data = integer(max_iterations)
  )
  
  # Iteration over seed percentages
  for (j in seq_len(max_iterations)) {
    num_seeds <- if (num_test_data > 100) percentages[j] else j
    auc_diff_values <- numeric(n_iterations)
    auc_neigh_values <- numeric(n_iterations)
    auc_neigh2_values <- numeric(n_iterations)
    aupr_diff_values <- numeric(n_iterations)
    aupr_neigh_values <- numeric(n_iterations)
    aupr_neigh2_values <- numeric(n_iterations)
    
    # Randomized iterations
    for (i in seq_len(n_iterations)) {
      seed_gene <- sample(namestrue, num_seeds)
      y <- as.numeric(V(g)$name %in% seed_gene)
      
      # Diffusion scores
      F_diffusion <- K %*% y
      F_dif <- F_diffusion[, 1]
      
      # Neighborhood lists
      neighbors_list <- unique(unlist(neighbors(g, seed_gene, mode = "all")))
      neighbors_list2 <- unique(unlist(neighbors(g, neighbors_list, mode = "all")))
      
      # Scores for different methods
      neigh_scores <- as.numeric(V(g)$name %in% c(seed_gene, neighbors_list))
      neigh_scores2 <- as.numeric(V(g)$name %in% c(seed_gene, neighbors_list, neighbors_list2))
      
      # Calculate AUC and AUPR
      roc_diff <- roc(test, F_dif)
      auc_diff_values[i] <- auc(roc_diff)
      aupr_diff_values[i] <- pr.curve(scores.class0 = F_dif[test == 1], scores.class1 = F_dif[test == 0], curve = TRUE)$auc.integral
      
      roc_neigh <- roc(test, neigh_scores)
      auc_neigh_values[i] <- auc(roc_neigh)
      aupr_neigh_values[i] <- pr.curve(scores.class0 = neigh_scores[test == 1], scores.class1 = neigh_scores[test == 0], curve = TRUE)$auc.integral
      
      roc_neigh2 <- roc(test, neigh_scores2)
      auc_neigh2_values[i] <- auc(roc_neigh2)
      aupr_neigh2_values[i] <- pr.curve(scores.class0 = neigh_scores2[test == 1], scores.class1 = neigh_scores2[test == 0], curve = TRUE)$auc.integral
    }
    
    # Store results
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
      percentage_of_seed = (num_seeds / num_test_data) * 100,
      num_test_data = num_test_data
    )
  }
  
  # Save results
  pathway_id_no_colon <- gsub(":", "_", pathway_id)
  write.csv(results, paste0(pathway_id_no_colon, "_kegg_new.csv"), row.names = FALSE)
}
