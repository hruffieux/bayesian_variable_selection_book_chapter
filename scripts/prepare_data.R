rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

data_dir <- file.path(CORE_DIR, "bayesian_variable_selection_book_chapter/data/")
main_dir <- file.path(CORE_DIR, "bayesian_variable_selection_book_chapter/scripts/")

setwd(main_dir)

require(echoseq)

load(file.path(data_dir, "replicated_data.RData"))

vec_type <- names(list_expr)

list_data <- NULL
for (type in vec_type) {

  replicated_snps <- list_repl_snps[[type]]
  original_expr <- list_expr[[type]]

  rownames(original_expr) <- rownames(replicated_snps)

  n <- nrow(replicated_snps)
  p <- ncol(replicated_snps)
  d <- ncol(original_expr)

  obj_snps <- convert_snps(replicated_snps)
  obj_expr <- convert_phenos(original_expr)

  ## simulates association pattern between the SNPs and the expression levels
  #
  p0 <- 10  # number of active SNPs
  ind_p0 <- sample(1:p, p0, replace = FALSE)

  d0 <- 250  # number of expression levels under genetic control
  ind_d0 <- sample(1:d, d0, replace = FALSE)

  # vec_prob_sh: vector of length p0 providing the probabilities with which each
  # active SNP will be associated with an additional active expression level
  vec_prob_sh <- rbeta(p0, shape1 = 1, shape2 = 2)

  max_tot_pve <- 0.5 # maximum total proportion of outcome variance explained by the SNPs.

  obj_data <- generate_dependence(obj_snps, obj_expr, ind_d0, ind_p0,
                                  vec_prob_sh, family = "gaussian",
                                  max_tot_pve = max_tot_pve)

  stopifnot(all.equal(rownames(obj_data$snps), rownames(obj_data$phenos)))

  list_data <- append(list_data,
                      list(list("snps" = obj_data$snps, "expr" = obj_data$phenos)))
}
names(list_data) <- vec_type

# Check that same set of SNPs and expression levels across conditions.
list_snp_names <- lapply(list_data, function(ll_type) names(ll_type$snps))
stopifnot(length(unique(list_snp_names)) == 1)

list_expr_names <- sapply(list_data, function(ll_type) names(ll_type$expr))
stopifnot(length(unique(list_expr_names)) == 1)

save(list_data, file = file.path(data_dir, "prepared_data.RData"))
