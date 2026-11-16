#' @export
print.summary.occu_gof <- function(x, ...) {

  cat("\nGoodness-of-fit tests for occupancy model\n")
  cat("------------------------------------------\n\n")

  cat("Number of sites:", x$n_sites, "\n")
  cat("Number of visits:", x$n_visits, "\n")
  cat("Bootstrap replicates:", x$n_boot, "\n\n")

  cat("GOF test results:\n\n")

  results <- x$results

  results$p_value <- format.pval(
    results$p_value,
    digits = 3,
    eps = 0.001
  )

  print(results, row.names = FALSE)

  invisible(x)
}
