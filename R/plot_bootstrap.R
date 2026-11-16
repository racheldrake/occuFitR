plot_bootstrap <- function(x, test) {

  test <- match.arg(
    test,
    choices = c("mb", "psi_p", "p")
  )

  ggplot2::ggplot(
    x$bootstrap,
    ggplot2::aes(x = .data[[test]])
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = "grey70",
      colour = "white"
    ) +
    ggplot2::geom_vline(
      xintercept = x$statistic[test],
      linetype = 2,
      linewidth = 1
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste(test, "GOF statistic"),
      y = "Frequency",
      title = paste0(
        "Bootstrap GOF test: p = ",
        round(x$p_value[test], 4)
      )
    )
}
