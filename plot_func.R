p1 <- function(data, start_date) {
 data %>% 
  dplyr::filter(date >= as.Date(start_date)) %>%
  ggplot2::ggplot(ggplot2::aes(x=date,y=med)) +
  ggplot2::geom_line(ggplot2::aes(col = model)) +
  ggplot2::geom_ribbon(ggplot2::aes(fill = model, ymin = LQ, ymax = UQ),alpha=0.2) + 
  cowplot::theme_cowplot() +
  ggplot2::geom_hline(yintercept = 1, lty = 2) +
  ggplot2::ylab("Effective reproduction\nnumber") +
  ggplot2::xlab("Date") +
  ggplot2::scale_x_date(date_breaks = "1 week", date_labels = "%b %d")
}

p2 <- function(data, start_date) {
 data %>% 
  dplyr::filter(dates >= as.Date(start_date)) %>%
  ggplot2::ggplot() +
  ggplot2::geom_bar(ggplot2::aes(x = dates, y = imported+local), stat = "identity") + 
  cowplot::theme_cowplot() +
  ggplot2::ylab("Cases") + 
  ggplot2::xlab("") +
  ggplot2::scale_x_date(date_breaks = "1 week", date_labels = "%b %d")
}

p3 <- function(data) {
 data %>%
  tidyr::gather(key = "model", value = "freq", -t) %>%
  dplyr::filter(t < 20) %>%
  ggplot2::ggplot(ggplot2::aes(x = t, y = freq, col = model)) +
  ggplot2::geom_line() +
  cowplot::theme_cowplot() +
  ggplot2::xlab("Time") +
  ggplot2::ylab("Frequency") +
  ggplot2::scale_color_discrete(guide = "none")
}