# harmonic

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

x <- seq(0, 1, length.out = 101)

# histograms
# greens 2 # https://www.color-hex.com/color-palette/5016
xpale <- "#F0F7DA"
xlight <- "#C9DF8A"
xmid <- "#77AB59"
xdark <- "#36802D"
xdata <- "#043927" # https://graf1x.com/shades-of-green-color-palette-html-hex-rgb-code/
xaxis <- "#999999"
xgrey <- "grey"

n <- 10
df <- vector("list", n)
set.seed(123)
for (i in 1:n){
  
  y0 <- i/5 - 0.1
  y1 <- y0 + 0.0 * x
  y2 <- y1 + 0.2 * sin(pi * x) * rnorm(1, 0, 0.4)
  y3 <- y2 + 0.2 * sin(2 * pi * x) * rnorm(1, 0, 0.2)
  
  df[[i]] <- tibble(i, x, Constant = y0, Linear = y1, First = y2, Second = y3) %>% 
    gather("Term", "y", -i, -x) %>% 
    mutate(Term = factor(Term, levels = c("Constant", "Linear", "First", "Second")))

}
df <- bind_rows(df) %>% 
  mutate(alpha = case_when(
    Term == "Constant" ~ 0,
    Term == "Linear" ~ 0,
    Term == "First" ~ 0.0,
    Term == "Second" ~ as.numeric(i)/as.numeric(n)
  )) %>% 
  filter(Term %in% c("Constant", "Second"))

ph <- ggplot(df) +
  labs(title = "Harmonic Curve Examples", y = expression(c[t])) +
  guides(colour = "none", linetype = "none") +
  theme_cowplot() +
  panel_border(colour = "black") +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c(xmid, xdark)) +
  geom_line(mapping = aes(x = x, y = y, colour = Term == "Constant", linetype = Term == "Constant", group = interaction(Term, i)))
  # geom_line(mapping = aes(x = x, y = y, linetype = Term, alpha = alpha, group = interaction(Term, i)))

print(ph)

save_plot("harmonic.png", ph, base_asp = 1, dpi = 600)
