library(dplyr)
library(ggplot2)
library(sjPlot)
library(patchwork)
colors_cb <- c("#0072B2", "#D55E00")  # two high-contrast colors
cb_colors_11 <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky Blue
  "#009E73",  # Green
  "#D55E00"   # Vermillion / Red-Orange
)

data <- read.csv("Figures and codes/DataTable.csv", header = T)
order_pca <- c("100% Pine", "70% Pine", "60% Pine", "50% Pine", "40% Pine", "30% Pine", "20% Pine", "10% Pine", "100% Oak", "70% Sweetgum","100% Sweetgum")
data$Category <- as.factor(data$Category)
data$Category <- factor(data$Category, 
                            levels = order_pca)
reference_values <- data %>%
  filter(Category %in% c("100% Pine", "100% Oak", "100% Sweetgum")) %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

reference_values


pine_meanMLR  <- reference_values$Mass_Loss_Rate[reference_values$Category == "100% Pine"]
oak_meanMLR  <- reference_values$Mass_Loss_Rate[reference_values$Category == "100% Oak"]
sweetgum_meanMLR  <- reference_values$Mass_Loss_Rate[reference_values$Category == "100% Sweetgum"]

data <- data %>%
  mutate(
    expected_mass_loss_rate =
      (Pine_Fraction      * pine_meanMLR) +
      (Oak_Fraction       * oak_meanMLR) +
      (Sweetgum_Fraction  * sweetgum_meanMLR)
  )

data <- data %>%
  mutate(
    NAE_MLR = Mass_Loss_Rate - expected_mass_loss_rate
  )

ggplot(data, aes(x = Category)) +
  geom_segment(aes(
    xend = Category,
    y = expected_mass_loss_rate,
    yend = Mass_Loss_Rate
  ), color = "gray50") +
  geom_point(aes(y = expected_mass_loss_rate), color = "blue", size = 3) +
  geom_point(aes(y = Mass_Loss_Rate), color = "red", size = 3) +
  labs(
    x = "Treatment",
    y = "Flame Height",
    title = "Observed (red) vs Expected (blue) Across Mixtures"
  ) +
  theme_bw()

ggplot(data, aes(x = expected_mass_loss_rate, y = Mass_Loss_Rate, color = Category)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 3) +
  theme_bw()


# Additive model
m_add <- lm(Bottom_Temperature ~ Pine_Fraction + Oak_Fraction, data = data)
summary(m_add)
# Interaction model (tests non-additivity)
m_int <- lm(Bottom_Temperature ~ Pine_Fraction * Oak_Fraction + Pine_Fraction + Oak_Fraction , data = data)
summary(m_int)

anova(m_add, m_int)


plot_model(m_int, type="int", mdrt.values = "all")


p <- plot_model(m_int, type = "int", mdrt.values = "quart", 
                axis.title = c("Pine Fraction", "Propagation Rate"), 
                title = ""
                ) + 
  theme_minimal(base_size = 14) +          # clean background
  theme(
    axis.title = element_text(),
    legend.title = element_text(),
    legend.position = "right"
  )

p




flamm_metrics <- c(
  "Mass_Loss_Rate",
  "Combustion_Fraction",
  "Bottom_Temperature",
  "Middle_Temperature",
  "Top_Temperature",
  "Propagation_Rate")


models_int <- list()
plots <- list()
model_summaries <- list()


for (metric in flamm_metrics) {
  
  # Interaction model: Pine x Oak (Sweetgum dropped if multicollinearity)
  formula_int <- as.formula(paste0(metric, " ~ Pine_Fraction * Oak_Fraction + Pine_Fraction + Oak_Fraction"))
  
  # Fit model
  mod <- lm(formula_int, data = data)
  models_int[[metric]] <- mod
  
  # Generate realistic fraction combinations
  pine_vals <- seq(0, 1, by = 0.1)
  oak_vals <- seq(0, 0.3, by = 0.3)
  grid <- expand.grid(Pine_Fraction = pine_vals, Oak_Fraction = oak_vals)
  grid <- grid[grid$Pine_Fraction + grid$Oak_Fraction <= 1, ]  # feasible combos only
  
  # Predict with confidence intervals
  pred <- predict(mod, newdata = grid, interval = "confidence", level = 0.95)
  grid$Predicted <- pred[, "fit"]
  grid$CI_lower <- pred[, "lwr"]
  grid$CI_upper <- pred[, "upr"]
  
  # Plot interaction with CI ribbon
  p <- ggplot(grid, aes(x = Pine_Fraction, y = Predicted, color = factor(Oak_Fraction), fill = factor(Oak_Fraction))) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = colors_cb[1:2], name = "Oak Fraction") +
    scale_fill_manual(values = colors_cb[1:2], name = "Oak Fraction") +
    labs(x = "Pine Fraction", y = metric) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(),
      legend.title = element_text(),
      legend.position = "right"
    )
  
  plots[[metric]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot
ggsave("PineOakInteraction.png", combined_plot, bg= "white", dpi = 600)


lapply(model_summaries, function(x) x$coefficients)

####################### SWEETGUM INTERACTION ##########################
models_int <- list()
plots <- list()
model_summaries <- list()


for (metric in flamm_metrics) {
  
  formula_int <- as.formula(paste0(metric, " ~ Pine_Fraction * Sweetgum_Fraction + Pine_Fraction + Sweetgum_Fraction"))
  
  # Fit model
  mod <- lm(formula_int, data = data)
  models_int[[metric]] <- mod
  model_summaries[[metric]] <- summary(mod)
  
  # Generate realistic fraction combinations
  pine_vals <- seq(0, 1, by = 0.1)
  swg_vals <- seq(0, 1, by = 0.3)
  grid <- expand.grid(Pine_Fraction = pine_vals, Sweetgum_Fraction = swg_vals)
  grid <- grid[grid$Pine_Fraction + grid$Sweetgum_Fraction <= 1, ]  # feasible combos only
  
  # Predict with confidence intervals
  pred <- predict(mod, newdata = grid, interval = "confidence", level = 0.95)
  grid$Predicted <- pred[, "fit"]
  grid$CI_lower <- pred[, "lwr"]
  grid$CI_upper <- pred[, "upr"]
  
  # Plot interaction with CI ribbon
  p <- ggplot(grid, aes(x = Pine_Fraction, y = Predicted, color = factor(Sweetgum_Fraction), fill = factor(Sweetgum_Fraction))) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +
    scale_color_manual(values = cb_colors_11 [1:4], name = "Sweetgum_Fraction") +
    scale_fill_manual(values = cb_colors_11 [1:4], name = "Sweetgum_Fraction") +
    labs(x = "Pine Fraction", y = metric) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(),
      legend.title = element_text(),
      legend.position = "right"
    )
  
  plots[[metric]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot
ggsave("PineSweetgumInteraction.png", combined_plot, bg= "white", dpi = 600)


lapply(model_summaries, function(x) x$coefficients)
