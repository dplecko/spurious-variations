

library(data.table)

# load census data
compas_data <- get(data("compas", package = "faircause"))

# get the SFM projection
mdata <- SFM_proj("compas")
mdata

# compute the decomposition
fcb <- fairness_cookbook(as.data.frame(data), X = mdata$X, Z = mdata$Z, 
                         W = mdata$W, Y = mdata$Y, x0 = mdata$x0, x1 = mdata$x1)

autoplot(fcb, decompose = "general")

compas_meas <- summary(fcb)$measures

# obtain the Exp-SE_{x_0} using logistic regression
nboot <- 50
res <- NULL
for (i in seq_len(nboot)) {
  
  data <- compas_data[sample.int(nrow(data), replace = TRUE), ]
  
  logreg <- glm(two_year_recid ~ sex + age + race, data = data, family = "binomial")
  idx0 <- data[[mdata$X]] == mdata$x0
  
  py_x0 <- mean(data[[mdata$Y]][idx0])
  data_x0 <- data
  data_x0[[mdata$X]] <- "White"
  yx0 <- predict(logreg, data_x0, type = "response")
  
  expse_x0 <- py_x0 - mean(yx0)
  
  # obtain the intermediate outcome; P(z1)P(z2 | x0)
  data_x0u1 <- data[idx0, ]
  
  # soft intervention to set Z1 to its marginal P(z1)
  data_x0u1$sex <- sample(data$sex, size = nrow(data_x0u1), replace = TRUE)
  
  y_x0u1 <- predict(logreg, data_x0u1, type = "response")
  
  expse_termI <- py_x0 - mean(y_x0u1)
  expse_termII <- mean(y_x0u1) - mean(yx0)
  
  res <- rbind(res, c(expse_x0, expse_termI, expse_termII, i))
}

res <- as.data.table(res)
res <- setnames(res, names(res), c("Exp-SE", "Exp-SE_I", "Exp-SE_II", "iter"))

res <- melt(res, id.vars = "iter")


ggplot(
  res[, list(mean(value), sd(value)), by = "variable"],
  aes(x = variable, y = V1, fill = variable) 
) +
  geom_col() + 
  geom_errorbar(aes(x = variable, ymin = V1 - 1.96 * V2, 
                    ymax = V1 + 1.96 * V2), width = 0.3) +
  cowplot::theme_minimal_hgrid() +
  xlab("Measure") + ylab("Value") +
  scale_x_discrete(
    labels = c(
      latex2exp::TeX("Exp-SE$_{x_0}(y)$"),
      latex2exp::TeX("Exp-SE$^{\\oslash, U_1}_{x_0}(y)$"),
      latex2exp::TeX("Exp-SE$^{U_1, \\{U_1, U_2\\}}_{x_0}(y)$")
    )
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0)
  )

ggsave("~/Desktop/compas-spurious-decomp.png", width = 6, height = 4)                  
