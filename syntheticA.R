
library(data.table)
library(ggplot2)
library(latex2exp)

# synthetic A


true_and_est <- function(l1, l2, l3) {
  
  a <- 0.2
  b <- 0.3
  c1 <- d1 <- l1 
  c2 <- d2 <- l2
  c3 <- d3 <- l3
  beta <- 0.2
  
  synth_A <- function(n) {
    
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rbinom(n, 1, 0.4 + a * z1)
    z3 <- rbinom(n, 1, 0.3 + b * z1 * z2)
    x <- rbinom(n, 1, 0.2 + c1 * z1 + c2 * z2 + c3 * z3)
    y <- rbinom(n, 1, 0.1 + beta * x + d1 * z1 + d2 * z2 + d3 * z3)
    
    data.table(z1, z2, z3, x, y)
  }
  
  {
    # compute the ground truth via rejection sampling
    n <- 10^6
    t_fr <- synth_A(n)
    py_x_t <- t_fr[x == 1, mean(y)]
    t_fr[, yx1 := rbinom(nrow(t_fr), 1, 0.1 + beta + d1 * z1 + d2 * z2 + d3 * z3)]
    pyx_t <- mean(t_fr$yx1)
    
    # use rejection sampling for intm1
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rbinom(n, 1, 0.4 + a * z1)
    z3 <- rbinom(n, 1, 0.3 + b * z1 * z2)
    x <- rbinom(n, 1, 0.2 + c1 * z1 + c2 * z2 + c3 * z3)
    while(any(x == 0)) {
      
      idx <- x == 0
      # cat(sum(idx), "X = 0 values remaining.\n")
      z2[idx] <- rbinom(sum(idx), 1, 0.4 + a * z1[idx])
      z3[idx] <- rbinom(sum(idx), 1, 0.3 + b * z1[idx] * z2[idx])
      x[idx] <- rbinom(sum(idx), 1, 0.2 + c1 * z1[idx] + c2 * z2[idx] + c3 * z3[idx])
    }
    y_intm1 <- rbinom(n, 1, 0.1 + beta * x + d1 * z1 + d2 * z2 + d3 * z3)
    intm1_t <- mean(y_intm1)
    
    # use rejection sampling for intm2
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rbinom(n, 1, 0.4 + a * z1)
    z3 <- rbinom(n, 1, 0.3 + b * z1 * z2)
    x <- rbinom(n, 1, 0.2 + c1 * z1 + c2 * z2 + c3 * z3)
    while(any(x == 0)) {
      
      idx <- x == 0
      # cat(sum(idx), "X = 0 values remaining.\n")
      z3[idx] <- rbinom(sum(idx), 1, 0.3 + b * z1[idx] * z2[idx])
      x[idx] <- rbinom(sum(idx), 1, 0.2 + c1 * z1[idx] + c2 * z2[idx] + c3 * z3[idx])
    }
    y_intm2 <- rbinom(n, 1, 0.1 + beta * x + d1 * z1 + d2 * z2 + d3 * z3)
    intm2_t <- mean(y_intm2)
    
    # compute the ground truth effects
    se_z1_t <- py_x_t - intm1_t
    se_z2_t <- intm1_t - intm2_t
    se_z3_t <- intm2_t - pyx_t
    
    gt <- data.frame(se1 = se_z1_t, se2 = se_z2_t, se3 = se_z3_t,
                     seed = "true")
    # gt <- data.frame(se1 = intm1_t, se2 = intm2_t, se3 = pyx_t,
    #                  seed = "true")
  }
  
  # get y_x and yx
  nrep <- 50
  res <- NULL
  n <- 10^4
  for (seed in seq_len(nrep)) {
    
    # set.seed(seed)
    dt <- synth_A(10^4)
    py_x <- dt[x == 1, mean(y)]
    sl1 <- merge(
      dt[, list(pz = .N / nrow(dt)), by = c("z1", "z2", "z3")],
      dt[x == 1, list(py_xz = mean(y)), by = c("x", "z1", "z2", "z3")],
      by = c("z1", "z2", "z3")
    )
    
    pyx <- sum(sl1$py_xz * sl1$pz)
    
    # intermediate quantity 1 \sum_z p(y | x,z)p(z1)p(z2,z3 | z1, x)
    intm1 <- 0
    for (z1s in c(0, 1)) {
      
      for (z2s in c(0, 1)) {
        
        for (z3s in c(0, 1)) {
          
          intm1 <- intm1 +
            mean(dt$z1 == z1s) * 
            dt[z1 == z1s & x == 1, mean(z2 == z2s & z3 == z3s)] *
            sl1[z1 == z1s & z2 == z2s & z3 == z3s]$py_xz
        }
      }
    }
    
    # intermediate quantity 2 \sum_z p(y | x,z)p(z1,z2)p(z3 | z1,z2,x)
    intm2 <- 0
    for (z1s in c(0, 1)) {
      
      for (z2s in c(0, 1)) {
        
        for (z3s in c(0, 1)) {
          
          intm2 <- intm2 +
            mean(dt$z1 == z1s & dt$z2 == z2s) * 
            dt[z1 == z1s & z2 == z2s & x == 1, mean(z3 == z3s)] *
            sl1[z1 == z1s & z2 == z2s & z3 == z3s]$py_xz
        }
      }
    }
    
    # compute the effects
    se_z1 <- py_x - intm1
    se_z2 <- intm1 - intm2
    se_z3 <- intm2 - pyx
    
    # save the run
    res <- rbind(res, data.frame(se1 = se_z1, se2 = se_z2, se3 = se_z3,
                                 seed = seed))
    # res <- rbind(res, data.frame(se1 = intm1, se2 = intm2, se3 = pyx, 
    #                              seed = seed))
    # cat("\r", seed)
  }
  
  # get a plot
  res <- as.data.table(res)
  est <- melt(res, id.vars = "seed")[, list(value = mean(value), sd = sd(value)), 
                                     by = "variable"]
  est[, seed := "estimate"]
  
  tr <- melt(as.data.table(gt), id.vars = "seed")
  tr[, sd := 0]
  rbind(est, tr)
}

exp_space <- rbind(
  expand.grid(l1 = seq(0, 0.2, 0.01), l2 = 0.2, l3 = 0.2, dir = "l1"),
  expand.grid(l1 = 0.2, l2 = seq(0, 0.2, 0.01), l3 = 0.2, dir = "l2"),
  expand.grid(l1 = 0.2, l2 = 0.2, l3 = seq(0, 0.2, 0.01), dir = "l3")
)
exp_space <- cbind(exp_space, x_ax = rep(seq(0, 0.2, 0.01), 3))

res <- NULL
set.seed(2023)
for (i in seq_len(nrow(exp_space))) {
  
  s <- true_and_est(l1 = exp_space$l1[i], l2 = exp_space$l2[i], 
                    l3 = exp_space$l3[i])
  s$x <- exp_space$x_ax[i]
  s$dir <- exp_space$dir[i]
  res <- rbind(res, s) 
  cat("\r", i)
}
res$dir <- factor(res$dir)

ggplot(
  res[seed == "estimate"],
  aes(x = x, color = variable, y = value)
) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin = value - 1.96 * sd, 
                  ymax = value + 1.96 * sd,
                  fill = variable),
              linewidth=0, alpha = 0.2,) +
  geom_point(data = res[seed == "true"], 
             aes(x = x, y = value, color = variable)) +
  theme_minimal() +
  ylab("Effect Value") +
  xlab("Coefficient Size") +
  scale_fill_discrete(
    name = "Effect",
    labels = c(TeX("$Exp-SE^{U_1}$"), TeX("$Exp-SE^{U_2}$"), 
               TeX("$Exp-SE^{U_3}$"))
  ) +
  scale_color_discrete(
    name = "Effect",
    labels = c(TeX("$Exp-SE^{U_1}$"), TeX("$Exp-SE^{U_2}$"), 
               TeX("$Exp-SE^{U_3}$"))
  ) +
  facet_grid(
    cols = vars(dir),
    labeller = labeller(dir = c("l1" = "Direction 1", 
                                "l2" = "Direction 2", 
                                "l3" = "Direction 3"))
  ) +
  # facet_grid(
  #   cols = vars(dir),
  #   labeller = labeller(dir = c("l1" = TeX("Direction $\\lambda_1$"), 
  #                               "l2" = TeX("Direction $\\lambda_2$"), 
  #                               "l3" = TeX("Direction $\\lambda_3$")))
  # ) +
  theme(legend.position = "bottom",
        legend.box.background = element_rect())

ggsave("syntheticA.png", width = 8, height = 3.5, bg = "white")
