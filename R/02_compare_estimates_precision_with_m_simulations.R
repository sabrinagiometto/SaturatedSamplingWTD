library(data.table)
library(tidyverse)
library(Deriv)
library(wtdr)
source("R/satwtdttt.R")

# ! remove if statement in ranwtdttt code when using numerical start and end

set.seed(123)

m <- 100

# 1) reverse, fixed

res_rev_fixed <- data.frame(replication = numeric(), 
                                logitp_est = numeric(), mu_est = numeric(), lnsigma_est = numeric(),
                                logitp_se = numeric(), mu_se = numeric(), lnsigma_se = numeric(),
                                disp_length_est = numeric(), disp_length_se = numeric())

# 2) reverse, random 5

res_rev_random5 <- data.frame(replication = numeric(), 
                              logitp_est = numeric(), mu_est = numeric(), lnsigma_est = numeric(),
                              logitp_se = numeric(), mu_se = numeric(), lnsigma_se = numeric(),
                              disp_length_est = numeric(), disp_length_se = numeric())

# 3) reverse, random 50

res_rev_random50 <- data.frame(replication = numeric(), 
                               logitp_est = numeric(), mu_est = numeric(), lnsigma_est = numeric(),
                               logitp_se = numeric(), mu_se = numeric(), lnsigma_se = numeric(),
                               disp_length_est = numeric(), disp_length_se = numeric())

# 4) reverse, saturated

res_rev_saturated <- data.frame(replication = numeric(), 
                                logitp_est = numeric(), mu_est = numeric(), lnsigma_est = numeric(),
                                logitp_se = numeric(), mu_se = numeric(), lnsigma_se = numeric(),
                                disp_length_est = numeric(), disp_length_se = numeric())


system.time({

for (replication in 1:m) {
  
  duration_generated_Rx_lengths <- 60
  
  N <- 1000
  
  log_sd <- 0.75
  
  # set.seed(123)
  
  s <- qnorm(runif(N*0.8))*log_sd + (log(duration_generated_Rx_lengths)+(log_sd^2))
  V <- exp(s)
  u <- runif(N*0.8, 0, V)
  y <- u
  
  results_p <- vector()
  tmp <- y
  i <- 1
  time <- 0
  
  # hist(y, breaks = 100)
  
  while(min(time) <= 730.5){
    
    x <- exp(qnorm(runif(N*0.8))*log_sd + log(duration_generated_Rx_lengths))
    tmp <- data.frame(cbind(tmp, x))
    time <- rowSums(tmp[1:i])
    results_p <- data.frame(cbind(results_p, time))
    
    i <- i + 1
    
  }
  
  # remove time values larger than 365, reshape,
  results_p <- results_p %>%
    mutate(across(everything(.), ~ ifelse(. >= 730.5, NA, .))) %>%
    mutate(pid = sample(1:dim(results_p)[1], N*0.8, replace = F)) %>%
    pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
    drop_na() %>%
    select(!count)
  
  # results_sel <- results_p %>% group_by(pid) %>% slice_min(disp_time)
  # hist(results_sel$disp_time, breaks = 100)
  
  ## WTD incident component (needed in case the ordinary version of model is used, reverse == F) ----
  
  # 20% incident users
  inc <- runif(N*0.2, 0, 730.5)
  # hist(inc, breaks = 100)
  
  results_i <- vector()
  tmp <- inc
  i <- 1
  time <- 0
  
  while(min(time) <= 730.5){
    
    x <- exp(qnorm(runif(N*0.2))*log_sd + log(duration_generated_Rx_lengths))
    tmp <- data.frame(cbind(tmp, x))
    time <- rowSums(tmp[1:i])
    results_i <- data.frame(cbind(results_i, time))
    
    i <- i + 1
    
  }
  
  # remove time values larger than 365, reshape,
  results_i <- results_i %>%
    mutate(across(everything(.), ~ ifelse(. >= 730.5, NA, .))) %>%
    mutate(pid = sample((length(unique(results_p$pid))+1):((length(unique(results_p$pid)))+dim(results_i)[1]), N*0.2, replace = F)) %>%
    pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
    drop_na() %>%
    select(!count)
  
  # results_sel <- results_i %>% group_by(pid) %>% slice_min(disp_time)
  # hist(results_sel$disp_time, breaks = 100)
  
  ## WTD stopping component (needed in case the reverse version of model is used, reverse == T)----
  
  # 20% stoppers
  inc <- runif(N*0.2, 0, 730.5)
  # hist(inc, breaks = 100)
  
  # results_i <- vector()
  results_s <- inc
  tmp <- inc
  i <- 1
  time <- 0
  
  while(max(time) >= 0){
    
    x <- exp(qnorm(runif(N*0.2))*log_sd + log(duration_generated_Rx_lengths))
    tmp <- data.frame(cbind(tmp, x))
    
    time <- apply(tmp, 1, function(tmp) {
      
      tmp[1] - sum(tmp[2:(1+i)])
      
    })
    
    results_s <- data.frame(cbind(results_s, time))
    
    i <- i + 1
    
  }
  
  # remove time values smaller than 0, reshape,
  results_s <- results_s %>% rename(first = results_s)
  
  results_s <- results_s %>%
    mutate(across(everything(.), ~ ifelse(. <= 0, NA, .))) %>%
    mutate(pid = base::sample(((length(unique(results_p$pid)))+(length(unique(results_i$pid)))+1):((length(unique(results_p$pid)))+(length(unique(results_i$pid)))+nrow(results_s)), nrow(results_s), replace = F)) %>%
    pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
    drop_na() %>%
    select(!count)
  
  results_sel <- results_s %>% group_by(pid) %>% slice_max(disp_time)
  # hist(results_sel$disp_time, breaks = 100)
  
  
  # bind prevalent, incident users and stoppers
  # df <- rbind(results_p, results_i, results_s)
  df <- rbind(results_p, results_i)
  # 
  # results_sel <- df %>% group_by(pid) %>% slice_min(disp_time)
  # hist(results_sel$disp_time, breaks = 100)
  
  # run analysis
  
  # 1) reverse, fixed
  
  rev1 <- wtdttt(data = df,
                  disp_time ~ dlnorm(logitp, mu, lnsigma),
                  id = "pid",
                  start = 365.25, 
                  end = 730.5, 
                  reverse = T
  )
  
  p1 <- predict(rev1, distrx = "disp_time", quantile = 0.8, se.fit = T)
  
  est1 <- rev1@coef
  se1 <- c(sqrt(rev1@vcov[1,1]), sqrt(rev1@vcov[2,2]), sqrt(rev1@vcov[3,3]))
  
  disp_est1 <- unique(p1$fit)
  disp_se1 <- unique(p1$se.fit)
  
  out1 <- data.frame(iter = replication, 
                    logitp_est = est1[1], mu_est = est1[2], lnsigma_est = est1[3],
                    logitp_se = se1[1], mu_se = se1[2], lnsigma_se = se1[3],
                    disp_length_est = disp_est1, disp_length_se = disp_se1)
  
  res_rev_fixed <- rbind(res_rev_fixed, out1)
  
  # 2) reverse, random (m=5)

  rev2 <- ranwtdttt(data = df,
                    disp_time ~ dlnorm(logitp, mu, lnsigma),
                    id = "pid",
                    start = 365.25,
                    end = 730.5,
                    reverse = T,
                    nsamp = 5,
                    robust = T
  )

  p2 <- predict(rev2, distrx = "disp_time", quantile = 0.8, se.fit = T)

  est2 <- rev2@coef
  se2 <- c(sqrt(rev2@vcov[1,1]), sqrt(rev2@vcov[2,2]), sqrt(rev2@vcov[3,3]))

  disp_est2 <- unique(p2$fit)
  disp_se2 <- unique(p2$se.fit)

  out2 <- data.frame(iter = replication,
                     logitp_est = est2[1], mu_est = est2[2], lnsigma_est = est2[3],
                     logitp_se = se2[1], mu_se = se2[2], lnsigma_se = se2[3],
                     disp_length_est = disp_est2, disp_length_se = disp_se2)

  res_rev_random5 <- rbind(res_rev_random5, out2)

  # 3) reverse, random (m=50)

  rev3 <- ranwtdttt(data = df,
                    disp_time ~ dlnorm(logitp, mu, lnsigma),
                    id = "pid",
                    start = 365.25,
                    end = 730.5,
                    reverse = T,
                    nsamp = 50,
                    robust = T
  )

  p3 <- predict(rev3, distrx = "disp_time", quantile = 0.8, se.fit = T)

  est3 <- rev3@coef
  se3 <- c(sqrt(rev3@vcov[1,1]), sqrt(rev3@vcov[2,2]), sqrt(rev3@vcov[3,3]))

  disp_est3 <- unique(p3$fit)
  disp_se3 <- unique(p3$se.fit)

  out3 <- data.frame(iter = replication,
                     logitp_est = est3[1], mu_est = est3[2], lnsigma_est = est3[3],
                     logitp_se = se3[1], mu_se = se3[2], lnsigma_se = se3[3],
                     disp_length_est = disp_est3, disp_length_se = disp_se3)

  res_rev_random50 <- rbind(res_rev_random50, out3)

  # 4) reverse, saturated

  rev4 <- satwtdttt(data = df,
                    disp_time ~ dlnorm(logitp, mu, lnsigma),
                    id = "pid",
                    start = 365.25,
                    end = 730.5,
                    robust = T,
                    reverse = F
  )

  p4 <- predict(rev4, distrx = "disp_time", quantile = 0.8, se.fit = T)

  est4 <- rev4@coef
  se4 <- c(sqrt(rev4@vcov[1,1]), sqrt(rev4@vcov[2,2]), sqrt(rev4@vcov[3,3]))

  disp_est4 <- unique(p4$fit)
  disp_se4 <- unique(p4$se.fit)

  out4 <- data.frame(iter = replication,
                     logitp_est = est4[1], mu_est = est4[2], lnsigma_est = est4[3],
                     logitp_se = se4[1], mu_se = se4[2], lnsigma_se = se4[3],
                     disp_length_est = disp_est4, disp_length_se = disp_se4)

  res_rev_saturated <- rbind(res_rev_saturated, out4)
  
 }
  
})

# save results
# save(res_rev_fixed, res_rev_random5, res_rev_random50, res_rev_saturated,
#      file = "res/res_rev_all.RData")


# 1) reverse, fixed
m1 <- apply(res_rev_fixed[-1], 2, median)
# 2) reverse, random (m=5)
m2 <- apply(res_rev_random5[-1], 2, median)
# 3) reverse, random (m=50)
m3 <- apply(res_rev_random50[-1], 2, median)
# 4) reverse, saturated
m4 <- apply(res_rev_saturated[-1], 2, median)

# compute relative reduction

# 2) logitp, mu, lnsigma, dispensation length
round((-(m2["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
round((-(m2["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
round((-(m2["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
round((-(m2["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

# 3) logitp, mu, lnsigma, dispensation length
round((-(m3["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
round((-(m3["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
round((-(m3["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
round((-(m3["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

# 4) logitp, mu, lnsigma, dispensation length
round((-(m4["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
round((-(m4["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
round((-(m4["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
round((-(m4["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

# plot estimates and SEs to comparare relative reductions visually

df <- data.frame(
  Model_variant = factor(c("rWTD, fixed (m = 1)","rWTD, random (m = 5)","rWTD, random (m = 50)","rWTD, saturated (m = 365)")),
  Estimate = c(m1[7], m2[7], m3[7], m4[7],
               m1[2], m2[2], m3[2], m4[2],
               m1[3], m2[3], m3[3], m4[3],
               m1[1], m2[1], m3[1], m4[1]),
  SE = c(m1[8], m2[8], m3[8], m4[8],
         m1[5], m2[5], m3[5], m4[5],
         m1[6], m2[6], m3[6], m4[6],
         m1[4], m2[4], m3[4], m4[4]),
  Parameter = c(rep("dispensation length", 4), rep("mu", 4), rep("lnsigma", 4), rep("logitp", 4))
)

y1min <- 80
y1max <- 120
# y1min <- 90
# y1max <- 130
deltay1 <- y1max-y1min

ratio2 <- df$SE[1]/df$SE[5]
ratio3 <- df$SE[1]/df$SE[9]
ratio4 <- df$SE[1]/df$SE[13]

deltay2 <- deltay1/ratio2
deltay3 <- deltay1/ratio3
deltay4 <- deltay1/ratio4

y2min <- df$Estimate[5]-deltay2/2
y2max <- df$Estimate[5]+deltay2/2

y3min <- df$Estimate[9]-deltay3/2
y3max <- df$Estimate[9]+deltay3/2

y4min <- df$Estimate[13]-deltay4/2
y4max <- df$Estimate[13]+deltay4/2

df_s <- df %>% filter(Parameter=="dispensation length")

# y_fix <- c(rep(112.7548, 4))
y_fix <- c(rep(91.36208, 4))

x_start <- levels(df_s$Model_variant)[1]
x_end   <- tail(levels(df_s$Model_variant), 1)

ggplot(df_s, aes(x = Model_variant, y = Estimate, group = 1)) +
  geom_line(aes(linetype = "Estimated value"), color = "black", linewidth = 0.4) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.1, color = "black") +
  geom_segment(
    aes(x = x_start,
        xend = x_end,
        y = y_fix,
        yend = y_fix,
        linetype = "True value"),
    linewidth = 0.4,
    color = "black"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "True value" = "solid",
      "Estimated value" = "dashed"
    )
  ) +
  labs(
    title = "Dispensation lengths computed as 90th percentile of the IAD",
    x = "Model variant",
    y = "Estimate ± SE"
  ) +
  ylim(c(y1min,y1max)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  )

a <- grid.grab()
# ggsave("report/Figure_2.png", width = 3000, height = 1500, dpi = 300, units = "px")

df_s <- df %>% filter(Parameter=="mu")

y_fix <- c(rep(4.094, 4))

x_start <- levels(df_s$Model_variant)[1]
x_end   <- tail(levels(df_s$Model_variant), 1)

ggplot(df_s, aes(x = Model_variant, y = Estimate, group = Parameter, colour = Parameter)) +
  geom_line(aes(linetype = "Estimated value"), color = "black", linewidth = 0.4) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.1, color = "black") +
  geom_segment(
    aes(x = x_start,
        xend = x_end,
        y = y_fix,
        yend = y_fix,
        linetype = "True value"),
    linewidth = 0.4,
    color = "black"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "True value" = "solid",
      "Estimated value" = "dashed"
    )
  ) +
  labs(
    title = "Estimated mu",
    x = "Model variant",
    y = "Estimate ± SE"
  ) +
  ylim(c(y2min,y2max)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  )

b <- grid.grab()

df_s <- df %>% filter(Parameter=="lnsigma")

# y_fix <- c(rep(-0.2876821, 4))
y_fix <- c(rep(-0.6931, 4))

x_start <- levels(df_s$Model_variant)[1]
x_end   <- tail(levels(df_s$Model_variant), 1)

ggplot(df_s, aes(x = Model_variant, y = Estimate, group = Parameter, colour = Parameter)) +
  geom_line(aes(linetype = "Estimated value"), color = "black", linewidth = 0.4) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.1, color = "black") +
  geom_segment(
    aes(x = x_start,
        xend = x_end,
        y = y_fix,
        yend = y_fix,
        linetype = "True value"),
    linewidth = 0.4,
    color = "black"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "True value" = "solid",
      "Estimated value" = "dashed"
    )
  ) +
  labs(
    title = "Estimated lnsigma",
    x = "Model variant",
    y = "Estimate ± SE"
  ) +
  ylim(c(y3min,y3max)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  )

c <- grid.grab()

df_s <- df %>% filter(Parameter=="logitp")

y_fix <- c(rep(2.197, 4))

x_start <- levels(df_s$Model_variant)[1]
x_end   <- tail(levels(df_s$Model_variant), 1)

ggplot(df_s, aes(x = Model_variant, y = Estimate, group = Parameter, colour = Parameter)) +
  geom_line(aes(linetype = "Estimated value"), color = "black", linewidth = 0.4) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.1, color = "black") +
  geom_segment(
    aes(x = x_start,
        xend = x_end,
        y = y_fix,
        yend = y_fix,
        linetype = "True value"),
    linewidth = 0.4,
    color = "black"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "True value" = "solid",
      "Estimated value" = "dashed"
    )
  ) +
  labs(
    title = "Estimated logitp",
    x = "Model variant",
    y = "Estimate ± SE"
  ) +
  ylim(c(y4min,y4max)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 12))
  )

d <- grid.grab()

p <- grid.arrange(a,b,c,d, ncol = 2)

ggsave("report/Figure_2.png", width = 6300, height = 3200, dpi = 300, units = "px", plot = p)

dev.off()

