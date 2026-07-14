library(data.table)
library(tidyverse)
library(Deriv)
library(wtdr)
library(grid)
library(gridExtra)
library(flextable)
library(officer)

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
  
  duration_generated_Rx_lengths <- 30
  
  N <- 1000
  
  log_sd <- 0.50
  
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
  
  # bind prevalent users and stoppers
  df <- rbind(results_p, results_s)
  
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
  
  p1 <- predict(rev1, distrx = "disp_time", quantile = 0.9, se.fit = T)
  
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

  p2 <- predict(rev2, distrx = "disp_time", quantile = 0.9, se.fit = T)

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

  p3 <- predict(rev3, distrx = "disp_time", quantile = 0.9, se.fit = T)

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

  p4 <- predict(rev4, distrx = "disp_time", quantile = 0.9, se.fit = T)

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

# Median

# 1) reverse, fixed
m1 <- round(apply(res_rev_fixed[-1], 2, median),3)
# 2) reverse, random (m=5)
m2 <- round(apply(res_rev_random5[-1], 2, median),3)
# 3) reverse, random (m=50)
m3 <- round(apply(res_rev_random50[-1], 2, median),3)
# 4) reverse, saturated
m4 <- round(apply(res_rev_saturated[-1], 2, median),3)

# Mean

# 1) reverse, fixed
mean1 <- round(apply(res_rev_fixed[-1], 2, mean),3)
# 2) reverse, random (m=5)
mean2 <- round(apply(res_rev_random5[-1], 2, mean),3)
# 3) reverse, random (m=50)
mean3 <- round(apply(res_rev_random50[-1], 2, mean),3)
# 4) reverse, saturated
mean4 <- round(apply(res_rev_saturated[-1], 2, mean),3)

# RMSE
data <- list(res_rev_fixed, res_rev_random5, res_rev_random50, res_rev_saturated)

rmse_f <- NULL

for (i in data) {
  
  rmse_vec <- vector()
  
  # logitp
  predicted <- i[, "logitp_est"]
  actual <- rep(2.197, length(predicted))
  rmse <- round(sqrt(mean((actual - predicted)^2)),3)
  
  rmse_vec <- c(rmse_vec, rmse)
  
  # mu 
  predicted <- i[, "mu_est"]
  actual <- rep(log(duration_generated_Rx_lengths), length(predicted))
  rmse <- round(sqrt(mean((actual - predicted)^2)),3)
  
  rmse_vec <- c(rmse_vec, rmse)
  
  # lnsigma
  predicted <- i[, "lnsigma_est"]
  actual <- rep(log(log_sd), length(predicted))
  rmse <- round(sqrt(mean((actual - predicted)^2)),3)
  
  rmse_vec <- c(rmse_vec, rmse)
  
  # dispensation length
  predicted <- i[, "disp_length_est"]
  actual <- rep(exp(log(duration_generated_Rx_lengths) + log_sd*qnorm(0.90)), length(predicted))
  rmse <- round(sqrt(mean((actual - predicted)^2)),3)
  
  rmse_vec <- c(rmse_vec, rmse)
  
  rmse_f <- rbind(rmse_f, rmse_vec)
  
}

# bias 

data <- list(res_rev_fixed, res_rev_random5, res_rev_random50, res_rev_saturated)

bias_f <- NULL

for (i in data) {
  
  bias_vec <- vector()
  
  # logitp
  predicted <- i[, "logitp_est"]
  actual <- rep(2.197, length(predicted))
  bias <- round(mean((predicted - actual)/actual),3)
  
  bias_vec <- c(bias_vec, bias)
  
  # mu 
  predicted <- i[, "mu_est"]
  actual <- rep(log(duration_generated_Rx_lengths), length(predicted))
  bias <- round(mean((predicted - actual)/actual),3)
  
  bias_vec <- c(bias_vec, bias)
  
  # lnsigma
  predicted <- i[, "lnsigma_est"]
  actual <- rep(log(log_sd), length(predicted))
  bias <- round(mean((predicted - actual)/actual),3)
  
  bias_vec <- c(bias_vec, bias)
  
  # dispensation length
  predicted <- i[, "disp_length_est"]
  actual <- rep(exp(log(duration_generated_Rx_lengths) + log_sd*qnorm(0.90)), length(predicted))
  bias <- round(mean((predicted - actual)/actual),3)
  
  bias_vec <- c(bias_vec, bias)
  
  bias_f <- rbind(bias_f, bias_vec)
  
}

# compute relative reduction

# 2) logitp, mu, lnsigma, dispensation length
rel_log_2 <- round((-(m2["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
rel_mu_2 <- round((-(m2["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
rel_ln_2 <- round((-(m2["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
rel_disp_2 <- round((-(m2["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

rel_21 <- c(rel_log_2, rel_mu_2, rel_ln_2, rel_disp_2)

# 3) logitp, mu, lnsigma, dispensation length
rel_log_3 <-round((-(m3["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
rel_mu_3 <- round((-(m3["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
rel_ln_3 <- round((-(m3["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
rel_disp_3 <- round((-(m3["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

rel_31 <- c(rel_log_3, rel_mu_3, rel_ln_3, rel_disp_3)

# 4) logitp, mu, lnsigma, dispensation length
rel_log_4 <-round((-(m4["logitp_se"] - m1["logitp_se"])/m1["logitp_se"])*100, 1)
rel_mu_4 <- round((-(m4["mu_se"] - m1["mu_se"])/m1["mu_se"])*100, 1)
rel_ln_4 <- round((-(m4["lnsigma_se"] - m1["lnsigma_se"])/m1["lnsigma_se"])*100, 1)
rel_disp_4 <- round((-(m4["disp_length_se"] - m1["disp_length_se"])/m1["disp_length_se"])*100, 1)

rel_41 <- c(rel_log_4, rel_mu_4, rel_ln_4, rel_disp_4)

rel_red <- cbind(rel_21, rel_31, rel_41)

# delta relative reduction
del_log_22 <- rel_log_2
del_log_32 <- rel_log_3 - rel_log_2
del_log_43 <- rel_log_4 - rel_log_3
del_mu_22 <- rel_mu_2
del_mu_32 <- rel_mu_3 - rel_mu_2
del_mu_43 <- rel_mu_4 - rel_mu_3
del_ln_22 <- rel_ln_2
del_ln_32 <- rel_ln_3 - rel_ln_2
del_ln_43 <- rel_ln_4 - rel_ln_3
del_disp_22 <- rel_disp_2
del_disp_32 <- rel_disp_3 - rel_disp_2
del_disp_43 <- rel_disp_4 - rel_disp_3

# Produce output table ----

## logitp ----
logitp_ref <- 2.197

# Estimate, median and (mean) 
est_log <- paste0(c(m1["logitp_est"], m2["logitp_est"], m3["logitp_est"], 
                    m4["logitp_est"])," (",
                  c(mean1["logitp_est"], mean2["logitp_est"], 
                    mean3["logitp_est"], mean4["logitp_est"]),")")

# Standard error, median (mean)
se_log <- paste0(c(m1["logitp_se"], m2["logitp_se"], m3["logitp_se"], 
                   m4["logitp_se"]), " (",
                 c(mean1["logitp_se"], mean2["logitp_se"], mean3["logitp_se"], 
                   mean4["logitp_se"]),")")

# SE relative reduction %
rel_red_log <- c("reference", rel_log_2, rel_log_3, rel_log_4)

# Delta relative reduction %
del_red_log <- c("-", del_log_22, del_log_32, del_log_43)

# RMSE and Bias
rmse_log <- rmse_f[, 1]
bias_log <- bias_f[, 1]

## mu ----
mu_ref <- log(duration_generated_Rx_lengths)

est_mu <- paste0(c(m1["mu_est"], m2["mu_est"], m3["mu_est"], m4["mu_est"])," (",
                 c(mean1["mu_est"], mean2["mu_est"], mean3["mu_est"], 
                   mean4["mu_est"]),")")

se_mu <- paste0(c(m1["mu_se"], m2["mu_se"], m3["mu_se"], m4["mu_se"])," (",
                c(mean1["mu_se"], mean2["mu_se"], mean3["mu_se"], 
                  mean4["mu_se"]),")")

rel_red_mu <- c("reference", rel_mu_2, rel_mu_3, rel_mu_4)
del_red_mu <- c("-", del_mu_22, del_mu_32, del_mu_43)

rmse_mu <- rmse_f[, 2]
bias_mu <- bias_f[, 2]

# lnsigma ----
lnsigma_ref <- log(log_sd)

est_ln <- paste0(c(m1["lnsigma_est"], m2["lnsigma_est"], m3["lnsigma_est"], 
                   m4["lnsigma_est"])," (",
                 c(mean1["lnsigma_est"], mean2["lnsigma_est"], 
                   mean3["lnsigma_est"], mean4["lnsigma_est"]),")")

se_ln <- paste0(c(m1["lnsigma_se"], m2["lnsigma_se"], m3["lnsigma_se"], 
                  m4["lnsigma_se"])," (",
                c(mean1["lnsigma_se"], mean2["lnsigma_se"], mean3["lnsigma_se"],
                  mean4["lnsigma_se"]),")")

rel_red_ln <- c("reference", rel_ln_2, rel_ln_3, rel_ln_4)
del_red_ln <- c("-", del_ln_22, del_ln_32, del_ln_43)

rmse_ln <- rmse_f[, 3]
bias_ln <- bias_f[, 3]

# dispensation length ----

disp_ref <- round(exp(log(duration_generated_Rx_lengths) + log_sd * qnorm(0.90)),
                  3)

est_disp <- paste0(c(m1["disp_length_est"], m2["disp_length_est"], 
                     m3["disp_length_est"], m4["disp_length_est"])," (",
                   c(mean1["disp_length_est"], mean2["disp_length_est"], 
                     mean3["disp_length_est"], mean4["disp_length_est"]),")")

se_disp <- paste0(c(m1["disp_length_se"], m2["disp_length_se"], 
                    m3["disp_length_se"], m4["disp_length_se"])," (",
                  c(mean1["disp_length_se"], mean2["disp_length_se"], 
                    mean3["disp_length_se"], mean4["disp_length_se"]),")")

rel_red_disp <- c("reference", rel_disp_2, rel_disp_3, rel_disp_4)
del_red_disp <- c("-", del_disp_22, del_disp_32, del_disp_43)
rmse_disp    <- rmse_f[, 4]
bias_disp    <- bias_f[, 4]

# bind info
col_names <- c(
  "Metric",
  "Reverse WTD\n(1 fixed date)",
  "Reverse WTD\n(5 random dates)",
  "Reverse WTD\n(50 random dates)",
  "Reverse WTD\n(saturated)"
)

param_labels <- c(
  paste0("logitp (ref=", logitp_ref, ")"),
  paste0("mu (ref=", round(mu_ref, 3), ")"),
  paste0("lnsigma (ref=", round(lnsigma_ref, 3), ")"),
  paste0("Dispensation length (ref=", disp_ref, ")")
)

make_block <- function(param_labels, est, se, rel_red, del_red, rmse, bias) {
  data.frame(
    Metric = c(
      param_labels,        
      "Estimate",
      "Standard error",
      "SE relative reduction, %",
      "\u0394 relative reduction, %",
      "RMSE",
      "Bias"
    ),
    C1 = c("", est[1], se[1], rel_red[1], del_red[1], rmse[1], bias[1]),
    C2 = c("", est[2], se[2], rel_red[2], del_red[2], rmse[2], bias[2]),
    C3 = c("", est[3], se[3], rel_red[3], del_red[3], rmse[3], bias[3]),
    C4 = c("", est[4], se[4], rel_red[4], del_red[4], rmse[4], bias[4]),
    stringsAsFactors = FALSE
  )
}

tab_data <- bind_rows(
  make_block(param_labels[1], est_log,  se_log,  rel_red_log,  del_red_log,  rmse_log,  bias_log),
  make_block(param_labels[2], est_mu,   se_mu,   rel_red_mu,   del_red_mu,   rmse_mu,   bias_mu),
  make_block(param_labels[3], est_ln,   se_ln,   rel_red_ln,   del_red_ln,   rmse_ln,   bias_ln),
  make_block(param_labels[4], est_disp, se_disp, rel_red_disp, del_red_disp, rmse_disp, bias_disp)
)

names(tab_data) <- col_names

ft <- flextable(tab_data) |>
  width(j = 1, width = 2.2) |>
  width(j = 2:5, width = 1.5) |>
  
  bold(i = which(tab_data$Metric %in% param_labels), bold = TRUE) |>
  bg(i = which(tab_data$Metric %in% param_labels), bg = "#D9D9D9") |>
  
  bold(part = "header") |>
  align(j = 2:5, align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) |>
  border_inner_h(part = "body", border = fp_border(color = "#AAAAAA", width = 0.5)) |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Arial", part = "all") |>
  set_table_properties(layout = "fixed")

# export
doc <- read_docx() |>
  body_add_par("Table.S2", style = "heading 1") |>
  body_add_flextable(ft)

print(doc, target = "Table.S2_v2.docx")


# plot estimates and SEs to compare relative reductions visually

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
# y1min <- 35
# y1max <- 75
y1min <- 55
y1max <- 95
# y1min <- 80
# y1max <- 120
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

y_fix <- c(rep(exp(log(duration_generated_Rx_lengths) + log_sd*qnorm(0.90)), 4))

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

y_fix <- c(rep(log(duration_generated_Rx_lengths), 4))

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

y_fix <- c(rep(log(log_sd), 4))

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

ggsave("report/Figure_2_v5.png", width = 6300, height = 3200, dpi = 300, units = "px", plot = p)

dev.off()

