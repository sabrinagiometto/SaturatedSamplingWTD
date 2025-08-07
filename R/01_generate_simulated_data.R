library(data.table)
library(tidyverse)
library(Deriv)


# First try: Generating a small dataset for one person with dispensations uniformly distributed over 5 years ----

set.seed(123)

# generate data from a uniform distribution over 5 y
dt <- data.table(disp_time = round(runif(30, 0, 365.25*5),0))

# keep only unique dispensing times
dt <- dt[!duplicated(disp_time)]

# add an artificial last observation to include index dates after the last observed dispensation until the end of the window
dt <- rbind(dt, data.table(disp_time = round(5*365.25)))

# order by dispensing time
dt <- dt[order(disp_time)]

# compute time from consecutive dispensations
dt <- dt[, time_to_next := shift(disp_time, type = "lead") - disp_time]

# remove the last row as it doesn't have time to the next one
dt <- dt[!is.na(time_to_next)]

# expand disp_time time_to_next times
dt_e <- dt[rep(1:.N, times = time_to_next), .(disp_time)]

# compute time since last dispensation, considering all the possible index dates until the next dispensation
dt_e <- dt_e[, days_since_last := (seq(.N)-1)+0.5, by = disp_time]

# count days_since_last
reps <- dt_e[days_since_last < 365.25, .(reps = .N), by = days_since_last]

# same histograms when plotting the expanded version of days_since_last (a) and the one with weights (b)
hist(dt_e$days_since_last, breaks = 50) # a

ggplot(reps, aes(x=days_since_last, weight = reps)) +
  geom_histogram(bins = 50) # b


# Increase the number of dispensings sampled per individual, and add more than one individual ----

# generate data from a uniform distribution over 5 y
set.seed(123)

ids <- sample(0:50, 50, replace = FALSE) # 2 individuals
n_per_id <- 50 # 50 dispensings per individual

dt_list <- lapply(ids, function(i) {
  data.table(
    id = i,
    disp_time = round(runif(n_per_id, 0, 365.25), 0)
  )
})

dt <- rbindlist(dt_list)

rm(dt_list)

# keep only unique dispensing times
dt <- dt[, .SD[!duplicated(disp_time)], by = id]

# add an artificial last observation to include index dates after the last observed dispensation until the end of the window
add_artif_last <- dt %>%
                   distinct(id) %>% 
                   mutate(disp_time = round(5*365.25))

dt <- rbind(dt, add_artif_last)

# order by dispensing time
dt <- dt[order(id,disp_time)]

# compute time from consecutive dispensations
dt <- dt[, time_to_next := shift(disp_time, type = "lead") - disp_time , by = id]

# remove last row as it doesn't have time to the next one
dt <- dt[!is.na(time_to_next)]

# expand disp_time time_to_next times
dt_e <- dt[rep(1:.N, times = time_to_next), .(id, disp_time)]

# compute time since last dispensation, considering all the possible index dates until the next dispensation
dt_e <- dt_e[, days_since_last := (seq(.N)-1)+0.5, by = disp_time]

# count days_since_last
reps <- dt_e[days_since_last < 365.25, .(reps = .N), by = days_since_last]

# same histograms when plotting the expanded version of days_since_last (a) and the one with weights (b)
hist(dt_e$days_since_last, breaks = 50) # a

ggplot(reps, aes(x=days_since_last, weight = reps)) +
  geom_histogram(bins = 300) # b

# Add potential index dates before the first dispensing (between start and end), as we are currently truncating the observations ----

set.seed(123)

start <- 1095.75
end <- 1461

# generate data from a uniform distribution over 5 y
dt <- data.table(disp_time = round(runif(30, 0, 365.25*5),0))

# keep only unique dispensing times
dt <- dt[!duplicated(disp_time)]


## 1) between first dispensation and last dispensation ----

# add an artificial last observation to include index dates after the last observed dispensation until the end of the window
dt <- rbind(dt, data.table(disp_time = end))

# order by dispensing time
dt <- dt[order(disp_time)]

dt1 <- dt %>% filter(disp_time >= start & disp_time <= end)

# compute time from consecutive dispensations
dt1 <- dt1 %>% 
  mutate(time_to_next = shift(disp_time, type = "lead") - disp_time)

# remove the last row as it doesn't have time to the next one
dt1 <- dt1[!is.na(time_to_next)]

# expand disp_time time_to_next times
dt_e <- dt1[rep(1:.N, times = time_to_next), .(disp_time)]

# compute time since last dispensation, considering all the possible index dates until the next dispensation
dt_e <- dt_e[, days_since_last := (seq(.N)-1)+0.5, by = disp_time]

dt_e_post <- dt_e %>% filter(disp_time >= start & disp_time <= end)

## 2) before first dispensation ----

# select first dispensation within the sampling window, i.e., start-end window
first_disp <- as.numeric(dt %>% filter(disp_time >= start & disp_time <= end) %>% summarise(first_disp = min(disp_time)))

# select dispensations in the moving window, i.e., 365 days before the first dispensation within the sampling window
# and keep the last dispensation within the moving window
dt_pre_first_disp <- dt %>% filter(disp_time < first_disp & disp_time >= first_disp-365.25) %>% slice_max(disp_time)

# bind this last dispensation within the moving window with the first dispensation within the sampling window
dt_pre_first_disp <- rbind(dt_pre_first_disp, data.table(disp_time = first_disp))

# compute time to next dispensation
dt_pre_first_disp <- dt_pre_first_disp %>% 
                           mutate(time_to_next = shift(disp_time, type = "lead") - disp_time)

# remove the last row as it doesn't have time to the next one
dt_pre_first_disp <- dt_pre_first_disp[!is.na(time_to_next)]

# expand disp_time time_to_next times
dt_e_pre_first_disp <- dt_pre_first_disp[rep(1:.N, times = time_to_next), .(disp_time)]

# compute time since last dispensation, considering all the potential index dates until the next dispensation
dt_e_pre_first_disp <- dt_e_pre_first_disp[, days_since_last := (seq(.N)-1)+0.5, by = disp_time]

# select last dispensation within the moving window
max_disp_last_year <- as.numeric(dt %>% filter(disp_time < first_disp & disp_time >= first_disp-365.25) %>% summarise(last_disp = max(disp_time)))

# compute the time interval between this last dispensation within the moving window and the start of sampling window
min_diff <- start-max_disp_last_year

# drop days_since_last smaller than the min_diff interval computed (it means to remove those potential index dates prior to the start of the sampling window)
dt_e_pre_first_disp <- dt_e_pre_first_disp %>% filter(days_since_last > min_diff) 

## 3) bind time window ----
# AFTER FIRST DISPENSING, i.e., between first dispensing (within the sampling window) and end, AND BEFORE FIRST DISPENSING, i.e., between first dispensing (within the sampling window) and last observation (within the moving window) 
dt_e <- rbind(dt_e_pre_first_disp, dt_e_post)


# Generate data from a lognormal instead of a uniform distribution ----

## WTD prevalent component ----

duration_generated_Rx_lengths <- 60

set.seed(123)

s <- qnorm(runif(1000))*0.5 + (log(duration_generated_Rx_lengths)+(0.5^2))
V <- exp(s)
u <- runif(1000, 0, V)
y <- u

results_p <- vector()
tmp <- y
i <- 1
time <- 0

hist(y, breaks = 100)

while(min(time) <= 730.5){
  
  x <- exp(qnorm(runif(1000))*0.5 + log(duration_generated_Rx_lengths))
  tmp <- data.frame(cbind(tmp, x))
  time <- rowSums(tmp[1:i])
  results_p <- data.frame(cbind(results_p, time))
  
  i <- i + 1
  
}

# remove time values larger than 365, reshape,
results_p <- results_p %>%
  mutate(across(everything(.), ~ ifelse(. >= 730.5, NA, .))) %>%
  mutate(pid = sample(1:dim(results_p)[1], 1000, replace = F)) %>%
  pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
  drop_na() %>%
  select(!count)

results_sel <- results_p %>% group_by(pid) %>% slice_min(disp_time)
hist(results_sel$disp_time, breaks = 100)

## WTD incident component (needed in case the ordinary version of model is used, reverse == F) ----

# 10% incident users
inc <- runif(1000, 0, 730.5)
hist(inc, breaks = 100)

results_i <- vector()
tmp <- inc
i <- 1
time <- 0

while(min(time) <= 730.5){
  
  x <- exp(qnorm(runif(1000))*0.5 + log(duration_generated_Rx_lengths))
  tmp <- data.frame(cbind(tmp, x))
  time <- rowSums(tmp[1:i])
  results_i <- data.frame(cbind(results_i, time))
  
  i <- i + 1
  
}

# remove time values larger than 365, reshape,
results_i <- results_i %>%
  mutate(across(everything(.), ~ ifelse(. >= 730.5, NA, .))) %>%
  mutate(pid = sample(dim(results_p)[1]+1:dim(results_p)[1]+1+dim(results_i)[1], 1000, replace = F)) %>%
  pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
  drop_na() %>%
  select(!count)

results_sel <- results_i %>% group_by(pid) %>% slice_min(disp_time)
hist(results_sel$disp_time, breaks = 100)

## WTD stopping component (needed in case the reverse version of model is used, reverse == T)----

# 10% stoppers
inc <- runif(1000, 0, 730.5)
hist(inc, breaks = 100)

# results_i <- vector()
results_s <- inc
tmp <- inc
i <- 1
time <- 0

while(max(time) >= 0){
  
  x <- exp(qnorm(runif(1000))*0.5 + log(duration_generated_Rx_lengths))
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
  mutate(pid = base::sample((nrow(results_p)+1):(nrow(results_p)+1+nrow(results_s)), nrow(results_s), replace = F)) %>%
  pivot_longer(cols = -c(pid), names_to = "count", values_to = "disp_time") %>%
  drop_na() %>%
  select(!count)

results_sel <- results_s %>% group_by(pid) %>% slice_max(disp_time)
hist(results_sel$disp_time, breaks = 100)


# bind prevalent, incident users and stoppers
results_a <- rbind(results_p, results_i, results_s)

results_sel <- results_a %>% group_by(pid) %>% slice_min(disp_time)
hist(results_sel$disp_time, breaks = 100)

# save(results_a, file = file.path("extdata", "results_a.rda"))

