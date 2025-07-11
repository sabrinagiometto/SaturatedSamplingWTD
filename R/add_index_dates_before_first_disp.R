set.seed(123)

start <- 1095.75
end <- 1461

# generate data from a uniform distribution over 5 y
dt <- data.table(disp_time = round(runif(30, 0, 365.25*5),0))

# keep only unique dispensing times
dt <- dt[!duplicated(disp_time)]

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

# bind time windows:
# AFTER first dispensing, i.e., between first dispensing (within the sampling window) and end, AND BEFORE first dispensing, i.e., between first dispensing (within the sampling window) and last observation (within the moving window) 
dt_e <- rbind(dt_e_pre_first_disp, dt_e_post)
