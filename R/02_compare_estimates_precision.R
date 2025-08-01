library(wtdr)
source("R/satwtdttt.R")

# Fit different models to compare estimates precision ----

## using df_sat_ex.csv (sample of JH synthetic data) ----

set.seed(345)

df <- read.csv(file.path("extdata", "df_sat_ex.csv"))
df <- as.data.table(df)

df <- df[, X:=NULL]
df <- df[, rxdate:=as.Date(rxdate)]

# 1) single index date - forward
forw_1 <- wtdttt(data = df,
                 rxdate ~ dlnorm(logitp, mu, lnsigma),
                 id = "pid",
                 start = as.Date('2014-01-01'), 
                 end = as.Date('2014-12-31'), 
                 reverse = F
)

summary(forw_1)

predict(forw_1, distrx = "rxdate", quantile = 0.8)

# 2) single index date - reverse
rev_1 <- wtdttt(data = df,
                rxdate ~ dlnorm(logitp, mu, lnsigma),
                id = "pid",
                start = as.Date('2014-01-01'), 
                end = as.Date('2014-12-31'), 
                reverse = T
)

summary(rev_1)

predict(rev_1, distrx = "rxdate", quantile = 0.8)

# 3) multiple random index date (m =5)
rev_ran_5 <- ranwtdttt(data = df,
                       rxdate ~ dlnorm(logitp, mu, lnsigma),
                       id = "pid",
                       start = as.Date('2014-01-01'),
                       end = as.Date('2014-12-31'),
                       reverse = T,
                       nsamp = 5, 
                       robust = T
)

summary(rev_ran_5)

predict(rev_ran_5, distrx = "rxdate", quantile = 0.8)

# 4) multiple random index date (m =50)
rev_ran_50 <- ranwtdttt(data = df,
                        rxdate ~ dlnorm(logitp, mu, lnsigma),
                        id = "pid",
                        start = as.Date('2014-01-01'),
                        end = as.Date('2014-12-31'),
                        reverse = T,
                        nsamp = 50, 
                        robust = T
)

summary(rev_ran_50)

predict(rev_ran_50, distrx = "rxdate", quantile = 0.8)

# 5) saturated sampling
rev_sat_e <- satwtdttt(data = df,
                       rxdate ~ dlnorm(logitp, mu, lnsigma),
                       id = "pid",
                       start = as.Date('2014-01-01'), 
                       end = as.Date('2014-12-31'),
                       robust = T,
                       reverse = F
)


summary(rev_sat_e)

predict(rev_sat_e, distrx = "rxdate", quantile = 0.8)

## using simulated data (TO BE FIXED)----

load(file.path("extdata", "results_a.rda"))

## 1) single index date ----
fit1 <- wtdttt(data = results_a,
               disp_time ~ dlnorm(logitp, mu, lnsigma),
               id = "id",
               start = 0, 
               end = 365, 
               reverse = T
)

summary(fit1)

wtdr::predict(fit1, quantile = 0.8) 

## 2) multiple random index date (m = 5) ----
fit2 <- ranwtdttt(data = results_a,
                   disp_time ~ dlnorm(logitp, mu, lnsigma),
                   id = "id",
                   start = 0,
                   end = 365,
                   reverse = T,
                   nsamp = 5, 
                   robust = T
)

summary(fit2)

predict(fit2, quantile = 0.8) 

## 3) satured sampling ----
fit3 <- satwtdttt(data = results_a,
                  disp_time ~ dlnorm(logitp, mu, lnsigma),
                  id = "id",
                  start = 365, 
                  end = 730,
                  reverse = T, # con reverse = F funziona
                  robust = T
)



summary(fit3)

predict(fit3, quantile = 0.8) 
