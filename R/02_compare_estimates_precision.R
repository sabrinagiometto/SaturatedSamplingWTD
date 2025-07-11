library(wtdr)

# Fit different models to compare estimates precision ----

## 1) single index date ----
fit1 <- wtdttt(data = results_a,
               disp_time ~ dlnorm(logitp, mu, lnsigma),
               id = "id",
               start = 0, 
               end = 365, 
               reverse = T
)

summary(fit1)

wtdr::predict(fit1, quantile = 0.9) # SE = 5.15, Estimate = 107

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

predict(fit2, quantile = 0.9) # SE = 2.02, Estimate = 116

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

predict(fit3, quantile = 0.9) # SE = 1.39, Estimate = 83
