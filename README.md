# SaturatedSamplingWTD
Implementation of saturated sampling in the use of parametric waiting time distribution for assigning dispensation lengths with secondary use of datasources.
The approach aims to leverage as much information as possible from available data by increasing the effective sample size and, consequently, improving statistical precision. This is achieved by treating each day within a prespecified time window as an index date.

There are two scripts available in the repository:
1) satwtdttt.R: it implements the saturated sampling based on the wtdttt() function from the wtdr package, which is available for download at https://github.com/mbg-unsw/wtdr
2) compare_estimates_precision_with_m_simulations.R: it includes the data generation process, the analysis according to four analytic methods (i) reverse WTD with a single fixed index date; ii) reverse WTD with a small number of random index dates (m=5); iii) reverse WTD with a larger number of random index dates (m=50); iv) reverse WTD with saturated sampling, including all possible index dates), and the summary results comparing the precision of estimates between analytic methods.
