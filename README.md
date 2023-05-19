# **Code for article A Quantile Shift Approach To Main Effects And Interactions In A 2-By-2 Design**
Rand R. Wilcox & Guillaume A. Rousselet  
(submitted)  
Preprint = <URL>  

Code to reproduce the figures, analyses and simulations in the article.
The data, knitted notebooks and figures in pdf format are available on [figshare](https://doi.org/10.6084/m9.figshare.22927520.v1).  

The main functions are illustrated in the notebook `examples.Rmd`. 

Method 1 is implemented in the function `decinter()` available in `./code/decinter.R`.

Method 2 is implemented in the function `iband()` available in `./code/iband.R`.

Both methods are implemented in the wrapper function `dec_apd_all_wrap` (available in `./code/examples.R`) that calls the Rcpp code `dec_apd_all.cpp` (also available in the `./code/` folder).  

Below is a breakdown of the different notebooks.  

## Method 1: compare quantiles of marginal distributions

**`sim_fp.Rmd`**  
Estimate FWER = probability of at least one false negative among 9 deciles. Use percentile bootstrap and bootstrap-t, with the same samples for all deciles. Use a percentile bootstrap approach as in Wilcox et al. (2014), but with the same samples for all deciles. 
Compare HD (Harrell & Davis, 1982) to quantile(type=7), the default in R (Hyndman & Fan, 1996; Wicklin, 2017).
Compare two methods to correct for multiple comparisons: Hochberg (1988), as used in Wilcox et al. (2014) and FDR (Benjamini & Hochberg, 1995). Compare to standard ANOVA on means and 20% trimmed means (no bootstrap, see chapter 7 in Wilcox, 2017).

*Conclusion:* All methods tested are a bit too conservative. FDR gives FWER closer to nominal level than Hochberg across all conditions. FDR dominates QT7 in all situations. Recommendation: use HD + FDR as default.

*Generated article figures:*  
- fig_art_fp_norm_lnorm.pdf  = Figure 2  
- fig_art_fp_pois_bbr9.pdf = Figure 3  
- fig_art_fp_ind_norm_bbr9.pdf  = Figure 4  

**`sim_fp_b1b9.Rmd`**  
Check whether using separate bootstrap samples for every decile improves coverage.

*Conclusion:* No evidence for better performance using separate bootstrap samples for every decile.
Keep using the same bootstrap samples for all deciles for computational efficiency.
Overall the approach is conservative, with FWER below the nominal level at all sample sizes, for all distributions, and for main effects and interaction.

**`sim_tp.Rmd`**  
Estimate family-wise power = probability of at least one true positive among 9 deciles.
Use percentile bootstrap with the same samples for all deciles. 
Compare Harrell-Davis quantile estimator to quantile(type=7).
Compare to standard ANOVA on means and 20% trimmed means (no bootstrap).

*Conclusion:* Overall, no method dominates: situations are found in which ANOVA on means outperforms quantile methods; in other situations power for ANOVA on means collapses while quantile methods retain high power. FDR performs better or as well as Hochberg in all situations. HD is always more powerful than QT7, especially when dealing with tied values.  
Recommend to use HD + FDR by default.  

*Generated article figures:*  
- fig_art_tp_norm_lnorm.pdf = Figure 5  
- fig_art_tp_pois_bb.pdf = Figure 6  

## Method 2: compare quantiles of distributions of all pairwise differences

**`sim_fp_apd.Rmd`**  
Strategy: for each level of A, compute all pairwise differences between B1 and B2. Then compare the quantiles (here deciles) of these two distributions. Estimate FWER = probability of at least one false negative among 9 deciles.
Use percentile bootstrap only and boot1 method = use the same bootstrap samples for all deciles. 

*Conclusion:* FWER is lower than the nominal level in all situations. Type I error rate is close to the expected level at individual quantiles, particularly for central quantiles and n>20. HD and QT7 perform similarly for continuous distributions. When sampling from distributions with tied values, HD is much closer to the nominal level than QT7. FDR outperforms Hochberg in all situations.
Recommend to use HD + FDR by default.  

*Generated article figures:*  
- fig_art_fp_apd.pdf = Figure 7  
- fig_art_fp_apd_ind_norm_clnorm.pdf = Figure 8  
- fig_art_fp_apd_ind_pois_bbr1.pdf = Figure 9    

**`sim_fp_apd_b1b9.Rmd`**  
Strategy: for each level of A, compute all pairwise differences between B1 and B2. Then compare the quantiles (here deciles) of these two distributions. Estimate FWER = probability of at least one false negative among 9 deciles. Use HD only; compare bootstrap with the same samples for all deciles vs. different samples. Do that for one sample size only. With n=40, the distribution of all pairwise differences = 1600 observations!

*Conclusion:* No evidence for better performance using separate bootstrap samples for every decile.
Keep using the same bootstrap samples for all deciles for computational efficiency.
Overall the approach is conservative for all distributions. Looking at individual deciles, we're closer to the nominal level nearer the centre of the distribution, a bit conservative at the extremes, but overall close to 0.05.

**`sim_tp_apd.Rmd`**  
Estimate familywise power = probability of at least one true positive among 9 deciles.
Use percentile bootstrap with the same samples for all deciles. 
Compare HD to QT7.  
Compare to standard ANOVA on means and 20% trimmed means (no bootstrap).

*Conclusion:* No method dominates. In some situations ANOVA on means outperforms the others, yet in other situations, power completely collapses for the ANOVA on means but not for the quantile method. HD and QT7 perform similarly for continuous distributions, but in the presence of tied values, HD dominates QT7.  

*Generated article figure:*  
- fig_art_tp_apd.pdf = Figure 10  

## Illustrations

**`hd.Rmd`**  
Illustrate the beta weights used to compute the Harrell-Davis quantile estimates.  

*Generated article figure:*  
- fig_ex_hd_beta_weights.pdf  = Figure 1  

**`apd_ex.Rmd`**  
Demonstrate that when dealing with distributions of all pairwise differences, interchanging the rows and columns can yield different interaction results.  

**`kurtosis_estimation.Rmd`**  
Estimate the kurtosis of samples from a lognormal distribution and a contaminated lognormal distribution.  

**`examples.Rmd`**  
Illustrate distributions used in the power simulations. Generate samples at maximum sample size. Plot marginals, shift functions for main effects and interaction, distributions of all pairwise differences, their quantiles, and their quantile differences (interaction).
Illustrate health example presented in the article. 

*Generated article figures:*  
- fig_ex_A1B1C_dec.pdf = Figure 11  
- fig_ex_A1B1C_apd.pdf = Figure 12  





