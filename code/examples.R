# Code used in examples.Rmd

# Compute p values
bootpval <- function(bootdist){
  # bootdist = nboot x quantile matrix
  temp <- colSums(bootdist < 0) / nboot + colSums(bootdist == 0) / (2*nboot)
  boot.pv <- 2*(apply(rbind(temp, 1-temp), 2, min))
  boot.pv
}

# Calculate everything
# x is a list
# w = pre-calculated beta weights for method 1
# w_apd = pre-calculated beta weights for method 2
# nboot = number of bootstrap samples
# probs = quantiles of the bootstrap distribution to return
dec_apd_all_wrap <- function(x, nboot, probs){
  
  # pre-compute HD beta weights once here
  na <- length(x[[1]]) # sample sizes
  nb <- length(x[[2]])
  nc <- length(x[[3]])
  nd <- length(x[[4]])
  hd.w <- hd.w.calc(c(na,nb,nc,nd), qseq)
  wa <- hd.w[[1]]
  wb <- hd.w[[2]]
  wc <- hd.w[[3]]
  wd <- hd.w[[4]]
  hd.w <- hd.w.calc(c(na*nb,nc*nd), qseq)
  wab_apd <- hd.w[[1]]
  wcd_apd <- hd.w[[2]]
  
  # compute bootstrap distributions
  res <- dec_apd_all_cpp(x, wa, wb, wc, wd, wab_apd, wcd_apd, nboot)
  # x = list with 4 groups
  # w... = HD beta weights
  # nboot = number of bootstrap samples
  
  # ------------------------------
  # get results
  # Decile estimates for the marginals
  ahd <- res[["ahd"]]
  bhd <- res[["bhd"]]
  chd <- res[["chd"]]
  dhd <- res[["dhd"]]
  
  # differences and interaction: bootstrap distributions
  boot_ahd <- res[["boot_ahd"]] # marginals
  boot_bhd <- res[["boot_bhd"]]
  boot_chd <- res[["boot_chd"]]
  boot_dhd <- res[["boot_dhd"]]
  boot_ab_diff <- boot_ahd - boot_bhd # shift functions at each level of A1
  boot_cd_diff <- boot_chd - boot_dhd
  boot_mainA <- boot_ahd + boot_bhd - boot_chd - boot_dhd # main effect of A
  boot_mainB <- boot_ahd - boot_bhd + boot_chd - boot_dhd # main effect of B
  boot_inter <- boot_ahd - boot_bhd - boot_chd + boot_dhd # interaction
  
  # all pairwise differences
  apd_ab <- res[["apd_ab"]] # distributions of all pairwise differences (apd)
  apd_cd <- res[["apd_cd"]]
  apd_ab_hd <- res[["apd_ab_hd"]] # deciles of apd
  apd_cd_hd <- res[["apd_cd_hd"]]
  apd_ab_boot <- res[["apd_ab_boot"]] # bootstrap distributions of deciles of apd
  apd_cd_boot <- res[["apd_cd_boot"]]
  apd_inter_boot <- apd_ab_boot - apd_cd_boot # bootstrap distributions of apd interaction
  # ------------------------------
  
  # Compute confidence intervals
  ab_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  cd_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_ab_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_cd_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  mainA_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  mainB_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  inter_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  ab_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  cd_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_ab_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_cd_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  mainA_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  mainB_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  inter_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_inter_ci <- matrix(data = NA, nrow = 9, ncol = 2)
  apd_inter_ci50 <- matrix(data = NA, nrow = 9, ncol = 2)
  
  for(q in 1:9){
    # individual shift functions
    ab_ci[q,] <- quantile(sort(boot_ab_diff[,q]), probs = probs)
    cd_ci[q,] <- quantile(sort(boot_cd_diff[,q]), probs = probs)
    ab_ci50[q,] <- quantile(sort(boot_ab_diff[,q]), probs = c(.25, .75))
    cd_ci50[q,] <- quantile(sort(boot_cd_diff[,q]), probs = c(.25, .75))
    apd_ab_ci[q,] <- quantile(sort(apd_ab_boot[,q]), probs = probs)
    apd_cd_ci[q,] <- quantile(sort(apd_cd_boot[,q]), probs = probs)
    apd_ab_ci50[q,] <- quantile(sort(apd_ab_boot[,q]), probs = c(.25, .75))
    apd_cd_ci50[q,] <- quantile(sort(apd_cd_boot[,q]), probs = c(.25, .75))
    
    # main effects
    mainA_ci[q,] <- quantile(sort(boot_mainA[,q]), probs = probs)
    mainB_ci[q,] <- quantile(sort(boot_mainB[,q]), probs = probs)
    mainA_ci50[q,] <- quantile(sort(boot_mainA[,q]), probs = c(.25, .75))
    mainB_ci50[q,] <- quantile(sort(boot_mainB[,q]), probs = c(.25, .75))
    
    # interactions
    inter_ci[q,] <- quantile(sort(boot_inter[,q]), probs = probs)
    apd_inter_ci[q,] <- quantile(sort(apd_inter_boot[,q]), probs = probs)
    inter_ci50[q,] <- quantile(sort(boot_inter[,q]), probs = c(.25, .75))
    apd_inter_ci50[q,] <- quantile(sort(apd_inter_boot[,q]), probs = c(.25, .75))
  }
  
  # Compute p values
  ab.pv <- bootpval(boot_ab_diff) # shift function for A1: A1B1 - A1B2
  cd.pv <- bootpval(boot_cd_diff) # shift function for A2: A2B1 - A2B2
  mainA.pv <- bootpval(boot_mainA)
  mainB.pv <- bootpval(boot_mainB)
  inter.pv <- bootpval(boot_inter)
  apd.inter.pv <- bootpval(apd_inter_boot)
  
  # export everything
  list(ahd = ahd, 
       bhd = bhd, 
       chd = chd, 
       dhd = dhd,
       ab = ahd - bhd,
       cd = chd - dhd,
       ab_ci = ab_ci, 
       cd_ci = cd_ci, 
       ab_ci50 = ab_ci50, 
       cd_ci50 = cd_ci50, 
       # ab_pv = ab_pv,
       # cd_pv = cd_pv,
       mainA = ahd + bhd - chd - dhd,
       mainB = ahd - bhd + chd - dhd,
       # ---------------------------------------------------
       # options to illustrate main effects and interactions
       # average deciles to plot along x axis
       A1_mean = (ahd + bhd) / 2, 
       A2_mean = (chd + dhd) / 2, 
       B1_mean = (ahd + chd) / 2, 
       B2_mean = (bhd + dhd) / 2, 
       # deciles of pooled values
       A1_pool = deciles(c(x[[1]], x[[2]])), 
       A2_pool = deciles(c(x[[3]], x[[4]])),
       B1_pool = deciles(c(x[[1]], x[[3]])),
       B2_pool = deciles(c(x[[2]], x[[4]])),
       # ---------------------------------------------------
       inter = ahd - bhd - chd + dhd,
       mainA_ci = mainA_ci, 
       mainB_ci = mainB_ci, 
       inter_ci = inter_ci, 
       mainA_ci50 = mainA_ci50, 
       mainB_ci50 = mainB_ci50, 
       inter_ci50 = inter_ci50, 
       apd_ab = apd_ab,
       apd_cd = apd_cd,
       apd_inter = apd_ab_hd - apd_cd_hd,
       apd_ab_hd = apd_ab_hd,
       apd_cd_hd = apd_cd_hd,
       apd_ab_ci = apd_ab_ci, 
       apd_cd_ci = apd_cd_ci, 
       apd_inter_ci = apd_inter_ci,
       apd_ab_ci50 = apd_ab_ci50, 
       apd_cd_ci50 = apd_cd_ci50, 
       apd_inter_ci50 = apd_inter_ci50,
       ab.pv = ab.pv,
       cd.pv = cd.pv,
       mainA.pv = mainA.pv,
       mainB.pv = mainB.pv,
       inter.pv = inter.pv,
       apd.inter.pv = apd.inter.pv)
}
# ------------------------------

# Plot functions
plot_marg <- function(x, res, col = "black"){
  
  # x <- matl(x)
  na <- length(x[[1]])
  nb <- length(x[[2]])
  nc <- length(x[[3]])
  nd <- length(x[[4]])
  
  df <- tibble(y = c(unlist(x)),
               # Conditions = factor(rep(c("A1B1", "A1B2", "A2B1", "A2B2"), each = n)),
               Conditions = c(rep(1, na), rep(2, nb), rep(3, nc), rep(4, nd)),
  )
  
  df.hd <- tibble(Conditions = rep(1:4, each = 9),
                  HD = c(res$ahd, res$bhd, res$chd, res$dhd))
  
  df.md <- tibble(Conditions = 1:4,
                  HD = c(res$ahd[5], res$bhd[5], res$chd[5], res$dhd[5]))
  
  p <- ggplot(data = df, aes(x = Conditions, y = y, fill = Conditions)) + theme_gar +
    geom_jitter(position=position_jitter(0.2), show.legend = NA, 
                colour = col, fill = col,
                shape = 21, size = 3, alpha = 0.15) +
    # add deciles
    geom_errorbar(data = df.hd,
                  aes(x = Conditions,
                      ymin = HD,
                      y = HD,
                      ymax = HD),
                  colour = col,
                  linewidth = 0.75,
                  width = 0.3,
                  lineend = "round") +
    # make median thicker than other quantiles
    geom_errorbar(data = df.md,
                  aes(x = Conditions,
                      ymin = HD,
                      y = HD,
                      ymax = HD),
                  colour = col,
                  linewidth = 1.5,
                  width = 0.4,
                  lineend = "round") +
    scale_x_continuous(breaks=c(1,2,3,4), labels=c("A1B1", "A1B2", "A2B1", "A2B2")) +
    theme(legend.position = "none") +
    labs(y = "Observations") +
    # horizontal reference line at zero
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75)
  p
}

plot_ind_sf_ab <- function(res, plot_param){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "A1B1"
  lab.y <- "A1: B1 - B2"
  
  # median reference line
  xintercept <- res$ahd[5]
  
  df <- tibble(x = res$ahd, 
               y = res$ab,
               ci_lower = res$ab_ci[,1],
               ci_upper = res$ab_ci[,2],
               ci_lower50 = res$ab_ci50[,1],
               ci_upper50 = res$ab_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   linewidth = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   linewidth = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              linewidth = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  p
}

plot_ind_sf_cd <- function(res, plot_param){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "A2B1"
  lab.y <- "A2: B1 - B2"
  # median reference line
  xintercept <- res$chd[5]
  
  df <- tibble(x = res$chd, 
               y = res$cd,
               ci_lower = res$cd_ci[,1],
               ci_upper = res$cd_ci[,2],
               ci_lower50 = res$cd_ci50[,1],
               ci_upper50 = res$cd_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   linewidth = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   linewidth = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              linewidth = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  p
}

plot_main_A <- function(res, plot_param, plotx = "pool"){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "A1"
  lab.y <- "A1 - A2"
  if(plotx == "mean"){  
    # x axis: plot mean of quantiles 
    xintercept <- res$A1_mean[5] # median reference line
    xplot <- res$A1_mean
    # or quantiles of observations pooled across 2 levels
  } else if(plotx == "pool"){
    xintercept <- res$A1_pool[5] # median reference line
    xplot <- res$A1_pool
  }
  
  df <- tibble(x = xplot, 
               y = res$mainA,
               ci_lower = res$mainA_ci[,1],
               ci_upper = res$mainA_ci[,2],
               ci_lower50 = res$mainA_ci50[,1],
               ci_upper50 = res$mainA_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   linewidth = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   linewidth = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              linewidth = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  
  p
}

plot_main_B <- function(res, plot_param, plotx = "pool"){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "B1"
  lab.y <- "B1 - B2"
  if(plotx == "mean"){  
    # x axis: plot mean of quantiles 
    xintercept <- res$B1_mean[5] # median reference line
    xplot <- res$B1_mean
    # or quantiles of observations pooled across 2 levels
  } else if(plotx == "pool"){
    xintercept <- res$B1_pool[5] # median reference line
    xplot <- res$B1_pool
  }
  
  df <- tibble(x = xplot, 
               y = res$mainB,
               ci_lower = res$mainB_ci[,1],
               ci_upper = res$mainB_ci[,2],
               ci_lower50 = res$mainB_ci50[,1],
               ci_upper50 = res$mainB_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   linewidth = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   linewidth = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              linewidth = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  
  p
}

plot_inter <- function(res, plot_param, plotx = "pool"){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "A1"
  lab.y <- "Interaction"

  if(plotx == "mean"){  
    # x axis: plot mean of quantiles 
    xintercept <- res$A1_mean[5] # median reference line
    xplot <- res$A1_mean
    # or quantiles of observations pooled across 2 levels
  } else if(plotx == "pool"){
    xintercept <- res$A1_pool[5] # median reference line
    xplot <- res$A1_pool
  }

  df <- tibble(x = xplot, 
               y = res$inter,
               ci_lower = res$inter_ci[,1],
               ci_upper = res$inter_ci[,2],
               ci_lower50 = res$inter_ci50[,1],
               ci_upper50 = res$inter_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   linewidth = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   linewidth = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              linewidth = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  
  p
}

plot_apd <- function(res, plot_param){
  
  df <- tibble(apd = c(res$apd_ab,res$apd_cd),
               Conditions = factor(c(rep("A1: B1-B2", length(res$apd_ab)), rep("A2: B1-B2", length(res$apd_cd)))),
  )
  
  df.hd <- tibble(Conditions = factor(rep(c("A1: B1-B2", "A2: B1-B2"), each = 9)),
                  HD = c(res$apd_ab_hd, res$apd_cd_hd))
  
  df.md <- tibble(Conditions = factor(c("A1: B1-B2", "A2: B1-B2")),
                  HD = c(res$apd_ab_hd[5], res$apd_cd_hd[5]))
  
  p <- ggplot(data = df, aes(x = apd)) + theme_gar +
    geom_density(fill = "grey95", trim = TRUE) +
    facet_wrap(vars(Conditions), nrow = 2, strip.position = "right") +
    # add deciles
    geom_vline(data = df.hd,
               mapping=aes(xintercept=HD)) +
    # make median thicker than other quantiles
    geom_vline(data = df.md,
               mapping=aes(xintercept=HD),
               size = 2) +
    labs(x = "All pairwise differences",
         y = "Density") +
    geom_vline(xintercept = 0, linetype = 2) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(size=16))
  
  p
}

plot_apd_inter <- function(res, plot_param){
  
  # vertical bars for uncertainty intervals
  diffint_col <- plot_param$diffint_col
  diffint_size <- plot_param$diffint_size
  diffint_col50 <-plot_param$diffint_col50
  diffint_size50 <- plot_param$diffint_size50
  # line joining the quantiles
  q_line_col <- plot_param$q_line_col
  q_line_alpha <- plot_param$q_line_alpha
  q_line_size <- plot_param$q_line_size
  # symbols marking the quantiles
  symb_col <- plot_param$symb_col
  symb_size <- plot_param$symb_size
  symb_shape <- plot_param$symb_shape
  symb_fill <- plot_param$symb_fill
  
  lab.x <- "A1: B1-B2"
  lab.y <- "Interaction: A1-A2"
  # median reference line
  xintercept <- res$apd_ab_hd[5]
  
  df <- tibble(x = res$apd_ab_hd, 
               y = res$apd_inter,
               ci_lower = res$apd_inter_ci[,1],
               ci_upper = res$apd_inter_ci[,2],
               ci_lower50 = res$apd_inter_ci50[,1],
               ci_upper50 = res$apd_inter_ci50[,2])
  
  p <- ggplot(df, aes(x = x, y = y)) + theme_gar +
    # x=0 reference line
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
    # y=median reference line
    geom_vline(xintercept = xintercept, linetype = 2, alpha = 0.5) +
    xlab(lab.x) +
    ylab(lab.y) +
    # scale_y_continuous(limits = ylim) +
    # vertical bars for uncertainty intervals
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), colour = diffint_col,
                   size = diffint_size) +
    geom_linerange(aes(ymin = ci_lower50, ymax = ci_upper50), colour = diffint_col50,
                   size = diffint_size50) +
    # line joining the quantiles
    geom_line(colour = q_line_col, alpha = q_line_alpha, linetype = "solid",
              size = q_line_size) +
    # symbols marking the quantiles
    geom_point(colour = symb_col, size = symb_size, shape = symb_shape, fill = symb_fill)
  
  p
}
# ------------------------------

