
# Plot FWER for HD/QT, FDR/hochberg + ANOVA
plot_fp_fwer <- function(fwer, simres, plot.title){
  
  df <- tibble(prob = c(apply(fwer$AB.qt.raw, 2, mean), 
                        apply(fwer$AB.qt.fdr, 2, mean),
                        apply(fwer$AB.qt.hoch, 2, mean),
                        apply(fwer$AB.hd.raw, 2, mean), 
                        apply(fwer$AB.hd.fdr, 2, mean),
                        apply(fwer$AB.hd.hoch, 2, mean),
                        apply(fwer$A.qt.raw, 2, mean), 
                        apply(fwer$A.qt.fdr, 2, mean),
                        apply(fwer$A.qt.hoch, 2, mean),
                        apply(fwer$A.hd.raw, 2, mean), 
                        apply(fwer$A.hd.fdr, 2, mean),
                        apply(fwer$A.hd.hoch, 2, mean),
                        apply(fwer$B.qt.raw, 2, mean), 
                        apply(fwer$B.qt.fdr, 2, mean),
                        apply(fwer$B.qt.hoch, 2, mean),
                        apply(fwer$B.hd.raw, 2, mean), 
                        apply(fwer$B.hd.fdr, 2, mean),
                        apply(fwer$B.hd.hoch, 2, mean)),
               n = rep(nseq, 6*3),
               method = factor(rep(rep(c("No corr.", "FDR", "Hochberg"), each = nn), 2*3)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*3), 3)),
               comp = factor(rep(c("Interaction", "Main A", "Main B"), each = nn*6)) 
  )
  
  df2 <- tibble(prob = c(apply(simres$ANOVA.m[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,3]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,3]<alpha.val, 2, mean)),
                n = rep(nseq, 3 * 2),
                method = factor(rep(c("ANOVA.M", "ANOVA.TM"), each = nn*3)),
                comp = factor(rep(rep(c("Main A", "Main B", "Interaction"), each = nn), 2))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    geom_point(data = df2, aes(x = n, y = prob, colour = method)) +
    geom_line(data = df2, aes(x = n, y = prob, colour = method)) +
    scale_x_continuous(breaks = nseq) +
    # scale_y_continuous(limits = c(0,0.15)) +
    scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00", "#56B4E9", "#999999")) +
    labs(x = "Sample size", y = "Familywise Error Rate", title = plot.title) +
    facet_grid(rows = vars(comp),
               cols = vars(estimator)) +
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_blank()) +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
  
  p
}

# Plot only FDR results: boot 1 vs. boot 9
plot_fp_fwer_fdr_b1b9 <- function(fwer, plot.title){
  
  df <- tibble(prob = c(apply(fwer$AB.b9.raw, 2, mean), 
                        apply(fwer$AB.b9.fdr, 2, mean),
                        apply(fwer$AB.b1.raw, 2, mean), 
                        apply(fwer$AB.b1.fdr, 2, mean),
                        apply(fwer$A.b9.raw, 2, mean), 
                        apply(fwer$A.b9.fdr, 2, mean),
                        apply(fwer$A.b1.raw, 2, mean), 
                        apply(fwer$A.b1.fdr, 2, mean),
                        apply(fwer$B.b9.raw, 2, mean), 
                        apply(fwer$B.b9.fdr, 2, mean),
                        apply(fwer$B.b1.raw, 2, mean), 
                        apply(fwer$B.b1.fdr, 2, mean)),
               n = rep(nseq, 4*3),
               Test = factor(rep(rep(c("No corr.", "FDR"), each = nn), 2*3)),
               Boot = factor(rep(rep(c("boot9", "boot1"), each = nn*2), 3)),
               Comp = factor(rep(c("Interaction", "Main A", "Main B"), each = nn*4)) 
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = Test)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line(aes(linetype=Boot)) +
    scale_x_continuous(breaks = nseq) +
    # scale_y_continuous(limits = c(0,0.15)) +
    scale_color_manual(values=c("#E69F00", "#999999")) +
    # "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
    labs(x = "Sample size", y = "Familywise Error Rate", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) +
    facet_grid(rows = vars(Comp))
  
  p
}

# Plot FWER for HD/QT, FDR only + ANOVA
plot_fp_fwer_fdr <- function(fwer, simres, plot.title){
  
  df <- tibble(prob = c(apply(fwer$AB.qt.raw, 2, mean), 
                        apply(fwer$AB.qt.fdr, 2, mean),
                        # apply(fwer$AB.qt.hoch, 2, mean),
                        apply(fwer$AB.hd.raw, 2, mean), 
                        apply(fwer$AB.hd.fdr, 2, mean),
                        # apply(fwer$AB.hd.hoch, 2, mean),
                        apply(fwer$A.qt.raw, 2, mean), 
                        apply(fwer$A.qt.fdr, 2, mean),
                        # apply(fwer$A.qt.hoch, 2, mean),
                        apply(fwer$A.hd.raw, 2, mean), 
                        apply(fwer$A.hd.fdr, 2, mean),
                        # apply(fwer$A.hd.hoch, 2, mean),
                        apply(fwer$B.qt.raw, 2, mean), 
                        apply(fwer$B.qt.fdr, 2, mean),
                        # apply(fwer$B.qt.hoch, 2, mean),
                        apply(fwer$B.hd.raw, 2, mean), 
                        apply(fwer$B.hd.fdr, 2, mean)),
               # apply(fwer$B.hd.hoch, 2, mean)
               n = rep(nseq, 12),
               method = factor(rep(rep(c("No corr.", "FDR"), each = nn), 6)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*2), 3)),
               comp = factor(rep(c("Interaction", "Main A", "Main B"), each = nn*4)) 
  )
  
  df2 <- tibble(prob = c(apply(simres$ANOVA.m[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,3]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,3]<alpha.val, 2, mean)),
                n = rep(nseq, 3 * 2),
                method = factor(rep(c("ANOVA.M", "ANOVA.TM"), each = nn*3)),
                comp = factor(rep(rep(c("Main A", "Main B", "Interaction"), each = nn), 2))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line(aes(linetype=estimator)) +
    geom_point(data = df2, aes(x = n, y = prob, colour = method)) +
    geom_line(data = df2, aes(x = n, y = prob, colour = method)) +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0,0.2)) +
    scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00", "#999999")) +
    labs(x = "Sample size", y = "Familywise Error Rate", title = plot.title) +
    facet_grid(rows = vars(comp)) +
    theme(legend.position = "right",
          panel.grid.minor.x = element_blank()) +
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  
  p
}

# --------------------------------------------
# Plot results at individual Quantiles for QT7 and HD
# Plot QT7 and HD separately, to contrast sample sizes.
plot_fp_ind <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*6)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq*3)),
               comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
               decile = rep(rep(qseq, each = nn), 6) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       labels = c(".1","",".3","",".5","",".7","",".9")) +
    scale_y_continuous(limits = c(0,0.15)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Type I Error Rate", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) + 
    facet_grid(rows = vars(comp),
               cols = vars(estimator))
  
  p
}

# Plot results at individual Quantiles for QT7 and HD
# Plot separately for each sample size to compare QT7 to HD
plot_fp_ind_n <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*6)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq*3)),
               comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
               decile = rep(rep(qseq, each = nn), 6) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, linetype = estimator)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       labels = c(".1","",".3","",".5","",".7","",".9")) +
    scale_y_continuous(limits = c(0,0.15)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Type I Error Rate", title = plot.title) +
    facet_grid(rows = vars(n),
               cols = vars(comp)) +
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_blank()) +
    labs(linetype='quantile method')
  
  p
}

# Plot family-wise power: at least one positive result among the deciles
# Used in sim_tp.Rmd
plot_tp_fwpr <- function(fwpr, simres, alpha.val, plot.title){
  
  df <- tibble(prob = c(apply(fwpr$AB.qt.raw, 2, mean), 
                        apply(fwpr$AB.qt.fdr, 2, mean),
                        apply(fwpr$AB.qt.hoch, 2, mean),
                        apply(fwpr$AB.hd.raw, 2, mean), 
                        apply(fwpr$AB.hd.fdr, 2, mean),
                        apply(fwpr$AB.hd.hoch, 2, mean),
                        apply(fwpr$A.qt.raw, 2, mean), 
                        apply(fwpr$A.qt.fdr, 2, mean),
                        apply(fwpr$A.qt.hoch, 2, mean),
                        apply(fwpr$A.hd.raw, 2, mean), 
                        apply(fwpr$A.hd.fdr, 2, mean),
                        apply(fwpr$A.hd.hoch, 2, mean),
                        apply(fwpr$B.qt.raw, 2, mean), 
                        apply(fwpr$B.qt.fdr, 2, mean),
                        apply(fwpr$B.qt.hoch, 2, mean),
                        apply(fwpr$B.hd.raw, 2, mean), 
                        apply(fwpr$B.hd.fdr, 2, mean),
                        apply(fwpr$B.hd.hoch, 2, mean)),
               n = rep(nseq, 6*3),
               method = factor(rep(rep(c("No corr.", "FDR", "Hochberg"), each = nn), 2*3)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*3), 3)),
               comp = factor(rep(c("Interaction", "Main A", "Main B"), each = nn*6)) 
  )
  
  df2 <- tibble(prob = c(apply(simres$ANOVA.m[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,3]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,3]<alpha.val, 2, mean)),
                n = rep(nseq, 3 * 2),
                method = factor(rep(c("ANOVA.M", "ANOVA.TM"), each = nn*3)),
                comp = factor(rep(rep(c("Main A", "Main B", "Interaction"), each = nn), 2))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    geom_point(data = df2, aes(x = n, y = prob, colour = method)) +
    geom_line(data = df2, aes(x = n, y = prob, colour = method)) +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00", "#56B4E9", "#999999")) +
    labs(x = "Sample size", y = "Familywise Power", title = plot.title) +
    facet_grid(rows = vars(comp),
               cols = vars(estimator)) +
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_blank()) +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
  
  p
}

# Plot family-wise power for HD/QT, FDR only + ANOVA
plot_tp_fwpr_fdr <- function(fwpr, simres, alpha.val, plot.title){
  
  df <- tibble(prob = c(apply(fwpr$AB.qt.raw, 2, mean), 
                        apply(fwpr$AB.qt.fdr, 2, mean),
                        # apply(fwpr$AB.qt.hoch, 2, mean),
                        apply(fwpr$AB.hd.raw, 2, mean), 
                        apply(fwpr$AB.hd.fdr, 2, mean),
                        # apply(fwpr$AB.hd.hoch, 2, mean),
                        apply(fwpr$A.qt.raw, 2, mean), 
                        apply(fwpr$A.qt.fdr, 2, mean),
                        # apply(fwpr$A.qt.hoch, 2, mean),
                        apply(fwpr$A.hd.raw, 2, mean), 
                        apply(fwpr$A.hd.fdr, 2, mean),
                        # apply(fwpr$A.hd.hoch, 2, mean),
                        apply(fwpr$B.qt.raw, 2, mean), 
                        apply(fwpr$B.qt.fdr, 2, mean),
                        # apply(fwpr$B.qt.hoch, 2, mean),
                        apply(fwpr$B.hd.raw, 2, mean), 
                        apply(fwpr$B.hd.fdr, 2, mean)),
               # apply(fwpr$B.hd.hoch, 2, mean)
               n = rep(nseq, 12),
               method = factor(rep(rep(c("No corr.", "FDR"), each = nn), 6)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*2), 3)),
               comp = factor(rep(c("Interaction", "Main A", "Main B"), each = nn*4)) 
  )
  
  df2 <- tibble(prob = c(apply(simres$ANOVA.m[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.m[,,3]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,1]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,2]<alpha.val, 2, mean),
                         apply(simres$ANOVA.tm[,,3]<alpha.val, 2, mean)),
                n = rep(nseq, 3 * 2),
                method = factor(rep(c("ANOVA.M", "ANOVA.TM"), each = nn*3)),
                comp = factor(rep(rep(c("Main A", "Main B", "Interaction"), each = nn), 2))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line(aes(linetype=estimator)) +
    geom_point(data = df2, aes(x = n, y = prob, colour = method)) +
    geom_line(data = df2, aes(x = n, y = prob, colour = method)) +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00", "#999999")) +
    labs(x = "Sample size", y = "Familywise Power", title = plot.title) +
    facet_grid(rows = vars(comp)) +
    theme(legend.position = "right",
          panel.grid.minor.x = element_blank()) +
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  
  p
}

# --------------------------------------------
# Plot power results at individual Quantiles for QT7 and HD
# Plot QT7 and HD separately, to contrast sample sizes.
plot_tp_ind <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*6)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq*3)),
               comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
               decile = rep(rep(qseq, each = nn), 6) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       labels = c(".1","",".3","",".5","",".7","",".9")) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Power", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) +
    facet_grid(rows = vars(comp),
               cols = vars(estimator))
  
  p
}

# Plot power results at individual Quantiles for QT7 and HD
# Plot separately for each sample size to compare QT7 to HD
plot_tp_ind_n <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.qt<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$A.hd<alpha.val, c(2,3), mean)), 
                        as.vector(apply(simres$B.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*6)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq*3)),
               comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
               decile = rep(rep(qseq, each = nn), 6) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, linetype = estimator)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       labels = c(".1","",".3","",".5","",".7","",".9")) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Power", title = plot.title) +
    facet_grid(rows = vars(n),
               cols = vars(comp)) +
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_blank()) +
    labs(linetype='quantile method')
  
  p
}

# ----------------------------------------------
# APDINTER method: boot1 vs. boot9
plot_fp_apd_ind_b1b9 <- function(simres, plot.title, alpha.val){
  
  b1.ci <- matrix(data = NA, nrow = 2, ncol = 9)
  b9.ci <- matrix(data = NA, nrow = 2, ncol = 9)
  for(Q in 1:9){
    todo <- simres$AB.b1 < alpha.val
    b1.ci[,Q] <- binom.test(sum(todo[,Q]), length(todo[,Q]), p = alpha.val)$conf.int
    todo <- simres$AB.b9 < alpha.val
    b9.ci[,Q] <- binom.test(sum(todo[,Q]), length(todo[,Q]), p = alpha.val)$conf.int
  }
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.b1<alpha.val, 2, mean)), 
                        as.vector(apply(simres$AB.b9<alpha.val, 2, mean))),
               Boot = factor(rep(c("boot1", "boot9"), each = nq)),
               decile = rep(qseq, 2),
               ci.lo = c(b1.ci[1,], b9.ci[1,]),
               ci.hi = c(b1.ci[2,], b9.ci[2,])
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, colour = Boot)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point(size = 3, position=position_dodge(0.02)) +
    geom_line(linewidth = 0.75, position=position_dodge(0.02)) +
    geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.05,
                  position=position_dodge(0.02),
                  linewidth = 0.75) +
    scale_x_continuous(breaks = qseq,
                       # labels = c(".1","",".3","",".5","",".7","",".9")
                       labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
    scale_y_continuous(limits = c(0,0.15)) +
    # scale_color_viridis_d(end = 0.9, option = "magma") +
    scale_color_manual(values=c("#E69F00", "black")) + ##999999
    labs(x = "Quantiles", y = "Type I Error Rate", title = plot.title)
  
  p
}

# ----------------------------------------------
# APDINTER method: boot1 only 
plot_fp_apd_fwer <- function(fwer, plot.title){
  
  df <- tibble(prob = c(apply(fwer$AB.qt.raw, 2, mean), 
                        apply(fwer$AB.qt.fdr, 2, mean),
                        apply(fwer$AB.qt.hoch, 2, mean),
                        apply(fwer$AB.hd.raw, 2, mean), 
                        apply(fwer$AB.hd.fdr, 2, mean),
                        apply(fwer$AB.hd.hoch, 2, mean)),
               n = rep(nseq, 6),
               method = factor(rep(rep(c("No corr.", "FDR", "Hochberg"), each = nn), 2)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*3)))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0, 0.17)) +
    scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
    # "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
    labs(x = "Sample size", y = "Familywise Error Rate", title = plot.title) +
    facet_grid(cols = vars(estimator)) +
    theme(legend.position = "bottom")  
    # guides(color=guide_legend(nrow=2, byrow=TRUE))
  
  p
}

# Plot only FDR results
plot_fp_apd_fwer_fdr <- function(fwer, plot.title){
  
  df <- tibble(prob = c(apply(fwer$AB.qt.raw, 2, mean), 
                        apply(fwer$AB.qt.fdr, 2, mean),
                        apply(fwer$AB.hd.raw, 2, mean), 
                        apply(fwer$AB.hd.fdr, 2, mean)),
               n = rep(nseq, 4),
               method = factor(rep(rep(c("No corr.", "FDR"), each = nn), 2)),
               estimator = factor(rep(rep(c("QT7", "HD"), each = nn*2)))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line(aes(linetype=estimator)) +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0, 0.17)) +
    scale_color_manual(values=c("#E69F00", "#999999")) +
    labs(x = "Sample size", y = "Familywise Error Rate", title = plot.title)
  
  p
}

# Plot QT7 and HD separately, to contrast sample sizes.
plot_fp_apd_ind <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)),
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*2)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq)),
               decile = rep(rep(qseq, each = nn), 2) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       # labels = c(".1","",".3","",".5","",".7","",".9")
                       labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
    scale_y_continuous(limits = c(0,0.15)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Type I Error Rate", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) +
    facet_grid(cols = vars(estimator))
  
  p
}

# Plot results at individual Quantiles for QT7 and HD
# Plot separately for each sample size to compare QT7 to HD
plot_fp_apd_ind_n <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)),
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*2)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq)),
               decile = rep(rep(qseq, each = nn),2) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, linetype = estimator)) + theme_gar +
    # Bradley's (1978) satisfactory range
    geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
                fill = "grey80", colour = "grey80", show.legend=FALSE) +
    # Bradley's (1978) ideal range
    geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
                fill = "grey90", colour = "grey90", show.legend=FALSE) +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       # labels = c(".1","",".3","",".5","",".7","",".9")
                       labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
    scale_y_continuous(limits = c(0,0.15)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Type I Error Rate", title = plot.title) +
    facet_wrap(vars(n)) +
    theme(legend.position = "bottom") +
    labs(linetype='quantile method')
  
  p
}

# Plot power results
plot_tp_apd_fwpr <- function(fwpr, simres, alpha.val, plot.title){
  
  df <- tibble(prob = c(apply(fwpr$AB.qt.raw, 2, mean), 
                        apply(fwpr$AB.qt.fdr, 2, mean),
                        apply(fwpr$AB.hd.raw, 2, mean), 
                        apply(fwpr$AB.hd.fdr, 2, mean)),
               n = rep(nseq, 4),
               method = factor(rep(rep(c("No corr.", "FDR"), each = nn), 2)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*2))
  )
  
  df2 <- tibble(prob = c(apply(simres$ANOVA.m < alpha.val, 2, mean),
                         apply(simres$ANOVA.tm < alpha.val, 2, mean)),
                n = rep(nseq, 2),
                method = factor(rep(c("ANOVA.M", "ANOVA.TM"), each = nn))
  )
  
  p <- ggplot(df, aes(x = n, y = prob, colour = method)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line(aes(linetype=estimator)) +
    geom_point(data = df2, aes(x = n, y = prob, colour = method)) +
    geom_line(data = df2, aes(x = n, y = prob, colour = method)) +
    scale_x_continuous(breaks = nseq) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_manual(values=c("#009E73", "#0072B2", "#E69F00", "#999999")) +
    #, "#56B4E9",, "#D55E00")) +
    labs(x = "Sample size", y = "Familywise Power", title = plot.title) +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  
  p
}

# Plot QT7 and HD separately, to contrast sample sizes.
plot_tp_apd_ind <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)),
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*2)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq)),
               decile = rep(rep(qseq, each = nn), 2) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       # labels = c(".1","",".3","",".5","",".7","",".9")
                       labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Power", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) +
    facet_grid(cols = vars(estimator))
  
  p
}

# Plot sample sizes separately, to contrast quantile estimators.
plot_tp_apd_ind_n <- function(simres, plot.title, alpha.val){
  
  df <- tibble(prob = c(as.vector(apply(simres$AB.qt<alpha.val, c(2,3), mean)),
                        as.vector(apply(simres$AB.hd<alpha.val, c(2,3), mean))),
               n = factor(rep(nseq, nq*2)),
               estimator = factor(rep(c("QT7", "HD"), each = nn*nq)),
               decile = rep(rep(qseq, each = nn), 2) 
  )
  
  p <- ggplot(df, aes(x = decile, y = prob, linetype = estimator)) + theme_gar +
    # 0.05 reference line
    geom_abline(intercept = 0.05, slope = 0, colour="black") +  
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = qseq,
                       labels = c(".1","",".3","",".5","",".7","",".9")
                       #labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
                       ) +
    scale_y_continuous(limits = c(0,1)) +
    scale_color_viridis_d(end = 0.9, option = "magma") +
    labs(x = "Quantiles", y = "Power", title = plot.title) +
    theme(panel.grid.minor.x = element_blank()) +
    facet_wrap(vars(n))
  
  p
}
# ----------------------------------------------

# # Plot results for Poisson simulation: 
# # contrast qt7 and hd results
# 
# plot_fp_fwer_pois <- function(fwer, plot.title){
#   
#   df <- tibble(prob = c(apply(fwer$AB.h0.raw.qt, 2, mean), # qt7 
#                         apply(fwer$AB.h0.fdr.qt, 2, mean),
#                         apply(fwer$AB.h1.raw.qt, 2, mean), 
#                         apply(fwer$AB.h1.fdr.qt, 2, mean),
#                         apply(fwer$A.h0.raw.qt, 2, mean), 
#                         apply(fwer$A.h0.fdr.qt, 2, mean),
#                         apply(fwer$A.h1.raw.qt, 2, mean), 
#                         apply(fwer$A.h1.fdr.qt, 2, mean),
#                         apply(fwer$B.h0.raw.qt, 2, mean), 
#                         apply(fwer$B.h0.fdr.qt, 2, mean),
#                         apply(fwer$B.h1.raw.qt, 2, mean), 
#                         apply(fwer$B.h1.fdr.qt, 2, mean),
#                         apply(fwer$AB.h0.raw.hd, 2, mean), # hd 
#                         apply(fwer$AB.h0.fdr.hd, 2, mean),
#                         apply(fwer$AB.h1.raw.hd, 2, mean), 
#                         apply(fwer$AB.h1.fdr.hd, 2, mean),
#                         apply(fwer$A.h0.raw.hd, 2, mean), 
#                         apply(fwer$A.h0.fdr.hd, 2, mean),
#                         apply(fwer$A.h1.raw.hd, 2, mean), 
#                         apply(fwer$A.h1.fdr.hd, 2, mean),
#                         apply(fwer$B.h0.raw.hd, 2, mean), 
#                         apply(fwer$B.h0.fdr.hd, 2, mean),
#                         apply(fwer$B.h1.raw.hd, 2, mean), 
#                         apply(fwer$B.h1.fdr.hd, 2, mean)),
#                n = rep(nseq, 24),
#                Test = factor(rep(rep(c("No corr.", "FDR"), each = nn), 12)),
#                Boot = factor(rep(rep(c("bootstrap-t", "perc. boot."), each = nn*2), 6)),
#                Comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*4), 2)),
#                estimator = factor(rep(c("qt7", "hd"), each = nn*12))
#   )
#   
#   p <- ggplot(df, aes(x = n, y = prob, colour = Test, linetype = estimator)) + theme_gar +
#     # Bradley's (1978) satisfactory range
#     geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
#                 fill = "grey80", colour = "grey80", show.legend=FALSE) +
#     # Bradley's (1978) ideal range
#     geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
#                 fill = "grey90", colour = "grey90", show.legend=FALSE) +
#     # 0.05 reference line
#     geom_abline(intercept = 0.05, slope = 0, colour="black") +  
#     geom_point() +
#     geom_line() +
#     scale_x_continuous(breaks = nseq) +
#     # scale_y_continuous(limits = c(0,0.15)) +
#     scale_color_manual(values=c("#E69F00", "#999999")) +
#     #, "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
#     labs(x = "Sample size", y = "False positives", title = plot.title) +
#     facet_grid(rows = vars(Comp),
#                cols = vars(Boot))
#   
#   p
# }
# 
# # FP: plot only FDR results
# plot_fp_fwer_fdr_pois <- function(fwer, plot.title){
#   
#   df <- tibble(prob = c(apply(fwer$AB.h0.fdr.hd, 2, mean),
#                         apply(fwer$AB.h1.fdr.hd, 2, mean),
#                         apply(fwer$A.h0.fdr.hd, 2, mean),
#                         apply(fwer$A.h1.fdr.hd, 2, mean),
#                         apply(fwer$B.h0.fdr.hd, 2, mean),
#                         apply(fwer$B.h1.fdr.hd, 2, mean),
#                         apply(fwer$AB.h0.fdr.qt, 2, mean),
#                         apply(fwer$AB.h1.fdr.qt, 2, mean),
#                         apply(fwer$A.h0.fdr.qt, 2, mean),
#                         apply(fwer$A.h1.fdr.qt, 2, mean),
#                         apply(fwer$B.h0.fdr.qt, 2, mean),
#                         apply(fwer$B.h1.fdr.qt, 2, mean)),
#                n = rep(nseq, 12),
#                Boot = factor(rep(rep(c("bootstrap-t", "perc. boot."), each = nn), 6)),
#                Comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*2), 2)),
#                estimator = factor(rep(c("hd", "qt7"), each = nn*6)) 
#   )
#   
#   p <- ggplot(df, aes(x = n, y = prob, colour = estimator)) + theme_gar +
#     # Bradley's (1978) satisfactory range
#     geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
#                 fill = "grey80", colour = "grey80", show.legend=FALSE) +
#     # Bradley's (1978) ideal range
#     geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
#                 fill = "grey90", colour = "grey90", show.legend=FALSE) +
#     # 0.05 reference line
#     geom_abline(intercept = 0.05, slope = 0, colour="black") +  
#     geom_point() +
#     geom_line(aes(linetype=Boot)) +
#     scale_x_continuous(breaks = nseq) +
#     # scale_y_continuous(limits = c(0,0.15)) +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     # "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
#     labs(x = "Sample size", y = "False positives", title = plot.title) +
#     facet_grid(rows = vars(Comp))
#   
#   p
# }
# 
# # TP: plot only FDR results
# plot_tp_fwer_fdr_pois <- function(fwer, plot.title){
#   
#   df <- tibble(prob = c(apply(fwer$AB.h0.fdr.hd, 2, mean),
#                         apply(fwer$AB.h1.fdr.hd, 2, mean),
#                         apply(fwer$A.h0.fdr.hd, 2, mean),
#                         apply(fwer$A.h1.fdr.hd, 2, mean),
#                         apply(fwer$B.h0.fdr.hd, 2, mean),
#                         apply(fwer$B.h1.fdr.hd, 2, mean),
#                         apply(fwer$AB.h0.fdr.qt, 2, mean),
#                         apply(fwer$AB.h1.fdr.qt, 2, mean),
#                         apply(fwer$A.h0.fdr.qt, 2, mean),
#                         apply(fwer$A.h1.fdr.qt, 2, mean),
#                         apply(fwer$B.h0.fdr.qt, 2, mean),
#                         apply(fwer$B.h1.fdr.qt, 2, mean)),
#                n = rep(nseq, 12),
#                Boot = factor(rep(rep(c("bootstrap-t", "perc. boot."), each = nn), 6)),
#                Comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*2), 2)),
#                estimator = factor(rep(c("hd", "qt7"), each = nn*6)) 
#   )
#   
#   p <- ggplot(df, aes(x = n, y = prob, colour = estimator)) + theme_gar +
#     geom_point() +
#     geom_line(aes(linetype=Boot)) +
#     scale_x_continuous(breaks = nseq) +
#     scale_color_manual(values=c("#E69F00", "#56B4E9")) +
#     # "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
#     labs(x = "Sample size", y = "True positives", title = plot.title) +
#     facet_grid(rows = vars(Comp))
#   
#   p
# }
# 
# 
# plot_tp_fwpr_pois <- function(fwer, plot.title){
#   
#   df <- tibble(prob = c(apply(fwer$AB.h0.raw.qt, 2, mean), # qt7 
#                         apply(fwer$AB.h0.fdr.qt, 2, mean),
#                         apply(fwer$AB.h1.raw.qt, 2, mean), 
#                         apply(fwer$AB.h1.fdr.qt, 2, mean),
#                         apply(fwer$A.h0.raw.qt, 2, mean), 
#                         apply(fwer$A.h0.fdr.qt, 2, mean),
#                         apply(fwer$A.h1.raw.qt, 2, mean), 
#                         apply(fwer$A.h1.fdr.qt, 2, mean),
#                         apply(fwer$B.h0.raw.qt, 2, mean), 
#                         apply(fwer$B.h0.fdr.qt, 2, mean),
#                         apply(fwer$B.h1.raw.qt, 2, mean), 
#                         apply(fwer$B.h1.fdr.qt, 2, mean),
#                         apply(fwer$AB.h0.raw.hd, 2, mean), # hd 
#                         apply(fwer$AB.h0.fdr.hd, 2, mean),
#                         apply(fwer$AB.h1.raw.hd, 2, mean), 
#                         apply(fwer$AB.h1.fdr.hd, 2, mean),
#                         apply(fwer$A.h0.raw.hd, 2, mean), 
#                         apply(fwer$A.h0.fdr.hd, 2, mean),
#                         apply(fwer$A.h1.raw.hd, 2, mean), 
#                         apply(fwer$A.h1.fdr.hd, 2, mean),
#                         apply(fwer$B.h0.raw.hd, 2, mean), 
#                         apply(fwer$B.h0.fdr.hd, 2, mean),
#                         apply(fwer$B.h1.raw.hd, 2, mean), 
#                         apply(fwer$B.h1.fdr.hd, 2, mean)),
#                n = rep(nseq, 24),
#                Test = factor(rep(rep(c("No corr.", "FDR"), each = nn), 12)),
#                Boot = factor(rep(rep(c("bootstrap-t", "perc. boot."), each = nn*2), 6)),
#                Comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*4), 2)),
#                estimator = factor(rep(c("qt7", "hd"), each = nn*12))
#   )
#   
#   p <- ggplot(df, aes(x = n, y = prob, colour = Test, linetype = estimator)) + theme_gar +
#     geom_point() +
#     geom_line() +
#     scale_x_continuous(breaks = nseq) +
#     scale_y_continuous(limits = c(0,1)) +
#     scale_color_manual(values=c("#E69F00", "#999999")) +
#     #, "#56B4E9", "#009E73", "#0072B2", "#D55E00")) +
#     labs(x = "Sample size", y = "True positives", title = plot.title) +
#     facet_grid(rows = vars(Comp),
#                cols = vars(Boot))
#   
#   p
# }

# plot_fp_ind_qthd <- function(simres, plot.title, alpha.val, qthd){
#   
#   if(qthd == "hd"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h1.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h1.hd<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*6)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq*3)),
#                  comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
#                  decile = rep(rep(qseq, each = nn), 6))
#   } else if(qthd == "qt"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h1.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h1.qt<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*6)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq*3)),
#                  comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
#                  decile = rep(rep(qseq, each = nn), 6))
#   }
#   
#   p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
#     # Bradley's (1978) satisfactory range
#     geom_ribbon(aes(ymin = 0.025, ymax = 0.075), 
#                 fill = "grey80", colour = "grey80", show.legend=FALSE) +
#     # Bradley's (1978) ideal range
#     geom_ribbon(aes(ymin = 0.045, ymax = 0.055), 
#                 fill = "grey90", colour = "grey90", show.legend=FALSE) +
#     # 0.05 reference line
#     geom_abline(intercept = 0.05, slope = 0, colour="black") +  
#     geom_point() +
#     geom_line() +
#     scale_x_continuous(breaks = qseq,
#                        labels = c(".1","",".3","",".5","",".7","",".9")) +
#     scale_y_continuous(limits = c(0,0.15)) +
#     scale_color_viridis_d(end = 0.9, option = "magma") +
#     labs(x = "Quantiles", y = "False positives", title = plot.title) +
#     facet_grid(rows = vars(comp),
#                cols = vars(bootmet))
#   
#   p
# }

# plot_tp_ind_qthd <- function(simres, plot.title, alpha.val, qthd){
#   
#   if(qthd == "hd"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h1.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h1.hd<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*6)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq*3)),
#                  comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
#                  decile = rep(rep(qseq, each = nn), 6))
#   } else if(qthd == "qt"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$A.h1.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$B.h1.qt<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*6)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq*3)),
#                  comp = factor(rep(rep(c("Interaction", "Main A", "Main B"), each = nn*nq),2)),
#                  decile = rep(rep(qseq, each = nn), 6))
#   }
#   
#   p <- ggplot(df, aes(x = decile, y = prob, colour = n)) + theme_gar +
#     geom_point() +
#     geom_line() +
#     scale_x_continuous(breaks = qseq,
#                        labels = c(".1","",".3","",".5","",".7","",".9")) +
#     scale_y_continuous(limits = c(0,1)) +
#     scale_color_viridis_d(end = 0.9, option = "magma") +
#     labs(x = "Quantiles", y = "True positives", title = plot.title) +
#     facet_grid(rows = vars(comp),
#                cols = vars(bootmet))
#   
#   p
# }

# plot_tp_ind_qthd_comp <- function(simres, plot.title, alpha.val, qthd){
#   
#   if(qthd == "hd"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.hd<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.hd<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*2)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq)),
#                  decile = rep(rep(qseq, each = nn), 2))
#   } else if(qthd == "qt"){
#     df <- tibble(prob = c(as.vector(apply(simres$AB.h0.qt<alpha.val, c(2,3), mean)), 
#                           as.vector(apply(simres$AB.h1.qt<alpha.val, c(2,3), mean))),
#                  n = factor(rep(nseq, nq*2)),
#                  bootmet = factor(rep(c("bootstrap-t", "perc. boot."), each = nn*nq)),
#                  decile = rep(rep(qseq, each = nn), 2))
#   }
#   
#   p <- ggplot(df, aes(x = decile, y = prob)) + theme_gar +
#     geom_point() +
#     geom_line(aes(linetype=bootmet)) +
#     scale_x_continuous(breaks = qseq,
#                        labels = c(".1","",".3","",".5","",".7","",".9")) +
#     scale_y_continuous(limits = c(0,1)) +
#     scale_color_viridis_d(end = 0.9, option = "magma") +
#     labs(x = "Quantiles", y = "True positives", title = plot.title) +
#     facet_wrap(vars(n))
#   
#   p
# }