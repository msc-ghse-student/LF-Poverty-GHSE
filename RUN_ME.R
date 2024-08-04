################################# LF & Poverty #################################
############# MSc in Global Health Science & Epidemiology 2023/24 ##############
############################# University of Oxford #############################
## Code adapted from Clark 2023 https://github.com/iamjessclark/non_access_LF ##

# Load in packages
source("Setup.R")

# Set up population characteristics
nPops <- 1 # number of populations (set to 1 if not using movement functionality)
popSize <- 1000 # mean population size
popData <- generatePops(nPops,popSize)
popData$prev <- NA
popData$nation <- 1
dMatrix <- matrix(0,1,1) # only one population, so just a zero matrix

# Generate populations & investigate biological drivers (Supplementary Material)
aggk <- c(0.06, 0.16, 0.27, 0.37, 0.48) # from Irvine, et al.
vths <- seq(0, 160, 20) # from Irvine, et al.
baseline_prevs <- c()
vth_bp <- c()
k_bp <- c()
pops <- c()
baseline_prevs_range <- c()

for (k in aggk) {
  print(k)
  for (vth in vths) {
    print(vth)
    for (i in 1:100) {
      newPop <- Population$new(popData,dMatrix)
      nHost <- newPop$nHost
      newPop$VtH <- rep(vth, nHost)
      newPop$k <- k
      newPop$biteRisk <- rgamma(nHost, shape = k, rate = k)
      
      newPop$burnin()
      prevalence <- length(which(newPop$Mf>0))/newPop$nHosts
      baseline_prevs <- c(baseline_prevs, prevalence)
      vth_bp <- c(vth_bp, vth)
      k_bp <- c(k_bp, k)
      
      if (prevalence < 0.25 & prevalence >= 0.05) {
        pops <- c(pops, newPop)
        baseline_prevs_range <- c(baseline_prevs_range, prevalence)
      }
    }
  }
}

prev_df <- tibble(vth=vth_bp, k=k_bp, mf=baseline_prevs*100)
prev_df$k <- as.factor(prev_df$k)

for_plot <- aggregate(mf ~ vth+k, data = prev_df, FUN = mean)
var <- aggregate(mf ~ vth+k, data = prev_df, FUN = sd)
for_plot$sd <- var$mf

# Figure S1
ggplot(for_plot, aes(x=vth, y=mf, color=k)) + geom_line(linewidth=2) + geom_point(size=3) +
  labs(x="vector-to-host ratio", y="baseline mf prevalence (%)") +
  theme_classic(base_size=20) + geom_errorbar(aes(ymin=mf-sd, ymax=mf+sd), position=position_dodge(0.05)) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=30)) +
  scale_x_continuous(breaks=seq(0, 160, 20)) + scale_y_continuous(breaks=seq(0, 70, 10))

# Impact of Established Factors on Elimination
# Separate into bins by baseline prevalence
bins <- rep(0, length(pops))
bins[baseline_prevs_range < 0.1] = 1
bins[baseline_prevs_range >= 0.1 & baseline_prevs_range < 0.15] = 2
bins[baseline_prevs_range >= 0.15 & baseline_prevs_range < 0.2] = 3
bins[baseline_prevs_range >= 0.2] = 4

# Randomly sample 100 populations from each baseline prevalence category
low_pops <- pops[bins == 1]
low_sample <- sample(1:length(low_pops), 100)
low_pops <- low_pops[low_sample]
med_pops <- pops[bins == 2]
med_sample <- sample(1:length(med_pops), 100)
med_pops <- med_pops[med_sample]
hi_pops <- pops[bins == 3]
hi_sample <- sample(1:length(hi_pops), 100)
hi_pops <- hi_pops[hi_sample]
vhi_pops <- pops[bins == 4]
vhi_sample <- sample(1:length(vhi_pops), 100)
vhi_pops <- vhi_pops[vhi_sample]

comps <- c(0.2, 0.4, 0.6)
coverage <- c(0.65, 0.65)

# Low
res_low <- matrix(nrow=length(low_pops)*length(comps), ncol=5)
low_prevs_0.2 <- list()
low_prevs_0.4 <- list()
low_prevs_0.6 <- list()
j <- 1

for (curr_pop in low_pops) {
  for (comp in comps) {
    # Run simulation on each population with low baseline prevalence for each systematic non-compliance value
    sample <- sample.R(1, curr_pop, comp, coverage)
    
    res_low[j, 1] <- comp
    res_low[j, 2] <- sample$rounds
    res_low[j, 3] <- sample$restart
    res_low[j, 4] <- sample$tas1
    res_low[j, 5] <- sample$elim
    if (comp == 0.2) {
      prevs <- unlist(sample$prevs)
      low_prevs_0.2 <- c(low_prevs_0.2, list(prevs))
    } else if (comp == 0.4) {
      prevs <- unlist(sample$prevs)
      low_prevs_0.4 <- c(low_prevs_0.4, list(prevs))
    } else {
      prevs <- unlist(sample$prevs)
      low_prevs_0.6 <- c(low_prevs_0.6, list(prevs))
    }
    j <- j + 1
  }
}

res_low_df <- data.frame(res_low)
colnames(res_low_df) <- c("comps", "rounds", "restart", "tas1", "elim")
res_low_df$comp <- as.numeric(res_low_df$comp)

# Make categories for first time passing TAS-1
res_low_df$tas1_cat <- as.character(res_low_df$tas1)
res_low_df$tas1_cat[res_low_df$tas1 > 11] <- "> 11"
res_low_df$tas1_cat[res_low_df$tas1 == 0] <- "time-out"

res_low_2 <- data.frame(res_low_df[res_low_df$comps == 0.2,])
res_low_4 <- data.frame(res_low_df[res_low_df$comps == 0.4,])
res_low_6 <- data.frame(res_low_df[res_low_df$comps == 0.6,])

# Make boxplot of number of rounds for each systematic non-compliance value
group1 <- res_low_df[, 1:2]
group1$comps <- as.factor(group1$comps)
syscomp_low <- ggplot(data = group1, aes(x = comps, y = rounds)) +
  geom_boxplot() + labs(x=expression(rho), y="number of treatment rounds to EPHP") +
  theme_classic(base_size=20) + scale_y_continuous(breaks=seq(0, 100, 10), limits = c(0, 101)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Make barplot of proportion of simulations for TAS-1 category
tab1 <- as.data.frame(table(res_low_2$tas1_cat))
colnames(tab1) <- c("rounds", "freq")
tab1$x <- rep("5-10%", nrow(tab1))
tab1$rounds <- factor(tab1$rounds, levels=c("5", "7", "9", "11", "> 11", "time-out"))
baseprev_low <- ggplot(tab1, aes(fill=rounds, y=freq, x=x)) + geom_bar(position="fill", stat="identity") + theme_classic(base_size=20) +
  labs(x="baseline prevalence", y="proportion of simulations", title="Baseline 5-10%") + guides(fill=guide_legend(title="number of treatment rounds\nto pass TAS-1"))

# Medium
res_med <- matrix(nrow=length(med_pops)*length(comps), ncol=5)
med_prevs_0.2 <- list()
med_prevs_0.4 <- list()
med_prevs_0.6 <- list()
j <- 1

for (curr_pop in med_pops) {
  for (comp in comps) {
    # Run simulation on each population with medium baseline prevalence for each systematic non-compliance value
    sample <- sample.R(1, curr_pop, comp, coverage)
    
    res_med[j, 1] <- comp
    res_med[j, 2] <- sample$rounds
    res_med[j, 3] <- sample$restart
    res_med[j, 4] <- sample$tas1
    res_med[j, 5] <- sample$elim
    if (comp == 0.2) {
      prevs <- unlist(sample$prevs)
      med_prevs_0.2 <- c(med_prevs_0.2, list(prevs))
    } else if (comp == 0.4) {
      prevs <- unlist(sample$prevs)
      med_prevs_0.4 <- c(med_prevs_0.4, list(prevs))
    } else {
      prevs <- unlist(sample$prevs)
      med_prevs_0.6 <- c(med_prevs_0.6, list(prevs))
    }
    j <- j + 1
  }
}

res_med_df <- data.frame(res_med)
colnames(res_med_df) <- c("comps", "rounds", "restart", "tas1", "elim")
res_med_df$rounds <- as.numeric(res_med_df$rounds)

# Make categories for first time passing TAS-1
res_med_df$tas1_cat <- as.character(res_med_df$tas1_cat)
res_med_df$tas1_cat[res_med_df$tas1 > 11] <- "> 11"
res_med_df$tas1_cat[res_med_df$tas1 == 0] <- "time-out"

res_med_2 <- data.frame(res_med_df[res_med_df$comps == 0.2,])
res_med_4 <- data.frame(res_med_df[res_med_df$comps == 0.4,])
res_med_6 <- data.frame(res_med_df[res_med_df$comps == 0.6,])

# Make boxplot of number of rounds for each systematic non-compliance value
group2 <- res_med_df[, 1:2]
group2$comps <- as.factor(group2$comps)
syscomp_med <- ggplot(data = group2, aes(x = comps, y = rounds)) +
  geom_boxplot() + labs(x=expression(rho), y="number of treatment rounds to EPHP") +
  theme_classic(base_size=20) + scale_y_continuous(breaks=seq(0, 100, 10), limits = c(0, 101)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Make barplot of proportion of simulations for TAS-1 category
tab2 <- as.data.frame(table(res_med_2$tas1_cat))
colnames(tab2) <- c("rounds", "freq")
tab2$x <- rep("10-15%", nrow(tab2))
tab2$rounds <- factor(tab2$rounds, levels=c("5", "7", "9", "11", "> 11", "restart", "time-out"))
baseprev_med <- ggplot(tab2, aes(fill=rounds, y=freq, x=x)) + geom_bar(position="fill", stat="identity") +
  labs(x="baseline prevalence", y="proportion of simulations", title="Baseline 10-15%") + guides(fill=guide_legend(title="number of treatment rounds\nto pass TAS-1")) + theme_classic(base_size=20)

# High
res_hi <- matrix(nrow=length(hi_pops)*length(comps), ncol=5)
hi_prevs_0.2 <- list()
hi_prevs_0.4 <- list()
hi_prevs_0.6 <- list()
j <- 1

for (curr_pop in hi_pops) {
  for (comp in comps) {
    # Run simulation on each population with high baseline prevalence for each systematic non-compliance value
    sample <- sample.R(1, curr_pop, comp, coverage)
    
    res_hi[j, 1] <- comp
    res_hi[j, 2] <- sample$rounds
    res_hi[j, 3] <- sample$restart
    res_hi[j, 4] <- sample$tas1
    res_hi[j, 5] <- sample$elim
    if (comp == 0.2) {
      prevs <- unlist(sample$prevs)
      hi_prevs_0.2 <- c(hi_prevs_0.2, list(prevs))
    } else if (comp == 0.4) {
      prevs <- unlist(sample$prevs)
      hi_prevs_0.4 <- c(hi_prevs_0.4, list(prevs))
    } else {
      prevs <- unlist(sample$prevs)
      hi_prevs_0.6 <- c(hi_prevs_0.6, list(prevs))
    }
    j <- j + 1
  }
}

res_hi_df <- data.frame(res_hi)
colnames(res_hi_df) <- c("comps", "rounds", "restart", "tas1", "elim")
res_hi_df$rounds <- as.numeric(res_hi_df$rounds)

# Make categories for first time passing TAS-1
res_hi_df$tas1_cat <- as.character(res_hi_df$tas1)
res_hi_df$tas1_cat[res_hi_df$tas1 > 11] <- "> 11"
res_hi_df$tas1_cat[res_hi_df$tas1 == 0] <- "time-out"

res_hi_2 <- data.frame(res_hi_df[res_hi_df$comps == 0.2,])
res_hi_4 <- data.frame(res_hi_df[res_hi_df$comps == 0.4,])
res_hi_6 <- data.frame(res_hi_df[res_hi_df$comps == 0.6,])

# Make boxplot of number of rounds for each systematic non-compliance value
group3 <- res_hi_df[, 1:2]
group3$comps <- as.factor(group3$comps)
syscomp_hi <- ggplot(data = group3, aes(x = comps, y = rounds)) +
  geom_boxplot() + labs(x=expression(rho), y="number of treatment rounds to EPHP") +
  theme_classic(base_size=20) + scale_y_continuous(breaks=seq(0, 100, 10), limits = c(0, 101)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Make barplot of proportion of simulations for TAS-1 category
tab3 <- as.data.frame(table(res_hi_2$tas1_cat))
colnames(tab3) <- c("rounds", "freq")
tab3$x <- rep("15-20%", nrow(tab3))
tab3$rounds <- factor(tab3$rounds, levels=c("5", "7", "9", "11", "> 11", "time-out"))
baseprev_hi <- ggplot(tab3, aes(fill=rounds, y=freq, x=x)) + geom_bar(position="fill", stat="identity") +
  labs(x="baseline prevalence", y="proportion of simulations", title="Baseline 15-20%") + guides(fill=guide_legend(title="number of treatment rounds\nto pass TAS-1")) + theme_classic(base_size=20)

# Very High
res_vhi <- matrix(nrow=length(vhi_pops)*length(comps), ncol=5)
vhi_prevs_0.2 <- list()
vhi_prevs_0.4 <- list()
vhi_prevs_0.6 <- list()
j <- 1

for (curr_pop in vhi_pops) {
  for (comp in comps) {
    # Run simulation on each population with very high baseline prevalence for each systematic non-compliance value
    sample <- sample.R(1, curr_pop, comp, coverage)
    
    res_vhi[j, 1] <- comp
    res_vhi[j, 2] <- sample$rounds
    res_vhi[j, 3] <- sample$restart
    res_vhi[j, 4] <- sample$tas1
    res_vhi[j, 5] <- sample$elim
    if (comp == 0.2) {
      prevs <- unlist(sample$prevs)
      vhi_prevs_0.2 <- c(vhi_prevs_0.2, list(prevs))
    } else if (comp == 0.4) {
      prevs <- unlist(sample$prevs)
      vhi_prevs_0.4 <- c(vhi_prevs_0.4, list(prevs))
    } else {
      prevs <- unlist(sample$prevs)
      vhi_prevs_0.6 <- c(vhi_prevs_0.6, list(prevs))
    }
    j <- j + 1
  }
}

res_vhi_df <- data.frame(res_vhi)
colnames(res_vhi_df) <- c("comps", "rounds", "restart", "tas1", "elim")
res_vhi_df$rounds <- as.numeric(res_vhi_df$rounds)

# Make categories for first time passing TAS-1
res_vhi_df$tas1_cat <- as.character(res_vhi_df$tas1)
res_vhi_df$tas1_cat[res_vhi_df$tas1 > 11] <- "> 11"
res_vhi_df$tas1_cat[res_vhi_df$tas1 == 0] <- "time-out"

res_vhi_2 <- data.frame(res_vhi_df[res_vhi_df$comps == 0.2,])
res_vhi_4 <- data.frame(res_vhi_df[res_vhi_df$comps == 0.4,])
res_vhi_6 <- data.frame(res_vhi_df[res_vhi_df$comps == 0.6,])

# Make boxplot of number of rounds for each systematic non-compliance value
group4 <- res_vhi_df[, 1:2]
group4$comps <- as.factor(group4$comps)
syscomp_vhi <- ggplot(data = group4, aes(x = comps, y = rounds)) +
  geom_boxplot() + labs(x=expression(rho), y="number of treatment rounds to EPHP") +
  theme_classic(base_size=20) + scale_y_continuous(breaks=seq(0, 100, 10), limits = c(0, 101)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Make barplot of proportion of simulations for TAS-1 category
tab4 <- as.data.frame(table(res_vhi_2$tas1_cat))
colnames(tab4) <- c("rounds", "freq")
tab4$x <- rep("20-25%", nrow(tab4))
tab4$rounds <- factor(tab4$rounds, levels=c("5", "7", "9", "11", "> 11", "time-out"))
baseprev_vhi <- ggplot(tab4, aes(fill=rounds, y=freq, x=x)) + geom_bar(position="fill", stat="identity") +
  labs(x="baseline prevalence", y="proportion of simulations", title="Baseline 20-25%") + guides(fill=guide_legend(title="number of treatment rounds\nto pass TAS-1")) + theme_classic(base_size=20)

# Figure 7 - barplots
tot <- rbind(tab1, tab2, tab3, tab4)
tot$x <- factor(tot$x, levels=c("5-10%", "10-15%", "15-20%", "20-25%"))
tot$rounds <- factor(tot$rounds, levels=c("5", "7", "9", "11", "> 11", "time-out"))

blues <- brewer.pal(6, "Blues")
blue_range <- colorRampPalette(blues)
ggplot(tot, aes(fill=rounds, y=freq, x=x)) + geom_bar(position="fill", stat="identity", colour="black") +
  xlab("baseline mf prevalence") + ylab("proportion of simulations") + guides(fill=guide_legend(title="number of treatment rounds\nto first time passing TAS-1")) +
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 25), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = blue_range(6)) + theme_classic(base_size=20) + theme(legend.direction="horizontal")

# Figure 8
bottom <- textGrob(expression(paste(rho, " - systematic non-compliance parameter")), gp = gpar(fontsize = 30))
left <- textGrob("number of treatment rounds to EPHP", gp = gpar(fontsize = 30), rot=90)
grid.arrange(syscomp_low, syscomp_med, syscomp_hi, syscomp_vhi, ncol=4, left=left, bottom=bottom)

# Figure 7 - mf prevalence plots
curr_prevs <- vhi_prevs_0.2 # low_prevs_0.2, med_prevs_0.2, hi_prevs_0.2, vhi_prevs_0.2
prevs_0.2_sample <- sample(1:100, 5) # Get random sample of 5 simulations to plot
prevs_0.2 <- curr_prevs[prevs_0.2_sample]

prevs_0.2_1 <- unlist(prevs_0.2[1])
prevs_0.2_2 <- unlist(prevs_0.2[2])
prevs_0.2_3 <- unlist(prevs_0.2[3])
prevs_0.2_4 <- unlist(prevs_0.2[4])
prevs_0.2_5 <- unlist(prevs_0.2[5])
# Pad the prevalences to get same axis
max_length <- 1340 # Maximum time is ~ 111 years (1332 months)
prevs_0.2_1 <- c(prevs_0.2_1, rep(NA, max_length-length(prevs_0.2_1)))
prevs_0.2_2 <- c(prevs_0.2_2, rep(NA, max_length-length(prevs_0.2_2)))
prevs_0.2_3 <- c(prevs_0.2_3, rep(NA, max_length-length(prevs_0.2_3)))
prevs_0.2_4 <- c(prevs_0.2_4, rep(NA, max_length-length(prevs_0.2_4)))
prevs_0.2_5 <- c(prevs_0.2_5, rep(NA, max_length-length(prevs_0.2_5)))

ggplot() +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_1))/12, y = prevs_0.2_1), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_1))/12, y = prevs_0.2_1), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_2))/12, y = prevs_0.2_2), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_3))/12, y = prevs_0.2_3), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_4))/12, y = prevs_0.2_4), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_0.2_5))/12, y = prevs_0.2_5), linewidth=0.8, color='black') +
  xlab("time (years)") + ylab("mf prevalence in adults aged 20+ (%)") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1, size=1.5, color="grey") + scale_y_continuous(breaks=seq(0, 30, 5), limits=c(0, 30)) + scale_x_continuous(breaks=seq(0, 100, 10), limits=c(0, 110))

# Figure S2
#prevs_comp_sample <- sample(1:100, 5) # Get random sample of 5 simulations to plot
curr_comp <- hi_prevs_0.6 # hi_prevs_0.2, hi_prevs_0.4, hi_prevs_0.6
prevs_comp <- curr_comp[prevs_comp_sample]

prevs_comp_1 <- unlist(prevs_comp[1])
prevs_comp_2 <- unlist(prevs_comp[2])
prevs_comp_3 <- unlist(prevs_comp[3])
prevs_comp_4 <- unlist(prevs_comp[4])
prevs_comp_5 <- unlist(prevs_comp[5])
# Pad the prevalences to get same axis
max_length <- 1340 # Maximum time is ~ 111 years (1332 months)
prevs_comp_1 <- c(prevs_comp_1, rep(NA, max_length-length(prevs_comp_1)))
prevs_comp_2 <- c(prevs_comp_2, rep(NA, max_length-length(prevs_comp_2)))
prevs_comp_3 <- c(prevs_comp_3, rep(NA, max_length-length(prevs_comp_3)))
prevs_comp_4 <- c(prevs_comp_4, rep(NA, max_length-length(prevs_comp_4)))
prevs_comp_5 <- c(prevs_comp_5, rep(NA, max_length-length(prevs_comp_5)))

ggplot() +
  geom_line(mapping = aes(x = (1:length(prevs_comp_1))/12, y = prevs_comp_1), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_comp_2))/12, y = prevs_comp_2), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_comp_3))/12, y = prevs_comp_3), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_comp_4))/12, y = prevs_comp_4), linewidth=0.8, color='black') +
  geom_line(mapping = aes(x = (1:length(prevs_comp_5))/12, y = prevs_comp_5), linewidth=0.8, color='black') +
  xlab("time (years)") + ylab("mf prevalence in adults aged 20+ (%)") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1, size=1.5, color="grey") + scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0, 25)) + scale_x_continuous(breaks=seq(0, 100, 10), limits=c(0, 110))

# Impact of Socioeconomic Inequalities on Elimination
# Tanga calibration: heterogeneous
# Step 1: Wide range
vths <- seq(10, 160, 5)
calibrate_prev_het <- c()
calibrate_low_het <- c()
calibrate_high_het <- c()
calibrate_vth1_het <- c()
calibrate_vth2_het <- c()

for (i in 1:(length(vths)-1)) {
  vth1 <- vths[i]
  print(vth1)
  for (j in (i+1):length(vths)) {
    vth2 <- vths[j]
    print(vth2)
    collect_prev <- c()
    collect_low <- c()
    collect_high <- c()
    for (m in 1:25) {
      newPop <- Population$new(popData,dMatrix)
      nHost <- newPop$nHost
      VtHprob <- rbinom(nHost, 1, 0.4)
      newPop$VtH[which(VtHprob == 0)] <- vth1
      newPop$VtH[which(VtHprob == 1)] <- vth2
      
      newPop$k <- 0.1
      newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
      
      newPop$burnin()
      prevalence <- length(which(newPop$Mf>0))/newPop$nHosts
      low_prev <- length(which(newPop$Mf[newPop$VtH == vth1]>0))/length(which(newPop$VtH == vth1))
      high_prev <- length(which(newPop$Mf[newPop$VtH == vth2]>0))/length(which(newPop$VtH == vth2))
      collect_prev <- c(collect_prev, prevalence)
      collect_low <- c(collect_low, low_prev)
      collect_high <- c(collect_high, high_prev)
    }
    calibrate_prev_het <- c(calibrate_prev_het, median(collect_prev))
    calibrate_low_het <- c(calibrate_low_het, median(collect_low))
    calibrate_high_het <- c(calibrate_high_het, median(collect_high))
    calibrate_vth1_het <- c(calibrate_vth1_het, vth1)
    calibrate_vth2_het <- c(calibrate_vth2_het, vth2)
  }
}

calibrate_agg <- data.frame(all=calibrate_prev_het, low=calibrate_low_het, high=calibrate_high_het, vth1=calibrate_vth1_het, vth2=calibrate_vth2_het)

calibrate_agg$all_bin = "outside target range"
calibrate_agg$all_bin[calibrate_agg$all < 0.25 & calibrate_agg$all >= 0.24] = "within target range"
calibrate_agg$low_bin = "outside target range"
calibrate_agg$low_bin[calibrate_agg$low < 0.21 & calibrate_agg$low >= 0.2] = "within target range"
calibrate_agg$high_bin = "outside target range"
calibrate_agg$high_bin[calibrate_agg$high < 0.31 & calibrate_agg$high >= 0.30] = "within target range"

calibrate_agg$all_bin <- as.factor(calibrate_agg$all_bin)
calibrate_agg$low_bin <- as.factor(calibrate_agg$low_bin)
calibrate_agg$high_bin <- as.factor(calibrate_agg$high_bin)

# Figure 9
ggplot(calibrate_agg, mapping=aes(x=vth1, y=vth2)) + geom_point(aes(colour=all_bin), size=6) +
  xlab("low-risk vector-to-host ratio") + ylab("high-risk vector-to-host ratio") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20),
        legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c('grey', '#619CFF')) + scale_y_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) +
  scale_x_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) + theme(legend.position="none")

ggplot(calibrate_agg, mapping=aes(x=vth1, y=vth2)) + geom_point(aes(colour=high_bin), size=6) +
  xlab("low-risk vector-to-host ratio") + ylab("high-risk vector-to-host ratio") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20),
        legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c('grey', '#F8766D')) + scale_y_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) +
  scale_x_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) + theme(legend.position="none")

ggplot(calibrate_agg, mapping=aes(x=vth1, y=vth2)) + geom_point(aes(colour=low_bin), size=6) +
  xlab("low-risk vector-to-host ratio") + ylab("high-risk vector-to-host ratio") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20),
        legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c('grey', '#00BA38')) + scale_y_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) +
  scale_x_continuous(breaks=seq(10, 160, 10), limits=c(10, 160)) + theme(legend.position="none")

# Step 2: Narrow range
vths_low <- seq(20, 30, 1)
vths_high <- seq(100, 120, 2)
calibrate_prev_fine <- c()
calibrate_low_fine <- c()
calibrate_high_fine <- c()
calibrate_vth1_fine <- c()
calibrate_vth2_fine <- c()

for (i in 1:length(vths_low)) {
  vth1 <- vths_low[i]
  print(vth1)
  for (j in 1:length(vths_high)) {
    vth2 <- vths_high[j]
    print(vth2)
    collect_prev <- c()
    collect_low <- c()
    collect_high <- c()
    for (m in 1:25) {
      newPop <- Population$new(popData,dMatrix)
      nHost <- newPop$nHost
      VtHprob <- rbinom(nHost, 1, 0.4)
      newPop$VtH[which(VtHprob == 0)] <- vth1
      newPop$VtH[which(VtHprob == 1)] <- vth2
      
      newPop$k <- 0.1
      newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
      
      newPop$burnin()
      prevalence <- length(which(newPop$Mf>0))/newPop$nHosts
      low_prev <- length(which(newPop$Mf[newPop$VtH == vth1]>0))/length(which(newPop$VtH == vth1))
      high_prev <- length(which(newPop$Mf[newPop$VtH == vth2]>0))/length(which(newPop$VtH == vth2))
      collect_prev <- c(collect_prev, prevalence)
      collect_low <- c(collect_low, low_prev)
      collect_high <- c(collect_high, high_prev)
    }
    calibrate_prev_fine <- c(calibrate_prev_fine, median(collect_prev))
    calibrate_low_fine <- c(calibrate_low_fine, median(collect_low))
    calibrate_high_fine <- c(calibrate_high_fine, median(collect_high))
    calibrate_vth1_fine <- c(calibrate_vth1_fine, vth1)
    calibrate_vth2_fine <- c(calibrate_vth2_fine, vth2)
  }
}

calibrate_agg_fine <- data.frame(all=calibrate_prev_fine, low=calibrate_low_fine, high=calibrate_high_fine, vth1=calibrate_vth1_fine, vth2=calibrate_vth2_fine)

calibrate_agg_fine$all_bin = "outside target range"
calibrate_agg_fine$all_bin[calibrate_agg_fine$all < 0.25 & calibrate_agg_fine$all >= 0.24] = "within target range"
calibrate_agg_fine$low_bin = "outside target range"
calibrate_agg_fine$low_bin[calibrate_agg_fine$low < 0.21 & calibrate_agg_fine$low >= 0.2] = "within target range"
calibrate_agg_fine$high_bin = "outside target range"
calibrate_agg_fine$high_bin[calibrate_agg_fine$high < 0.31 & calibrate_agg_fine$high >= 0.30] = "within target range"

# Step 3: Pinpoint
finalists <- calibrate_agg_fine[which(calibrate_agg_fine$all_bin == "within target range" & calibrate_agg_fine$low_bin == "within target range" & calibrate_agg_fine$high_bin == "within target range"),]
finalists$diffs <- abs(finalists$all-0.245) + abs(finalists$low-0.207) + abs(finalists$high-0.302)
min(finalists$diffs)

# Table 2
vth1 <- 23
vth2 <- 100
newPop <- Population$new(popData,dMatrix)
nHost <- newPop$nHost
VtHprob <- rbinom(nHost, 1, 0.4)
newPop$VtH[which(VtHprob == 0)] <- vth1
newPop$VtH[which(VtHprob == 1)] <- vth2
newPop$k <- 0.1
newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
newPop$burnin()
length(which(newPop$Mf>0))/newPop$nHosts
length(which(newPop$Mf[newPop$VtH == vth1]>0))/length(which(newPop$VtH == vth1))
length(which(newPop$Mf[newPop$VtH == vth2]>0))/length(which(newPop$VtH == vth2))

tanga_het <- newPop

# Tanga calibration: homogeneous
# Step 1: Narrow range
vths <- seq(30, 50, 1) # range based on Figure 9
calibrate_prev_hom <- c()
calibrate_vth_hom <- c()

for (i in 1:length(vths)) {
  vth1 <- vths[i]
  vth2 <- vth1
  print(vth1)
  collect_prev <- c()
  collect_pop <- c()
  for (m in 1:25) {
    newPop <- Population$new(popData,dMatrix)
    nHost <- newPop$nHost
    VtHprob <- rbinom(nHost, 1, 0.4)
    newPop$VtH[which(VtHprob == 0)] <- vth1
    newPop$VtH[which(VtHprob == 1)] <- vth2
    
    newPop$k <- 0.1
    newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
    
    newPop$burnin()
    prevalence <- length(which(newPop$Mf>0))/newPop$nHosts
    collect_prev <- c(collect_prev, prevalence)
    collect_pop <- c(collect_pop, newPop)
  }
  calibrate_prev_hom <- c(calibrate_prev_hom, median(collect_prev))
  calibrate_vth_hom <- c(calibrate_vth_hom, vth1)
}

calibrate_hom <- data.frame(all=calibrate_prev_hom, vth=calibrate_vth_hom)

calibrate_hom$all_bin = "outside target range"
calibrate_hom$all_bin[calibrate_hom$all < 0.25 & calibrate_homo$all >= 0.24] = "within target range"

calibrate_hom$all_bin <- as.factor(calibrate_hom$all_bin)

# Figure 10
ggplot(calibrate_homo, mapping=aes(x=vth, y=all*100)) + geom_point(aes(colour=all_bin), size=6) +
  xlab("vector-to-host ratio") + ylab("mf baseline prevalence") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20),
        legend.text = element_text(size=20), legend.title = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c('grey', '#619CFF')) + scale_y_continuous(breaks=seq(16, 30, 2), limits=c(16, 30)) +
  scale_x_continuous(breaks=seq(30, 50, 5), limits=c(30, 50)) + theme(legend.position="none")

# Table 2
finalists_hom <- calibrate_hom[which(calibrate_hom$all_bin == "within target range"),]
finalists_hom$diffs <- abs(finalists_hom$all-0.245)
min(finalists_hom$diffs)

newPop <- Population$new(popData,dMatrix)
newPop$k <- 0.1
nHost <- newPop$nHost
newPop$VtH <- rep(42, nHost)
newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
newPop$burnin()
length(which(newPop$Mf>0))/newPop$nHosts

tanga_hom <- newPop

# Tanga analysis
tanga_hom_res <- sample.R(100, tanga_hom, 0.2, c(0.65, 0.65))
tanga_het_res <- sample.R(100, tanga_het, 0.2, c(0.65, 0.65))
summary(tanga_hom_res$rounds)
summary(tanga_het_res$rounds)

# Get populations with median number of rounds to EPHP (for prototypical illustration)
hom_sample <- sample(which(tanga_hom_res$rounds == 21), 1)
het_sample <- sample(which(tanga_het_res$rounds == 47 | tanga_het_res$rounds == 49), 1)
temp <- tanga_hom_res$prevs
homo_prevs <- data.frame(all=unlist(temp[hom_sample]))
temp <- tanga_het_res$prevs
temp_low <- tanga_het_res$low
temp_hi <- tanga_het_res$hi
het_prevs <- data.frame(all=unlist(temp[het_sample]), low=unlist(temp_low[het_sample]), hi=unlist(temp_hi[het_sample]))

# Figure 11
ggplot() +
  geom_line(data = het_prevs, mapping = aes(x = (1:nrow(het_prevs))/12, y = all, color="All"), linewidth=1.5) +
  geom_line(data = het_prevs, mapping = aes(x = (1:nrow(het_prevs))/12, y = low, color="Low risk"), linewidth=1.5) +
  geom_line(data = het_prevs, mapping = aes(x = (1:nrow(het_prevs))/12, y = hi, color="High risk"), linewidth=1.5) +
  xlab("time (years)") + ylab("mf prevalence in adults aged 20+ (%)") + guides(color=guide_legend(title="Group")) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1, size=1.5, color="grey") + scale_color_manual(values=c('#619CFF', '#F8766D', '#00BA38')) +
  geom_vline(xintercept=(nrow(hetero_prevs)-120-24)/12, size=1.5, color="orange") + scale_x_continuous(breaks=seq(0, 65, 10), limits=c(0, 65))

ggplot() +
  geom_line(data = hom_prevs, mapping = aes(x = (1:nrow(hom_prevs))/12, y = all, color="All"), size=2) +
  xlab("time (years)") + ylab("mf prevalence in adults aged 20+ (%)") +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1, size=1.5, color="grey") + scale_color_manual(values=c('#619CFF')) + guides(color=guide_legend(title="Group")) +
  geom_vline(xintercept=(nrow(homo_prevs)-120-24)/12, size=1.5, color="orange") + scale_x_continuous(breaks=seq(0, 65, 10), limits=c(0, 65))

# Equal coverage in both groups
coverage50_eq <- sample.R(100, tanga_het, 0.2, c(0.5, 0.5))
coverage65_eq <- sample.R(100, tanga_het, 0.2, c(0.65, 0.65))
coverage80_eq <- sample.R(100, tanga_het, 0.2, c(0.8, 0.8))
coverage95_eq <- sample.R(100, tanga_het, 0.2, c(0.95, 0.95))

# Higher coverage in low-risk group
coverage65_low <- sample.R(100, tanga_het, 0.2, c(0.65, 0.5))
coverage80_low <- sample.R(100, tanga_het, 0.2, c(0.8, 0.5))
coverage95_low <- sample.R(100, tanga_het, 0.2, c(0.95, 0.5))

# Higher coverage in high-risk group
coverage65_hi <- sample.R(100, tanga_het, 0.2, c(0.5, 0.65))
coverage80_hi <- sample.R(100, tanga_het, 0.2, c(0.5, 0.8))
coverage95_hi <- sample.R(100, tanga_het, 0.2, c(0.5, 0.95))

# Table 3
equals <- c(median(coverage50_eq$rounds), median(coverage65_eq$rounds), median(coverage80_eq$rounds), median(coverage95_eq$rounds))
equals_cov <- c(50, 65, 80, 95)
lows <- c(median(coverage50_eq$rounds), median(coverage65_low$rounds), median(coverage80_low$rounds), median(coverage95_low$rounds))
lows_cov <- c(50, 59, 68, 77)
highs <- c(median(coverage50_eq$rounds), median(coverage65_hi$rounds), median(coverage80_hi$rounds), median(coverage95_hi$rounds))
highs_cov <- c(50, 56, 62, 68)
altogether <- data.frame(equals=equals, lows=lows, highs=highs, equals_cov=equals_cov, lows_cov=lows_cov, highs_cov=highs_cov)

# Figure S3
ggplot() +
  geom_line(data=altogether, mapping=aes(x=equals_cov, y=equals, color="equal coverage"), linewidth=2) +
  geom_point(data=altogether, aes(x=equals_cov, y=equals, color="equal coverage"), size=6) +
  geom_line(data=altogether, mapping=aes(x=lows_cov, y=lows, color="higher coverage in low-risk group"), linewidth=2) +
  geom_point(data=altogether, aes(x=lows_cov, y=lows, color="higher coverage in low-risk group"), size=6) +
  geom_line(data=altogether, mapping=aes(x=highs_cov, y=highs, color="higher coverage in high-risk group"), linewidth=2) +
  geom_point(data=altogether, aes(x=highs_cov, y=highs, color="higher coverage in high-risk group"), size=6) +
  xlab("overall population treatment coverage (%)") + ylab("number of treatmenr rounds to EPHP") + guides(color=guide_legend(title="treatment coverage scheme")) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 20), legend.text = element_text(size=20), legend.title = element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c('#619CFF', '#F8766D', '#00BA38')) + scale_y_continuous(breaks=seq(10, 70, 10), limits=c(10, 75))

# General analysis
# Generate heterogeneous populations, aiming to get >= 100 in each baseline prevalence category
vths <- seq(10, 160, 5)
het_pops <- c()
het_prevs <- c()

for (i in 1:(length(vths)-1)) {
  vth1 <- vths[i]
  print(vth1)
  for (j in (i+1):length(vths)) {
    vth2 <- vths[j]
    print(vth2)
    for (m in 1:40) {
      newPop <- Population$new(popData,dMatrix)
      nHost <- newPop$nHost
      VtHprob <- rbinom(nHost, 1, 0.4)
      newPop$VtH[which(VtHprob == 0)] <- vth1
      newPop$VtH[which(VtHprob == 1)] <- vth2
      
      newPop$k <- 0.1
      newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)
      
      newPop$burnin()
      prevalence <- length(which(newPop$Mf>0))/newPop$nHosts
      low_prev <- length(which(newPop$Mf[newPop$VtH == vth1]>0))/length(which(newPop$VtH == vth1))
      high_prev <- length(which(newPop$Mf[newPop$VtH == vth2]>0))/length(which(newPop$VtH == vth2))
      ratio <- high_prev / low_prev
      if (!is.nan(ratio)) {
        if (ratio >= 1.5) {
          if (prevalence >= 0.05 & prevalence < 0.25) {
            het_pops <- c(het_pops, newPop)
            het_prevs <- c(het_prevs, prevalence)
          }
        }
      }
    }
  }
}

# Separate into bins
bins <- rep(0, length(het_pops))
bins[het_prevs < 0.1] = 1
bins[het_prevs >= 0.1 & het_prevs < 0.15] = 2
bins[het_prevs >= 0.15 & het_prevs < 0.2] = 3
bins[het_prevs >= 0.2] = 4

low_pops_het <- het_pops[bins == 1]
low_sample_het <- sample(1:length(low_pops_het), 100)
low_pops_het <- low_pops_het[low_sample_het]
med_pops_het <- het_pops[bins == 2]
med_sample_het <- sample(1:length(med_pops_het), 100)
med_pops_het <- med_pops_het[med_sample_het]
hi_pops_het <- het_pops[bins == 3]
hi_sample_het <- sample(1:length(hi_pops_het), 100)
hi_pops_het <- hi_pops_het[hi_sample_het]
vhi_pops_het <- het_pops[bins == 4]
vhi_sample_het <- sample(1:length(vhi_pops_het), 100)
vhi_pops_het <- vhi_pops_het[vhi_sample_het]

# Run simulations
res_low_het <- matrix(nrow=length(low_pops_het), ncol=2)
j <- 1
for (curr_pop in low_pops_het) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_low_het[j, 1] <- sample$rounds
  res_low_het[j, 2] <- sample$elim
  j <- j + 1
}
res_low_df_het <- as.data.frame(res_low_het)

res_med_het <- matrix(nrow=length(med_pops_het), ncol=2)
j <- 1
for (curr_pop in med_pops_het) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_med_het[j, 1] <- sample$rounds
  res_med_het[j, 2] <- sample$elim
  j <- j + 1
}
res_med_df_het <- as.data.frame(res_med_het)

res_hi_het <- matrix(nrow=length(hi_pops_het), ncol=2)
j <- 1
for (curr_pop in hi_pops_het) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_hi_het[j, 1] <- sample$rounds
  res_hi_het[j, 2] <- sample$elim
  j <- j + 1
}
res_hi_df_het <- as.data.frame(res_hi_het)

res_vhi_het <- matrix(nrow=length(vhi_pops_het), ncol=2)
j <- 1
for (curr_pop in vhi_pops_het) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_vhi_het[j, 1] <- sample$rounds
  res_vhi_het[j, 2] <- sample$elim
  j <- j + 1
}
res_vhi_df_het <- as.data.frame(res_vhi_het)

# Generate homogeneous populations, aiming to get >= 100 in each baseline prevalence category
vths <- seq(20, 60, 1)
hom_pops <- c()
hom_prevs <- c()

for (i in 1:length(vths)) {
  vth1 <- vths[i]
  print(vth1)
  vth2 <- vth1
  for (m in 1:100) {
    newPop <- Population$new(popData,dMatrix)
    nHost <- newPop$nHost
    VtHprob <- rbinom(nHost, 1, 0.4)
    newPop$VtH[which(VtHprob == 0)] <- vth1
    newPop$VtH[which(VtHprob == 1)] <- vth2

    newPop$k <- 0.1
    newPop$biteRisk <- rgamma(nHost, shape = newPop$k, rate = newPop$k)

    newPop$burnin()
    prevalence <- length(which(newPop$Mf>0))/newPop$nHosts

    if (prevalence < 0.25 & prevalence >= 0.05) {
      hom_pops <- c(hom_pops, newPop)
      hom_prevs <- c(hom_prevs, prevalence)
    }
  }
}

# Separate into bins
bins <- rep(0, length(hom_pops))
bins[hom_prevs < 0.1] = 1
bins[hom_prevs >= 0.1 & hom_prevs < 0.15] = 2
bins[hom_prevs >= 0.15 & hom_prevs < 0.2] = 3
bins[hom_prevs >= 0.2] = 4

low_pops_hom <- hom_pops[bins == 1]
low_sample_hom <- sample(1:length(low_pops_hom), 100)
low_pops_hom <- low_pops_hom[low_sample_hom]
med_pops_hom <- hom_pops[bins == 2]
med_sample_hom <- sample(1:length(med_pops_hom), 100)
med_pops_hom <- med_pops_hom[med_sample_hom]
hi_pops_hom <- hom_pops[bins == 3]
hi_sample_hom <- sample(1:length(hi_pops_hom), 100)
hi_pops_hom <- hi_pops_hom[hi_sample_hom]
vhi_pops_hom <- hom_pops[bins == 4]
vhi_sample_hom <- sample(1:length(vhi_pops_hom), 100)
vhi_pops_hom <- vhi_pops_hom[vhi_sample_hom]

# Run simulations
res_low_hom <- matrix(nrow=length(low_pops_hom), ncol=2)
j <- 1
for (curr_pop in low_pops_hom) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_low_hom[j, 1] <- sample$rounds
  res_low_hom[j, 2] <- sample$elim
  j <- j + 1
}
res_low_df_hom <- as.data.frame(res_low_hom)

res_med_hom <- matrix(nrow=length(med_pops_hom), ncol=2)
j <- 1
for (curr_pop in med_pops_hom) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_med_hom[j, 1] <- sample$rounds
  res_med_hom[j, 2] <- sample$elim
  j <- j + 1
}
res_med_df_hom <- as.data.frame(res_med_hom)

res_hi_hom <- matrix(nrow=length(hi_pops_hom), ncol=2)
j <- 1
for (curr_pop in hi_pops_hom) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_hi_hom[j, 1] <- sample$rounds
  res_hi_hom[j, 2] <- sample$elim
  j <- j + 1
}
res_hi_df_hom <- as.data.frame(res_hi_hom)

res_vhi_hom <- matrix(nrow=length(vhi_pops_hom), ncol=2)
j <- 1
for (curr_pop in vhi_pops_hom) {
  sample <- sample.R(1, curr_pop, 0.2, c(0.65, 0.65))
  
  res_vhi_hom[j, 1] <- sample$rounds
  res_vhi_hom[j, 2] <- sample$elim
  j <- j + 1
}
res_vhi_df_hom <- as.data.frame(res_vhi_hom)

# Table 4
summary(res_low_df_het$V1)
summary(res_low_df_hom$V1)
summary(res_med_df_het$V1)
summary(res_med_df_hom$V1)
summary(res_hi_df_het$V1)
summary(res_hi_df_hom$V1)
summary(res_vhi_df_het$V1)
summary(res_vhi_df_hom$V1)

# Perform Mann-Whitney U test, comparing number of rounds between heterogeneous and homogeneous populations in each baseline prevalence category
wilcox.test(res_low_df_het$V1, res_low_df_hom$V1)
wilcox.test(res_med_df_het$V1, res_med_df_hom$V1)
wilcox.test(res_hi_df_het$V1, res_hi_df_hom$V1)
wilcox.test(res_vhi_df_het$V1, res_vhi_df_hom$V1)

# Sensitivity analysis of model stochasticity
# Table S7
for (i in 1:10) {
  tanga_hom_res_test <- sample.R(100, tanga_hom, 0.2, c(0.65, 0.65))
  tanga_het_res_test <- sample.R(100, tanga_het, 0.2, c(0.65, 0.65))
  print(i)
  print(summary(tanga_hom_res_test$rounds))
  print(summary(tanga_het_res_test$rounds))
}

############################### END OF ANALYSIS ################################