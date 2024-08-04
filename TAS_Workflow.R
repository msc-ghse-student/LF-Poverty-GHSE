# Main function to implement treatment and TAS
callrunsMDA <- function(newPop, compliance, coverage){
  # Initialise variables
  HALT <- FALSE
  MDArounds <- 0
  maxTime <- 100
  month_total <- 1
  prevs <- c()
  low_prevs <- c()
  high_prevs <- c()
  simPop <- newPop$clone(deep=T)
  simPop$bnCov <- 0
  restart <- 0
  tas1 <- 0
  tas1_flag <- 0
  
  # Run initial FIVE rounds of MDA and TAS-1
  res <- TAS1(simPop, MDArounds, month_total, coverage, compliance)
  prevs <- c(prevs, res$prevs)
  low_prevs <- c(low_prevs, res$low_prevs)
  high_prevs <- c(high_prevs, res$high_prevs)
  if (res$EPHP == FALSE) {
    count <- 0
  } else {
    count <- 1
    tas1 <- 5
    tas1_flag <- 1
  }
  
  while (HALT == FALSE) {
    if (res$EPHP == FALSE) {
      # Run TWO more rounds of MDA if failed TAS
      if (res$time <= maxTime*12){
        res <- goForTwo(res$pop, res$rounds, res$time, coverage, compliance)
        prevs <- c(prevs, res$prevs)
        low_prevs <- c(low_prevs, res$low_prevs)
        high_prevs <- c(high_prevs, res$high_prevs)
      }
      if (res$EPHP == FALSE) {
        count <- 0
      } else {
        count <- 1
      }
      # If previously passed TAS-1 but failed most recent TAS, then MDA had to be restarted
      if (tas1_flag == 1) {
        restart <- 1
      }
    } else {
      # Record first time passing TAS-1
      if (tas1_flag == 0) {
        tas1 <- res$rounds
        tas1_flag <- 1
      }
      # Stop MDA for TWO years if passed TAS
      if (res$time <= maxTime*12){
        res <- stopForTwo(res$pop, res$rounds, res$time, coverage, compliance)
        prevs <- c(prevs, res$prevs)
        low_prevs <- c(low_prevs, res$low_prevs)
        high_prevs <- c(high_prevs, res$high_prevs)
      }
      if (res$EPHP == FALSE) {
        count <- 0
        restart <- 1
      } else {
        count <- count + 1
      }
    }
    
    # If passed TAS-1, TAS-2, and TAS-3, then EPHP!
    if (count == 3) {
      rounds <- res$rounds
      # Observe what happens over the next 10 years without MDA or TAS
      res <- runOut(simPop, res$time)
      prevs <- c(prevs, res$prevs)
      low_prevs <- c(low_prevs, res$low_prevs)
      high_prevs <- c(high_prevs, res$high_prevs)
      elim <- res$elim
      HALT = TRUE
    }
    # If maxed out on time, end simulation
    if (HALT == FALSE & res$time > maxTime*12) {
      rounds <- res$rounds
      elim <- 0
      HALT = TRUE
    }
  }
  return(list(pop=res$pop, prevs=prevs, rounds=rounds, time=res$time, restart=restart, tas1=tas1, elim=elim, low=low_prevs, high=high_prevs))
}

# Helper function to implement MDA for two more years
goForTwo <- function(simPop, MDArounds, month_total, coverage, compliance) {
  # Initialise variables
  prevalence <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
  prevs <- prevalence
  month <- 6
  treatment = "IA"
  EPHP <- FALSE
  HALT <- FALSE
  vth_low <- min(simPop$VtH)
  vth_high <- max(simPop$VtH)
  low_prevs <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
  high_prevs <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
  
  # Execute dynamics every month and MDA every year
  while(HALT == FALSE){
    simPop$runTimestep()
    
    if(month==12|month==24){
      simPop$runMDA(coverage=coverage, drug=treatment, compliance=compliance)
      MDArounds <- MDArounds + 1
    }
    
    mfprev <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
    month <- month + 1
    month_total <- month_total + 1
    prevs[month-5] <- mfprev
    low_prevs[month-5] <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
    high_prevs[month-5] <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
    
    # Wait for 6 months after last MDA to check prevalence
    if(mfprev <= 1 && month==30){
      EPHP <- TRUE
    }
    if(month==30){
      HALT <- TRUE
    }
  }
  
  return(list(pop=simPop, EPHP=EPHP, prevs=prevs[2:length(prevs)], rounds=MDArounds, time=month_total, low_prevs=low_prevs[2:length(low_prevs)], high_prevs=high_prevs[2:length(high_prevs)]))
}

# Helper function to stop MDA for two years
stopForTwo <- function(simPop, MDArounds, month_total, coverage, compliance) {
  # Initialise variables
  prevalence <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
  prevs <- prevalence
  month <- 0
  treatment = "IA"
  EPHP <- FALSE
  HALT <- FALSE
  vth_low <- min(simPop$VtH)
  vth_high <- max(simPop$VtH)
  low_prevs <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
  high_prevs <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
  
  # Execute dynamics every month
  while(HALT == FALSE){
    simPop$runTimestep()
    
    mfprev <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
    month <- month + 1
    month_total <- month_total + 1
    prevs[month+1] <- mfprev
    low_prevs[month+1] <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
    high_prevs[month+1] <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
    
    if(mfprev <= 1 && month==24){
      EPHP <- TRUE
    }
    if(month==24){
      HALT <- TRUE
    }
  }
  
  return(list(pop=simPop, EPHP=EPHP, prevs=prevs[2:length(prevs)], rounds=MDArounds, time=month_total, low_prevs=low_prevs[2:length(low_prevs)], high_prevs=high_prevs[2:length(high_prevs)]))
}

# Helper function to implement MDA for initial five years
TAS1 <- function(simPop, MDArounds, month_total, coverage, compliance) {
  prevalence <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
  prevs <- prevalence
  month <- 1
  treatment = "IA"
  vth_low <- min(simPop$VtH)
  vth_high <- max(simPop$VtH)
  low_prevs <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
  high_prevs <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
  
  EPHP <- FALSE
  HALT <- FALSE
  
  # Execute dynamics every month and MDA every year
  while(HALT == FALSE){
    simPop$runTimestep()
    
    if(month==12|month==24|month==36|month==48|month==60){
      simPop$runMDA(coverage=coverage, drug=treatment, compliance=compliance)
      MDArounds <- MDArounds + 1
    }
    
    mfprev <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
    month <- month + 1
    month_total <- month_total + 1
    prevs[month] <- mfprev
    low_prevs[month] <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
    high_prevs[month] <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
    
    # Wait for 6 months after last MDA to check prevalence
    if(mfprev <= 1 && month==66){
      EPHP <- TRUE
    }
    if(month == 66){
      HALT <- TRUE
    }
  }
  return(list(pop=simPop, EPHP=EPHP, prevs=prevs, rounds=MDArounds, time=month_total, low_prevs=low_prevs, high_prevs=high_prevs))
}

# Helper function to observe next 10 years after EPHP
runOut <- function(simPop, month_total) {
  HALT <- FALSE
  rebounce_time <- month_total + 60
  elim_time <- month_total + 120
  rebounce <- 0
  elim <- 0
  month <- 1
  prevalence <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
  prevs <- prevalence
  vth_low <- min(simPop$VtH)
  vth_high <- max(simPop$VtH)
  low_prevs <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
  high_prevs <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
  
  while (HALT == FALSE) {
    simPop$runTimestep()
    mfprev <- length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20))*100
    month <- month + 1
    month_total <- month_total + 1
    prevs[month] <- mfprev
    low_prevs[month] <- length(which((simPop$age[simPop$VtH==vth_low]/12)>=20 & simPop$Mf[simPop$VtH==vth_low]>0))/length(which((simPop$age[simPop$VtH==vth_low]/12)>=20))*100
    high_prevs[month] <- length(which((simPop$age[simPop$VtH==vth_high]/12)>=20 & simPop$Mf[simPop$VtH==vth_high]>0))/length(which((simPop$age[simPop$VtH==vth_high]/12)>=20))*100
    
    # Check for rebounce = prevalence exceeds 1% in any month 5-10 years after EPHP
    if (month_total >= rebounce_time) {
      if (mfprev > 1) {
        rebounce <- 1
      }
    }
    # Check for probability of elimination = 1 if prevalence exceeds 1% exactly 10 years after EPHP, 0 otherwise
    if (month_total == elim_time) {
      if (mfprev <= 1) {
        elim <- 1
      }
      HALT <- TRUE
    }
  }
  return(list(pop=simPop, time=month_total, prevs=prevs[2:length(prevs)], rebounce=rebounce, elim=elim, low_prevs=low_prevs[2:length(low_prevs)], high_prevs=high_prevs[2:length(high_prevs)]))
}

# Call function to run multiple simulations on the same population and collate results
sample.R <- function(N, newPop, compliance, coverage){
  all_prevs <- list()
  all_rounds <- c()
  all_time <- c()
  all_pops <- c()
  all_restarts <- c()
  all_tas1 <- c()
  all_elim <- c()
  all_low <- list()
  all_high <- list()
  
  # For the number of iterations specified, run the simulation and collect results
  for (j in 1:N){
    ans=callrunsMDA(newPop, compliance, coverage)
    
    all_prevs <- c(all_prevs, list(ans$prevs))
    all_rounds <- c(all_rounds, ans$rounds)
    all_time <- c(all_time, ans$time)
    all_pops <- c(all_pops, ans$pop)
    all_restarts <- c(all_restarts, ans$restart)
    all_tas1 <- c(all_tas1, ans$tas1)
    all_elim <- c(all_elim, ans$elim)
    all_low <- c(all_low, list(ans$low))
    all_high <- c(all_high, list(ans$high))
  }
  
  return(list(pops=all_pops, prevs=all_prevs, rounds=all_rounds, time=all_time, restart=all_restarts, tas1=all_tas1, elim=all_elim, low=all_low, high=all_high))
}