###############################
### Functions Spatial model ###
###############################
library(R6)

generatePops <- function(nPops=10,sizePop=100){
  x <- runif(nPops,0,sizePop)
  y <- runif(nPops,0,sizePop)
  name <- LETTERS[1:nPops] #Only nPops < 27
  pop <- round(rnorm(nPops,sizePop,sd=10))
  
  return(data.frame(name,pop,x,y))
}

# NOT IN USE
findBorder <- function(data,nPops){
  #define border as midway between two central populations
  #(i.e. same number of populations on each side)
  midpoint_x <- (data$x[order(data$x)][c(ceiling((nPops+1)/2))]+data$x[order(data$x)][c(floor((nPops+1)/2))])/2
  return(midpoint_x)
}

# NOT IN USE
splitBorder <- function(data,nPops,midpoint){
  data$nation <- ifelse(data$x<midpoint,1,2) #classify Nation 1 or 2
  data <- data[order(data$x),]
  row.names(data)<-c()
  data$name <- LETTERS[1:nPops]
  return(data)
}

Parameters <- R6Class(
  "Parameters",
  portable = F,
  public = list(
    ## Host ##
    tau = 0.00167, #death rate hosts
    mImpRate = .00006, #mean importation rate
    maxImpRate = 0.0001, #max importation rate
    commRate = 0.5, #commuting rate per timestep (month)
    mink = 0.01, #min k (overdispersion gamma bite risk) 0.01
    maxk = 0.1, #max k (overdispersion gamma bite risk) 0.1
    gl = 0.02, # lymphedema
    gh = 0.71, # hydrocele
    
    ## Vector ##
    species = 1, #0 culex, 1 anopheles
    initialL3 = 3, #initial larval density
    lbda = 5, #bite rate mosq/monthx [5-15] range
    infecMosq = 0.37, #Proportion of mosquitos that get infected
    kappas1 = 4.395, #L3 uptake and development
    r1 = 0.055, #L3 uptake and development
    deathMosq = 5, #death rate of mosquitos
    minvth = 10, #min V/H ratio = 0  
    maxvth = 50, #max V/H ratio = 70
    highVtH = 0,
    lowVtH = 0,
    
    ## Worm ##
    nu = 0, #poly-monogamy parameter
    alpha = 1, #mf birth rate per FW
    psi1 = .414, #Proportion L3 leaving mosquito per bite
    psi2 = .32, #Proportion of L3 leaving mosquito that enter host
    s2 = .00275, #Proportion of L3 entering host that develop into adults
    mu = .0104, #death rate worms
    gamma = .1, #mf death rate
    
    ## Intervention ##
    sysCompN = .99, #systematic non-compliance bednets
    rhoBU = 0, #correlation non-compliance MDA & bite risk
    rhoCN = 0, #correlation non-compliance MDA & bednets
    sysComp = 0.2, #systematic non-compliance - the correlation between rounds 
    fecreds = c(6,9,9),#6, #months of reduction in fecundity due to MDA (DA) - developmental period in the host is 6-12 months
    sysEx = 0, # set this using the model running function, this is the proportion of the pop that are excluded
    #IDA: 100% mf kill and 100% remove adult worms (to represent sterlisation of adult worms)
    #IA: 99% mf kill, 35% adult kill, 9 month sterilisation
    mfKillMDAs =c(.95,.99,1), #DA, IA, IDA
    wKillMDAs = c(.55,.35,1), #DA, IA, IDA
    drugRef = c('DA','IA','IDA') #diethylcarbamazine&Abdnzl, Ivrmctn&Abdnzl, or all three?
  )
)

Population <- R6Class(
  "Population",
  portable = F,
  public = list (
    p = Parameters$new(), #Load parameters
    nHost = NULL, #Total number host in pop
    nHosts = NULL, #Vector with number of hosts in each pop/town
    nPops = NULL, #Number of populations/towns
    nationIndex = NULL, #Index of which population/town is in which nation
    treatCov = c(0,0), #Treatment coverage
    vectorSpec = 1, #0 for culex, 1 for anopheles
    larvae = 50, #average mean free stage larvae
    L3=array(), # Number of L3
    k = NULL, #host overdispersion parameter
    VtH = array(), #Vector to host ratio
    town = array(),
    gl=0.02,
    gh=0.71,
    nation = array(),
    WF = array(), #Female number worms
    WM = array(), #Male number worms
    numberWormsDiedF = array(),
    numberNewWormsF = array(),
    numberWormsDiedM = array(),
    numberNewWormsM = array(),
    Wormyears = array(),
    timestep=array(),
    newinfections =array(),
    Mf = array(), #Number microfilarae
    age = array(), #Age host (in months)
    gender=array(),
    susceptibility=array(),
    susceptibilityh=array(),
    weights = matrix(), #weights according to population commuter model chosen
    borderWeights = matrix(), #border adjusted weights
    Tij = matrix(),  #Average number of individuals moving from i to j in one timestep
    treatprob = array(), # the probability someone will be treated given the coverage and degree of correlation between treatment rounds. 
    sysExProb = array(), # the probablity of a person being someone who is systematicaly excluded from surveillance
    treat = F, #Treatment received
    ttreat = array(),
    bednet = F, #Bednet used
    bnCov = 0, #Bednet coverage
    biteRisk = NULL, #bi - bite risk (heterogeneity)
    biteRate = NULL, #h(a) - age dependent biting rate
    #border = 0.5, #border strength 0=never cross border, 1=always cross border, 0.5=no effect
    dMatrix = matrix(), #distances between towns
    moveMatrix = matrix(), #who's moving, where to and where from
    
    # Initialize function - specify individual ages/bite risks and sample population parameters (e.g. k, VtH)#
    initialize = function(data,mydMatrix,nationalVtH=F){
      self$nHosts <- data$pop
      self$nationIndex <- data$nation
      self$dMatrix <- mydMatrix
      nHost <<- sum(nHosts)
      nPops <<- length(nHosts)
      larvae <<- rep(larvae,nPops)
      L3 <<-     rep(0,nHost)
      Tij <<- matrix(0,nPops,nPops)
      town <<- unlist(sapply(seq_along(nHosts),function(x){rep(LETTERS[x],nHosts[x])}), use.names=F)
      WF <<- rep(0,nHost) #Number of female adult worms
      WM <<- rep(0,nHost) #Number of male adult worms
      numberWormsDiedF<<- rep(0,nHost)
      numberNewWormsF<<- rep(0,nHost)
      numberWormsDiedM<<- rep(0,nHost)
      numberNewWormsM<<- rep(0,nHost)
      Wormyears<<- rep(0,nHost)
      timestep<<- rep(0,nHost)
      
      newinfections<<- rep(0,nHost)
      Mf <<- rep(0,nHost) #Number of microfilaria
      age <<- runif(nHost,0,100)*12 #in months
      gender <<- sample(0:1,size=nHost,replace=TRUE)
      susceptibility <<- rgamma(nHost, shape=gl,scale=1/gl)
      susceptibilityh <<- ifelse(gender!=0,rgamma(nHost, shape=gh,scale=1/gh),rgamma(nHost, shape=0,scale=0))
      treat <<- rep(F,nHost) #Boolean treated or not
      ttreat <<- rep(12*15,nHost) #Time since treatment (default to very high - 15 years, so that there is no influence)
      #treatprob <<- rbeta(nHost-(nHost*p$sysEx),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp) #Prob each individual received MDA
      sysExProb <<- rbinom(nHost, 1, p$sysEx)
      
      whichDrug <<- NA
      wKillMDA <<- NA
      mfKillMDA <<- NA
      drugID <<- NA
      
      #bednet <<- rep(F,nHost) #Boolean, using bednets or not
      bednet <<- rbinom(nHost,1,bnCov)
      # if(nationalVtH) {
      #   VtH <<- runif(max(self$nationIndex),p$minvth,p$maxvth)
      #   print(VtH)
      # } # sample one VtH for each nation from a uniform
      VtHprob <<- rbinom(nHost, 1, 0.4) # assigned specified proportion of the population to high-risk group
      VtH[which(VtHprob == 1)] <<- p$highVtH # assign high-risk VtH
      VtH[which(VtHprob == 0)] <<- p$lowVtH # assign low-risk VtH
      #if(!nationalVtH) { VtH <<- rep(runif(1,p$minvth,p$maxvth),length(nationIndex)) }
      #print(VtH)
      #print(p$minvth)
      #print(p$maxvth)
      k <<- runif(1,min = p$mink, max = p$maxk) # sample one k for the population from a uniform
      biteRisk <<- rgamma(nHost,shape = k, rate = k) #individual bite risk shape is k, mean=1 so rate=shape
      biteRate <<- pmin(age/(9*12),1) # h(a) age in months - linear increase max out at 9 years - age scale factor for bite risk
      
      # If heterogeneous (two risk groups), assign treatment probability based on specified coverage in each group
      if (p$lowVtH != p$highVtH) {
        treatprob[which(sysExProb==0 & age>=5*12 & VtH == p$lowVtH)] <<- rbeta(length(which(sysExProb==0 & age>=5*12 & VtH == p$lowVtH)),shape1=treatCov[1]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[1])*(1-p$sysComp)/p$sysComp)
        treatprob[which(sysExProb==0 & age>=5*12 & VtH == p$highVtH)] <<- rbeta(length(which(sysExProb==0 & age>=5*12 & VtH == p$highVtH)),shape1=treatCov[2]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[2])*(1-p$sysComp)/p$sysComp) # for those who are NOT part of that group this is their treatment probability 
      } else {
        treatprob[which(sysExProb==0 & age>=5*12)] <<- rbeta(length(which(sysExProb==0 & age>=5*12)),shape1=treatCov[1]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[1])*(1-p$sysComp)/p$sysComp)
      }
      treatprob[which(sysExProb==1 | age<5*12)] <<- 0 # for those who are part of that group, their treatment probability is 0 - they will never receive treatment
    },
    
    ## Do worm dynamics
      wormDynamics = function(){
      names <- LETTERS[1:nPops] #list of population names
      for (i in seq_along(names)){
        inds <- which(town==names[i]) #indexes of who is in population i
        townHost <- length(inds)
        bnRed <- ifelse(bednet[inds],.03,1) #bednet reduction, if individual has bednets reduce by 97%, if false no reduction
        numberWormsDiedF[inds] <- rpois(townHost,p$mu*WF[inds])+numberWormsDiedF[inds]
        numberNewWormsF[inds] <- rpois(townHost,.5 * p$lbda * biteRisk[inds] * VtH[inds] * p$psi1 * p$psi2 * p$s2 * biteRate[inds] * larvae[i] * bnRed)+numberNewWormsF[inds]
        numberWormsDiedM[inds]<- rpois(townHost,p$mu*WM[inds])+numberWormsDiedM[inds]
        numberNewWormsM[inds]<- rpois(townHost,.5 * p$lbda * biteRisk[inds] * VtH[inds] * p$psi1 * p$psi2 * p$s2 * biteRate[inds] * larvae[i] * bnRed)+numberNewWormsM[inds]
        WF[inds] <<- pmax(WF[inds] - numberWormsDiedF[inds] +
                            numberNewWormsF[inds],0)
        WM[inds] <<- pmax(WM[inds] - numberWormsDiedM[inds]+
                            numberNewWormsM[inds],0)

        newinfections[inds]=ifelse(numberNewWormsF[inds]>0 | numberNewWormsM[inds]>0,1,0)
      }
      # calculating the worm years
      timestep=which(newinfections==1) # check the index for each new infection
      wormyrs=c()
      timestep1=c()
      for (i in (1:length(timestep))){
        
        timestep1[i]=timestep[i+1]-timestep[i] # timestep= new worms gained at i+1 - new worm gained at i
        wormyrs[i]=(WF[i]+WM[i])*timestep1[i]
      }
      
      return(list(first=numberWormsDiedF,second=numberWormsDiedM,third=numberNewWormsF,fourth=numberNewWormsM,fifth=newinfections,sixth=wormyrs))
    },
      
      
    ## Do mf dynamics
    mfDynamics = function(){
      I = ifelse(WM>0,1,0) #Need males for reproduction
      Tr = ifelse(treat==T,1,0) #If treatment, no new Mf
      Mf <<- pmax(Mf - p$gamma * Mf + p$alpha * WF * I * (1-Tr),0)
    },
    
    ## Do larval dynamics
    larvaeDynamics = function(){
      names <- LETTERS[1:nPops]
      for(i in seq_along(names)){
        inds <- which(town==names[i]) #indexes of who is in population i
        if(p$species==0){ #Culex
          L3<<-  p$kappas1 * (1 - exp(-p$r1 * Mf[inds] / p$kappas1))
        } else { #Anopheles
          L3<<- p$kappas1 * (1 - exp(-p$r1 * Mf[inds] / p$kappas1))^2
        }
        L3 <<- L3* biteRisk[inds]
        bnRed <- ifelse(bnCov==0,1,1-bnCov*(1-.03)) #bednet reduction
        meanlarvae <- sum(L3)/sum(biteRisk[inds]) * bnRed #Sum across town population rather than whole population
        larvae[i] <<- pmax(p$lbda * p$infecMosq * meanlarvae / (p$deathMosq + p$lbda * p$psi1),0) 
      }
    },
    
    ## Run one timestep
    runTimestep = function(burnin=F){
      if (burnin==F){
        runCommute()
      }
      wormDynamics()
      mfDynamics()
      #rngImportation(x=which(rnorm(nHost)<(1-exp(-.0006))))
      larvaeDynamics()
      aging()
      if(burnin==F){
        returnHome()
        ttreat <<- ttreat + 1 #update how long it's been since treatment
      }
      treat[which(ttreat >= rbinom(length(ttreat),2*p$fecreds[drugID],0.5))] <<- F #stochastic change
      #treat[which(ttreat == p$fecred)] <<- F #step change
    },
    
    # NOT IN USE
    assignexclusion = function(exclusion){
      p$sysEx <<- exclusion
      sysExProb <<- rbinom(nHost, 1, p$sysEx)
    },
    
    # NOT IN USE
    AssignHighVectorVars = function(){
      p$lbda <<- 611
      p$minvth <<- 60 # deafault:60
      p$maxvth <<- 150 # default:150
    },
    
    # NOT IN USE
    AssignLowVectorVars = function(){
      p$lbda <<- 3
      p$minvth <<- 20# deafault:20
      p$maxvth <<- 40# deafault:40
    },
    
    runMDA = function(pop=1,towns=NA,coverage,drug, compliance){
      if(all(coverage != treatCov) | compliance!= p$sysComp ){ #check syst non-compliance probabilities been set for this coverage level already, if not then set them
        #| exclusion!= p$sysEx
        treatCov <<- coverage
        p$sysComp <<- compliance 
        #p$sysEx <<- exclusion
        #sysExProb <<- rbinom(nHost, 1, p$sysEx)
        sysExProb <<- rep(0, nHost)
        lowVtH <- min(VtH)
        highVtH <- max(VtH)
        # If heterogeneous (two risk groups), assigned treatment probability based on specified coverage in each group
        if (lowVtH != highVtH) {
          treatprob[which(sysExProb==0 & age>=5*12 & VtH == lowVtH)] <<- rbeta(length(which(sysExProb==0 & age>=5*12 & VtH == lowVtH)),shape1=treatCov[1]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[1])*(1-p$sysComp)/p$sysComp)
          treatprob[which(sysExProb==0 & age>=5*12 & VtH == highVtH)] <<- rbeta(length(which(sysExProb==0 & age>=5*12 & VtH == highVtH)),shape1=treatCov[2]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[2])*(1-p$sysComp)/p$sysComp)
        } else {
          treatprob[which(sysExProb==0 & age>=5*12)] <<- rbeta(length(which(sysExProb==0 & age>=5*12)),shape1=treatCov[1]*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov[1])*(1-p$sysComp)/p$sysComp) # for those who are NOT part of that group this is their treatment probability 
        }
          treatprob[which(sysExProb==1 | age<5*12)] <<- 0 # for those who are part of that group, their treatment probability is 0 - they will never receive treatment 
        #treatprob <<- rbeta(nHost-(nHost*p$sysEx),shape1=treatCov*(1-p$sysComp)/p$sysComp,shape2=(1-treatCov)*(1-p$sysComp)/p$sysComp)
      }
      
      drugID <<- which(p$drugRef == drug)
      wKillMDA <<- p$wKillMDAs[drugID]
      mfKillMDA <<- p$mfKillMDAs[drugID]
      if(is.na(towns)){
        names = LETTERS[1:nPops]
        towns = names[which(nationIndex == pop)]
      }
      
      treatPops = which(town %in% towns) # this selects the whole population associated with the right town
      treated = treatPops[which(rbinom(length(treatPops),1,treatprob)==1)] # this then randomly assigns threatment to each individual in the population with a probablity equal to louise's term 
      # insert a check to make sure that the treated group is the right proportion of the population 
      WM[treated] <<- rbinom(length(treated),WM[treated],1-wKillMDA)
      WF[treated] <<- rbinom(length(treated),WF[treated],1-wKillMDA)
      Mf[treated] <<- Mf[treated]*rbinom(Mf[treated],1,1-mfKillMDA)
      treat[treated] <<- T #this will cause no mf production
      ttreat[treated] <<- rep(0,length(treated))
      whichDrug <- drugID
    },
    
    # NOT IN USE
    runsurvey=function(){
      val=FALSE
    },

    # NOT IN USE
    #Random importation
    rngImportation = function(x,xi=10 * 9.2 * .414 * .32 * .00275,mu=.0104){
       WF[x] <<- rep(round(.5 * xi * biteRate * 10/mu),length=length(x))
       WM[x] <<- rep(round(.5 * xi * biteRate * 10/mu),length=length(x))
     },
    
    ##Aging & births
    aging = function(){
      age <<- age+1
      resetHosts(which(runif(nHost) < 1-exp(-p$tau)| age>1200)) #simulates deaths/births
      biteRate <<- pmin(age/(9*12),1) #Update biteRate - age in months - linear increase max out at 9 years
    },
    
    ## Do burnin period - lb length burnin (in months)
    burnin = function(lb=1200){
      i=0
      repeat{
        runTimestep(burnin=T)
        i = i+1
        if (i==lb){ break}
      }
    },
    
    ## reset hosts
    resetHosts = function(ID){
      sysExProb[ID] <<- rbinom(1, 1, p$sysEx)
      age[ID] <<- 0 #birth age set to 0
      WF[ID] <<- 0 #female worms set to 0
      WM[ID] <<- 0 #male worms set to 0

      Mf[ID] <<- 0 #Mf set to 0
      L3[ID] <<- 0
      treat[ID] <<- F #no treatment
      ttreat[ID] <<- 15*12
      biteRisk[ID] <<- rgamma(length(ID),shape = k, rate = k) #Get new biterisk - shape is k, mean=1 so rate=shape
      biteRate[ID] <<- pmin(age[ID]/(9*12),1)
      bednet[ID] <<- F #no bednet
    },
    
    # NOT IN USE
    #Gravity model of trade, mass*mass/distance (all variables to the power of alphas close to 1)
    weightsTrade = function(){
      weights_Grav1 <- matrix(0,nPops,nPops)
      alpha <- rnorm(3,1,0.01) #powers close to 1
      i = 1
      repeat{
        weights_Grav1[i,] <- (nHosts[i]^alpha[1])*(nHosts)^alpha[2]/(dMatrix[i,]^alpha[3])
        if(i==nPops){break}
        i = i+1
      }
      weights_Grav1[which(weights_Grav1==Inf)] <- NA
      weights <<- weights_Grav1
    },
    
    # NOT IN USE
    #Einstein's gravity model, mass*mass/distance^2
    weightsGravity = function(){
      weights_Grav2 <- matrix(0,nPops,nPops)
      i = 1
      repeat{
        weights_Grav2[i,] <- nHosts[i]*nHosts/(dMatrix[i,]^2)
        if(i==nPops){break}
        i = i+1
      }
      weights_Grav2[which(weights_Grav2==Inf)] <- NA
      weights <<- weights_Grav2
    },
    
    # NOT IN USE
    #Radiation model
    #https://en.wikipedia.org/wiki/Radiation_law_for_human_mobility
    weightsRadiation = function(){
      weights_Rad <- matrix(0,nPops,nPops)
      s <- matrix(0,nPops,nPops)
      #check <- matrix(0,nPops,nPops)
      m <- nHosts
      i = 1
      repeat{
        j = 1
        repeat{
          #check[i,j] <- check[i,j]+1
          dist <- dMatrix[i,j]
          s[i,j] <- sum(nHosts[which((dMatrix[i,]<dist) & (dMatrix[i,]!=0))])
          weights_Rad[i,j] <- ifelse(i==j,NA,(m[i]*m[j])/((m[i]+s[i,j])*(m[i]+m[j]+s[i,j])))
          if(j==nPops){break}
          j = j+1
        }
        if(i==nPops){break}
        i = i+1
      }
      weights <<- weights_Rad
    },
    
    # NOT IN USE
    # Modify weights to incorporate border effect
    # border = some fraction that defines relative probability of crossing border during migration
    # 0.5 = no effect
    # 1 = always cross border
    # 0 = never cross border
    commuterCalc = function(border=0){
      borderWeights <- matrix(NA,nPops,nPops)
      i = 1
      repeat{
        j = 1
        repeat{
          if (i!=j) {
            borderWeights[i,j] <- ifelse(nationIndex[i]==nationIndex[j],weights[i,j]*(1-border),weights[i,j]*border)
          }
          if(j==nPops){break}
          j = j+1
        }
        borderWeights[i,] <- borderWeights[i,]/sum(borderWeights[i,],na.rm=T)
        if(i==nPops){break}
        i = i+1
      }
      borderWeights <<- borderWeights
      
      commuters <- rbinom(nPops,nHosts,p$commRate)
      Tij <<- matrix(0,nPops,nPops)
      Tij <<- commuters*borderWeights #Average number of individuals moving from i to j in one timestep
      if(nPops==1) {Tij <<- matrix(0,1,1)}
    },
    
    # NOT IN USE
    runCommute = function(){ #runs one time step commute, simultaneously across all towns
      names <- LETTERS[1:nPops]
      commuters <- matrix(rpois(matrix(1,nPops,nPops),Tij),nPops,nPops)
      
      totalMoves = data.frame(whoMove=c(),whereTo=c(),whereFrom=c(),stringsAsFactors = F)
      i=1
      repeat {
        numMove <- sum(commuters[i,],na.rm=T)
        numMove <- min(numMove,length(which(town==names[i])))
        if(numMove > 0){
          whoMove <- sample(which(town==names[i]),numMove)
          probs <- replace(borderWeights[i,],which(is.na(borderWeights[i,])),0)
          whereTo <- sample(names,length(whoMove),replace=T,probs)
          whereFrom <- rep(names[i],numMove)
          moves <- data.frame(whoMove,whereTo,whereFrom,stringsAsFactors = F)
          totalMoves <- rbind(totalMoves,moves,stringsAsFactors=F)
        }
        if(i==nPops){break}
        i = i+1
      }
      
      moveMatrix <<- totalMoves
      town[moveMatrix$whoMove] <<- moveMatrix$whereTo
      nHosts <<- as.vector(table(town)) #recalculate population sizes
      #Needed if permanently moving, but if just commuting then home, weights etc. don't change:
      #weightsTrade() #recalculate commuting weights
      #commuterCalc() #recalculate commuter numbers/border weights
    },
    
    # NOT IN USE
    returnHome = function(){
      town[moveMatrix$whoMove] <<- moveMatrix$whereFrom
      nHosts <<- as.vector(table(town)) #recalculate population sizes
      #weightsTrade() #recalculate commuting weights
      #commuterCalc() #recalculate commuter numbers/border weights
    }
  )
)