## This script is for simulating the future characteristics of coral populations following differing future annual thermal stress regimes.
## This script requires constructed IPMs capturing the differing dynamics of coral populations exposed to both thermal stress, and non-bleaching conditions - see see the IPM cotruction script (https://github.com/CantJ/Subtropical-Coral-Demographics) for further details on the construction of IPM models 
## This script also requires 
## Finally this script also requires a document detailing predicted recurrent disturbance frequencies. See see the Predicting DHWs from CMIP5 data script (https://github.com/CantJ/Subtropical-Coral-Demographics) for further details on predicting future thermal distrubance frequencies.
# Author: James Cant

# set working directory
setwd("Directory/pathway")

# load required packages
library(gmodels)
library(reshape)
library(Rmisc)

# load data detailing the projected frequencies of thermal stress events - see the 'Predicting DHWs from CMIP5 data' script for how to source this information from open source climate model data.
# this data needs to consist of columns detailing the annual probabilities of thermal stress events under differing future climate scenarios, with each row respresenting subsequent years. 
p.bleaching <- read.csv("Bleaching regimes.csv", stringsAsFactors = FALSE) 

# load formatted data file - see Data formating script
ipm <- read.csv("File name3.csv")

######################################################
# Step 1: Extract initial population vectors and create necessary storage.
######################################################

# Colate relevant bleaching and non-bleaching models from each coral population.
bleaching.models <- list(acro_2016, turb_2016, poci_2016, encrust_2016) 
normal.models <- list(acro_2017, turb_2017, poci_2017, encrust_2017)
# (throughout this script components representing the different coral populations will be stored in a standard order - 1 = Acropora, 2 = Turbinaria, 3 = Pocillopora, 4 = Encrusting)

# store the model variables for each coral population (1 = Acropora, 2 = Turbinaria, 3 = Pocillopora, 4 = Encrusting)
L.store = list(0.9*min(acro[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(turb[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(poci[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(robust[,c("Size.t","Size.t1")], na.rm = T)) 
Umax.store = list(1.1*max(acro[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(turb[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(poci[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(robust[,c("Size.t","Size.t1")], na.rm = T))
m.store = list(110, 160, 110, 263)
Ut.store = list(discrete_threshold_acro, discrete_threshold_turbinaria, discrete_threshold_poci, discrete_threshold_r)

# Extract colony sizes from most recent surveys for each coral population. These will be retained to operate as initial population size distributions for the following simulations.
acro.2018 = subset(ipm, t.year == "2017" & !is.na(Size.t1) & Genus == "Acropora") 
turb.2018 = subset(ipm, t.year == "2017" & !is.na(Size.t1) & Genus == "Turbinaria") 
poci.2018 = subset(ipm, t.year == "2017" & !is.na(Size.t1) & Genus == "Pocillopora") 
encrust.2018 = subset(ipm, t.year == "2017" & !is.na(Size.t1) & morphology == "Encrusting_Massive")
initial.pop = list(acro.2018, turb.2018, poci.2018, encrust.2018) #store intial size distributions in one location

# create storage locations for simulation outputs, and define simulation length.
tmax <- (2100-2018) + 1 #length of projection
Nt.list <- list(NULL) # total population sizes
nt.list <- list(NULL) # size distributions
# these storage lists need repeating so that multipe loop outputs can be stored into one single array.
Nt.total <- list(NULL)
nt.total <- list(NULL)
# create storage for estimated mean sizes of colonies within the discrete size class for each coral population
mean_size.group <- list(NULL)

# define subsetting functions that will be used to place corals within the continuous or discrete size classes based on colony size.
Index.continous <- function(pop, x){ 
  index <- which(pop$Size.t1 < x)
  return(index)
}  #this will seperate out corals within the continous class

Index.discrete <- function(pop, x){ 
  index <- which(pop$Size.t1 >= x)
  return(index)
} # this will seperate out corals within the discrete class.

######################################################
# Step 2: Conduct simulations
######################################################

for (i in 1:1000) { #this will run the simulation 1000 times to allow for calculating variance in the projections
  # Within each run through a large loop will run through a stochastic projection for each of the RCP pathways, for each coral group storing the population vectors at each ioteration of the projection cycle.
  # it is run this way so that the timing of thermal stress events are identical for each coral group.
  
  for (z in 1:4) { #selects the RCP pathway: 1 = 2.6, 2 = 4.5, 3 = 6.0, 4 = 8.5
    
    # create sub storage areas
    Nt.RCP <- data.frame("Acropora" = numeric(tmax),"Turbinaria" = numeric(tmax), "Pocillopora" = numeric(tmax),"Encrusting" = numeric(tmax)) #one column for each coral group per RCP pathway
    nt.RCP <- list(NULL) #this will contain 4 matrices - 1 for each coral group
    
    #select the annual bleaching frequency sequence for the selected RCP pathway
    bleaching.frequency <- p.bleaching[z+2] 
    
    # define a vector of random numbers - this will be used to create a random probability occurence function in the loop below.
    prob <- runif(100) #this way each coral groups will be subjected to corresponding bleaching periods (if they occur)
    
    for (g in 1:4) { # now within each RCP pathway run a loop for each coral group
      #1 = Acropora, 2 = Turbinaria, 3 = Pocillopora, 4 = Encrusting 
      
      # extract the relevant model variables
      U <- Umax.store[[g]]
      L <- L.store[[g]]
      m <- m.store[[g]]
      Ut <- Ut.store[[g]]
      
      # -----------------------------------------------------------
      # Set up the initial population vector ready for projection. 
      # Currently the stored intial population vector is a simple vector. However due to the two stage nature of the IPMs these vectors need formatting to consist of a cotinuous scale of colony sizes upto the populations size threshold.
      # After this the vector will consist of one value indicating the number of colonies existing within the discrete size class.
      
      #define the bin widths
      h <- (Ut-L)/m 
      
      # storage matrix
      nt.group <- matrix(NA, m+1, tmax) #the +1 adds a section for the discrete catagory. This will be the size distribution of each during each iteration
      
      #determine the bin boundaries specific to the selected coral group.
      bin.bound <- L + c(0:m) * (Ut - L)/m
      threshold <- max(bin.bound) #this is the largest colony size that fits into the continous class (will be equal to the size threshold)
      
      # format the initial population vector.
      start.pop <- initial.pop[[g]] # select the corresponding initial population vector (colony sizes in 2018) 
      # identify the colonies within the starting population that fit into the continous class and determine their size distribution.
      n <- hist(start.pop[Index.continous(start.pop, threshold),]$Size.t1,   
                breaks = bin.bound, plot = FALSE)$counts # split the population down by the precalculated bin bundaries - so that our population vector matches the dimensions of the projection matrix.
      limit <- length(n) #this is the size of the continous class (for subsetting later in the loop)
      # IPM simulations use probability density functions rather than count vectors, so the continuous population vector needs to be converted into an integration kernel format.
      n <- length(start.pop[Index.continous(start.pop, threshold),]$Size.t1)* n/(sum(h * n))
      
      # now to determine the number of individuals in the discrete class - i.e the number of colonies larger than the size threshold.
      n_large <- dim(start.pop[Index.discrete(start.pop, threshold),])[1]
      # and determine the mean size of colonies within this size class
      large.class.size.group <- mean(start.pop[Index.discrete(start.pop, threshold),]$Size.t1, na.rm = T)
      # and store for later
      mean_size.group[[g]] <- large.class.size.group
      
      # Now combine the continous and discrete classes to form the initial population for projection
      n0 <- c(n, n_large) #this way all but the last entry in this vector is a probability density kernel.
      
      # a little check to ensure formatting the population worked
      check <- sum(h * n0[1:limit]) + n0[limit+1] #convert the probability density functions back to counts and add the discrete class size.
      if (check == length(start.pop$Size.t1)) {print("Yeah!")} else {print("Error")}
      # the population is now ready for projection.
      
      #--------------------------------------------------------------
      
      #store the inital population states
      nt.group[,1] <- n0
      Nt.RCP[1,g] <- (h * sum(nt.group[1:limit,1])) + nt.group[limit+1,1] #this is the appropriate way to determine the total size of a population from its combined density function and discrete class
      
      # -------------------------------------------------------------
      # now project the populations one year at a time.
      
      for(t in 2:tmax){
        if (prob[t] < bleaching.frequency[t,]){ # this randomises the chance of bleaching but ensures that as the probability of thrmal stress events increases it becomes more likely.
          Kt = bleaching.models[[g]]$K #select the matrix
          nt.group[,t] <- Kt%*%nt.group[,t-1] #project the population integral
          Nt.RCP[t, g] <- (h * sum(nt.group[1:limit,t])) + nt.group[limit+1,t] # store the population size (combining the PDF and discrete stage count)
        } else { # if bleaching doesn't occur.
          Kt = normal.models[[g]]$K 
          nt.group[,t] <- Kt%*%nt.group[,t-1]
          Nt.RCP[t, g] <- (h * sum(nt.group[1:limit,t])) + nt.group[limit+1,t]
        } 
      } #this completes the RCP projection for each of the individual coral groups and 
      # has already stored Nt estimates in the appropriate location but still needs to store the nt matrix 
      nt.RCP[[g]] <- nt.group
    } # this completes the coral group loop
    # rename list items to make it easier to find at later date
    names(nt.RCP) <- c("Acropora", "Turbinaria","Pocillopora","Encrusting")
    # Now store the outputs for each RCP pathway
    Nt.list[[z]] <- Nt.RCP
    nt.list[[z]] <- nt.RCP
  } # this completes the RCP pathway loop
  # rename list items to make it easier to find at later date
  names(nt.list) <- c("RCP2.6", "RCP4.5","RCP6.0","RCP8.5")
  names(Nt.list) <- c("RCP2.6", "RCP4.5","RCP6.0","RCP8.5")
  # store the multiple runs for each 
  Nt.total[[i]] <- Nt.list
  nt.total[[i]] <- nt.list
} #this completes the variance loop

######################################################
# Step 3: Format storage lists ready for calculating projected cover and group-specific lambdaS values under the different RCP pathways
######################################################
# the way the simulation loop has run has meant that the way population size and vector estimates have been store in a way that is difficult to use in any further analyses.
# The following section of code is designed to re-format the simulation outputs to allow for analyses of lambdaS or coverage estimates.

# a list of names for naming extracted data in a loop
RCP.names <- c("RCP2.6","RCP4.5","RCP6.0","RCP8.5")

# extract Nt values from storage list
for (d in 1:4) { #this will run once for each RCP pathway
  #create an extraction matrix
  mat <- matrix(NA,83,4000)
  
  #then for each variance run, extract the coral population sizes for the associated RCP.
  for (x in 1:1000) { # run through each variance run
    mat[,x] <- Nt.total[[x]][[RCP.names[d]]][["Acropora"]]
    mat[,x + 1000] <- Nt.total[[x]][[RCP.names[d]]][["Turbinaria"]]
    mat[,x + 2000] <- Nt.total[[x]][[RCP.names[d]]][["Pocillopora"]]
    mat[,x + 3000] <- Nt.total[[x]][[RCP.names[d]]][["Encrusting"]] #this places everything in the same matrix
  }
  
  assign(paste0("Nt.", RCP.names[d]), mat) #rename it so that it is saved based on the associated pathway
  # rm(mat) # it may be necessary to remove the original extraction to save processor space
}

# extract the nt vectors 
for (d in 1:4) { #this will run once for each RCP pathway
  #create an initial matrix of the first entry
  mat1 <- nt.total[[1]][[RCP.names[d]]][["Acropora"]]
  mat2 <- nt.total[[1]][[RCP.names[d]]][["Turbinaria"]]
  mat3 <- nt.total[[1]][[RCP.names[d]]][["Pocillopora"]]
  mat4 <- nt.total[[1]][[RCP.names[d]]][["Encrusting"]]
  #then for each variance run extract the coral populations vecotrs for the associated RCP.
  for (x in 2:1000) { # run through each variance run
    mat1 <- cbind(mat1, nt.total[[x]][[RCP.names[d]]][["Acropora"]])
    mat2 <- cbind(mat2, nt.total[[x]][[RCP.names[d]]][["Turbinaria"]])
    mat3 <- cbind(mat3, nt.total[[x]][[RCP.names[d]]][["Pocillopora"]])
    mat4 <- cbind(mat4, nt.total[[x]][[RCP.names[d]]][["Encrusting"]]) #this places everything in the same matrix
  }
  
  assign(paste0("nt.Acropora.", RCP.names[d]), mat1) #rename it so that it is saved based on the associated pathway
  # rm(mat1) # again it may be necessary to remove the original extraction to save processor space
  assign(paste0("nt.Turbinaria.", RCP.names[d]), mat2) 
  # rm(mat2)
  assign(paste0("nt.Pocillopora.", RCP.names[d]), mat3) 
  # rm(mat3)
  assign(paste0("nt.Encrusting.", RCP.names[d]), mat4) 
  # rm(mat4)
}

######################################################
# Step 4: Calculate lambdaS and its variance
######################################################
# Now to workout the lambdaS values for each population. Lambd is the asympotic growth rate of a population. LambdaS is the stochastic format, and easy essentially the mean annual growth rate expressed by a population throughout a simulation.

# Use the population size to determine the stepwise lambda values for each coral group for every run of the simulation (all 4000 of them), for each RCP pathway
df <- list(Nt.RCP2.6, Nt.RCP4.5, Nt.RCP6.0, Nt.RCP8.5)

# create a matrix to store results
lambda.s <- matrix(NA, 4, 4000)

for(v in 1:4){ #do it once for each RCP
  for(i in 1:4000) { # then once for every simulation in that RCP pathway work down each column calculating lambda for each iteration.
    # create a blank vector to use in the calculation
    lambda.step <- numeric(tmax-1) 
    for(t in 1:tmax-1) {
      lambda.step[t] <- df[[v]][t+1,i]/df[[v]][t,i] #this will go through calculating the ratio between Nt and Nt+1 to give the lambda value associated with that single time step.
    }
    # calculate and store the lambda s value 
    lambda.s[v,i] <- exp(mean(log(lambda.step), na.rm = T))
    # clean up 
    rm(lambda.step)
  } 
}


# Calculate the mean LambdaS value across the simulations, for each coral group under each RCP pathway and the associated variance.
# RCP2.6
mean(lambda.s[1,1:1000], na.rm = T); sd(lambda.s[1,1:1000], na.rm = T)  # Acropora 
mean(lambda.s[1,1001:2000], na.rm = T); sd(lambda.s[1,1001:2000], na.rm = T) # Turbinaria
mean(lambda.s[1,2001:3000], na.rm = T); sd(lambda.s[1,2001:3000], na.rm = T) # Pocillopora
mean(lambda.s[1,3001:4000], na.rm = T); sd(lambda.s[1,3001:4000], na.rm = T) # Encrusting
CI(na.omit(lambda.s[1,1:1000])) # Acropora
CI(na.omit(lambda.s[1,1001:2000])) # Turbinaria
CI(na.omit(lambda.s[1,2001:3000])) # Pocillopora
CI(na.omit(lambda.s[1,3001:4000])) # Encrusting

# RCP4.5
mean(lambda.s[2,1:1000], na.rm = T); sd(lambda.s[2,1:1000], na.rm = T)  # Acropora 
mean(lambda.s[2,1001:2000], na.rm = T); sd(lambda.s[2,1001:2000], na.rm = T) # Turbinaria
mean(lambda.s[2,2001:3000], na.rm = T); sd(lambda.s[2,2001:3000], na.rm = T) # Pocillopora
mean(lambda.s[2,3001:4000], na.rm = T); sd(lambda.s[2,3001:4000], na.rm = T) # Encrusting
CI(na.omit(lambda.s[2,1:1000])) # Acropora
CI(na.omit(lambda.s[2,1001:2000])) # Turbinaria
CI(na.omit(lambda.s[2,2001:3000])) # Pocillopora
CI(na.omit(lambda.s[2,3001:4000])) # Encrusting

# RCP6.0
mean(lambda.s[3,1:1000], na.rm = T); sd(lambda.s[3,1:1000], na.rm = T)  # Acropora 
mean(lambda.s[3,1001:2000], na.rm = T); sd(lambda.s[3,1001:2000], na.rm = T) # Turbinaria
mean(lambda.s[3,2001:3000], na.rm = T); sd(lambda.s[3,2001:3000], na.rm = T) # Pocillopora
mean(lambda.s[3,3001:4000], na.rm = T); sd(lambda.s[3,3001:4000], na.rm = T) # Encrusting
CI(na.omit(lambda.s[3,1:1000])) # Acropora
CI(na.omit(lambda.s[3,1001:2000])) # Turbinaria
CI(na.omit(lambda.s[3,2001:3000])) # Pocillopora
CI(na.omit(lambda.s[3,3001:4000])) # Encrusting

# RCP8.5
mean(lambda.s[4,1:1000], na.rm = T); sd(lambda.s[4,1:1000], na.rm = T)  # Acropora 
mean(lambda.s[4,1001:2000], na.rm = T); sd(lambda.s[4,1001:2000], na.rm = T) # Turbinaria
mean(lambda.s[4,2001:3000], na.rm = T); sd(lambda.s[4,2001:3000], na.rm = T) # Pocillopora
mean(lambda.s[4,3001:4000], na.rm = T); sd(lambda.s[4,3001:4000], na.rm = T) # Encrusting
CI(na.omit(lambda.s[4,1:1000])) # Acropora
CI(na.omit(lambda.s[4,1001:2000])) # Turbinaria
CI(na.omit(lambda.s[4,2001:3000])) # Pocillopora
CI(na.omit(lambda.s[4,3001:4000])) # Encrusting

# Run an anova to test the significance of any differences between population growth rates under the different RCP pathways for the different coral groups.
# First bring together the different RCP lambdaS estimates for each coral group.
lambda.s.acropora <- list(RCP_2.6 = lambda.s[1,1:1000], RCP_4.5 = lambda.s[2,1:1000], RCP_6.0 = lambda.s[3,1:1000], RCP_8.5 = lambda.s[4,1:1000]) #extract the first model to begin the dataframe
lambda.s.acropora <- plyr::ldply(lambda.s.acropora)
lambda.s.acropora <- reshape::melt(lambda.s.acropora); lambda.s.acropora$.id <- as.factor(lambda.s.acropora$.id) #this combines the different RCP rows so that I can compare the effect of RCP.
lambda.s.turbinaria <- list(RCP_2.6 = lambda.s[1,1001:2000], RCP_4.5 = lambda.s[2,1001:2000], RCP_6.0 = lambda.s[3,1001:2000], RCP_8.5 = lambda.s[4,1001:2000]) #extract the first model to begin the dataframe
lambda.s.turbinaria <- plyr::ldply(lambda.s.turbinaria)
lambda.s.turbinaria <- reshape::melt(lambda.s.turbinaria); lambda.s.turbinaria$.id <- as.factor(lambda.s.turbinaria$.id)
lambda.s.pocillopora <- list(RCP_2.6 = lambda.s[1,2001:3000], RCP_4.5 = lambda.s[2,2001:3000], RCP_6.0 = lambda.s[3,2001:3000], RCP_8.5 = lambda.s[4,2001:3000]) #extract the first model to begin the dataframe
lambda.s.pocillopora <- plyr::ldply(lambda.s.pocillopora)
lambda.s.pocillopora <- reshape::melt(lambda.s.pocillopora); lambda.s.pocillopora$.id <- as.factor(lambda.s.pocillopora$.id)
lambda.s.encrusting <- list(RCP_2.6 = lambda.s[1,3001:4000], RCP_4.5 = lambda.s[2,3001:4000], RCP_6.0 = lambda.s[3,3001:4000], RCP_8.5 = lambda.s[4,3001:4000]) #extract the first model to begin the dataframe
lambda.s.encrusting <- plyr::ldply(lambda.s.encrusting)
lambda.s.encrusting <- reshape::melt(lambda.s.encrusting); lambda.s.encrusting$.id <- as.factor(lambda.s.encrusting$.id)

# Run the anovas
lambdas.s.mod.a <- aov(value~.id, data = lambda.s.acropora)
summary(lambdas.s.mod.a); TukeyHSD(lambdas.s.mod.a)
lambdas.s.mod.t <- aov(value~.id, data = lambda.s.turbinaria)
summary(lambdas.s.mod.t); TukeyHSD(lambdas.s.mod.t)
lambdas.s.mod.p <- aov(value~.id, data = lambda.s.pocillopora)
summary(lambdas.s.mod.p); TukeyHSD(lambdas.s.mod.p)
lambdas.s.mod.r <- aov(value~.id, data = lambda.s.encrusting)
summary(lambdas.s.mod.r); TukeyHSD(lambdas.s.mod.r)

######################################################
# Step 5: Calculate coverage projections & variance
######################################################

# create a list to store everything for each coral group
acropora <- list(nt.Acropora.RCP2.6, nt.Acropora.RCP4.5, nt.Acropora.RCP6.0, nt.Acropora.RCP8.5, NA)
turbinaria <- list(nt.Turbinaria.RCP2.6, nt.Turbinaria.RCP4.5, nt.Turbinaria.RCP6.0, nt.Turbinaria.RCP8.5, NA)
pocillopora <- list(nt.Pocillopora.RCP2.6, nt.Pocillopora.RCP4.5, nt.Pocillopora.RCP6.0, nt.Pocillopora.RCP8.5, NA)
encrusting <- list(nt.Encrusting.RCP2.6, nt.Encrusting.RCP4.5, nt.Encrusting.RCP6.0, nt.Encrusting.RCP8.5, NA)
# store them all together
all.groups <- list(acropora,turbinaria,pocillopora,encrusting)

for (coral in 1:4) { # for each coral group
  # Create a dataframe to store results.
  coverage.results <- data.frame(RCP26_mean.coverage = numeric(), RCP45_mean.coverage = numeric(), RCP60_mean.coverage = numeric(), RCP85_mean.coverage = numeric(),
                                 RCP26_CI_upper = numeric(), RCP45_CI_upper = numeric(), RCP60_CI_upper = numeric(), RCP85_CI_upper = numeric(),
                                 RCP26_CI_lower = numeric(), RCP45_CI_lower = numeric(), RCP60_CI_lower = numeric(), RCP85_CI_lower = numeric())
  
  for (RCP in 1:4) { #for each RCP pathway
    
    # select the matrix for analysis
    matrix.use <- all.groups[[coral]][[RCP]]
    
    # select relevant model variables
    Ut <- Ut.store[[coral]]
    L <- L.store[[coral]]
    m <- m.store[[coral]]
    h <- (Ut-L)/m
    
    # the mesh sizes determine the average size of each bin in the continous matrix but not the discrete. This will be the mean size of the discrete class estimated earlier.
    # this is needed for calculaing coverage later in the loop
    discrete.size <- mean_size.group[[coral]]
    # determine dimensions being worked with
    dimensions <- dim(matrix.use)
    limit2 <- dim(matrix.use)[1]-1 
    
    # matrices for storing the mean and upper and lower limits of the different projections
    mat.mean <- matrix(NA, dimensions[1], tmax)
    mat.upper <- matrix(NA, dimensions[1], tmax)
    mat.lower <- matrix(NA, dimensions[1], tmax)
    
    # work through each element of the selected matrix and determine the mean number of individuals in each size class during each iteration
    for (r in 1:dimensions[1]) {
      for (c in 1:tmax) {
        mat.mean[r,c] <- mean(matrix.use[r,seq(from = c, by = tmax, length.out = dimensions[2]/tmax)], na.rm = T)
        mat.upper[r,c] <- as.numeric(CI(na.omit(matrix.use[r,seq(from = c, by = tmax, length.out = dimensions[2]/tmax)]))[1])
        mat.lower[r,c] <- as.numeric(CI(na.omit(matrix.use[r,seq(from = c, by = tmax, length.out = dimensions[2]/tmax)]))[3])
      }} # completes mean/CI loop for each RCP pathway
    
    # Now work out the coverage of the different coral groups over time (mean/upper/lower)
        for (t in 1:tmax) {         
      coverage.results[t,RCP] <- (sum(exp(normal.models[[coral]]$meshpts) * (h * mat.mean[1:limit2,t]))) + (mat.mean[limit2+1,t] * exp(discrete.size)) #this determines the average size of each meshpt and multiplies it by the total number in each bin.
      coverage.results[t,RCP + 4] <- (sum(exp(normal.models[[coral]]$meshpts) * (h * mat.upper[1:limit2,t]))) + (mat.upper[limit2+1,t] * exp(discrete.size)) 
      coverage.results[t,RCP + 8] <- (sum(exp(normal.models[[coral]]$meshpts) * (h * mat.lower[1:limit2,t]))) + (mat.lower[limit2+1,t] * exp(discrete.size)) 
    }
  } # completes RCP pathway loop
  
  # Now store the results with the corresponding coral group
  all.groups[[coral]][[5]] <- coverage.results
}

############################################################################## End of Code ####################################################################################