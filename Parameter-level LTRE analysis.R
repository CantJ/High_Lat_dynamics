## Script for conducting a parameter-level LTRE contribution analysis to determine the driving mechanisms underpinning the differences between the bleaching and non-bleaching dyanmics for each morphological coral groups.
## This script requires constructed bleaching and non-bleaching IPM models for each coral group, along with their associated estimates of Lambda and model parameters.
## This script will also estimate the variance in parameter-level contributions and so requires the list of model parameters produced during Jacknife analyses for each IPM model.
# Author: James Cant
# Date: August 2019

#load required packages
library(popdemo)
library(Hmisc)

# Ensure required models are constructed using relevant scripts available at https://github.com/CantJ/Subtropical-Coral-Demographics

############################################
# Step 1: Store the required models and variables
############################################

# collate the required IPM kernels.
# First bring together the bleaching models for each morphological coral group following the order Acropora, Turbinaria, Pocillopora and Encursting (This ordering will be maintained throughout this script)
bleaching.models <- list(acro_2016, turb_2016, poci_2016, encrust_2016) 
# repeat for non-bleaching models
nonbleaching.models <- list(acro_2017, turb_2017, poci_2017, encrust_2017)

# store required model variables relevant to each morphological coral group (L, Ut, Umax, m) - see IPM construction script for the necessity for each of these parameters.
L.store = list(0.9*min(acro[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(turb[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(poci[,c("Size.t","Size.t1")], na.rm = T), 1.1*min(robust[,c("Size.t","Size.t1")], na.rm = T)) 
Umax.store = list(1.1*max(acro[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(turb[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(poci[,c("Size.t","Size.t1")], na.rm = T), 1.1*max(robust[,c("Size.t","Size.t1")], na.rm = T))
m.store = list(110, 160, 110, 263)
Ut.store = list(discrete_threshold_acro, discrete_threshold_turbinaria, discrete_threshold_poci, discrete_threshold_encrust)
# as before the order remains Acropora[1], Turbinaria [2], Pocillopora [3], Encrusting [4].

# store lambda estimates for the original models - both bleaching and non-bleaching.
bleaching.lambda <- c(lam.a16, lam.t16, lam.p16, lam.e16)
nonbleaching.lambda <- c(lam.a17, lam.t17, lam.p17, lam.e17)

# store the model parameters used to construct the orginal IPMs for each morphological coral group
store.par.16 <- list(m.par.a16, m.par.t16, m.par.p16, m.par.e16)
store.par.17 <- list(m.par.a17, m.par.t17, m.par.p17, m.par.e17) # these will be used to estimate the initial parameter-level LTRE contributions

# store parameter lists created during the jacknife analysis for each morphological coral group. This will create a list within a list
store.par.var.16 <- list(a.par.16, t.par.16, p.par.16, e.par.16)
store.par.var.17 <- list(a.par.17, t.par.17, p.par.17, e.par.17) # these will allow for the estimation of the variance in parameter-level LTRE contributions.

# In the manuscript associated with this R code the structures of the different IPMs models differed between each morphological groups.
# As such the the LTRE analyses needed conducting seperately for each coral group before re-defining the structure of IPMs and therefore the arrangement of model parameters before moving to the next coral group.

############################################
# Step 2: Parameter-level LTRE analysis
############################################

###### Acropora (Postion [[1]] in storage lists)
# First calculate the sensitivity of lambda to small changes in mean parameter values positioned halfway between the actual parameters estimated when constructing the bleaching and non-bleaching IPM models.
# create a blank vector for storing parameter sensitivites
sParMeanMat <- numeric(length(store.par.16[[1]])) 
# calculate the mean parameter values positioned halfway between original parameter estimates for both bleaching and non-bleaching models.
m.par.halfway <- (store.par.16[[1]]+store.par.17[[1]])/2
# extract relevant model variables (L, Ut, Umax, m)
L <- L.store[[1]]
U <- Umax.store[[1]]
Ut <- Ut.store[[1]]
m <- m.store[[1]]   

# sequentially apply a perturbation to each mean parameter estimate and estimate the new value of lambda following the perturbation.
# this can then be used to estimate the sensitivity of lambda to each individual parameter
for(j in 1:length(store.par.16[[1]])){ # repeat the loop once for each model parameter 
  try(m.par <- m.par.halfway) # store the mean parameters so that the loop doesn't affect the lead vector
  try(m.par[j] <- m.par.halfway[j] - 0.001) # apply a negative perturbation to the selected parameter value
  try(IPM.down <- mk_K.n(m = m, m.par = m.par, L = L, Umax = U, Ut = Ut)) # generate a new IPM model using the peturbed parameter
  try(lambda.down <- Re(eigen(IPM.down$K)$values[1])) # calculate a new estimate for lambda
  try(m.par[j] <- m.par.halfway[j] + 0.001) # now apply a positive perturbation   
  try(IPM.up <- mk_K.n(m = m, m.par = m.par, L = L, Umax = U, Ut = Ut)) # generate a new IPM using the perturbed parameter
  try(lambda.up <- Re(eigen(IPM.up$K)$values[1])) #reestimate a new lambda
  try(sj <- (lambda.up-lambda.down)/(2*0.001)) # the sensitivity of lambda is equal to the change in lambda over the change in the paramter value
  try(sParMeanMat[j] <- sj) # store the sensitity estimate 
}    # the output here is a vector of the relative sensitivites of lambda to independant changes in each parameter

## To complete the etimates of parameter-level contributions multiply the sensitivity estimates by the differences in each parameter between the orginal bleaching and non-bleaching models.
LTREcontributionsPars <- (store.par.16[[1]] - store.par.17[[1]])*sParMeanMat

# Now run a similar loop to estimate the variance in these contributions estimates.
# This will use a loop that will work through all the stored parameter variants produced during the bootstrap analysis of the acropora IPMs.
# extract required bootstrapped parameter list
variance.use.a16 <- store.par.var.16[[1]] 
variance.use.a17 <- store.par.var.17[[1]]
# create a storage matrix
acro.ltre.variance <- matrix(NA, 1000, length(variance.use.a16[[1]]))
# run a loop that will go through each of the 1000 parameter set variants and run an LTRE analysis as before.
for (x in 1:1000){
  jsParMeanMat <- numeric(length(variance.use.a16[[x]])) 
  jm.par.halfway <- (variance.use.a16[[x]]+variance.use.a17[[x]])/2 
  jL <- L.store[[1]]
  jU <- Umax.store[[1]]
  jUt <- Ut.store[[1]]
  jm <- m.store[[1]]   
  
  # sequentially apply a perturbation to each mean parameter estimate and estimate the new value of lambda following the perturbation.
  # this can then be used to estimate the sensitivity of lambda to each individual parameter
  for(z in 1:length(jsParMeanMat)){ 
    try(jm.par <- jm.par.halfway) 
    try(jm.par[z] <- jm.par.halfway[z] - 0.001)
    try(jIPM.down <- mk_K.n(m = jm, m.par = jm.par, L = jL, Umax = jU, Ut = jUt)) 
    try(jlambda.down <- Re(eigen(jIPM.down$K)$values[1])) 
    try(jm.par[z] <- jm.par.halfway[z] + 0.001)   
    try(jIPM.up <- mk_K.n(m = jm, m.par = jm.par, L = jL, Umax = jU, Ut = jUt)) 
    try(jlambda.up <- Re(eigen(jIPM.up$K)$values[1])) 
    try(jsj <- (jlambda.up-jlambda.down)/(2*0.001)) 
    try(jsParMeanMat[z] <- jsj) 
  }  
  
  # multiply the sensitivity estimates by the change in the main parameter estimates between bleaching and non-bleaching models, before storing the vector of parameter contributions.
  acro.ltre.variance[x,] <- (variance.use.a16[[x]] - variance.use.a17[[x]]) * jsParMeanMat
}

# use the LTRE contribution variants to estimate the standard deviation for each paramters contribution estimates 
# create a blank storage vector witha position for each parameter.
LTREerror <- numeric(dim(acro.ltre.variance)[2])
# work through each  parameter (corresponding with each column of the variance matrix) and determine the standard deviation of contribution estimates
for (i in 1:dim(acro.ltre.variance)[2]){
  stan.dev <- sd(acro.ltre.variance[,i], na.rm = T)
  LTREerror[i] <- stan.dev
} # this will be used for displaying error bars on LTRE contribution plots

############################################
# Step 3: Plot parameter-level LTRE contributions
############################################

# set plotting window
par(mfrow = c(1,1))

# create the base plot
acro.plot <- barplot(LTREcontributionsPars, 
                     ylim = c(-0.6,0.6), # these values can be altered depending on the level of resolution equired. 
                     yaxs = "i", xaxs = "i", # forces the y and x axis to coverge at zero
                     xaxt = "n", 
                     ylab = NULL,  main = NULL, las = 2, col = "lightgoldenrod4", cex.axis = 1.5)

# in the manuscript associated with this code the plot colour was used to define a scale representing the magnitude of the difference in lambda between bleaching and non-bleaching models for each coral group.
# to achieve this a coloured background is added to the base plot
rect(par("usr")[1], par("usr")[3], par("usr")[2], 
     par("usr")[4],col = "lightgoldenrod2") 
# re-draw the base plot over the top of the coloured rectangle to accound for some of the details masked by adding the colour.
acro.plot <- barplot(LTREcontributionsPars, ylim = c(-0.6,0.6), yaxs = "i", xaxs = "i", xaxt = "n", ylab = NULL,  main = NULL, las = 2, col = "lightgoldenrod4", cex.axis = 1.5, add = TRUE)
# and final add error bars to the plot
arrows(acro.plot, LTREcontributionsPars - LTREerror, acro.plot, LTREcontributionsPars + LTREerror, code = 3, angle = 90, length = 0.05)


####### The above code can now be repeated for each coral group by simply replacing the place holder value [[...]] when extracting stored parameter lists, or model and variables - Acropora[1], Turbinaria [2], Pocillopora [3], Encrusting [4].

####################################################################### End of code #####################################################################################