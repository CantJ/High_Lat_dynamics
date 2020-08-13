# This script will use the formated demographic data frame to create Intergral Projection models for the four morphological coral groups
# Here we document the subsetting of the datafile the create separate datasets for each morphpological group before detailing the construction of an IPM framework.
# This script requires the formated demographic data frame and a dataframe detailing the Genus and species ID of each new colony fragment found during surveys, survey year, the size of each fragment and the original size of the parent colonyprior to fragmentation  
# Author: James Cant
# Date: July 2019

# Clear work space memory
rm(list=ls(all=TRUE))

# load in required packages.
library(popdemo)
library(Hmisc)
library(gam)
library(mgcv)
library(arm)
library(fields)
library(lme4)
library(ggeffects)
library(reshape2)
library(sn)
library(Rmisc)
library(reconPlots)

##################################################
# STEP 1: Load and format data
##################################################

# set working directory
setwd("insert directory pathway")

# load formatted data file - see Data formating script
ipm <- read.csv("File name3.csv")

# load data file contining fragment colony sizes
frag.size.data <- read.csv("File name4.csv")
# and sort data to match main dataframe.
frag.size.data$Size.t <- frag.size.data$Size.t/100 #convert Size to cm2
frag.size.data$frag.size <- frag.size.data$frag.size/100 #convert fragment size to cm2 also
# and convert parent colony sizes measuresured in April to reflect the size adjustment needed - this will require the mean monthly growth increment estimates obtained during data formating.
frag.size.data$Size.t[which(frag.size.data$t.year == "2016")] <- (a - (a - frag.size.data$Size.t[which(frag.size.data$t.year == "2016")])) + (b * 4)
# finally log transform all size data
frag.size.data$Size.t <- log(frag.size.data$Size.t) 
frag.size.data$frag.size <- log(frag.size.data$frag.size)

##################################################
# Step 2: Subset datafile to represent morphological groupings used in associated manuscript
##################################################

# Three of the four morphological coral groups align with the three coral genera Acropora, Turbinaria and Pocillopora.
# Therefore genus ID can be used to subset those three groups.
turb <- subset(ipm, Genus == "Turbinaria")
acro <- subset(ipm, Genus == "Acropora")
poci <- subset(ipm, Genus == "Pocillopora")

# The final morphological group is Encrusting/massive. This group can be subset out using the morphology variable inserted during data formatting.
encrust <- subset(ipm, morphology == "Encrusting_Massive")
  # Now remove the Porites genus as too many colonies died during the bleaching so this genus is displaying demographic charactersitics unrepresentatve of the overall encrusting group.
  encurst <- subset(encrust, Genus != "Porites")

# Check the species included within each of the groups
summary(turb$Species) 
summary(poci$Species) 
summary(acro$Species) 
summary(encrust$Genus)

# Subsetting the data in this way means that the only corals being ignored are Stylophora pistillata and Porites heronensis - which are excluded due to conflicting distributions (S. pistillata) and unrepresentative mortaility schedules (P. heronensis)

# subset fragment data to reflect the four morphological coral groups above.
acro.frag.size <- subset(frag.size.data, Genus == "Acropora")
poci.frag.size <- subset(frag.size.data, Genus == "Pocillopora")
turb.frag.size <- subset(frag.size.data, Genus == "Turbinaria")
encrust.frag.size <- subset(frag.size, Genus %in% c("Acanthastrea","Astrea","Micromussa","Paragoniastrea"))

##################################################
# Step 3: Vital Rate analysis 
##################################################

# IMPORTANT!!!! The following lines of code will need to be repeated separately (or built into a loop) for each of the different morphological coral groups.
# This is because the regression formats used to estimate the different vital rates for each group may need altering to accurately capture the specific survival, growth and fragmentation patterns of the different corals
# Here the code documents estimating the vital rates of one coral group - Turbinaria.

### 1. Colony growth ---------------------------------------------------------------------------------------------------------------------------------------------------------
# Whilst quantifying the growth patterns of the different coral groups its is crucial to also determin the size at which the growth trends during bleaching and non-bleaching phases converge.
# In our sample there was only a small number of very large colonies tagged. Therefore this prevented the accurate modelling of the dynamics of very large corals on a continuous scale.
# Subsequently our IPMs included a discrete size stage detailing the dynamics of large colonies. The threshold size at which corals transistion between the continous and discrete model stages was 
# determined as the size at which growth patterns during both bleaching and non-bleaching converged.

# Plot size data for visual reference (Size at t+1  ~  size at t)
plot(Size.t1 ~ Size.t, col = t.year, data = turb); abline(0,1)
# Using multiple regression formats model the growth patterns keeping year as a fixed effect allowing for the evaluation of changes between bleaching and non-bleaching periods.
# Here we test the accuracy of both linear and non-linear regression structures, both excluding and including different survey location variables as random effects.
grow.t1 <- lm(Size.t1 ~ Size.t * t.year, data = turb) # basic linear format
grow.t2 <- lmer(Size.t1 ~ Size.t * t.year + (Size.t|Site), data = turb) # Site as a random effect 
grow.t3 <- lmer(Size.t1 ~ Size.t * t.year + (Size.t|location), data = turb) # Inshore/Offshore locality as a random effect
grow.t4 <- lmer(Size.t1 ~ Size.t * t.year + (Size.t|Bleaching.State.t) + (Size.t|Site), data = turb) # Site and observed bleaching state as a random effect 
grow.t5 <- glm(Size.t1 ~ Size.t * t.year + I(Size.t^2) * t.year, data = turb) # non-linear model format
grow.t6 <- glm(Size.t1 ~ Size.t * t.year + I(Size.t^2), data = turb)

# Check model AIC values
AIC(grow.t1,grow.t2,grow.t3, grow.t4, grow.t5, grow.t6) 

# Now check the visual fits of the models with the lowest AIC values - this is because as a mathematical entity the AIC value doesn't determine how well a regression model reflect biological reality so it is important to check the models display patterns that are indeed possible.
# Use the selected models to predict colony size estimates for plotting a regression line.
grow.turb <- ggpredict(grow.t1, terms = c("Size.t[0:10]", "t.year"))
grow.turb <- split(grow.turb, grow.turb$group) # Split the predictions for the different years (fixed effects in the model)
grow.turb2 <- ggpredict(grow.t5, terms = c("Size.t[0:10]", "t.year"))
grow.turb2 <- split(grow.turb2, grow.turb2$group) 
grow.turb3 <- ggpredict(grow.t6, terms = c("Size.t[0:10]", "t.year"))
grow.turb3 <- split(grow.turb3, grow.turb3$group)

# Plot the raw size data with point colour differentiating between size transitions reported during bleaching and non-bleaching - this will aid with determining the size threshold at which the growth patterns meet.
plot(Size.t1 ~ Size.t, data = turb, col=c('red', 'black')[as.numeric(turb$t.year)],
     pch=20, cex=0.9,
     xlab=expression("Size at t"), 
     ylab=expression("Size at t+1"),
     ylim = c(0,10),
     xlim = c(0,10), xaxs = "i", yaxs = "i") 
abline(0,1,lty = "solid") #a reference line for the null growth threshold.
# add the predicted regression lines for the best fitting models
lines(grow.turb$'2016'$x, grow.turb$'2016'$predicted, col = "Red", lwd = 2)
lines(grow.turb$'2017'$x, grow.turb$'2017'$predicted, col = "Black", lwd = 2) #model 1
lines(grow.turb2$'2016'$x, grow.turb2$'2016'$predicted, col = "Red", lty = "dashed")
lines(grow.turb2$'2017'$x, grow.turb2$'2017'$predicted, col = "Blue", lty = "dashed") #model 5
lines(grow.turb3$'2016'$x, grow.turb3$'2016'$predicted, col = "Red", lty = "dashed")
lines(grow.turb3$'2017'$x, grow.turb3$'2017'$predicted, col = "Black", lty = "dashed") #model 6
legend(6, 2.5, c("2016", "2017"), fill = c("Red", "Blue"), bg = NULL, bty = "l")
# Both AIC values and visual checks confirm non-linear model 6 as the most suitable regression model.

# To incoperate a discrete size class witin our IPM there is a need to determine the point at which the growth patterns during both bleaching and non-bleaching intersect.
# This function involves the use of the curve_intersect function from the reconplots package. Further details and function code can be found on the reconplots Github page (https://github.com/andrewheiss/reconPlots).
# Create intersecting function to determine the threshold size at t at which the predicted growth regression lines intercept each other.
intersectpoint <- function(model, max, min) {
  x <- seq(min, max, 0.0001) # define the x axis
  coeffs1 <- c(model$coefficient[1], model$coefficient[2], model$coefficient[4])  # coefficients of bleaching pattern
  coeffs2 <- c(model$coefficient[1] + model$coefficient[3], model$coefficient[2] + model$coefficient[5], model$coefficient[4])  # coefficients of non-bleaching pattern
  
  # for each model estimate the value of y for each value of x
  y1 <- coeffs1[1] + (coeffs1[2]*x) + (coeffs1[3]*(x^2)) 
  y2 <- coeffs2[1] + (coeffs2[2]*x) + (coeffs2[3]*(x^2))
  
  # bring everything together
  trend1 <- data.frame(x, y1)
  trend2 <- data.frame(x, y2)
  # determine the intersect of the 2 curves.
  intersect <- curve_intersect(trend1, trend2, empirical = TRUE)$x
  return(intersect)
}  

# The function can now be used to define the intersect size from the chosen regression model.
# this size value will then be used throuhgout the rest of the IPM to define the threshold point at which corals transition between continuous and discrete size classes.
# determine size threshold - note the use of maximum and minimum recorded colony sizes.
discrete_threshold_turbinaria <- intersectpoint(grow.t6, max = max(turb$Size.t, na.rm = T), min = min(turb$Size.t, na.rm = T))

# the data set being used can now be split based on this threshold size - however the full data set will still be used to estimate the overall vital rates for the model so that our artificial threshold point has no bearing on reported dynamics within the continuous size class.
turb_discrete <- turb[which(turb$Size.t >= discrete_threshold_turbinaria),] #corals larger or equal to the threshold size are in this discrete class.
turb_continous <- turb[which(turb$Size.t < discrete_threshold_turbinaria),] #smaller corals remain in the continous class.
# recruits/fragments and any corals with NA as Size.t will remain in the original dataframe

### 2. Progression into the discrete size class -------------------------------------------------------------------------------------------------------------------------------
# Estimate the size-specific probability of progression into the discrete stage out of the continuous stage (This will only involve colonies identified as being in the continuous class at time t)
turb_continous$becomelarge <- 0 # all corals initially assigned a progression probability of 0%
turb_continous$becomelarge[which(turb_continous$Size.t1 >= discrete_threshold_turbinaria)] = 1 # colonies larger than the estimated threshold size, during time t+1, assigned a progression probability of 100%
turb_continous$becomelarge[which(turb_continous$surv == 0)] = NA # colonies that didn't survive can not be included in this progression analysis.
# Using binomial regression the sized based probability of transferring into the discrete stage can be estimated
plot(becomelarge ~ Size.t, data = turb_continous); growth.progression.turb <- glm(becomelarge ~ Size.t, family = "binomial", data = turb_continous)
# predict trend line for plotting
progression <- ggpredict(growth.progression.turb, terms = c("Size.t"), full.data = T) 
#plot regression lines
lines(spline(progression$predicted~progression$x, method = "n", n =250), col = "Red", lwd = 2) 

### 3. Retrogression out of the discrete size class --------------------------------------------------------------------------------------------------------------------------
# Estimate the probability of retrogression out of discrete stage into the continuous stage (This will only involve colonies identified as being in the discrete class at time t)
turb_discrete$becomesmall <- 0 # probability intially set as 0
turb_discrete$becomesmall[which(turb_discrete$Size.t1 < discrete_threshold_turbinaria)] = 1 # large colonies that end time t+1 smaller than the estimated threshold size are assigned a retrogression probability of 100%
turb_discrete$becomesmall[which(turb_discrete$surv == 0)] = NA # only colonies that survive are included
# Now these can be used to estimate the probability of transferring out off the discrete stage and the likely resultant size.
no_shrinking <- dim(turb_discrete[which(turb_discrete$becomesmall == 1),])[1] #how many colonies shrunk?
stasis <- dim(turb_discrete)[1] #how many colonies started out in the discrete size class?
retrogression_turb <- no_shrinking/stasis #so the ratio of colonies shirnking out of the discrete class is?
# chek the size distribution of colonies following shrinkage out of the discrete stage
hist(turb_discrete[which(turb_discrete$becomesmall == 1),]$Size.t1)
shrink_size_turb <- glm(Size.t1~1, data = turb_discrete[which(turb_discrete$becomesmall == 1),]) # forcing this regression model to a slope of provides an estimate of the mean colony size following shrinkage (model intercept) and size variance (SD)

### 4. Growth variation -----------------------------------------------------------------------------------------------------------------------------------------------------
# From the growth plots above it does appear that growth variability varies with colony size.
# this can be observed by plotting growth regression residuals
plot(grow.t6$model$Size.t, abs(resid(grow.t6)), xlab='size', ylab='residual')
# determine the relationship between model residuals and colony size.
growth.sd.turb = glm(abs(resid(grow.t6))~grow.t6$model$Size.t, family = Gamma(link = "log")) # non-linear - here a Gamma distribution is applied to prevent the the model from allowing negatuve variances
growth.sd.turb2 = lm(abs(resid(grow.t6))~grow.t6$model$Size.t) # linear
# How do these models compare visually and through AIC values
AIC(growth.sd.turb, growth.sd.turb2)
lines(grow.t6$model$Size.t, exp(predict(growth.sd.turb, list(grow.t6$model$Size.t))))
lines(grow.t6$model$Size.t, (predict(growth.sd.turb2, list(grow.t6$model$Size.t))))
# the non-linear trend fits better here 

### 5. Survival -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Continuous class
# Vital rates for the continuous class will be determine using the full data set (ignoring the threshold size) to ensure the model is not effected by artificial groupings.
# Survival is estimated using a binomial regression format - again multiple structures should be tested for accuracy with year as a fixed effect
surv.t1 <- glm(surv ~ Size.t * t.year, family = "binomial", data = turb)
surv.t2 <- glmer(surv ~ Size.t * t.year + (Size.t|Site), family = "binomial", data = turb) # Site as a random effect
surv.t3 <- glmer(surv ~ Size.t * t.year + (Size.t|location), family = "binomial", data = turb) # Inshore/Offshore locality as a random effect.

# Which model fits the best?
# Check AIC scores
AIC(surv.t1, surv.t2, surv.t3)
# Check the visual accuracy
surv.turb <- ggpredict(surv.t1, terms = c("Size.t[0:10]", "t.year"))
surv.turb <- split(surv.turb, surv.turb$group)
# plot the data
plot(surv ~Size.t, type = "n", data = turb_continous,
     # the x axis needs to match the one used when predicting values.
     xlab=expression("Size at t"), 
     ylab=expression("Survival Probability"),
     ylim = c(0,1),
     xlim = c(0,10),
     xaxs = "i", yaxs = "i") #this forces the plot to meet at (0,0)
# add the predicted lines
lines(spline(surv.turb$'2016'$predicted~surv.turb$'2016'$x, method='n', n=250), col = "Red", lwd = 2) #splines makes the plot fit a smooth curve which is a more natural visualisation of this type of data.
lines(spline(surv.turb$'2017'$predicted~surv.turb$'2017'$x, method='n', n=250), col = "Black", lwd = 2)

# Survival in the discrete size class (this is a single probability and not size dependant)
# Survival during bleaching period
turb.mortality_bleaching <- dim(turb_discrete[which(turb_discrete$surv == 0 & turb_discrete$t.year == "2016"),])[1] #how many corals in the discrete group died
turb.group_size_bleaching <- dim(turb_discrete[which(turb_discrete$t.year == "2016"),])[1] #total size of group
turb.survival_bleaching <- 1-(turb.mortality_bleaching/turb.group_size_bleaching) # probability of dying within the discrete class converted to reflect survival probability.
# Repeat the same estimates for survival during non-bleaching
turb.mortality <- dim(turb_discrete[which(turb_discrete$surv == 0 & turb_discrete$t.year == "2017"),])[1]
turb.group_size <- dim(turb_discrete[which(turb_discrete$t.year == "2017"),])[1]
turb.survival_n.bleaching <- 1-(turb.mortality/turb.group_size)

### 6. Fragmentation --------------------------------------------------------------------------------------------------------------------------------------------------------
# Continuous class
# As with survival fragmentation is estimated using a binomial regression format - again multiple structures should be tested for accuracy with year as a fixed effect
frag.t1 <- glm(frag ~ Size.t * t.year, family = "binomial", data = turb)
frag.t2 <- glmer(frag ~ Size.t * t.year + (Size.t|Site), family = "binomial", data = turb)

# Check model fits using visual checks and AIC scores
AIC(frag.t1,frag.t2)
# Predicted fragmentation probabilities
fragmentation.t <- ggpredict(frag.t1, terms = c("Size.t[0:10]", "t.year"))
fragmentation.t <- split(fragmentation.t, fragmentation.t$group)
# Plot data
plot(frag~Size.t, type = "p", data=turb_continous,
     # the x axis needs to match the one used when predicting values.
     xlab=expression("Size at t"), 
     ylab=expression("Fragmentation Probability"),
     ylim = c(0,1),
     xlim = c(0,10),
     xaxs = "i", yaxs = "i")
# Add predicted lines
lines(spline(fragmentation.t$'2016'$predicted ~ fragmentation.t$'2016'$x, method = "n", n = 250), col = "Red", lwd = 2)
lines(fragmentation.t$'2017'$predicted ~ fragmentation.t$'2017'$x, col = "Black", lwd = 2)

# Fragmentation prbability within the Discrete class
# Fragmentation during bleaching
turb.fragmenting_bleaching <- dim(turb_discrete[which(turb_discrete$frag == 1 & turb_discrete$t.year == "2016"),])[1]
turb.fragmentation_bleaching <- turb.fragmenting_bleaching/turb.group_size_bleaching # the group size used in the survival estimates can be re-used here.
# Fragmentation during non-bleaching
turb.fragmenting <- dim(turb_discrete[which(turb_discrete$frag == 1 & turb_discrete$t.year == "2017"),])[1]
turb.fragmentation_n.bleaching <- turb.fragmenting/turb.group_size

### 7. Number of fragments produced as a function of parent colony size ------------------------------------------------------------------------------------------------------
# Continuous class
# For this we are only interested in the corals that did fragment.
t.frag.corals <- subset(turb, No.fragments > 0)
# fragmentation was not common enough to accuratly model the number of fragments produced with year as a fixed effect
# The number of fragments produced should be modelled as a poisson distribution as the number of fragments produced can not be negative and is a form of discrete count data.
t.no.frag <- glm(No.fragments ~ Size.t, family = "poisson", data = t.frag.corals)

# predicted fragment number produced as a function of size.
t.frag.no <- ggpredict(t.no.frag, terms = c("Size.t[-2:10]"))
# plot raw data
plot(No.fragments ~ Size.t, type = "p", data = t.frag.corals,
     # the x axis needs to match the one used when predicting values.
     xlab=expression("Size at t"), 
     ylab=expression("Number of fragments"),
     ylim = c(0,5),
     xlim = c(-2,9.5),
     xaxs = "i", yaxs = "i")
# add in predicted regression lines
lines(t.frag.no$predicted ~ t.frag.no$x, col = "Blue")

# This sized based function needs extrapolating for the discrete class - as no fragmentation was acctually observed in large colonies. Again this will be done ignoring the effect of year (Bleaching/non-bleaching)
mean_turb_size <- mean(turb_discrete$Size.t, na.rm = T) # what is the mean colony size of corals in the discrete size class 
# and now extrapolate the number of fragments predicted by the continuous model adove for colonies of mean discrete colony size. 
turb_frag_no <- as.numeric(exp(predict(t.no.frag, list(Size.t = mean_turb_size))))

### 8a. Fragment size --------------------------------------------------------------------------------------------------------------------------------------------------------
# Fragment size needs to be modelled as a function of parent colony size. However again due to the lack of observed fragmentation this will be carried out without using year as a fixed effect.
# Plot the size of fragments against initial colony size
plot(frag.size~Size.t, data = turb.frag.size, col = t.year, xlim = c(-2,11), ylim = c(-2,11))
# Again multiple regression formats should be used to test how to best represent the relationship. 
t.frag.size.mod <- nls(frag.size ~ exp(a + b * Size.t), data = turb.frag.size, start = list(a = 1, b = 1))
t.frag.size.mod2 <- lm(frag.size ~ Size.t, data = turb.frag.size)
# Check AIC scores
AIC(t.frag.size.mod, t.frag.size.mod2)

# plot model predictions to also evaluate model fit
lines(turb.frag.size$Size.t, predict(t.frag.size.mod, list(Size.t = turb.frag.size$Size.t)), col = "blue") # non-linear
abline(coefficients(t.frag.size.mod2)[1], coefficients(t.frag.size.mod2)[2]) # linear model
# Here it is important to ensure that the model selected does not result in smaller corals producing fragments larger than themselves which is biologically unfeasible
# In our case this was the case in the non-linear model so we worked with the linear model format.

# Fragment size produced by the fragmentation of discrete class colonies
# this again involves the use of an estimated mean size of colonies within this discrete size class.
# and now extrapolate the size of fragments predicted by the continuous model adove. 
turb_frag_size <- as.numeric(predict(t.frag.size.mod2, list(Size.t = mean_turb_size)))

### 8b. Fragment size variability -------------------------------------------------------------------------------------------------------------------------------------------
# as with the main growth patterns variance in the size of fragements varies depending on the size of the initial parent colony.
# This again can be seen by plotting the residuals of the fragment size model.
plot(turb.frag.size$Size.t,abs(resid(t.frag.size.mod2)),xlab='Size',ylab='residual')
# again regression models can be used to accurately reflect this pattern
t.frag.sd.mod = glm(abs(resid(t.frag.size.mod2))~turb.frag.size$Size.t, family = Gamma(link = "log")) # again the gam,a distribution serves to prevent negative variance.
t.frag.sd.mod2 <- lm(abs(resid(t.frag.size.mod2))~turb.frag.size$Size.t) # compare with a linear model format
# check AIC scores
AIC(t.frag.sd.mod,t.frag.sd.mod2) 
# And check which model fits better visually
lines(turb.frag.size$Size.t, exp(predict(t.frag.sd.mod, list(turb.frag.size$Size.t)))) 
abline(coefficients(t.frag.sd.mod2)[1], coefficients(t.frag.sd.mod2)[2]) 
# Now the Standard deviation of fragment size can vary in colony size.

### 9. Fecundity ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Continuous class
# Quantifying the size specific fecundity patterns and how this changes between years work better if the two years are access separately.
fec.t16 <- subset(turb_continous, t.year == 2016); fec.t17 <- subset(turb_continous, t.year == 2017) 
# basic none linear model - the fecundity data was estimated artificially so there is no need to include random effects here when outlining the relationship specific to our populations
# equally including fecundity in our IPMs only serves to numerically complete the loop between exisiting colonies and new recruits - but below we calculate a recruit settlement parameter that ensures our IPMs are not sensitive to colony fecundity estimates.
fec.turb.16 <- lm(log(fec) ~ Size.t, data = fec.t16)
fec.turb.17 <- lm(log(fec) ~ Size.t, data = fec.t17)

# Fecundity of the discrete class colonies
# This requires determining the per capita fecundity of colonies of mean size within this discrete size class.
mean_turb_size.b <- mean(turb_discrete[which(turb_discrete$t.year == "2016"),]$Size.t, na.rm = T) #during bleaching
mean_turb_size.nb <- mean(turb_discrete[which(turb_discrete$t.year == "2017"),]$Size.t, na.rm = T) #during non-bleaching
# and now extrapolate the fecundity predicted by the continuous models adove. 
turb_fec_bleaching <- as.numeric(exp(predict(fec.turb.16, list(Size.t = mean_turb_size.b))))
turb_fec_n.bleaching <- as.numeric(exp(predict(fec.turb.17, list(Size.t = mean_turb_size.nb))))
# this is the mean per captita fecundity of colonies within this size class for bleaching and none bleaching

### 10. Recruit Settlement probability (Recruit settlement factor) ----------------------------------------------------------------------------------------------------------
# The recruit settlement factor serves to convert colony fecundity (which measures the larval volumne produced as a function of colony size) into a size-based proportional contribution of the actual number of recruits observed in reality.
# This ratio is therefore designed to constrain coral recruitment within our IPMs to reflect recruitment trends observed within the field. 
# As with fecundity this easier to determine by keeping survey data from the different years separate.
turb.16 <- subset(turb, t.year == 2016); turb.17 <- subset(turb, t.year == 2017)
# first estimate the total larval volume expected from our tagged population during 2016
turb.fec.16 <- sum(fec.t16$fec, na.rm = T) + #the total fecundity predicted for corals within the continous class during 2016
  (turb_fec_bleaching * dim(turb_discrete[which(turb_discrete$t.year == "2016" & turb_discrete$surv == 1),])[1]) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
# Now determine the actual number of recruits iobserved during 2016 surveys 
turb.rec.16 <- subset(turb.16, Recruit.t1 == "Yes") 
turb.no.recruit.16 <- dim(turb.rec.16)[1] #total number of reported recruits.
# The recruit settlement factor is then the ratio between number of recruits and expected larval volume
turb.prob.rec.16 <- turb.no.recruit.16/turb.fec.16 

# This can then be repeated to estimate the pattern for 2017
turb.fec.17 <- sum(fec.t17$fec, na.rm = T) + #the total fecundity predicted for corals within the continous class during 2017
  (turb_fec_n.bleaching * dim(turb_discrete[which(turb_discrete$t.year == "2017" & turb_discrete$surv == 1),])[1]) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
turb.rec.17 <- subset(turb.17, Recruit.t1 == "Yes") 
turb.no.recruit.17 <- dim(turb.rec.17)[1] #total number of reported recruits.
turb.prob.rec.17 <- turb.no.recruit.17/turb.fec.17 # 2017 recruit settlement factor/ratio

### 11. Recruit size distribution -------------------------------------------------------------------------------------------------------------------------------------------
# Tagged recruits for each year have already been subset whilst estimating the recruit settlement factor above.
# Start by checking the size distribution of the recruits for 2016.
hist(turb.rec.16$Size.t1, breaks = 20,
     xlab = expression("Recruit Size"),
     main = NULL)
# and again during 2017
hist(turb.rec.17$Size.t1, breaks = 20,
     xlab = expression("Recruit Size"),
     main = NULL) 
# Because the number of recruits was poor during 2016 (for Turbinaria) the recruit size distribution will be modelled with years combined. Though it may be possible to capture the interannual variation in the other coral groups.
# subset all recruits
recruits <- subset(turb, Recruit.t1 == "Yes")
#  check distribution
hist(recruits$Size.t1, breaks = 20,
     xlab = expression("Recruit Size"),
     main = NULL) 
# Because this is not a normal distribution we will use a skewed normal distribution to account for negative log(size) values.
# Importantly recruit size is being modelled independantly of adult colony sizes.
t.rec.size <- selm(Size.t1~1, data = recruits) # this provides distribution parameters that can be used to reconstructed a distribution of expected recruit colony sizes in a probability density function.

######################################
# STEP 4: Store the relevant model coefficients.
######################################
# Here it is simply a case of bringing together all the key model coefficients that will be used to construct our IPMs.
# The coefficients relating to the vital rates in both 2016 and 2017 will be stored seperately allowing our models to capture interannual variation in the demographics of this coral population

# 2016 coefficients
m.par.t16 <- c(
  ##### size-structured vital rates for the continous class portion of our IPMs
  # survival
  surv.int        =  coefficients(surv.t1)[1],
  surv.slope      =  coefficients(surv.t1)[2],
  # growth 
  grow.int        =  coefficients(grow.t6)[1],
  grow.slope      =  coefficients(grow.t6)[2],
  grow.slope2     =  coefficients(grow.t6)[4],
  # class progression
  prog.int        =  coefficients(growth.progression.turb)[1],
  prog.slope      =  coefficients(growth.progression.turb)[2],
  # growth variation
  grow.sd.int     =  coefficients(growth.sd.turb)[1], 
  grow.sd.slope   =  coefficients(growth.sd.turb)[2],
  # fragmentation probability
  frag.int        =  coefficients(frag.t1)[1],
  frag.slope      =  coefficients(frag.t1)[2],
  # Number of fragments
  frag.no.int     =  coefficients(t.no.frag)[1],
  frag.no.slope   =  coefficients(t.no.frag)[2],
  # fragment size
  frag.size.int   =  coefficients(t.frag.size.mod2)[1],
  frag.size.slope =  coefficients(t.frag.size.mod2)[2],
  frag.sd.int     =  coefficients(t.frag.sd.mod)[1],
  frag.sd.slope   =  coefficients(t.frag.sd.mod)[2],
  ##### Discrete probability reflecting the vital rates for the discrete class portion of our IPMs
  # retrogression
  retrogression   =  retrogression_turb,
  shrink.size.int =  coefficients(shrink_size_turb)[1],
  shrink.size.sd  =  sigma(shrink_size_turb),
  #survival
  survival        =  turb.survival_bleaching,
  #fragmentation
  fragmentation   =  turb.fragmentation_bleaching,
  frag.number     =  turb_frag_no,
  frag.size       =  turb_frag_size,
  frag.size.var   =  sqrt(pi/2) * exp((coefficients(t.frag.sd.mod)[1] + coefficients(t.frag.sd.mod)[2] * mean_turb_size.b)), #this is the extrapolated fragment size variance for the discrete class.
  #fecundity
  fec             =  turb_fec_bleaching,
  # recruitment into discrete class
  rec             =  0, #recruitment of new colonies into the discrete class is not possible
  #unrelated to classes
  # fecundity as a function of size
  fec.int         =  coefficients(fec.turb.16)[1],
  fec.slope       =  coefficients(fec.turb.16)[2],
  # recruitment probability
  prob.rec        =  turb.prob.rec.16,
  # recruit size
  rcsz.xi         =  t.rec.size@param$dp[1],
  rcsz.omega      =  t.rec.size@param$dp[2],
  rcsz.alpha      =  t.rec.size@param$dp[3])

# 2017 coefficients
m.par.t17 <- c(
  ##### Continous class
  # survival
  surv.int        =  coefficients(surv.t1)[1] + coefficients(surv.t1)[3],
  surv.slope      =  coefficients(surv.t1)[2] + coefficients(surv.t1)[4],
  # growth 
  grow.int        =  coefficients(grow.t6)[1] + coefficients(grow.t6)[3],
  grow.slope      =  coefficients(grow.t6)[2] + coefficients(grow.t6)[5],
  grow.slope2     =  coefficients(grow.t6)[4],
  # class progression
  prog.int        =  coefficients(growth.progression.turb)[1],
  prog.slope      =  coefficients(growth.progression.turb)[2],
  # growth variation
  grow.sd.int     =  coefficients(growth.sd.turb)[1], 
  grow.sd.slope   =  coefficients(growth.sd.turb)[2],
  # fragmentation probability 
  frag.int        =  coefficients(frag.t1)[1] + coefficients(frag.t1)[3],
  frag.slope      =  coefficients(frag.t1)[2] + coefficients(frag.t1)[4],
  # Number of fragments
  frag.no.int     =  coefficients(t.no.frag)[1],
  frag.no.slope   =  coefficients(t.no.frag)[2],
  # fragment size
  frag.size.int   =  coefficients(t.frag.size.mod2)[1],
  frag.size.slope =  coefficients(t.frag.size.mod2)[2],
  frag.sd.int     =  coefficients(t.frag.sd.mod)[1],
  frag.sd.slope   =  coefficients(t.frag.sd.mod)[2],
  ##### Discrete class
  # retrogression
  retrogression   =  retrogression_turb,
  shrink.size.int =  coefficients(shrink_size_turb)[1],
  shrink.size.sd  =  sigma(shrink_size_turb),
  #survival
  survival        =  turb.survival_n.bleaching,
  #fragmentation
  fragmentation   =  turb.fragmentation_n.bleaching,
  frag.number     =  turb_frag_no,
  frag.size       =  turb_frag_size,
  frag.size.var   =  sqrt(pi/2) * exp((coefficients(t.frag.sd.mod)[1] + coefficients(t.frag.sd.mod)[2] * mean_turb_size.nb)),  
  #fecundity
  fec             =  turb_fec_n.bleaching,
  # recruitment into discrete class
  rec             =  0,
  #unrelated to classes
  # fecundity as a function of size
  fec.int         =  coefficients(fec.turb.17)[1],
  fec.slope       =  coefficients(fec.turb.17)[2],
  # recruitment probability
  prob.rec        =  turb.prob.rec.17,
  # recruit size
  rcsz.xi         =  t.rec.size@param$dp[1],
  rcsz.omega      =  t.rec.size@param$dp[2],
  rcsz.alpha      =  t.rec.size@param$dp[3])

# Rename the parameter elements for ease when calling them during the construction of our IPMs.
names(m.par.t16) <- names(m.par.t17) <- c("surv.int","surv.slope", #survival
                                          "grow.int","grow.slope", "grow.slope2", # growth
                                          "prog.int", "prog.slope", #progression to discrete
                                          "grow.sd.int","grow.sd.slope", #growth variance
                                          "frag.int","frag.slope", #fragmentation with size
                                          "frag.no.int", "frag.no.slope", #number of fragments
                                          "frag.s.int", "frag.s.slope", "frag.sd.int", "frag.sd.slope", #fragment size as a function of adult size.
                                          "retrogression", #retrogression from discrete
                                          "size.mean", "size.var", #shrinkage size
                                          "survival", #discrete survival
                                          "frag.prob", #discrete fragmentation
                                          "frag.number", # mean number of frags produced by large coral class
                                          "frag.size", "frag.var", # mean size and variance of frags produced by large coral class
                                          "fecundity", #discrete fecundity
                                          "recruitment", #recruitment into discrete class
                                          "fec.int", "fec.slope", #fecundity with size
                                          "rec.prob", #recruitment liklihood (success)
                                          "rcsz.xi", "rcsz.omega", "rcsz.alpha")  #recruit size

# and store them together in a list
m.par.year <- list(m.par.t16, m.par.t17)

######################################
# STEP 5: Build IPM functions.
######################################
# Our IPMs were built to consist of 3 distinct sub-Kernels. P: reflecting the survival and growth of none fragmenting colonies. H: Refelecting the survival and fragmentation of fragmenting colonies and F: Reflecting the recruitment and settlement of new colonies.
# However because our models consisted of a continuous and discrete size class each kernel require 4 panels to appropriately capture the dynamics of corals to reflect this.
# The first a main panel (1) describes the size-based dynamics of corals remaining within the continuous class. The second panel (2) reflects the size based transistion of colonies from the continuous class into the discrete group.
# The third panel (3) represents the discrete dynamics of colonies starting and remaining within the discrete size class. The final panel (4) describes the probabilities of corals regressing/recruiting from the discrete class into the continuous class.
# These four panels are stitched together internally as part of the final construction of the IPMs.
# The functions below will call upon the store coefficients to describe the underlying patterns being used to build each sub-Kernel and subsequently the overall model.

###### Vital rate functions
# Survival probability function
s_z <- function(Size.t, m.par){
  linear_p <- m.par["surv.int"] + m.par["surv.slope"] * Size.t  
  p <- 1/(1+exp(-linear_p)) #this adjusts for the fact that the relationship involves logits (binomial distribution) 
  return(p) 
}

# Growth function, given you are size z now returns the population density function of size t+1
G_z1z <- function(Size.t1, Size.t, m.par){
  mean <- m.par["grow.int"] + (m.par["grow.slope"] * Size.t)
  sd <- sqrt(pi/2) * exp((m.par["grow.sd.int"] + m.par["grow.sd.slope"] * Size.t)) #allows for the variance to change with size.
  p_den_grow <- dnorm(Size.t1, mean = mean, sd = sd)
  return(p_den_grow)
} 

# Fragmentation probability function
h_z <- function(Size.t, m.par){ 
  linear_p <- m.par["frag.int"] + m.par["frag.slope"] * Size.t
  p <- 1/(1+exp(-linear_p))
  return(p)
}

# Number of fragments produced
ph_z <- function(Size.t, m.par){
  p <- exp(m.par["frag.no.int"] + m.par["frag.no.slope"] * Size.t)
  return(p)
}

# Fragment size function, returns the population density function
f_z1 <- function(Size.t1, Size.t, m.par){
  a <- m.par["frag.s.int"]
  b <- m.par["frag.s.slope"]
  mean <- a + b * Size.t
  sd <- sqrt(pi/2) * exp((m.par["frag.sd.int"] + m.par["frag.sd.slope"] * Size.t)) #allows for the variance in fragment size to change with adult colony size.
  p_den_fsz <- dnorm(Size.t1, mean = mean, sd = sd) #<- calculate the probability of your size given any starting size. The Dnorm function creates a normal size distribution from this data - because remember IPMs create population density distributions.
  return(p_den_fsz)
}

# size based probability of progressing into the discrete stage.
prog_z <- function(Size.t, m.par){
  linear_prog <- m.par["prog.int"] + m.par["prog.slope"] * Size.t #this is the progression likelihood given a certain starting size. 
  prog <- 1/(1+exp(-linear_prog)) 
  return(prog)
}

# probability of shrinking out of discrete class
shrink_discrete <- function(m.par){
  prob <- m.par["retrogression"]
  return(prob)
}

# colony size following retrogression
shrink_size_z1 <- function(Size.t1, m.par){
  mean <- m.par["size.mean"]
  sd <- m.par["size.var"]
  p_den_sz <- dnorm(Size.t1, mean = mean, sd = sd) #<- calculate the probability of your size given any starting size. The Dnorm function creates a normal size distribution from this data - because remember IPMs create population density distributions.
  return(p_den_sz)
}

# probability of surviving in the discrete class
surv_discrete <- function(m.par){
  prob <- m.par["survival"]
  return(prob)
}

# probability of fragmenting in the discrete class
frag_discrete <- function(m.par){
  prob <- m.par["frag.prob"]
  return(prob)
}

# number of fragments produced by fragmentation in the discrete class
ph_discrete <- function(m.par){
  number <- m.par["frag.number"]
  return(number)
}

# the size of fragements produced by fragmentation in the discrete class
F_z1 <- function(Size.t1, m.par){
  mean <- m.par["frag.size"]
  sd <- m.par["frag.var"]
  p_den_Fsz <- dnorm(Size.t1, mean = mean, sd = sd) 
  return(p_den_Fsz)
}

#Recruit Size
# the pdf of size z1 recruits next summer - 
c0_z1 <- function(Size.t1, m.par){
  xi <- m.par["rcsz.xi"]
  omega <- m.par["rcsz.omega"]
  alpha <- m.par["rcsz.alpha"]                          
  p_den_rcsz <- dsn(Size.t1, xi = xi, omega = omega, alpha = alpha) # this returns a skewed normal probability density function relfecting the skewed nature of recruit sizes.
  return(p_den_rcsz)
}

#probability of settlement success
rec.prob <- function(m.par){
  success.prob <- m.par["rec.prob"]
  return(success.prob)
}

# larval density output - assuming exponential fecundity 
# continuous class
pb_z <- function(Size.t, m.par){
  a <- m.par["fec.int"]
  b <- m.par["fec.slope"]
  p <- exp(a + b * Size.t) #produces an exponential relationship between size and predicted fecundity.
  return(p)
}

# discrete class
pb_discrete <- function(m.par){
  pb <- m.par["fecundity"]
  return(pb)
}

#### Sub-kernel functions. The function bring together the relevant vital rates reflecting the different P, H and F dynamics
# Panel 1
# These function represent the dynamics of corals that REMAIN within the continuous size class.
# Under normal circumstances the internal model functions wouldn't change between models - just the coeffcients being used to parameterise them. However because our vital rates of survival and fragementation probabilities for colonies during 2016 still need adjusting 
# to account for the mismathing survey time frames this will be carried out here as part of the functions. Therefore seperate model functions are needed for 2016 and 2017 to reflect this. 

# 2016 (Bleaching) Functions  
# P kernel for annual transitions for corals that DONT fragment
P.16_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - prog_z(Size.t, m.par)) * ((1 - (h_z(Size.t, m.par)^(12/16))) * (s_z(Size.t, m.par)^(12/16)) * G_z1z(Size.t1, Size.t, m.par)) ) # corals that don't fragment survive and grow 
} # to adjust the transition duration the probability of surviving increases and the probability of fragmenting goes down

# fragmentation kernel for corals that do fragment
H.16_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - prog_z(Size.t, m.par)) * ((h_z(Size.t, m.par)^(12/16)) * (s_z(Size.t, m.par)^(12/16)) * ph_z(Size.t, m.par) * f_z1(Size.t1, Size.t, m.par)) ) # corals that do fragment produce a certain number of fragments of a certain size.
}

# Bleaching F kernels ).
F.16_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - prog_z(Size.t, m.par)) * ((s_z(Size.t, m.par)^(12/16)) * pb_z(Size.t, m.par) * # survival and fecundity of 'adult' colonies
                                           rec.prob(m.par) * c0_z1(Size.t1, m.par)) ) # recruitment success and size of recruits.
} 

# 2017 (Non-bleaching) Functions
# P kernel for annual transitions for corals that DONT fragment
P_z1z <- function (Size.t1, Size.t, m.par) {#<- the kernals are created as straightforward functions.
  return( (1 - prog_z(Size.t, m.par)) * ((1 - h_z(Size.t, m.par)) * s_z(Size.t, m.par) * G_z1z(Size.t1, Size.t, m.par)) ) # corals that don't fragment survive and grow 
} 

# fragmentation kernel for corals that do fragment (offshore)
H_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - prog_z(Size.t, m.par)) * (h_z(Size.t, m.par) * s_z(Size.t, m.par) * ph_z(Size.t, m.par) * f_z1(Size.t1, Size.t, m.par)) ) # corals that do fragment produce a certain number of fragments of a certain size.
}

# Define the F kernels (same between assemblages).
F_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - prog_z(Size.t, m.par)) * (s_z(Size.t, m.par) * pb_z(Size.t, m.par) * # survival and fecundity of 'adult' colonies
                                           rec.prob(m.par) * c0_z1(Size.t1, m.par)) ) # recruitment success and size of recruits.
} 

# Panel 2 (The transion of corals from the continuous size class to the discrete class)
# Again any mention of survival and fragmentation requires adjusting in the 2016 model construction
tran_P16 <- function (Size.t, m.par){
  return( prog_z(Size.t, m.par) * (s_z(Size.t, m.par)^(12/16)) )
}

tran_P17 <- function(Size.t, m.par){
  return( prog_z(Size.t, m.par) * s_z(Size.t, m.par) )
} 

# the two functions below apply to both bleaching and none bleaching as neither fragmentation nore recruitment can result in the transition of individuals from the continuous class into the discrete class   
tran_H <- function (Size.t, m.par){
  return( 0 * prog_z(Size.t, m.par) * s_z(Size.t, m.par) )
} # there is no fragmentation that results in the progression of colonies to the discrete stage

tran_F <- function (Size.t, m.par){
  return( 0 * prog_z(Size.t, m.par) * s_z(Size.t, m.par) )
} # there is no recruitment that results in the progression of colonies to the discrete stage

# Panel 3 (Corals within the discrete size class)
large_class_P16 <- function(m.par){
  return( (1 - (frag_discrete(m.par)^(12/16))) * (1 - (shrink_discrete(m.par)^(12/16))) * (surv_discrete(m.par)^(12/16)) )
}

large_class_H16 <- function(m.par){
  return(0 * (frag_discrete(m.par)^(12/16)) * (surv_discrete(m.par)^(12/16)) )
}

large_class_F <- function(m.par){
  recruitment <- m.par["recruitment"]
  return(recruitment)
} #there is no recruitment into this class (the 'recruitment' coefficient is stored as 0)

large_class_P17 <- function(m.par){
  return( (1 - frag_discrete(m.par)) * (1 - shrink_discrete(m.par)) * surv_discrete(m.par) )
}

large_class_H17 <- function(m.par){
  return(0 * frag_discrete(m.par) * surv_discrete(m.par) )
}

# Panel 4 (Transitions of corals from the discrete class into the continuous class)
regress_P16 <- function(Size.t1, m.par){
  return( (shrink_discrete(m.par)^(12/16)) * (surv_discrete(m.par)^(12/16)) * shrink_size_z1(Size.t1, m.par) )
}

regress_H16 <- function(Size.t1, m.par){
  return( (frag_discrete(m.par)^(12/16)) * (surv_discrete(m.par)^(12/16)) * ph_discrete(m.par) * F_z1(Size.t1, m.par) )
}

regress_F16 <- function(Size.t1, m.par){
  return( (surv_discrete(m.par)^(12/16)) * pb_discrete(m.par) * rec.prob(m.par) * c0_z1(Size.t1, m.par) )
}

regress_P17 <- function(Size.t1, m.par){
  return( shrink_discrete(m.par) * surv_discrete(m.par) * shrink_size_z1(Size.t1, m.par) )
}

regress_H17 <- function(Size.t1, m.par){
  return( frag_discrete(m.par) * surv_discrete(m.par) * ph_discrete(m.par) * F_z1(Size.t1, m.par) )
}

regress_F17 <- function(Size.t1, m.par){
  return( surv_discrete(m.par) * pb_discrete(m.par) * rec.prob(m.par) * c0_z1(Size.t1, m.par) )
}

#####################################
### STEP 7: Construct the IPM by stitching together the sub-kernel panels and combining each sub-kernel.
#####################################
# Again this could normally be handled by one function (mK_k), however the different years required different functions to be called.
# Non-bleaching model function
mk_K.n <- function(m, m.par, L, Umax, Ut) {
  
  # mesh points 
  h <- (Ut - L)/m #this is the number of bins for the full population
  meshpts <- L + ((1:m) - 1/2) * h #these are the corresponding mesh points for the full model
  #S <- which(meshpts < Ut) # this is an index for the meshpoints corresponding to sizes smaller than the discrete threshold size
  # P subkernel  (survival & growth)
  P1 <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par)) #Panel 1 (continuous matix element) (REMAIN CONTINUOUS)
  P2 <- as.matrix(h * (as.numeric(sapply(meshpts, tran_P17, m.par = m.par)))) # Panel 2 (continuous vector element - with 1 row) (TRANSITION TO DISCRETE)
  P3 <- as.matrix(as.numeric(large_class_P17(m.par = m.par))) #panel 3 (single probability) (REMAIN DISCRETE)
  P4 <- as.matrix(as.numeric(sapply(meshpts, regress_P17, m.par = m.par))) #Panel 4 (a reverse vector element) (TRANSITION TO CONTINUOUS)
  # H subkernel  (fragmentation)
  H1 <- h * (outer(meshpts, meshpts, H_z1z, m.par = m.par)) 
  H2 <- as.matrix(h * (as.numeric(sapply(meshpts, tran_H, m.par = m.par)))) 
  H3 <- as.matrix(as.numeric(large_class_H17(m.par = m.par)))
  H4 <- as.matrix(as.numeric(sapply(meshpts, regress_H17, m.par = m.par)))  
  # keeping fragmentation seperate in this way prevents the model from mistaking the addition of new fragments as hypermortality at larger size classes.
  # F recruitment  
  F1 <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par)) 
  F2 <- as.matrix(h * as.numeric(sapply(meshpts, tran_F, m.par = m.par)))
  F3 <- as.matrix(as.numeric(large_class_F(m.par = m.par)))
  F4 <- as.matrix(as.numeric(sapply(meshpts, regress_F17, m.par = m.par))) 
  # Now sew the panels together into a complete continuos and discrete IPM
  P <- cbind(rbind(P1,t(P2)),rbind(P4,P3)) #this stitches the relevant panels together for each of the sub-kernels (and ensures the correct alignment of the panels 2 & 4)
  H <- cbind(rbind(H1,t(H2)),rbind(H4,H3))
  F <- cbind(rbind(F1,t(F2)),rbind(F4,F3))
  K <- P + F + H #and builds the full model.
  return(list(K = K, meshpts = meshpts, P = P, F = F, H = H))
} 

mk_K.b <- function(m, m.par, L, Umax, Ut) {
  
  # mesh points 
  h <- (Ut - L)/m #this is the number of bins for the full population
  meshpts <- L + ((1:m) - 1/2) * h #these are the corresponding mesh points for the full model
  #S <- which(meshpts < Ut) # this is an index for the meshpoints corresponding to sizes smaller than the discrete threshold size
  # P subkernel  (survival & growth)
  P1 <- h * (outer(meshpts, meshpts, P.16_z1z, m.par = m.par)) #Panel 1 (continuos matix element) (REMAIN CONTINUOUS)
  P2 <- as.matrix(h * (as.numeric(sapply(meshpts, tran_P16, m.par = m.par)))) # Panel 2 (continuous vector element - with 1 row) (TRANSITION TO DISCRETE)
  P3 <- as.matrix(as.numeric(large_class_P16(m.par = m.par))) #panel 3 (single probability) (REMAIN DISCRETE)
  P4 <- as.matrix(as.numeric(sapply(meshpts, regress_P16, m.par = m.par))) #Panel 4 (a reverse vector element) (TRANSITION TO CONTINUOUS)
  # H subkernel  (fragmentation)
  H1 <- h * (outer(meshpts, meshpts, H.16_z1z, m.par = m.par)) 
  H2 <- as.matrix(h * (as.numeric(sapply(meshpts, tran_H, m.par = m.par)))) 
  H3 <- as.matrix(as.numeric(large_class_H16(m.par = m.par)))
  H4 <- as.matrix(as.numeric(sapply(meshpts, regress_H16, m.par = m.par)))  
  # keeping fragmentation seperate in this way prevents the model from mistaking the addition of new fragments as hypermortality at larger size classes.
  # F recruitment  
  F1 <- h * (outer(meshpts, meshpts, F.16_z1z, m.par = m.par)) 
  F2 <- as.matrix(h * as.numeric(sapply(meshpts, tran_F, m.par = m.par)))
  F3 <- as.matrix(as.numeric(large_class_F(m.par = m.par)))
  F4 <- as.matrix(as.numeric(sapply(meshpts, regress_F16, m.par = m.par))) 
  # Now sew the panels together into a complete continuos and discrete IPM
  P <- cbind(rbind(P1,t(P2)),rbind(P4,P3)) #this puts on the panels together for each of the sub-kernels (and ensures the correct alignment of the panels 2 & 4)
  H <- cbind(rbind(H1,t(H2)),rbind(H4,H3))
  F <- cbind(rbind(F1,t(F2)),rbind(F4,F3))
  K <- P + F + H #and builds the full model.
  return(list(K = K, meshpts = meshpts, P = P, F = F, H = H))
} 

#####################################
### STEP 8: Implement the IPM.
#####################################
# To implement an IPM you need to define the boundary sizes and the number of bins in the kernal matrix.
# First you need to define the minimum and maximum colony sizes to be considered within the model - these are typically 10% smaller and larger than the smallest and largest colonies actually observed.
# Lowest possible colony size for the model to consider
L = 1.1*min(turb[,c("Size.t","Size.t1")], na.rm = T) # in this dataset the smallest colony size is a negative value (on the log scale) and so multiplying by 1.1 produces a lower balue
U = 1.1*max(turb[,c("Size.t","Size.t1")], na.rm = T) #maximum size ever 
# the model also needs to know the continuous to discrete size class threshold size.
Ut <- discrete_threshold_turbinaria 

# Now for the number of bins (m). A histogram of the relevant populations should tell us the appropriate (minimum) bin sizes for the population.
# the appropriate number of bins is the point at which adding further bins to the histogram no longer affects the size frequency distribution of the populations
matrix.16 <- subset(turb, t.year == "2016")
hist(matrix.16$Size.t,breaks = 85)
matrix.17 <- subset(turb, t.year == "2017")
hist(matrix.17$Size.t,breaks = 150) # the largest optimum number of bins is taken as the optimum number of bins across all years. The exact number of bins used is then slightly higher than this optimum number. 
# This ensures that all models are of equal dimension and accounts for extended size range being used for integration (U and L).
m = 160

#####################################
# STEP 9: Run the IPMs and estimate lambda
#####################################

#2016
# Build model
turb_2016 <- mk_K.b(m = m, m.par = m.par.year[[1]], L = L, Umax = U, Ut = Ut) # this model selects model coefficients representing the dynamics observed during 2016/17
# estimate and store lambda
lam.t16 <- Re(eigen(turb_2016$K)$values)[1]; lam.t16

# and again for 2017
turb_2017 <- mk_K.n(m = m.year, m.par = m.par.year[[2]], L = L, Umax = U, Ut = Ut) # this model then selects model coefficients representing the dynamics observed during 2017/18
lam.t17 <- Re(eigen(turb_2017$K)$values)[1]; lam.t17

#########################
#Step 10: Jacknife analysis to determine 95% confidence intervals for lambda estimates
#########################
# This involves creating a loop that will repeat the IPM construction above numerous times each time using a sample equal to 95% of the original sample. 
# During each run the loop will construct an IPM - forcing the data to fit the same functions outlined above and estimate a lambda value and store it.
# The mean and confidence intervals of the resultant lambda distributions can then be estimated to provide an indication of our bleachng and no-bleaching model estimates, and whether they are indeed capturing different or similar patterns. 

# reassign the required data to ensure the raw data isn't effected during the loop.
d = turb #restore the data so that it is not changed during the loops.
x = ipm
z = turb.frag.size
# a blank storage matrix for lambda estimates
lambdas.turb <- matrix(NA, 1000, 2) #to store the lambda values of each of the different models.
# also create a list for storing kernals produced during each loop to enable variance estimates to be obtained for later analyses
turb.16 = turb.17 = list(NULL)
t.par.16 = t.par.17 = list(NULL)

# run a loop repeating only the essential steps from above - all regression formats and model functions used will be the same in each loop - only the data using to parameterise each model changes.
for(i in 1:1000){
  
  # determine the new Jacknife sample
  # the use of try should mean the loop will return NA's rather than breaking for errors but mean that = signs can't be used.
  try(d.1 <- d[sample(nrow(d), 0.95*dim(d)[1], replace = F),]) #this randomly takes a 95% subset of the data.
  try(d3 <- x[sample(nrow(x), 0.95*dim(x)[1], replace = F),])
  try(z1 <- z[sample(nrow(z), 0.95*dim(z)[1], replace = F),]) #and the same for the fragment size data.
  
  #now break the subsamples down into the discrete and continuous classes.
  try(d1 <- d.1[which(d.1$Size.t >= discrete_threshold_turbinaria),])
  try(d1a <- d.1[which(d.1$Size.t < discrete_threshold_turbinaria),])
  
  # Now force the models to follow the same regression formats as the 'original' models.
  
  # 1. Growth
  # growth in continous class
  try(try(jgrow.t6 <- glm(Size.t1 ~ Size.t * t.year + I(Size.t^2), data = d.1)))
  # progression into discrete stage
  try(d1a$becomelarge <- 0)
  try(d1a$becomelarge[which(d1a$Size.t1 >= discrete_threshold_turbinaria)] <- 1)
  try(d1a$becomelarge[which(d1a$surv == 0)] <- NA) 
  try(jgrowth.progression.turb <- glm(becomelarge ~ Size.t, family = "binomial", data = d1a))
  # Retrogression out of the discrete stage.
  try(d1$becomesmall <- 0)
  try(d1$becomesmall[which(d1$Size.t1 < discrete_threshold_turbinaria)] <- 1)
  try(d1$becomesmall[which(d1$surv == 0)] <- NA) 
  try(jno_shrinking <- dim(d1[which(d1$becomesmall == 1),])[1])
  try(jstasis <- dim(d1)[1])
  # probability of retrogression is
  try(jretrogression_turb <- jno_shrinking/jstasis)
  # sizes produced by retrogression
  try(jshrink_size_turb <- glm(Size.t1~1, data = d1[which(d1$becomesmall == 1),]))
  
  ## 2. Growth variation
  try(jgrowth.sd.turb <- glm(abs(resid(jgrow.t6))~jgrow.t6$model$Size.t, family = Gamma(link = "log")))
  
  # 3. Survival
  # Continuous class
  try(jsurv.t1 <- glm(surv ~ Size.t * t.year, family = "binomial", data = d.1))
  # discrete survival during bleaching
  try(jturb.mortality_bleaching <- dim(d1[which(d1$surv == 0 & d1$t.year == "2016"),])[1])
  try(jturb.group_size_bleaching <- dim(d1[which(d1$t.year == "2016"),])[1])
  try(jturb.survival_bleaching <- 1-(jturb.mortality_bleaching/jturb.group_size_bleaching))
  # discrete survival during non-bleaching
  try(jturb.mortality <- dim(d1[which(d1$surv == 0 & d1$t.year == "2017"),])[1])
  try(jturb.group_size <- dim(d1[which(d1$t.year == "2017"),])[1])
  try(jturb.survival_n.bleaching <- 1-(jturb.mortality/jturb.group_size))
  # interestingly mortaility in larger corals was worse during the good year - possible lag effect?
  
  ## 4. Fragmentation
  # Continuous class
  try(jfrag.t1 <- glm(frag ~ Size.t * t.year, family = "binomial", data = d.1))
  # discrete fragmentation during bleaching
  try(jturb.fragmenting_bleaching <- dim(d1[which(d1$frag == 1 & d1$t.year == "2016"),])[1])
  try(jturb.fragmentation_bleaching <- jturb.fragmenting_bleaching/jturb.group_size_bleaching)
  # discrete fragmentation during none bleaching
  try(jturb.fragmenting <- dim(d1[which(d1$frag == 1 & d1$t.year == "2017"),])[1])
  try(jturb.fragmentation_n.bleaching <- jturb.fragmenting/jturb.group_size)
  
  # 5. Number of fragments as a function of size
  # for this I am only interested in the corals that did fragment.
  try(jt.frag.corals <- subset(d.1, No.fragments > 0))
  # model independantly of year
  try(jt.no.frag <- glm(No.fragments ~ Size.t, family = "poisson", data = jt.frag.corals)) 
  # now extrapolate for the discrete class
  try(jmean_turb_size <- mean(d1$Size.t, na.rm = T)) 
  # and now extrapolate the size of fragments predicted by the continuous model adove. 
  try(jturb_frag_no <- as.numeric(exp(predict(jt.no.frag, list(Size.t = jmean_turb_size)))))
  
  # 6. Fragment size
  try(jt.frag.size.mod2 <- lm(frag.size ~ Size.t, data = z1))
  # and now extrapolate the size of fragments predicted by the continuous model adove. 
  try(jturb_frag_size <- as.numeric(predict(jt.frag.size.mod2, list(Size.t = jmean_turb_size))))
  
  # 6b. fragment size variability
  try(jt.frag.sd.mod <- glm(abs(resid(jt.frag.size.mod2))~z1$Size.t, family = Gamma(link = "log")))
  
  # 7. Fecundity 
  # this analysis is also easier working with the two year transitions seperately.
  try(jfec.t16 <- subset(d1a, t.year == 2016))
  try(jfec.t17 <- subset(d1a, t.year == 2017))
  # basic none linear model - with this data added artificially there is no need to include random effects of size as these weren't considered when estimating colony fecundity.
  try(jfec.turb.16 <- lm(log(fec) ~ Size.t, data = jfec.t16))
  try(jfec.turb.17 <- lm(log(fec) ~ Size.t, data = jfec.t17))
  # for this I need to determine the per capita fecundity of the mean size of colonies within this discrete size class.
  try(jmean_turb_size.b <- mean(d1[which(d1$t.year == "2016"),]$Size.t, na.rm = T)) 
  try(jmean_turb_size.nb <- mean(d1[which(d1$t.year == "2017"),]$Size.t, na.rm = T)) 
  # and now extrapolate the fecundity predicted by the continuous models adove. 
  try(jturb_fec_bleaching <- as.numeric(exp(predict(jfec.turb.16, list(Size.t = jmean_turb_size.b)))))
  try(jturb_fec_n.bleaching <- as.numeric(exp(predict(jfec.turb.17, list(Size.t = jmean_turb_size.nb)))))
  # this is the mean per captita fecundity of colonies within this size class for bleaching and none bleaching
  
  # 8. Recruitment probability
  # this will be used to constrain the coral fecundity and force it to be similar to the solitarys observations
  try(jturb.16 <- subset(d.1, t.year == 2016))
  try(jturb.17 <- subset(d.1, t.year == 2017))
  # a. 2016
  try(jturb.fec.16 <- sum(jfec.t16$fec, na.rm = T) + #the total fecundity predicted for the continous class
        (jturb_fec_bleaching * dim(d1[which(d1$t.year == "2016" & d1$surv == 1),])[1])) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
  # subset out recruits 
  try(jturb.rec.16 <- subset(jturb.16, Recruit.t1 == "Yes")) 
  try(jturb.no.recruit.16 <- dim(jturb.rec.16)[1]) #total number of reported recruits.
  try(jturb.prob.rec.16 <- jturb.no.recruit.16/jturb.fec.16) #recruitment liklihood during bleaching years
  
  # a. 2017
  try(jturb.fec.17 <- sum(jfec.t17$fec, na.rm = T) + #the total fecundity predicted for the continous class
        (jturb_fec_n.bleaching * dim(d1[which(d1$t.year == "2017" & d1$surv == 1),])[1])) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
  # subset out recruits 
  try(jturb.rec.17 <- subset(jturb.17, Recruit.t1 == "Yes")) 
  try(jturb.no.recruit.17 <- dim(jturb.rec.17)[1]) #total number of reported recruits.
  try(jturb.prob.rec.17 <- jturb.no.recruit.17/jturb.fec.17) #recruitment liklihood during bleaching years
  
  ## 9. Recruit size
  try(jrecruits <- subset(d.1, Recruit.t1 == "Yes"))
  try(jt.rec.size <- selm(Size.t1~1, data = jrecruits))
  
  # 2016 coefficients
  try(jm.par.t16 <- c(
    #continous class
    # survival
    surv.int        =  coefficients(jsurv.t1)[1],
    surv.slope      =  coefficients(jsurv.t1)[2],
    # growth 
    grow.int        =  coefficients(jgrow.t6)[1],
    grow.slope      =  coefficients(jgrow.t6)[2],
    grow.slope2     =  coefficients(jgrow.t6)[4],
    # class progression
    prog.int        =  coefficients(jgrowth.progression.turb)[1],
    prog.slope      =  coefficients(jgrowth.progression.turb)[2],
    # growth variation
    grow.sd.int     =  coefficients(jgrowth.sd.turb)[1], 
    grow.sd.slope   =  coefficients(jgrowth.sd.turb)[2],
    # fragmentation probability as a function of size
    frag.int        =  coefficients(jfrag.t1)[1],
    frag.slope      =  coefficients(jfrag.t1)[2],
    # Number of fragments
    frag.no.int     =  coefficients(jt.no.frag)[1],
    frag.no.slope   =  coefficients(jt.no.frag)[2],
    # fragment size
    frag.size.int   =  coefficients(jt.frag.size.mod2)[1],
    frag.size.slope =  coefficients(jt.frag.size.mod2)[2],
    frag.sd.int     =  coefficients(jt.frag.sd.mod)[1],
    frag.sd.slope   =  coefficients(jt.frag.sd.mod)[2],
    # Discrete class
    # retrogression
    retrogression   =  jretrogression_turb,
    shrink.size.int =  coefficients(jshrink_size_turb)[1],
    shrink.size.sd  =  sigma(jshrink_size_turb),
    #survival
    survival        =  jturb.survival_bleaching,
    #fragmentation
    fragmentation   =  jturb.fragmentation_bleaching,
    frag.number     =  jturb_frag_no,
    frag.size       =  jturb_frag_size,
    frag.size.var   =  sqrt(pi/2) * exp((coefficients(jt.frag.sd.mod)[1] + coefficients(jt.frag.sd.mod)[2] * jmean_turb_size.b)), #this is the extrapolated fragment size variance for the discrete class. 
    #fecundity
    fec             =  jturb_fec_bleaching,
    # recruitment into discrete class
    rec             =  0,
    #unrelated to classes
    # fecundity as a function of size
    fec.int         =  coefficients(jfec.turb.16)[1],
    fec.slope       =  coefficients(jfec.turb.16)[2],
    # recruitment probability
    prob.rec        =  jturb.prob.rec.16,
    # recruit size
    rcsz.xi         =  jt.rec.size@param$dp[1],
    rcsz.omega      =  jt.rec.size@param$dp[2],
    rcsz.alpha      =  jt.rec.size@param$dp[3]))
  
  # 2017 coefficients
  try(jm.par.t17 <- c(
    #continous class
    # survival
    surv.int        =  coefficients(jsurv.t1)[1] + coefficients(jsurv.t1)[3],
    surv.slope      =  coefficients(jsurv.t1)[2] + coefficients(jsurv.t1)[4],
    # growth 
    grow.int        =  coefficients(jgrow.t6)[1] + coefficients(jgrow.t6)[3],
    grow.slope      =  coefficients(jgrow.t6)[2] + coefficients(jgrow.t6)[5],
    grow.slope2     =  coefficients(jgrow.t6)[4],
    # class progression
    prog.int        =  coefficients(jgrowth.progression.turb)[1],
    prog.slope      =  coefficients(jgrowth.progression.turb)[2],
    # growth variation
    grow.sd.int     =  coefficients(jgrowth.sd.turb)[1], 
    grow.sd.slope   =  coefficients(jgrowth.sd.turb)[2],
    # fragmentation probability as a function of size
    frag.int        =  coefficients(jfrag.t1)[1] + coefficients(jfrag.t1)[3],
    frag.slope      =  coefficients(jfrag.t1)[2] + coefficients(jfrag.t1)[4],
    # Number of fragments
    frag.no.int     =  coefficients(jt.no.frag)[1],
    frag.no.slope   =  coefficients(jt.no.frag)[2],
    # fragment size
    frag.size.int   =  coefficients(jt.frag.size.mod2)[1],
    frag.size.slope =  coefficients(jt.frag.size.mod2)[2],
    frag.sd.int     =  coefficients(jt.frag.sd.mod)[1],
    frag.sd.slope   =  coefficients(jt.frag.sd.mod)[2],
    # Discrete class
    # retrogression
    retrogression   =  jretrogression_turb,
    shrink.size.int =  coefficients(jshrink_size_turb)[1],
    shrink.size.sd  =  sigma(jshrink_size_turb),
    #survival
    survival        =  jturb.survival_n.bleaching,
    #fragmentation
    fragmentation   =  jturb.fragmentation_n.bleaching,
    frag.number     =  jturb_frag_no,
    frag.size       =  jturb_frag_size,
    frag.size.var   =  sqrt(pi/2) * exp((coefficients(jt.frag.sd.mod)[1] + coefficients(jt.frag.sd.mod)[2] * jmean_turb_size.nb)), #this is the extrapolated fragment size variance for the discrete class. 
    #fecundity
    fec             =  jturb_fec_n.bleaching,
    # recruitment into discrete class
    rec             =  0,
    #unrelated to classes
    # fecundity as a function of size
    fec.int         =  coefficients(jfec.turb.17)[1],
    fec.slope       =  coefficients(jfec.turb.17)[2],
    # recruitment probability
    prob.rec        =  jturb.prob.rec.17,
    # recruit size
    rcsz.xi         =  jt.rec.size@param$dp[1],
    rcsz.omega      =  jt.rec.size@param$dp[2],
    rcsz.alpha      =  jt.rec.size@param$dp[3]))
  
  
  try(names(jm.par.t16) <- names(jm.par.t17) <- c("surv.int","surv.slope", #survival
                                                  "grow.int","grow.slope", "grow.slope2", # growth
                                                  "prog.int", "prog.slope", #progression to discrete
                                                  "grow.sd.int","grow.sd.slope", #growth variance
                                                  "frag.int","frag.slope", #fragmentation with size
                                                  "frag.no.int", "frag.no.slope", #number of fragments
                                                  "frag.s.int", "frag.s.slope", "frag.sd.int", "frag.sd.slope", #fragment size as a function of adult size.
                                                  "retrogression", #retrogression from discrete
                                                  "size.mean", "size.var", #shrinkage size
                                                  "survival", #discrete survival
                                                  "frag.prob", #discrete fragmentation
                                                  "frag.number", # mean number of frags produced by large coral class
                                                  "frag.size", "frag.var", # mean size and variance of frags produced by large coral class
                                                  "fecundity", #discrete fecundity
                                                  "recruitment", #recruitment into discrete class
                                                  "fec.int", "fec.slope", #fecundity with size
                                                  "rec.prob", #recruitment liklihood (success)
                                                  "rcsz.xi", "rcsz.omega", "rcsz.alpha"))  #recruit size
  
  #this just renames parameters names to make it clearer when putting them into the Kernal functions
  try(jm.par.year <- list(jm.par.t16, jm.par.t17))
  
  # define the parameters
  try(jL <- 1.1*min(d.1[,c("Size.t","Size.t1")], na.rm = T)) #minimum size
  try(jU <- 1.1*max(d.1[,c("Size.t","Size.t1")], na.rm = T)) #maximum size ever 
  
  # create the models
  try(jturb_2016 <- mk_K.b(m = m.year, m.par = jm.par.year[[1]], L = jL, U = jU, Ut = Ut))
  try(jturb_2017 <- mk_K.n(m = m.year, m.par = jm.par.year[[2]], L = jL, U = jU, Ut = Ut))
  
  # Store kernels
  try(turb.16[[i]] <- jturb_2016$K)
  try(turb.17[[i]] <- jturb_2017$K)
  
  #store parameters
  try(t.par.16[[i]] <- jm.par.t16)
  try(t.par.17[[i]] <- jm.par.t17)
  
  #store lambdas
  try(lambdas.turb[i,1] <- Re(eigen(jturb_2016$K)$value)[1])
  try(lambdas.turb[i,2] <- Re(eigen(jturb_2017$K)$value)[1])
}

# this now produced a dataframe of possible lambda values for the bleaching and non-bleaching models.
# now to calculate the variation metrics

# standard deviation
sd(lambdas.turb[,1], na.rm = T)
sd(lambdas.turb[,2], na.rm = T)

# the confidence intervals and mean lambda values for each model
CI(na.omit(lambdas.turb[,1]))
CI(na.omit(lambdas.turb[,2]))

# test the difference between the lambda values
t.test(lambdas.turb[,1],lambdas.turb[,2])

################################################################################# End of Code ################################################################################
