# This script is for sorting raw demographic data in preperation for construction of Integral Projection Models representing the dynamics of morphological coral assemblages in the Solitary Islands Marine Park.
# Data sorting will include the adjustment of coral sizes recorded in April 2016 to estimate colony sizes in August 2016 - this is to correct for a mismatch in census intervals with surveys in 2016-2017 representing an interval of 16 months, whereas surveys in 2017-2018 represent a 12 month interval.
# Colonies will also be assigned fecundity estimates based on a size-specific fecundity data obtained from the Coral Trait database (Madin et al 2016), and orginally recorded by Hall & Hughes (1996).
# The raw data file used for the code below must at least contain columns denoting the taxonomy of each tagged colony (Family, Genus, and species), colony tag reference information (Tag number), census location, census year (year t), colony size recorded during year t, the colony size the following year (year t+1), the number of fragments produced by each colony between year t and t+1, and an indication of whether each colony represents either a new fragment, a new recruit or neither.
# During each resurvey colonies are entered as a new row, with the colony reference number and year t information used to distinguish between data collected during subsequent surveys.  
# NA must be used in the event of missing entries for any colony.
# Author: James Cant
# Date: July 2019

# Clear work space memory
rm(list=ls(all=TRUE))

# load in required packages.
library(Hmisc)
library(gam)
library(mgcv)
library(arm)
library(fields)
library(lme4)
library(reshape2)

##################################################
# STEP 1: Load and format data
##################################################

# set working directory
setwd("insert directory pathway")

# raw demographic data file
ipm <- read.csv("File name1.csv")

# load in colony fecundity data file extracted from the Coral trait database (https://coraltraits.org/traits?utf8=%E2%9C%93&trait_class=Reproductive) - select the colony fecundity dataset for download.
fec <- read.csv("File name2.csv")
  # the fecundity dataset needs subsetting based on colony ID so that the fecundity of survey colonies can be estimated using data from similar genera reported by  Hall & Hughes (1996)
  Branch.fec <- subset(fec, specie_name %in% c("Acropora gemmifera", "Acropora hyacinthus", "Acropora millepora", "Acropora nana"))
  Encrusting.fec <- subset(fec, specie_name == "Goniastrea retiformis")
  Small.fec <- subset(fec, specie_name == "Stylophora pistillata")
  # this will allow for the estimation of colony fecundity based on their taxonomy

# Remove all none hard corals form the demographic data 
ipm <- subset(ipm, Group == "Scleractinia")

# before continuing some of the data variables may need reformating.
# firstly remove any corals lost during each projection interval (these individuals should have 'Lost' entered instead of a size as either size at t or t+1) 
ipm <- subset(ipm, Size.t != "Lost" | is.na(Size.t)) # <- this removes all rows containing 'Lost' as the size at t, but importantly retains any NA entries!
ipm <- subset(ipm, Size.t1 != "Lost" | is.na(Size.t1)) #<- this does the same for Size at t+1
# Now remove colong entries lost after they fused with other colonies (These individual will have 'Fusion' entered instead of a size at t+1)
ipm <- subset(ipm, Size.t1 != "Fusion" | is.na(Size.t1))

# Size entries now need to be converted to a numercial format,
ipm$Size.t <- as.numeric(paste(ipm$Size.t))
ipm$Size.t1 <- as.numeric(paste(ipm$Size.t1))
# and t.year needs to be assigned factor level for use as a random effect variable in regression analyses
ipm$t.year <- as.factor(ipm$t.year) 

# Ensure size data is in cm2 (in the orginal analyses data was measured in mm2 and so required converting)
ipm$Size.t <- ipm$Size.t/100
ipm$Size.t1 <- ipm$Size.t1/100

##################################################
# STEP 2: Adjust size t entries for colonies measured in April 2016
##################################################
# The census interval between April 2016 and August 2017 is 16 months - to make this comparable to the annual interval between August 2017 and Aug/Sept 2018 this needs standardising to 12 months
# Therefore Size.t estimates for April 2016 need advancing by four months to estimate size.t for August 2016.
# This will have to be done using a mean estimate of monthly growth between April 2016 and Auguast 2017 to allow for sizes to adjusted for colonies that died before August 2017 and so weren't resurveyed.

# Extract relevant size data.
grow.trend <- matrix(NA, length(ipm$Size.t[which(ipm$t.year == "2016")]),3) #blank matrix to store sizes over time with each coral as a row.
grow.trend[,1]  <- ipm$Size.t[which(ipm$t.year == "2016")] #extract size data measured in April 2016
grow.trend[,2]  <- ipm$Size.t[which(ipm$t.month == "October")] # extract size data measured in October 2016
grow.trend[,3]  <- ipm$Size.t1[which(ipm$t.year == "2016")] # extract size data measured in August 2017.
# By using the different size estimates available across the April 2016 - August 2017 interval the estimate of monthly growth will account for chaging growth patterns throughout the year. 

# rename the extracted columns reflecting the number of months after April 2016 that those size measurements were taken - April(0), October (6), August (16).
month <- c(0,6,16)
colnames(grow.trend) <- month

# Remove any corals that weren't tagged during April 2016 (i.e size in Aprik 2016 = NA).
grow.trend <- subset(grow.trend, !is.na(grow.trend[,1]))
# store extract growth data for plotting
for.plot <- grow.trend

#rename the rows to relflect the individual identity of individual colonies. 
colony <- seq(1,dim(grow.trend)[1],1)
rownames(grow.trend) <- colony

# now the change in size of all colonies, tagged in April 2016, can be visualised through to August 2017
matplot(t(grow.trend), type = "l")
# This can be used to estimate a general growth trend.
grow.trend <- melt(grow.trend, value.name = "Size")
colnames(grow.trend) <- c("ColonyID","Month","Size")
grow.trend$Month <- as.numeric(paste(grow.trend$Month)) # ensure month is in numerical format

# run regression comparing colony size with number of months growing
grow.mod <- lm(Size ~ Month, data = grow.trend) # linear model
# how well does the model reflect the data?
abline(coefficients(grow.mod)[1],coefficients(grow.mod)[2], lwd = 3, col = "red") 
# check that a linear model is indeed the best fit
grow.mod2 <- glm(Size ~ Month + I(Month^2), data = grow.trend)
lines(grow.trend$Month, predict(grow.mod2, list(Month = grow.trend$Month)), lwd = 5)
# both models result in a more or less linear trend so whichever model provides the best fit is the better model?
AIC(grow.mod,grow.mod2) #marginally the original model.
# store the model coefficients
a <- coefficients(grow.mod)[1]
b <- coefficients(grow.mod)[2] #here the slope coefficient represents the mean monthly growth increment for the tagged sample 

# plot the growth trend for the supplamentary material
matplot(t(for.plot), type = "l", xaxt = "n", xaxs = "i", yaxs = "i", ylim = c(0,3500), cex.axis = 1.2)
axis(side = 1, at = c(1,2,3), labels = c(0,6,16), cex.axis = 1.2)
lines(grow.trend$Month, predict(grow.mod, list(Month = grow.trend$Month)), lwd = 5)

# The stored model coefficients can now be used to determine the size in August of all corals measured in April (whether they survived or not).
ipm$Size.t[which(ipm$t.year == "2016")] <- (a - (a - ipm$Size.t[which(ipm$t.year == "2016")])) + (b * 4) # here insert the actual size of each coral as the intercept and projecting it forward along the regression slope by four month steps.
# Adjustments are also required for the vital rates of fragmentation and colony survival but these  carried out latter within the Integral projection model framework.

##################################################
# STEP 3: Check and transform size data
##################################################

# Check the distributions of size data (generally size data has a skewed distribution with lots of little and few large individuals)
hist(ipm$Size.t)
hist(ipm$Size.t1)
hist(log(ipm$Size.t))
hist(log(ipm$Size.t1)) # a log transformation achieves a normal distribution when this is the case

# carry out transformation
ipm$Size.t <- log(ipm$Size.t)
ipm$Size.t1 <- log(ipm$Size.t1)

##################################################
# STEP 4: Estimate colony survival and fragmentation probabilities
##################################################
# Survival and fragmentation need to be added in a binary format denoting whether they occured (1) or not (0)

# Colonies are considered to have survived if they have size estimates for both time t and t+1.
# Colonies with size estimates only for time t are recorded as dead.
# Colonies without size estimates for time t but with size estimates for time t+1 are recruits or new fragements and so their survival is included by their presence in the estimates of recruitment or fragmentation as so is ignored here.

# Add a blank survival colony
ipm$surv = NA
ipm$surv[which(!is.na(ipm$Size.t) & !is.na(ipm$Size.t1))]=1 #<- when both Size t and Size t+1 are both NOT NA then surv = 1
ipm$surv[which(!is.na(ipm$Size.t) & is.na(ipm$Size.t1))]=0 #<- when Size t is NOT NA but Size t+1 is then surv = 0
table(ipm$surv) #<- good way to check the data 

# Add a blank fragmentation column
ipm$frag = NA
# Fragmentation is added in a similar way reflecting if colonies produced any fragments (1) or not (0)
ipm$frag[which(ipm$No.fragments >= 0.5)]=1 #<- when the number of fragments produced is greater than 0 then fragmentation = 1
ipm$frag[which(ipm$No.fragments == 0)]=0 #<- when the number of fragments = 0 then fragmentation = 0
ipm$frag[which(is.na(ipm$Size.t))] = NA # this ensures that any coral with no size at time t is not assigned a fragmentation level.
table(ipm$frag)

##################################################
# STEP 5: Assign morpholoical clusters ready for estimating fecundity
##################################################
summary(ipm$Genus) #the spread of diversity (and the importance of structure in coral communties) may mean its better using morphology for subsetting the population rather than genus

# Here we are clustering corals into one of four morphological groups: Encrusting/Massive, Small Branching (Pocilloporidae), Branching, Laminar.
# as outlined in our manuscript these morphological groups largely reflected coral genera.
# Add blank morphological group column
ipm$morphology = NA
ipm$morphology[which(ipm$Genus %in% c("Acanthastrea", "Astrea", "Goniopora", "Micromussa", "Montipora", "Paragoniastrea", "Porites", "Dispsastraea"))] = "Encrusting_Massive"
ipm$morphology[which(ipm$Genus %in% c("Pocillopora", "Stylophora"))] = "Small_branch"
ipm$morphology[which(ipm$Genus == "Turbinaria")] = "Laminar"
ipm$morphology[which(ipm$Genus == "Acropora")] = "Branching_plates"
ipm$morphology <- as.factor(paste(ipm$morphology)) #convert morphology type to a factor

##################################################
# STEP 6: Estimate sizespecific fecundity relationship and use to estimate the fecundity of observed colonies
##################################################
# The following section of code is for estimating the fecundity of observed colonies reflecting the size-specific fecundity relationship reported by Hall & Hughes (1996)
# In the Hall & Hughes (1996) dataset, total colony fecundity was estimated as the sum of egg and testes volume per polyp multiplied by the total number of polyps.
# This is therefore the relationship between colony size and fecundity (which will need dividing by two - to account for the inclusion of testes in the measure). This fecundity relationship will be estimate seperately corals of differing morphological types

# condense the different size-specific fecundity data sets for use in a for loop
fec.index <- list(Branch.fec, Encrusting.fec, Small.fec)
# create blank storage frames.
fec.data <- list()
fec.mod <- list()
fec.mod2 <- list()
fec.mod3 <- list()

# This for loop will work through each fecundity data set to estimate a size-specific fecundity relationship for each morphological group.
for(i in  1:3){
  fec <- fec.index[[i]] #which morphology type is being worked with
  #subset out the fecundity data
  indexfecundity = which(fec$trait_name == "Colony fecundity") 
  fecundity <- as.vector(fec$value[indexfecundity]) 
  fecundity <- as.integer(fecundity) # these are our gonad densities but need adjusting to cm3
  fecundity <- fecundity/100
  
  #subset out the colony size data
  indexsize=which(fec$trait_name == "Colony area") 
  size <- as.vector(fec$value[indexsize])
  size <- as.numeric(size) #this is our size data - but it doesn't need adjusting as it is already in cm2
  size <- log(size) # log transform to ensure scales match across datasets
  
  # store the extracted data
  fec.data[[i]] <- data.frame(fecundity,size)
  # Now to determine the relationship between size and fecundity for the coral trait database data
  # Regression analysis will be used to determine the relationship - however because the relationship will likely change between different coral groups the loop will run three different models
  # and store the outputs so that model fit can be evaluated ensuring we retain the most accurate relationship.
  fec.mod[[i]] <- nls(fecundity ~ exp(a + b * size), data = fec.data[[i]], start = list(a = 1, b = 1))
  fec.mod2[[i]] <- gam(fecundity ~ s(size), data = fec.data[[i]])
  fec.mod3[[i]] <- glm(fecundity ~ size + I(size^2) + I(size^2), data = fec.data[[i]])
}

# check the model fits and data for the different morpholigcal groups. 
#Branching
plot(fecundity~size, data = fec.data[[1]], pch=20, cex=0.9, xlim = c(0,10), ylim = c(0,8000))
lines(spline(fec.data[[1]]$size, 
             predict(fec.mod[[1]], list(Size = fec.data[[1]]$size)),
             method = "n", n =250), col = "black") 
lines(spline(fec.data[[1]]$size, 
             predict(fec.mod2[[1]], list(size = fec.data[[1]]$size)),
             method = "n", n =250), col = "blue")
lines(spline(fec.data[[1]]$size, 
             predict(fec.mod3[[1]], list(size = fec.data[[1]]$size)),
             method = "n", n =250), col = "red")

#Encrusting
plot(fecundity~size, data = fec.data[[2]], pch=20, cex=0.9, xlim = c(0,7), ylim = c(0,300), xaxs = "i", yaxs = "i") 
lines(spline(fec.data[[2]]$size, 
             predict(fec.mod[[2]], list(size = fec.data[[2]]$size)),
             method = "n", n =250), col = "black")
lines(spline(fec.data[[2]]$size, 
             predict(fec.mod2[[2]], list(size = fec.data[[2]]$size)),
             method = "n", n =250), col = "blue")
lines(spline(fec.data[[2]]$size, 
             predict(fec.mod3[[2]], list(size = fec.data[[2]]$size)),
             method = "n", n =250), col = "red") 

#small Branching
plot(fecundity~size, data = fec.data[[3]], pch=20, cex=0.9, xlim = c(0,6), ylim = c(0,300), xaxs = "i", yaxs = "i") 
lines(spline(fec.data[[3]]$size, 
             predict(fec.mod[[3]], list(size = fec.data[[3]]$size)),
             method = "n", n =250), col = "black") 
lines(spline(fec.data[[1]]$size, 
             predict(fec.mod2[[3]], list(size = fec.data[[1]]$size)),
             method = "n", n =250), col = "blue")
lines(spline(fec.data[[1]]$size, 
             predict(fec.mod3[[3]], list(size = fec.data[[1]]$size)),
             method = "n", n =250), col = "red") #gam again best at intermediate sizes.
# From these plots the most accurate relationship between colony size and fecundity can be determined and used below to estimat the funcdity of our observed colonies given their size and morphology.

# add fecundity column to main demographic data file
ipm$fec = NA
# The morphological groups Large Branching and Laminar (Acropora and Turbinaria) will have their fecundities predicted using the best fitting model from the Branching.fec dataset.
ipm$fec[which(ipm$morphology %in% c("Branching_plates", "Laminar"))] = predict(fec.mod[[1]], list(size = ipm$Size.t[which(ipm$morphology %in% c("Branching_plates", "Laminar"))])) # given their size this model will predict our colonies fecundity
# The morphological group Encrusting will have their fecundities predicted using the best fitting model from the Encrusting.fec dataset.
ipm$fec[which(ipm$morphology == "Encrusting_Massive")] = predict(fec.mod[[2]], list(size = ipm$Size.t[which(ipm$morphology == "Encrusting_Massive")]))
# The morphological group Small_branch will have their fecundities predicted using the best fitting model from the Small.fec data.
ipm$fec[which(ipm$morphology == "Small_branch")] = predict(fec.mod[[3]], list(size = ipm$Size.t[which(ipm$morphology == "Small_branch")]))

# Just a couple of fixes
ipm$fec[which(ipm$surv == 0)]= NA #No measure of fecundity if survival = 0
ipm$fec[which(ipm$Recruit.t1 == "Yes")]= NA #prevents assigning fecundity to new recruits.
ipm$fec[which(ipm$fec < 0)] = 0 #this ensures no corals have been predicted negative fecundity.
range(ipm$fec, na.rm = T) # check that has worked

# the Hall & Hughes (1996) data set determined fecundity as combined teste and egg density, we are assuming larval density output is half this so these values need dividing by 2.
ipm$fec <- ipm$fec/2

##################################################
# STEP 7: Final tidy up incase structures have been changed during formatting
##################################################
# The following variables are needed in factor format for use as random variables in vital rate regression analyses.

ipm$location <- as.factor(ipm$location)
ipm$Site <- as.factor(ipm$Site)
ipm$t.year <- as.factor(ipm$t.year)

# the dataframes are now ready for analysis.
# Store formatted data file
write.csv(ipm, file = "File name3.csv") 

 ####################################################################### End of Code #####################################################################################