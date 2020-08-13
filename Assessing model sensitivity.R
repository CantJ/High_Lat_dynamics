# In the manuscript associated with this R code we used a two-stage IPM model framework. 
# Our model consisted of two size stages, one continuous and the other discrete, to account for the census data regarding the dynamics of large colonies.
# As such colonies transition along a continuous size scale until they achieved a certain colony size, after which they progress into a discrete size stage subjected to fixed vital rate processes
# This script is designed to test the 
# Date last modified: January 2019
# Author: James Cant

# Testing the sensitivity of a two stage IPM to the positioning of a threshold size between continuous and discrete size classes requires a loop that sequentially changes the position of the size threshold and recalculates lambda.

# In this script the example provided was used for the testing the sensitivites of the Acropora IPMs.

####################################
# Step 1: Define the necessary data and generate storage outputs to be used
####################################

# store the nessecary data 
d = acro # this is the formatted demographic data file used to construct the original IPM models.
z1 = acro.frag.size # formatted data file used in constructing the fragmentation components of the original data files. 
# create a blank matrix for storing the lambda estimates and size threshold values for each model re run.
# the loop will run 1000 times - each time selecting a new size threshold value. This size threshold will be used to create new IPMs for both bleaching and non-bleaching periods.
# The loop will subsequently estimate lambda for the new models and store them with the corresponding size threshold value. 
outputs.acro <- matrix(NA, 1000, 3)
colnames(outputs.acro) <- c("size.threshold", "lambda.b", "lambda.nb") # rename matrix columns for ease
# create a sequence of threshold sizes to be used (each value will differ by 0.001 from the previous)
threshold.seq <- c(seq(to = discrete_threshold_acro, by = 0.001, length.out = 500),
                   seq(from = (discrete_threshold_acro + 0.001), by = 0.001, length.out = 500)) #this way the original threshold size is positioned centrally.
#store the size threshold sequence.
outputs.acro[,1] <- threshold.seq

####################################
# Step 2: Run the sensitivity analysis loop
####################################

# This loop will select a size threshold value from the available sequence.
# Bleaching and non-bleaching models will then be constructed following the same format as used in the construction of the original models - see the IPM cotruction script (https://github.com/CantJ/Subtropical-Coral-Demographics) for further details on the construction of IPM models
# the loop will then estimate and store the values of lambda from each new model.
for(i in 1:1000){ # the loop will run 1000 times.
  
  #select the discrete size threshold from the sequence.
  jdiscrete_threshold_acro <- outputs.acro[i,1]
  
  # now break the main data file down into discrete and continuous classes using this threshold size value.
  try(d1 <- d[which(d$Size.t >= jdiscrete_threshold_acro),])
  try(d1a <- d[which(d$Size.t < jdiscrete_threshold_acro),])
  
  # Now force the models to follow the same regression formats as the 'original' models.
  # 1. Growth
  try(jgrow.a5 <- glm(Size.t1 ~ Size.t * t.year + I(Size.t^2), data = d))
  
  # Now work out progression probabilities into the discrete stage for the inshore population.
  try(d1a$becomelarge <- 0)
  try(d1a$becomelarge[which(d1a$Size.t1 >= jdiscrete_threshold_acro)] <- 1)
  try(d1a$becomelarge[which(d1a$surv == 0)] <- NA) 
  try(jgrowth.progression.a1 <- glm(becomelarge ~ Size.t, family = "binomial", data = d1a))
  # Now to work out the retrogression out of the discrete stage.
  try(d1$becomesmall <- 0)
  try(d1$becomesmall[which(d1$Size.t1 < jdiscrete_threshold_acro)] <- 1)
  try(d1$becomesmall[which(d1$surv == 0)] <- NA) 
  # Now determine the probability of transferring out off the discrete stage and the likely resultant size of these colonies.
  try(jno_shrinking.a <- dim(d1[which(d1$becomesmall == 1),])[1])
  try(jstasis.a <- dim(d1)[1])
  # so the probability of retrogression is
  try(jretrogression_acro <- jno_shrinking.a/jstasis.a)
  # sizes produced
  try(if (jretrogression_acro == 0) {jshrink_size_acro <- 0} else {jshrink_size_acro <- glm(Size.t1~1, data = d1[which(d1$becomesmall == 1),])})
  
  # 2. Growth variability 
  try(jgrowth.sd.acro <- glm(abs(resid(jgrow.a5))~jgrow.a5$model$Size.t, family = Gamma(link = "log")))
  
  # 3. Survival
  try(jsurv.a1 <- glm(surv ~ Size.t * t.year, family = "binomial", data = d)) 
  
  # discrete stage survival during bleaching
  try(jacro.mortality_bleaching <- dim(d1[which(d1$surv == 0 & d1$t.year == "2016"),])[1])
  try(jacro.group_size_bleaching <- dim(d1[which(d1$t.year == "2016"),])[1])
  try(jacro.survival_bleaching <- 1-(jacro.mortality_bleaching/jacro.group_size_bleaching))
  # discrete stage survival during non-bleaching
  try(jacro.mortality <- dim(d1[which(d1$surv == 0 & d1$t.year == "2017"),])[1])
  try(jacro.group_size <- dim(d1[which(d1$t.year == "2017"),])[1])
  try(jacro.survival_n.bleaching <- 1-(jacro.mortality/jacro.group_size))
  
  # 4. Fragmentation
  try(jfrag.a1 <- glm(frag ~ Size.t * t.year, family = "binomial", data = d))
  # Discrete class
  # fragmentation during bleaching
  try(jacro.fragmenting_bleaching <- dim(d1[which(d1$frag == 1 & d1$t.year == "2016"),])[1])
  try(jacro.fragmentation_bleaching <- jacro.fragmenting_bleaching/jacro.group_size_bleaching)
  # fragmentation during none bleaching
  try(jacro.fragmenting <- dim(d1[which(d1$frag == 1 & d1$t.year == "2017"),])[1])
  try(jacro.fragmentation_n.bleaching <- jacro.fragmenting/jacro.group_size)
  
  # 5. Number of fragments as a function of size
  try(jfrag.corals <- subset(d, No.fragments > 0))
  try(ja.no.frag <- glm(No.fragments ~ Size.t, family = "poisson", data = jfrag.corals)) #this will be poisson distribution as number of fragments can't be negative and is a form of discrete count data.
  # now extrapolate for the discrete class
  try(jmean_acro_size <- mean(d1$Size.t, na.rm = T))
  # and now extrapolate the size of fragments predicted by the continuous model adove. 
  try(jacro_frag_no <- as.numeric(exp(predict(ja.no.frag, list(Size.t = jmean_acro_size)))))
  
  # 6. Fragment size
  try(ja.frag.size.mod <- lm(frag.size ~ Size.t, data = z1))
  try(jratio <- as.numeric(coefficients(ja.frag.size.mod)[1]/z1$Size.t[1]))
  try(jacro_frag_size <- jratio * jmean_acro_size)
  
  ########## Recruitment ########
  # 7. Fecundity
  # this analysis is also easier working with the two year transitions seperately.
  try(jfec.a16 <- subset(d1a, t.year == 2016))
  try(jfec.a17 <- subset(d1a, t.year == 2017)) 
  # basic none linear model - with this data added artificially there is no need to include random effects of size as these weren't considered when estimating colony fecundity.
  try(jfec.acro.16 <- lm(log(fec) ~ Size.t, data = jfec.a16))
  try(jfec.acro.17 <- lm(log(fec) ~ Size.t, data = jfec.a17))
  # Discrete class
  # determine the per capita fecundity of mean sized colonies within the discrete size class.
  try(jmean_acro_size.b <- mean(d1[which(d1$t.year == "2016"),]$Size.t, na.rm = T)) #during bleaching
  try(jmean_acro_size.nb <- mean(d1[which(d1$t.year == "2017"),]$Size.t, na.rm = T)) #during non-bleaching
  # and now extrapolate the fecundity predicted by the continuous models adove. 
  try(jacro_fec_bleaching <- exp(predict(jfec.acro.16, list(Size.t = jmean_acro_size.b))))
  try(jacro_fec_n.bleaching <- exp(predict(jfec.acro.17, list(Size.t = jmean_acro_size.nb))))
  # this is the mean per captita fecundity of colonies within the discrete size class for bleaching and none bleaching
  
  # 8. Recruitment probability
  # this will be used to constrain the coral fecundity and force it to be similar to field observations
  try(jacro.16 <- subset(d, t.year == 2016))
  try(jacro.17 <- subset(d, t.year == 2017))
  # a. 2016
  try(jacro.fec.16 <- sum(jfec.a16$fec, na.rm = T) + #the total fecundity predicted for the continous class
        (jacro_fec_bleaching * dim(d1[which(d1$t.year == "2016" & d1$surv == 1),])[1])) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
  # subset out recruits 
  try(jacro.rec.16 <- subset(jacro.16, Recruit.t1 == "Yes"))
  try(jacro.no.recruit.16 <- dim(jacro.rec.16)[1]) #total number of reported recruits.
  try(jacro.prob.rec.16 <- jacro.no.recruit.16/jacro.fec.16) #recruitment liklihood during bleaching years
  
  # b. 2017
  try(jacro.fec.17 <- sum(jfec.a17$fec, na.rm = T) + #the total fecundity predicted for the continous class
        (jacro_fec_n.bleaching * dim(d1[which(d1$t.year == "2017" & d1$surv == 1),])[1])) # plus the estimated fecundity given the number of individuals in the descrete stage and their mean fecundity
  # subset out recruits 
  try(jacro.rec.17 <- subset(jacro.17, Recruit.t1 == "Yes")) 
  try(jacro.no.recruit.17 <- dim(jacro.rec.17)[1]) #total number of reported recruits.
  try(jacro.prob.rec.17 <- jacro.no.recruit.17/jacro.fec.17) #recruitment liklihood during non-bleaching years
  
  ## 9. Recruit size
  try(jrecruits <- subset(d, Recruit.t1 == "Yes")) 
  try(ja.rec.size <- lm(Size.t1~1, data = jrecruits))
  
  #---------------------------------------- Store the relevant model coefficients.
  
  # 2016 coefficients
  try(jm.par.a16 <- c(
    #continous class
    # survival
    surv.int        =  coefficients(jsurv.a1)[1],
    surv.slope      =  coefficients(jsurv.a1)[2],
    # growth 
    grow.int        =  coefficients(jgrow.a5)[1],
    grow.slope      =  coefficients(jgrow.a5)[2],
    grow.slope2     =  coefficients(jgrow.a5)[4],
    # class progression
    prog.int        =  coefficients(jgrowth.progression.a1)[1],
    prog.slope      =  coefficients(jgrowth.progression.a1)[2],
    # growth variation
    grow.sd.int     =  coefficients(jgrowth.sd.acro)[1], 
    grow.sd.slope   =  coefficients(jgrowth.sd.acro)[2],
    # fragmentation probability as a function of size
    frag.int        =  coefficients(jfrag.a1)[1],
    frag.slope      =  coefficients(jfrag.a1)[2],
    # Number of fragments
    frag.no.int     =  coefficients(ja.no.frag)[1],
    # fragment size
    frag.size.ratio =  jratio,
    frag.sd         =  summary(ja.frag.size.mod)$sigma,
    # Discrete class
    # retrogression
    retrogression   =  jretrogression_acro,
    shrink.size.int =  coefficients(jshrink_size_acro)[1],
    shrink.size.sd  =  sqrt(pi/2) * exp((coefficients(jgrowth.sd.acro)[1] + coefficients(jgrowth.sd.acro)[2] * jmean_acro_size.b)),
    #survival
    survival        =  jacro.survival_bleaching,
    #fragmentation
    fragmentation   =  jacro.fragmentation_bleaching,
    frag.number     =  jacro_frag_no,
    frag.size       =  jacro_frag_size,
    frag.size.var   =  summary(ja.frag.size.mod)$sigma, #this is the extrapolated fragment size variance for the discrete class. 
    #fecundity
    fec             =  jacro_fec_bleaching,
    # recruitment into discrete class
    rec             =  0,
    #unrelated to classes
    # fecundity as a function of size
    fec.int         =  coefficients(jfec.acro.16)[1],
    fec.slope       =  coefficients(jfec.acro.16)[2],
    # recruitment probability
    prob.rec        =  jacro.prob.rec.16,
    # recruit size
    rcsz.mean       =  coefficients(ja.rec.size)[1],
    rcsz.sd         =  summary(ja.rec.size)$sigma))
  
  # 2017 coefficients
  try(jm.par.a17 <- c(
    #continous class
    # survival
    surv.int        =  coefficients(jsurv.a1)[1] + coefficients(jsurv.a1)[3],
    surv.slope      =  coefficients(jsurv.a1)[2] + coefficients(jsurv.a1)[4],
    # growth 
    grow.int        =  coefficients(jgrow.a5)[1] + coefficients(jgrow.a5)[3],
    grow.slope      =  coefficients(jgrow.a5)[2] + coefficients(jgrow.a5)[5],
    grow.slope2     =  coefficients(jgrow.a5)[4],
    # class progression
    prog.int        =  coefficients(jgrowth.progression.a1)[1],
    prog.slope      =  coefficients(jgrowth.progression.a1)[2],
    # growth variation
    grow.sd.int     =  coefficients(jgrowth.sd.acro)[1], 
    grow.sd.slope   =  coefficients(jgrowth.sd.acro)[2],
    # fragmentation probability as a function of size
    frag.int        =  coefficients(jfrag.a1)[1] + coefficients(jfrag.a1)[3],
    frag.slope      =  coefficients(jfrag.a1)[2] + coefficients(jfrag.a1)[4],
    # Number of fragments
    frag.no.int     =  coefficients(ja.no.frag)[1],
    # fragment size
    frag.size.ratio =  jratio,
    frag.sd         =  summary(ja.frag.size.mod)$sigma,
    # Discrete class
    # retrogression
    retrogression   =  jretrogression_acro,
    shrink.size.int =  coefficients(jshrink_size_acro)[1],
    shrink.size.sd  =  sqrt(pi/2) * exp((coefficients(jgrowth.sd.acro)[1] + coefficients(jgrowth.sd.acro)[2] * jmean_acro_size.nb)),
    #survival
    survival        =  jacro.survival_n.bleaching,
    #fragmentation
    fragmentation   =  jacro.fragmentation_n.bleaching,
    frag.number     =  jacro_frag_no,
    frag.size       =  jacro_frag_size,
    frag.size.var   =  summary(ja.frag.size.mod)$sigma, #this is the extrapolated fragment size variance for the discrete class. 
    #fecundity
    fec             =  jacro_fec_n.bleaching,
    # recruitment into discrete class
    rec             =  0,
    #unrelated to classes
    # fecundity as a function of size
    fec.int         =  coefficients(jfec.acro.17)[1],
    fec.slope       =  coefficients(jfec.acro.17)[2],
    # recruitment probability
    prob.rec        =  jacro.prob.rec.17,
    # recruit size
    rcsz.mean       =  coefficients(ja.rec.size)[1],
    rcsz.sd         =  summary(ja.rec.size)$sigma))
  
  try(names(jm.par.a16) <- names(jm.par.a17) <- c("surv.int","surv.slope", #survival
                                                  "grow.int","grow.slope", "grow.slope2", # growth
                                                  "prog.int", "prog.slope", #progression to discrete
                                                  "grow.sd.int","grow.sd.slope", #growth variance
                                                  "frag.int","frag.slope", #fragmentation with size
                                                  "frag.no", #number of fragments
                                                  "frag.s.ratio","frag.sd", #fragment size as a function of adult size.
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
                                                  "rcsz.mean", "rcsz.sd"))  #recruit size
  
  #this just renames parameters names to make it clearer when putting them into the Kernal functions
  try(jm.par.year <- list(jm.par.a16, jm.par.a17))
  
  # calculate model variables (U, L, Ut and m)
  try(L <- 0.9*min(d[,c("Size.t","Size.t1")], na.rm = T)) 
  try(U <- 1.1*max(d[,c("Size.t","Size.t1")], na.rm = T)) 
  try(Ut <- jdiscrete_threshold_acro)
  try(m <- 110)
  
  #build the IPMs
  try(jacro_2016 <- mk_K.b(m = m, m.par = jm.par.year[[1]], L = L, Umax = U, Ut = Ut))
  try(jacro_2017 <- mk_K.n(m = m, m.par = jm.par.year[[2]], L = L, Umax = U, Ut = Ut))
  
  # calculate and store lambda
  try(outputs.acro[i,2] <- Re(eigen(jacro_2016$K)$value)[1])
  try(outputs.acro[i,3] <- Re(eigen(jacro_2017$K)$value)[1])
}

####################################
# Step 3: Display the effect of the size threshold on model lambda estimates
####################################

# set plotting window.
par(mfrow = c(2,2))    

# use regression analysis to determine the relationship between lambda estimates and the positioning of size thresholds.
# Typically the sentivity of lambda to changes in model variables follows a non-linear shape. Does the sensitivity of the model differ between bleaching and non-bleaching models?
acro.lambda1 <- glm(outputs.acro[,2]~outputs.acro[,1] + I(outputs.acro[,1]^2)) #bleaching model
acro.lambda2 <- glm(outputs.acro[,3]~outputs.acro[,1] + I(outputs.acro[,1]^2)) #non-bleaching model

# plot the estimates of lambda against the value of the size threshold to show the relationship, which is clearer to see if only a few points are plotted from the full sample.
# select the points to be plotted
plot.acro <- outputs.acro[seq(1, 1000, by = 30),]
# create a blank plot
plot(outputs.acro[,2]~outputs.acro[,1], typ = "n", ylim = c(0.15,1.05), xlim = c(6,7), xlab = "", ylab = "", yaxs = "i", xaxs = "i", cex.axis = 1.5)
# add the selected points from both the bleaching and non-bleaching models
points(plot.acro[,2]~plot.acro[,1], typ = "p", cex = 0.8, col = "red", pch = 16)
points(plot.acro[,3]~plot.acro[,1], typ = "p", cex = 0.8, col = "blue", pch = 16)
# add regression lines to the plot to show the overall trend
lines(outputs.acro[,1], predict(acro.lambda1, type='response'), col = "red", lwd = 3)
lines(outputs.acro[,1], predict(acro.lambda2, type='response'), col = "blue", lwd = 3)
# add a vertical line to show the positioning of the actual threshold value used.
abline(v = discrete_threshold_acro, lty = "dashed", lwd = 2)

# Calculate the numerical sensitivity of the models to the postioning of the size threshold. Sensitivity = delta(lambda)/delta(size_threshold)
#Bleaching model sensitivity
(outputs.acro[1000,2]-outputs.acro[1,2])/(outputs.acro[1000,1]-outputs.acro[1,1])
range(outputs.acro[,2]) #range of lambda
# non-bleaching model sensitivity
(outputs.acro[1000,3]-outputs.acro[1,3])/(outputs.acro[1000,1]-outputs.acro[1,1])
range(outputs.acro[,3])

#################################################################################### End of Code ##############################################################################