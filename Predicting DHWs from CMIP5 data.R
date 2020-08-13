# Script is for extracting a temperature projection for a given location from different models projecting climates under different RCP pathways
# It is advisable to use multiple models for each projection and taking a mean to account for model variance.
# Author: James Cant


#load required packages 
library(ncdf4) #for working with nc files
library(dplyr)
library(lme4)
library(pscl) #for calculating pseudo-R2 values for glm models
library(ggeffects)
library(lsmeans)

##########################################
# Extract projected temperature data
##########################################

# first you will need to download a series of different climate projection files covering the RCP pathways you are interested in. 
# freely available climate models can be found at http://data.ceda.ac.uk/badc/cmip5/data/cmip5/output1
# each folder on the webpage represents a different climate projection institution. 
# within each of these subfolders there are then folders for each form of model they use. In these models you can then find data for a whole host of paramters - you'll want the rcp25 to rcp85 parameters (though not all models have all the RCPs)
# select day -> ocean -> day -> select a model run through -> files -> tos
# this gives you a lost of downloadable files of daily SST projections under the selected RCP pathway for the whole globe (we subset by location later)
# Download the files for the dates you want from as many different model types as you can.

# now use the code below
files = list.files('this is the folder pathway to allow R to find your downloaded climate model files',pattern='*.nc',full.names = TRUE) #this extracts all the nc files for each model
# now split the models down into the different RCP pathways they represent.
try(RCP.85_subset <- files[grep(pattern = "rcp85", files)])
try(RCP.60_subset <- files[grep(pattern = "rcp60", files)]) #finds the given pattern in the file names
try(RCP.45_subset <- files[grep(pattern = "rcp45", files)])
try(RCP.26_subset <- files[grep(pattern = "rcp26", files)])
# store these subsets in a list
RCP.store <- list(RCP.85_subset, RCP.60_subset, RCP.45_subset, RCP.26_subset)
# insert your GPSD points of interest
lat = XXXXXX
long = XXXXXX

# set up storage list for each of the different pathways (this will house each of the storage lists returned by each of the RCP pathways - 4 in total))
Pathways <- list(NULL)

# For each pathway loop over each nc file to extract the future projected temperature for the given GPS location predicted by each model.
try(for(RCP in 1:4) { #this loop will run once for each RCP pathway
  
  # define which subset is being worked with
  RCP.models <- RCP.store[[RCP]]
  
  # generate a storage list in which to store each model output
  projections <- list(NULL)
  
  # Now loop through each of the models for this RCP pathway and extract the projected temperature for the Solitaries and store in the projections list.
  for (x in 1:length(RCP.models)) { #this loop will process each avaliable model in turn.
    
    # open the nc file for the selected model
    nc <- nc_open(RCP.models[x])
    
    # determine the time frame for the selected model and the starting date information for the loop below
    year <- as.numeric(substr(nc[["filename"]], nchar(nc[["filename"]])-19, nchar(nc[["filename"]])-16)) #extracts the start year of the model
    day <- as.numeric(substr(nc[["filename"]], nchar(nc[["filename"]])-13, nchar(nc[["filename"]])-12)) #extracts the start day from the model
    
    # how many days does this model project over 
    duration <- dim(ncvar_get(nc, "time"))
    
    #create a storage dataframe for this model
    sst <- data.frame(Temp = rep(NA, duration),
                      Day = rep(NA, duration),
                      Year = rep(NA, duration))
    
    # determine the closest GPS coordinates in the model that corresponds with the GPS of the solitaries. 
    latitude <- ncvar_get(nc, "lat") #vector of latitudes within the nc file
    longitude <- ncvar_get(nc, "lon") #vector of longitudes within the nc file
    
    # different climate model use different GPS formats this next block of codes ensures R feeds the model GPS info in the correct format.
    if (c(nc[["var"]][["lon"]][["units"]], nc[["dim"]][["lon"]][["units"]]) == "degrees_east") {long.use = 360 - long
    } else {
      long.use = long
    }
    
    # now in the appropriate format the closest GPS location (to the one we care about) in the model can be determined. 
    if (is.matrix(latitude) == FALSE) {closest.lat <- which.min(abs(latitude - lat)) #vector location of closest latitude
    } else {
      closest.lat <- abs(latitude - lat)
      closest.lat <- as.numeric(which(closest.lat == min(closest.lat), arr.ind = TRUE)[1,2])
    }
    
    if (is.matrix(longitude) == FALSE) {closest.long <- which.min(abs(longitude - long.use)) #vector location of closest longitude
    } else {
      closest.long <- abs(longitude - long.use)
      closest.long <- as.numeric(which(closest.long == min(closest.long), arr.ind = TRUE)[1,1])
    }
    
    # Now for each day at the closest GPS point we will loop through the nc file and extract the daily projected temperatures
    # the if else block here is to deal with the fact that different climate models operate on different calender formats.
    if(nc[["dim"]][["time"]][["calendar"]] %in% c('noleap','365_day')) { #this works for the 365 day calenders
      
      for (i in 1:duration) {
        # first store the date information
        sst$Day[i] <- day
        sst$Year[i] <- year
        # next advance the date to the next day ready for the next loop
        day <- day + 1
        if (day > 365) {year = year + 1} # advances the year count on after 365 days
        if (day > 365) {day = 1} # after advancing the counts as nessecary this then resets the date at the end of the year
        # now the date is stored collect the temperature for that day
        kelvin <- ncvar_get(nc, "tos", start = c(closest.long,closest.lat, i) # this will start reading temperatures on the day in question
                            , count = c(1, 1, 1)) #and will just read one number
        # the model provides temperature in kelvin and so this needs converting and storing
        sst$Temp[i] <- kelvin - 273.15
      } #end of temperature extraction loop for selected model
      
    } else { #this works for the 360 day calenders
      
      for (i in 1:duration) {
        # first store the date information
        sst$Day[i] <- day
        sst$Year[i] <- year
        # next advance the date to the next day ready for the next loop
        day <- day + 1
        if (day > 360) {year = year + 1} # advances the year count on after 360 days
        if (day > 360) {day = 1} # after advancing the counts as nessecary this then resets the date at the end of the year
        # now the date is stored collect the temperature for that day
        kelvin <- ncvar_get(nc, "tos", start = c(closest.long,closest.lat, i) # this will start reading temperatures on the day in question
                            , count = c(1, 1, 1)) #and will just read one number
        # the model provides temperature in kelvin and so this needs converting and storing
        sst$Temp[i] <- kelvin - 273.15
      }
    }
    
    # store the model outputs in the projections list for the selected RCP
    projections[[x]] <- sst
    # and close the nc file ready for the next model
    nc_close(nc) #now the data has been extracted close the nc file to save space (or R crashes big time)
    
    try(rm(nc, latitude, longitude, sst, closest.lat, closest.long, day, year, duration)) #cleans ready for the next run.
  } #end of individual model loop
  
  # store all the extracted temperatures for each RCP pathway
  Pathways[[RCP]] <- projections
}) #end of full loop

# Now We have a metalist containing four separate lists.
# Each of these smaller lists represents the different RCP pathways and contains simulated future daily sea surface temperatures for the given GPS location from the different model outputs

# Now combine temperature chronologically and workout weekly averages across the different models for each RCP.
# Doing this requires extracting everything from the 'Pathway' list into dataframes representing each RCP scenario.
# first a vector of the RCP names (for tidying up later)
RCP.names <- c(8.5,6.0,4.5,2.6)
# and then an extraction loop
for (j in 1:4) { #this will run once for each RCP sublist
  # create a dataframe to start with using the first model output
  combine_temp <- data.frame(Pathways[[j]][[1]]) #extract the first model to begin the dataframe
  model_no <- length(Pathways[[j]]) #how many different models is there data for
  #now for each of these models
  for (d in 2:model_no) {
    combine_temp <- merge(combine_temp, data.frame(Pathways[[j]][[d]]), by = c("Year","Day"), all.x = TRUE, all.y = TRUE)
    # the inclusion of all.x and y means r puts them together even if there is no match.
  }
  # remove the days for which there is only one model output (means won't work)
  combine_temp <- combine_temp[which(combine_temp$Day <= 360),]
  
  # work out the mean projected daily temperatures.
  combine_temp$meanTemp <- rowMeans(combine_temp[,3:model_no], na.rm = TRUE)
  
  # now convert the daily mean temperatures in weekly means.
  # set up new datafile
  proj.temp <- data.frame(Year = numeric(),
                          MeanTemp = numeric())
  
  # now this loop will calculate the mean temperature every 7 days.
  day.select <- 1
  for(s in 1:(dim(combine_temp)[1]/7)) {
    proj.temp[s,1] <- combine_temp$Year[day.select]
    proj.temp[s,2] <- mean(combine_temp$meanTemp[day.select:day.select+6])
    day.select <- day.select + 7 #moves on a week
  }
  
  # and rename the dataframe so that it is easy to work out its associated RCP pathway
  assign(paste0("RCP", RCP.names[j]), proj.temp)
}

#########################################
# Estimate Mean Maximum Monthly Temperature for selected location
#########################################

# nc files documenting temperatures recorded in the selected location between 1985 and 1995 (Prior to any major bleaching) can be downloaded from the CoralTemp temperature data despository, and used to determine the Maximum Mean Monthly temperature for the region.
# the required files can be found at the NOAA page https://coralreefwatch.noaa.gov/product/5km/index_5km_sst.php#:~:text=The%20NOAA%20Coral%20Reef%20Watch,2%20to%2035%20%C2%B0C, downloaded on stored on a local drive.
# Source these downloaded files to extract required sst data
CW.files = list.files('Pathway to find CoralTemp files',pattern='*.nc',full.names = TRUE)
# each of these files represents one days worth of temperature data.
# define the time period needed
time.period <- c(1985:1995)
# subset the downloaded files based on the years actually needed.
try(CW.files <- CW.files[which(substr(CW.files, 55,58) %in% time.period)])

#setup a dataframe for storing extracted temperature variables
MMM.sst <- data.frame(Temp = rep(NA, length(CW.files)),
                      Month = rep(NA, length(CW.files)),
                      Year = rep(NA, length(CW.files)))

# Loop over each nc file to extract daily temperatures for the selected time frame.
try(for(x in 1:length(CW.files)) {
  
  #select correct nc file
  nc <- nc_open(CW.files[x])
  
  # Read the whole nc file and extract the important details
  latitude <- ncvar_get(nc, "lat") #vector of latitudes
  longitude <- ncvar_get(nc, "lon") #vector of longitudes
  sst_global <- ncvar_get(nc, "analysed_sst") #global matrix of sea surface temperatures for that particular day.
  closest.lat <- which.min(abs(latitude - lat)) #vector location of closest latitude 
  closest.long <- which.min(abs(longitude - long)) #and closest longitude
  date <- ncatt_get(nc,0,"time_coverage_start")$value # this is the year, month and day of the reading
  nc_close(nc) #now the data has been extracted close the nc file to save space (or R crashes big time)
  
  #extract the sst for the required location and store it 
  MMM.sst[x,1] <- sst_global[closest.long,closest.lat]
  MMM.sst[x,2] <- as.numeric(substr(date,5,6)) # to allow for the calculation of within year variation.
  MMM.sst[x,3] <- as.numeric(substr(date,1,4)) # this is to allow for subsetting the data by the years the population was being studied (1985 is year 1)
  
  try(rm(nc, latitude, longitude, sst_global, closest.lat, closest.long, date)) #cleans ready for the next run.
}) #end extraction loop

# The resulting dataframe contains a series of historical temperatures between 1985-1995 experienced in the selected location.
# Now calculate the average temperature for each month throughout the selected time period, and then use this to determine the highest monthly temperatures from each year
try(MMM.sst <- group_by(MMM.sst, Year, Month) %>%  
      summarise(MeanMonTemp = mean(Temp, na.rm = TRUE))) #calculate the variable for each month - keeping each year seperate
# Now work out the maximum from each year
try(MMM.sst <- group_by(MMM.sst, Year) %>%  
      summarise(MaxMonTemp = max(MeanMonTemp, na.rm = TRUE)) %>%
      summarise(MMMtemp = mean(MaxMonTemp, na.rm = TRUE))) #and finally take the mean of these maximum monthly temperatures. 
# MMM.sst$MMMtemp is now a singular value defining the maximum monthly temperature experienced in the selected location between 1985-1995

#########################################
# Determine daily coral hotspots and DHWs for a selected location
#########################################
# The MMM sea surface temperature obtained for the selected location can be used to determine the daily coral hotspots projected to occur within this region under scenarios expected by the differing RCP pathways
# From the daily coral hotspots we can estimate degree heating weeks.

# Store the different chronologically order RCP temperature projections in a single list
RCP.temp.store <- list(RCP8.5, RCP6, RCP4.5, RCP2.6)
# define the time period being used for the simulation
min.year <- XXXX
max.year <- XXXX

# the following loop will work through each daily projected temperature and determine if it exceeds the MMM temperature for the selected region and subsequently determine the associated coral hotspot magnitude.
# This is then used to estimate DHWs. 
for (i in 1:4) { #the loop works once for each RCP scenario
  data.select <- RCP.temp.store[[i]] #select each dataframe in turn
  #remove data from after max.year
  data.select <- data.select[which(data.select$Year <= max.year),]
  data.select <- data.select[which(is.finite(data.select$MeanTemp)),] # remove any inf estimates.
  data.select$hotspots <- NA # add a new column
  data.select$hotspots <- data.select$MeanTemp - MMM.sst$MMMtemp #this will for each day determine if the temperature exceeds the MMMtemp.
  # now use these hotspots to calculate the degree heating weeks.
  # DHWs are a running summation of the previous 3 months (12 weeks) of hotspot data (excluding negative values) 
  data.select$DHW <- NA #create a blank column for storing DHWs
  # work through each hotspot and sum the 84 values before it (12 weeks = ~84 days).
  for (z in 1:dim(data.select)[1]){
    # for the DHW analysis we are only concerned by values above the MMMtemp so negative values are removed prior to running anything
    if (data.select$hotspots[z] < 0) {data.select$hotspots[z] = 0}
    # now carry out a running summation of 84 hotspots to determine the DHWs for each day.
    # typically only coral hotspots greater than a 1 degree threshold are used here to estimate DHWs however for this analysis no threshold was applied and so all hotspots were used.
    if (z < 12) { #this prevents R breaking on the first 84 iterations where there isnt enough prior data to determine a cumalative sum.
      data.select$DHW[z] = NA
    } else {
      data.select$DHW[z] = sum(data.select$hotspots[(z-11):z]) 
    }
  } # end of DHW calculation
  
  # now remove past data
  data.select <- data.select[which(data.select$Year >= min.year),]
  
  # reassign the dataframe name so that it is easy to work out its associated pathway
  assign(paste0("DHW.RCP", RCP.names[i]), data.select)
  
} #end of complete loop

####################### 
# Use DHWs to estimate the frequncy of recurrent bleaching under each RCP simulation.
#######################

# create loop list
RCP.store.2 <- list(DHW.RCP2.6, DHW.RCP4.5, DHW.RCP6, DHW.RCP8.5)

# create storage data.frame
RCP_bleaching_trends <- data.frame(Year = c(min.year:max.year),
                                   bleaching2.6 = rep(NA, length.out = (max.year-min.year)+1), # equal to the number of years being simulated over.
                                   bleaching4.5 = rep(NA, length.out = (max.year-min.year)+1),
                                   bleaching6.0 = rep(NA, length.out = (max.year-min.year)+1),
                                   bleaching8.5 = rep(NA, length.out = (max.year-min.year)+1))

# determine the years in each simulation during which DHWs are expected to exceed a pre determined bleaching threshold at any point (typically 4 DHWs).  
for (RCP in 1:4) {
  RCP.use <- RCP.store.2[[RCP]]
  year.store <- unique(RCP.use$Year[which(RCP.use$DHW >= 4)]) # select all the years (removing duplicates) in which DHW exceeds 4 - even just once
  RCP_bleaching_trends[which(RCP_bleaching_trends$Year %in% year.store),RCP+1] = 1 #if 4DHW was reached bleaching happened
  RCP_bleaching_trends[is.na(RCP_bleaching_trends[,RCP+1]),RCP+1] = 0 # in all other years bleaching didn't happen
}

# Now these annual bleaching regimes can be used to determine the temporal bleaching probabilities under the different RCP pathways.

# stack together the different annual bleaching liklihood regimes
bleaching_trends <- data.frame(Year = rep(min.year:max.year, by = 4), bleaching = stack(RCP_bleaching_trends[,c(2:5)]))
colnames(bleaching_trends) <- c("Year", "Bleaching", "RCP")
bleaching_trends$RCP <- as.factor(bleaching_trends$RCP) #convert to factor for regression analysis

# Use a fixed effect logistic regression to investigate the differences between the pathways - and to quantify trends in annual probabilities of thermal stress events .
bleaching.mod <- glm(Bleaching ~ Year * RCP, family = "binomial", data = bleaching_trends)
summary(bleaching.mod) # very significant increase in bleaching with year and the trend differs significantly between RCP scenarios.
# the following code defines wether there is a significant difference between year and RCP interactions 
anova(bleaching.mod, test = "Chisq")
# now determine which of the interaction terms is driving this significance.
lsmeans(bleaching.mod,pairwise ~ RCP, adjust = "Tukey")$contrasts
# How well does the model fit
pR2(bleaching.mod)[6] #pseudo R2

# predict the trends for plotting
bleaching <- ggpredict(bleaching.mod, terms = c("Year[all]", "RCP"))
bleaching <- split(bleaching, bleaching$group) #seperate the different RCP trends
# and plot them.
plot(Bleaching~Year, data = bleaching_trends, typ = "n",
     xaxs = "i",
     yaxs = "i",
     xlab=expression("Years"), 
     ylab=expression("Probability of bleaching"),
     ylim = c(0,1))
# add predicted regression lines
lines(bleaching$bleaching2.6$predicted ~ bleaching$bleaching2.6$x, col = "Blue", lwd = 2)
lines(bleaching$bleaching4.5$predicted ~ bleaching$bleaching4.5$x, col = "Gray", lwd = 2)
lines(bleaching$bleaching6.0$predicted ~ bleaching$bleaching6.0$x, col = "Black", lwd = 2)
lines(bleaching$bleaching8.5$predicted ~ bleaching$bleaching8.5$x, col = "Red", lwd = 2)
legend(2020,0.9, c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), fill = c("Blue", "Gray", "Black", "Red"), bg = "White", bty = "l")

##############################
# Create a dataframe of projected bleaching probabilities for each RCP pathway based on the binomial model above.
RCP_bleaching <- data.frame(Years = rep(seq(2018,2100,1),times = 4),
                            RCP = rep(c("bleaching8.5","bleaching6.0", "bleaching4.5","bleaching2.6"), each = 83),
                            bleaching = rep(NA, length(332)))

# fill with the predictions
RCP_bleaching$bleaching <- predict(bleaching.mod, list(Year = RCP_bleaching$Year, RCP = RCP_bleaching$RCP), type = "response")
RCP_bleaching <- as.data.frame(split(RCP_bleaching, RCP_bleaching$RCP))
RCP_bleaching <- RCP_bleaching[,c(1,3,6,9,12)]
colnames(RCP_bleaching) <- c("Years","RCP2.6","RCP4.5","RCP6.0", "RCP8.5") #Just a little tidy up!

#####
#save the file
setwd("Directory/pathway")
write.csv(RCP_bleaching, file = "Bleaching regimes.csv") 

############################################### End of code ###################

################################################################################### End of code ############################################################################