# --- Script Setup ---
# Set the working directory. This line specifies the default location where R
# will look for input files and save output files.
#setwd("C:/Users/Leif Rauhoeft/Desktop/Desktop/bni/Publication/Process-based/Model 5.0/minimal example")
directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(directory)
# Record the current system time at the start of the script execution.
# This allows for measuring the total run time of the simulation.
start.time <- Sys.time()

# Load necessary R libraries.
library(zoo)
library(readxl)
library(geosphere)
library(ggplot2)
library(scales)         
require(lubridate)
library(zoo)
library(raster)
library(tidyverse)
library(readxl)
library(tidyr)
library(geosphere)

# --- Model Parameters ---
# These are key biological and environmental constants that define the
# behavior and interactions within the mosquito population model.

nP <- 80       #   Number of eggs laid by parous females
nN=141         #   Number of eggs laid by nulliparous females
KL <- 8*10^8   #   Standard environment carrying capacity for larvae (larvae ha???1)	
KP <- 10^7     #   Standard environment carrying capacity for pupae (pupae ha???1)
muem <- 0.1    #   Mortality rate during adult emergence (day???1)	
gam <- 0.003   #   Impact of density dependence
mur <- 0.08    #   Adult mortality rate related to seeking behavior (day???1)	
TAg <-9.8      #   Minimal temperature needed for egg maturation (?C)	
TDDAg <-64.4   #   Total number of degree-days necessary for egg maturation (?C)	
phi4 =0.143    #   Transition from emerging to host seeking
phi5=0.885     #   Transition from host seeking to resting
phi7 =0.2      #   Egg depositioning rate 
MuOv<- 0.02    #   Mortality during overwintering
MuE <- 0.0262  #   Egg Mortality
n0=5           #   Number of females each male can mate with
alpha=0.5      #   Proportion of Eggs which are female


# parameters for logistic function
asym<-0.5
xmid<-10
scal<--0.1



# # --- Geographic Locations and Time Series Data Setup ---
# 
# # Define the geographic locations (cities) for which the mosquito model will be run.
# # Each row contains the location name, its longitude, and its latitude.
locs_all<-as.data.frame(rbind(
  c("Hamburg", "9.993682", "53.551086"),
  c("Berlin", "13.4050", "52.5200"),
  c("Freiburg", "7.8522", "47.9990"),
  c("München", "11.5820", "48.1351")
))
# 
 locations<-locs_all

# Assign column names to the `locs_all` data frame for clarity.
names(locs_all)<- c("Location","lon","lat")

# Create a copy of the `locs_all` data frame for further manipulation.
locations<-locs_all

# Convert the `lat` and `lon` columns from character strings to numeric values.

 locations$lat<-as.numeric(locations$lat)
 locations$lon<-as.numeric(locations$lon)

# Define the date range for the simulation. The model will run daily from
# January 1, 2015, to December 31, 2023.
 dates <- seq.Date(as.Date("2015-01-01"), as.Date("2023-12-31"), by = "day")
 
# Mean daily temperature and daily rainfall data (European re-analysis and
# observations for monitoring, E-OBS, v29.0e) were obtained from the ECA&D
# project (http://www.ecad.eu, Cornes et al. 2018) and extracted for each
# sampling site for 2015 - 2023.
 

tg_Year<-readRDS("tg_Year.rds")
rr_Year<-readRDS("rr_Year.rds")
# Prepare latitude coordinates for the `daylength_function`.
# This loop creates a list where each element is a vector of repeated latitude values
# for a given location, matching the total number of days in the simulation period.
lat_coordinates<-list()
for(i in 1:nrow(locs_all)) {
  # Repeat the coordinates
  lat_coordinates[[i]] <- rep(locations$lat[i], each = length(seq_along(tg_Year[[1]])))
}
lat_coordinates <- lapply(lat_coordinates, function(vec) {
  matrix(vec, ncol = 1)})

lon_coordinates<-list()
for(i in 1:nrow(locs_all)) {
  # Repeat the coordinates
  lon_coordinates[[i]] <- rep(locations$lon[i], each = length(seq_along(tg_Year[[1]])))
}
lon_coordinates <- lapply(lon_coordinates, function(vec) {
  matrix(vec, ncol = 1)})

dates3<-list()

for (i in 1:nrow(locs_all)) {
  dates3[[i]] <- rep(dates, each = 1) 
}
dates3 <- lapply(dates3, function(vec) {
  matrix(vec, ncol = 1)})


daylength_function <- function(lat, dates) {
  daylength(lat, dates) 
}


daylength_results <- Map(daylength_function, lat_coordinates, dates3)


#transition rate Egg to larvae (temperature dependent) 
phi1 <-function(x){ asym/(1+exp(-(scal)*(xmid-x)))}#
#empty list for egg to larvae transition
phi1_res<-list()

# transition rate from Resting to ovipositioning (temperature dependent)
phi6 <-function(x) {ifelse(x>TAg,(x-TAg)/TDDAg,0)}
#empty list for Resting to ovipositioning transition rate
phi6_res<-list()

# Larvae mortality rate (temperature dependent)
MuL<-   function (x){ifelse(x>20, 0.86/(1+exp(-(-0.2)*(45-x))),
                            0.86/(1+exp((-5-x)/-5)))
}
# Adult mortality rate (temperature dependent)
MuAm<-   function (x){ifelse(x>20, 1/(1+exp(0.2*(40-x))),
                             1/(1+exp((-x)/-5)))
}
#empty list for larvae mortality
MuL_res<-list()

phi8<-function(x,y){ ifelse(y<15,1/(1+exp((12.5-x)/-0.5)),0)}
phi8_res<-list()

phi10<-function(x,y){ ifelse(y>6 ,exp((x - 11) / 0.5) / (1 + exp((x - 11) / 0.5)),0)} 
phi10_res<-list()
# Adult mortality rate
MuA <-MuL 
# empty list for Adult mosquito mortality rates    
muA_res<-list()
MuAm_res<-list()

# Renaming rain data
NIEDERSCHLAGSHOEHE<-rr_Year

#function calculating a rolling mean of 2 weeks
NIEDERSCHLAGSHOEHE2<-function(x)rollapply(abs(x), 14, sum)

# applying function to rain data
NIEDERSCHLAGSHOEHE2_res<-lapply(NIEDERSCHLAGSHOEHE,NIEDERSCHLAGSHOEHE2)

# normalising the rolling mean
newvalue<- function(x)(x-min(x, na.rm=T))/
  (max(x, na.rm=T)-min(x, na.rm=T))
newvalue_res<-lapply(NIEDERSCHLAGSHOEHE2_res, newvalue)

# Calculating environmental caryying capacity of larvae and pupae
Klt = function(x) KL * (x+1)
Klt_res= lapply(newvalue_res ,Klt)
KPt = function(x) KP * (x+1)
KPt_res= lapply(newvalue_res ,KPt)

# pre calculating all daily mortality and transition rates for each sampling site
for(i in 1:length(tg_Year)){
  
  phi1_res[[i]]<-sapply(tg_Year[i], phi1)
  phi6_res[[i]]<-sapply(tg_Year[i], phi6)
  MuL_res[[i]]<-sapply(tg_Year[i], MuL)
  muA_res[[i]]<-sapply(tg_Year[i], MuA)
  phi8_res[[i]] <- mapply(phi8, daylength_results[i], tg_Year[i])
  phi10_res[[i]]<-mapply(phi10, daylength_results[i], tg_Year[i])
  MuAm_res[[i]]<-sapply(tg_Year[i], MuAm)
  
}

# Non− diapause maintenance rate
# calculating phi9 to be the complement of phi8

phi9<- function(x){
  return(1-x)
}

phi9_res<-lapply(phi8_res,function(vec){
  sapply(vec,phi9)
}
)

#transition rate pupa to emerging	same as Egg to larvae
phi3_res <- phi1_res

# transition from larvae to pupae being one quarter of pupae to emerging
phi2<- function(x){
  return(x/4)
}

phi2_res<-lapply(phi3_res,function(vec){
  sapply(vec,phi2)
}
)


# Pupa mortality same as larvae mortality
MuP_res <-MuL_res 


# make empty lists for all parameters
empty_site <- matrix(rep(0,length(tg_Year[[1]])))
mylist <- rep(list(empty_site),length(tg_Year))


E_f<-E_m<-L_m<-L_f<-E_f<-P_f<-P_m<-P_D<-M<-E_F<-E_D<-F_mD<-HUN<-RUN<-OUN<-HUP<-RUP<-OUP<- density<-density_debug<-max_matings<-demand_EF<-demand_ED<-matings_EF<-matings_ED<-total_demand<-Location<-mylist

# defining starting population of 10^6 Eggs=
E_f <- lapply(E_f, function(x) {
  x[1] <- (10^6)*alpha
  return(x)
})
E_m <- lapply(E_m, function(x) {
  x[1] <- (10^6)*(1-alpha)
  return(x)
})
# defining starting pupulation for emerging femlaes and diapausing females
# since these are dependent on one another
E_F <- lapply(E_F, function(x) {
  x[1:10] <- 1000
  return(x)
})
E_D <- lapply(E_F, function(x) {
  x[1:10] <- 1000
  return(x)
})
p<-7
dates2<-list()
rain<-list()
temp<-list()


# morality rate of adult females is same as larvae morality rate
muA_res<-MuL_res


# Main loop for population dynamics by site and day
for (j in seq_along(tg_Year)) {
  for (i in 2:length(seq_along(tg_Year[[1]]))) {
    
    Location[[j]][i - 1] <- j
    
    # Egg female update
    E_f[[j]][i] <- max(
      E_f[[j]][i - 1] +
        (phi7 * alpha) * (nP * OUP[[j]][i - 1] + nN * OUN[[j]][i - 1]) -
        ((phi1_res[[j]][i - 1] + MuE) * E_f[[j]][i - 1]),
      0
    )
    
    # Egg male update
    E_m[[j]][i] <- max(
      E_m[[j]][i - 1] +
        (phi7 * (1 - alpha)) * (nP * OUP[[j]][i - 1] + nN * OUN[[j]][i - 1]) -
        (phi1_res[[j]][i - 1] + MuE) * E_m[[j]][i - 1],
      0
    )
    
    # Larvae female update
    L_f[[j]][i] <- max(
      L_f[[j]][i - 1] +
        (phi1_res[[j]][i - 1] * E_f[[j]][i - 1]) -
        (phi2_res[[j]][i - 1] + MuL_res[[j]][i - 1] *
           (1 + (L_f[[j]][i - 1] + L_m[[j]][i - 1]) / Klt_res[[j]][i - 1])) * L_f[[j]][i - 1],
      0
    )
    
    # Larvae male update (similar to female)
    L_m[[j]][i] <- max(
      L_m[[j]][i - 1] +
        (phi1_res[[j]][i - 1] * E_m[[j]][i - 1]) -
        (phi2_res[[j]][i - 1] + MuL_res[[j]][i - 1] *
           (1 + (L_f[[j]][i - 1] + L_m[[j]][i - 1]) / Klt_res[[j]][i - 1])) * L_m[[j]][i - 1],
      0
    )
    
    # Pupa female update
    P_f[[j]][i] <- max(
      P_f[[j]][i - 1] +
        (phi2_res[[j]][i - 1] * phi9_res[[j]][i - 1]) * L_f[[j]][i - 1] -
        (phi3_res[[j]][i - 1] + MuP_res[[j]][i - 1]) * P_f[[j]][i - 1],
      0
    )
    
    # Pupa male update
    P_m[[j]][i] <- max(
      P_m[[j]][i - 1] +
        phi2_res[[j]][i - 1] * L_m[[j]][i - 1] -
        (phi3_res[[j]][i - 1] + MuP_res[[j]][i - 1]) * P_m[[j]][i - 1],
      0
    )
    
    # Diapause induced pupae female update
    P_D[[j]][i] <- max(
      P_D[[j]][i - 1] +
        (phi2_res[[j]][i - 1] * phi8_res[[j]][i - 1]) * L_f[[j]][i - 1] -
        (MuP_res[[j]][i - 1] + phi3_res[[j]][i - 1]) * P_D[[j]][i - 1],
      0
    )
    
    # Adult males
    M[[j]][i] <- max(
      M[[j]][i - 1] +
        phi3_res[[j]][i - 1] * P_m[[j]][i - 1] * (1 - muem) *
        exp(-gam * (1 + (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1])) -
        MuAm_res[[j]][i - 1] * M[[j]][i - 1],
      0
    )
    
    # Calculate max matings and demands for matings
    max_matings[[j]][i] <- n0 * M[[j]][i - 1]
    
    total_demand[[j]][i] <- E_F[[j]][i - 1] + E_D[[j]][i - 1]
    if (is.na(total_demand[[j]][i])) total_demand[[j]][i] <- 0
    
    if (total_demand[[j]][i] > 0) {
      matings_EF[[j]][i] <- min(
        E_F[[j]][i - 1],
        (E_F[[j]][i - 1] / total_demand[[j]][i]) * max_matings[[j]][i - 1]
      )
      matings_ED[[j]][i] <- min(
        E_D[[j]][i - 1],
        (E_D[[j]][i - 1] / total_demand[[j]][i]) * max_matings[[j]][i - 1]
      )
    } else {
      matings_EF[[j]][i] <- 0
      matings_ED[[j]][i] <- 0
    }
    
    # Emerging females update
    E_F[[j]][i] <- max(
      E_F[[j]][i - 1] +
        phi3_res[[j]][i - 1] * P_f[[j]][i - 1] *(1-muem)*
        exp(-gam* (1 + (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1])) -
        (phi4 * matings_EF[[j]][i - 1]) -
        muA_res[[j]][i - 1] * E_F[[j]][i - 1],
      0
    )
    
    # Density calculations (for debugging or feedback)
    density[[j]][i] <- (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1]
    density_debug[[j]][i] <- exp(-gam * (1 + (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) /
                                           KPt_res[[j]][i - 1]))
    
    ##########################################################
    # Diapause conditioned emerging females
    E_D[[j]][i] <- max(
      E_D[[j]][i - 1] +
        phi3_res[[j]][i - 1] * P_D[[j]][i - 1] *(1-muem)*
        exp(-gam* (1 + (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1])) -
        phi4 * matings_ED[[j]][i - 1] -
        muA_res[[j]][i - 1] * E_D[[j]][i - 1],
      0
    )
    
    # Diapausing females update
    F_mD[[j]][i] <- max(
      F_mD[[j]][i - 1] +
        matings_ED[[j]][i - 1] -
        (phi10_res[[j]][i - 1] * phi4 + MuOv) * F_mD[[j]][i - 1],
      0
    )
    
    # Host-unfed nulliparous females update
    HUN[[j]][i] <- max(
      HUN[[j]][i - 1] +
        phi4 * matings_EF[[j]][i - 1] +
        phi10_res[[j]][i - 1] * phi4 * F_mD[[j]][i - 1] -
        (phi5 + mur + muA_res[[j]][i - 1]) * HUN[[j]][i - 1],
      0
    )
    
    # Resting unfed nulliparous females update
    RUN[[j]][i] <- max(
      RUN[[j]][i - 1] +
        phi5 * HUN[[j]][i - 1] -
        (phi6_res[[j]][i - 1] + muA_res[[j]][i - 1]) * RUN[[j]][i - 1],
      0
    )
    
    # Ovipositing unfed nulliparous females update
    OUN[[j]][i] <- max(
      OUN[[j]][i - 1] +
        phi6_res[[j]][i - 1] * RUN[[j]][i - 1] -
        (phi7 + muA_res[[j]][i - 1]) * OUN[[j]][i - 1],
      0
    )
    
    # Host-seeking parous females update
    HUP[[j]][i] <- max(
      HUP[[j]][i - 1] +
        phi7 * (OUP[[j]][i - 1] + OUN[[j]][i - 1]) -
        (phi5 + mur + muA_res[[j]][i - 1]) * HUP[[j]][i - 1],
      0
    )
    
    # Resting parous females update
    RUP[[j]][i] <- max(
      RUP[[j]][i - 1] +
        phi5 * HUP[[j]][i - 1] -
        (phi6_res[[j]][i - 1] + muA_res[[j]][i - 1]) * RUP[[j]][i - 1],
      0
    )
    
    # Ovipositing parous females update
    OUP[[j]][i] <- max(
      OUP[[j]][i - 1] +
        phi6_res[[j]][i - 1] * RUP[[j]][i - 1] -
        (phi7 + muA_res[[j]][i - 1]) * OUP[[j]][i - 1],
      0
    )
  }
  
  dates2[j] <- list(as.character(dates))
}

# --- Data Aggregation and Export ---

# Convert the `tg_Year` (temperature) data lists into single-column matrices.
# This prepares them for consistent binding into the final data frame.
tg_Year2 <- lapply(tg_Year, function(vec) {
  matrix(vec, ncol = 1)
})
# Remove any names from the elements of `tg_Year2`. This is often done to ensure
# `do.call(rbind.data.frame, ...)` works without unexpected column naming issues.
names(tg_Year2)<-NULL

# Convert the `rr_Year` (rainfall) data lists into single-column matrices,
# similar to `tg_Year2`.
rr_Year2 <- lapply(rr_Year, function(vec) {
  matrix(vec, ncol = 1)
})
# Remove any names from the elements of `rr_Year2`.
names(rr_Year2)<-NULL

# Combine all generated compartment time series data, dates, and environmental data
# into one large list called `all`. Each element of `all` is itself a list of time series
# (one time series for each location for that specific variable).
all<- list(Location,dates2,E_f,E_m,L_m,L_f,P_f,P_m,P_D,M,E_F,E_D,F_mD, HUN,RUN,OUN,HUP,RUP,OUP,tg_Year2, rr_Year2, lon_coordinates,lat_coordinates)



# Combine data for each location into a list of temporary data structures.
# This loop iterates through the elements of one of the compartment lists (e.g., `E_f`),
# treating each element's index `k` as a specific location. For each location, it then
# extracts the `k`-th element (which is the time series for that location) from *every*
# list within the `all` list. This effectively groups all data pertinent to one location.
all_list<-list()
for (k in seq_along(E_f)){
  all2<- sapply(all, "[", k) # Extracts the k-th time series (for location k) for each compartment.
  all_list<- append(all_list,list(all2)) # Appends this collection of time series (for location k) to `all_list`.
}

# Use `do.call(rbind.data.frame, ...)` to combine all the elements (each representing a location's
# full simulation output) from `all_list` into a single, large data frame `a`.
# `rbind.data.frame` stacks these data structures row by row.
a<- data.frame(do.call(rbind.data.frame, all_list))

# Assign clear and descriptive column names to the newly created data frame `a`.
names(a)<-c("Location","Date","E_f","E_m","L_m","L_f","P_f","P_m","P_D","M","E_F","E_D","F_mD", "HUN","RUN","OUN","HUP","RUP","OUP","temp", "rain","lon","lat")

# Subset the data frame to remove the very last row if its date is "2023-12-31".
# This is a common practice because the population values on the last day might represent
# the state *before* the final day's transitions, or might be incomplete.
a<- subset(a,!a$Date=="2023-12-31")

# Assign human-readable city names (from the `locations` data frame) to the 'Location'
# column in the final data frame `a`. This column is converted to a factor variable.
cities<-locations$Location
a$Location<-as.factor(a$Location)
levels(a$Location)<-c(cities)

# Export the final simulated data to a CSV file named "model_output_pub2.csv".
# `row.names=FALSE` prevents R from writing an additional column of row numbers into the CSV.
write.csv(a,"model_output_minimal_example.csv")


# --- Script End and Timing ---
# Record the current system time at the end of the script execution.
end.time <- Sys.time()
# Calculate the total time taken for the script to run, rounded to 2 decimal places.
time.taken <- round(end.time - start.time,2)
time.taken


### some quick plots
for (o in 1:4) {
  ##–– set a 2-panel layout: 2 rows, 1 column ––
  par(mfrow = c(2, 1),          # two plots, one below the other
      mar   = c(4, 4, 3, 2))    # bottom, left, top, right margins
  
  ## –– 1) Host-seeking females (HUN + HUP) ––
  plot(
    HUN[[o]] + HUP[[o]],
    type = "l", col = "black", lwd = 4,
    main = paste("Site", o, " – Host seeking"),
    xlab = "Day of year", ylab = "Normalised value"
  )
  abline(h = 0.5, col = "grey", lty = 2)
  abline(v = seq(   1, 3285, by = 365), col = "red"  , lty = 5) # 1 Jan
  abline(v = seq( 152, 3285, by = 365), col = "green", lty = 5) # 1 Jun
  abline(v = seq( 213, 3285, by = 365), col = "blue" , lty = 5) # 1 Aug
  
  ## –– 2) Temperature (tg_Year2) ––
  plot(
    tg_Year2[[o]],
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of year", ylab = "°C",
    main = paste("Site", o, " – Daily temperature")
  )
  abline(h = 0, col = "grey80")
  abline(v = seq(   1, 3285, by = 365), col = "red"  , lty = 5) # 1 Jan
  abline(v = seq( 152, 3285, by = 365), col = "green", lty = 5) # 1 Jun
  abline(v = seq( 213, 3285, by = 365), col = "blue" , lty = 5)
  
  ##–– reset layout for next loop iteration ––
  par(mfrow = c(1, 1))
}

par(mfrow = c(1, 1))  # Reset layout

par(mfrow = c(2, 2))  # 2 rows, 2 columns of plots

for (o in 1:4) {
  plot(
    ((E_f[[o]] )) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Eggs (black=E_f, red=E_m)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  lines( E_m[[o]], col="red")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1))  # Reset layout

par(mfrow = c(2, 2))  # 2 rows, 2 columns of plots

for (o in 1:4) {
  plot(
    ((L_f[[o]] )) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Larvae (black=L_f, red=L_m)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  lines( L_m[[o]], type="l",col="red")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1

}

par(mfrow = c(1, 1))  # Reset layout

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    ((P_f[[o]] )) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Pupae (red=P_D, black=P_f, blue=P_m) "),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  lines(P_D[[o]], col="red",lwd=3)
  lines(P_m[[o]], col="blue")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1))  # Reset layout



par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (M[[o]] ) ,
    col = "red",
    type = "l",
    main = paste("Plot for o =", o,"Emerging (black =E_F, red=M, blue= E_D, green=F_mD)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  lines( E_F[[o]], col="black",lwd=3)
  lines( E_D[[o]], col="blue",lwd=3)
  lines( F_mD[[o]], col="green",lwd=3)
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1

}

par(mfrow = c(1, 1))  # Reset layout

for (o in 1:4) {
  ##–– set a 2-panel layout: 2 rows, 1 column ––
  par(mfrow = c(2, 1),          # two plots, one below the other
      mar   = c(4, 4, 3, 2))    # bottom, left, top, right margins
  
  ## –– 1) Host-seeking females (HUN + HUP) ––
  plot(
    HUN[[o]] + HUP[[o]],
    type = "l", col = "black", lwd = 4,
    main = paste("Site", o, " – Host seeking"),
    xlab = "Day of year", ylab = "Normalised value"
  )
  abline(h = 0.5, col = "grey", lty = 2)
  #lines(E_f[[o]] , col = "red")
  abline(v = seq(   1, 3285, by = 365), col = "red"  , lty = 5) # 1 Jan
  abline(v = seq( 152, 3285, by = 365), col = "green", lty = 5) # 1 Jun
  abline(v = seq( 213, 3285, by = 365), col = "blue" , lty = 5) # 1 Aug
  
  ## –– 2) Temperature (tg_Year2) ––
  plot(
    tg_Year2[[o]],
    type = "l", col = "steelblue", lwd = 2,
    xlab = "Day of year", ylab = "°C",
    main = paste("Site", o, " – Daily temperature")
  )
  abline(h = 0, col = "grey80")
  abline(v = seq(   1, 3285, by = 365), col = "red"  , lty = 5) # 1 Jan
  abline(v = seq( 152, 3285, by = 365), col = "green", lty = 5) # 1 Jun
  abline(v = seq( 213, 3285, by = 365), col = "blue" , lty = 5)
  
  ##–– reset layout for next loop iteration ––
  par(mfrow = c(1, 1))
}

par(mfrow = c(1, 1))  # Reset layout


par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (RUN[[o]]+RUP[[o]] ) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Resting"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1

}

par(mfrow = c(1, 1)) 

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (RUN[[o]]+RUP[[o]] ) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Resting"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1)) 

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (OUN[[o]]+OUP[[o]] ) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Ovipositioning"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 4
  )
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1)) 

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (phi8_res[[o]]) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Transition (black = phi8, red=phi10)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    ylim=c(0,1)
    , lwd = 1
  )
  lines(phi10_res[[o]], type="l", col= "red")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1)) 
par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (phi6_res[[o]]) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Transition (black = phi6, red=phi1,phi3, phi2*4)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    ylim=c(0,1)
    , lwd = 1
  )
  lines(phi1_res[[o]], type="l", col= "red")
  #lines(phi8_res[[o]], type="l", col="green")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1)) 

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (MuL_res[[o]]) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o,"Mortalities (black = MuL,MuA,MuP red=MuAm)"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    ylim=c(0,1)
    , lwd = 1
  )
  lines(MuAm_res[[o]], type="l", col= "red")
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (density_debug[[o]]) ,
    col = "black",
    type = "l",
    main = paste("density_debug[[j]][i] <- exp(-muem * (1 + (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1]))"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    ylim=c(0,1)
    , lwd = 1
  )
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
  
}

par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
for (o in 1:4) {
  plot(
    (density[[o]]) ,
    col = "black",
    type = "l",
    main = paste("Plot for o =", o," (density[[j]][i] <- (P_m[[j]][i - 1] + P_f[[j]][i - 1] + P_D[[j]][i - 1]) / KPt_res[[j]][i - 1])"),
    ylab = "Normalized Value",
    xlab = "Time (days)",
    , lwd = 1
  )
  abline(h = 0.5, col = "gray", lty = 2)
  abline(v = seq(213, 3285, by = 365), col = "blue", lty = 5)#august 1
  abline(v = seq(152, 3285, by = 365), col = "green", lty = 5)#june 1
  abline(v = seq(1, 3285, by = 365), col = "red", lty = 5)#jan 1
}

par(mfrow = c(1, 1))