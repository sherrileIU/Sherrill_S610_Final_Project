## Elizabeth Sherrill
## STAT-S610
## Final Project Functions

## library(readr)
## library(lubridate)
## must have steps.txt file in working directory

## will produce the following files if write_file = 1:
## 1. raw data - data_stationname.csv
## 2. EQ times - eq_times_stationname.csv
## 3. final data and fits - fits_stationname_tR_alpha.csv
## 4. mhat - mhat_stationname_tR_alpha.csv

# --- data download --- #

# extract_data downloads GNSS time series data from the UNR 
# repository for the station name provided.

# INPUT: station name
# OUTPUT: raw GNSS data - times, east displacements, and north displacements
# diplacement units = meters

extract_data = function(station_name) {
	URL = 'http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/X.tenv3'
	data_path = gsub("X", station_name, URL)

	raw_data = read.table(data_path)

	last = dim(raw_data)[1]
	times = as.numeric(raw_data[2:last,3])
	east = as.numeric(raw_data[2:last,9])
	north = as.numeric(raw_data[2:last,11])
	datafile = data.frame(cbind(times, east, north))
	
	filename = gsub("X", station_name, 'data_X.csv')
	write.csv(datafile, filename, row.names = FALSE)
	
	return(datafile)
}

# find_EQ_times finds the times of previously identified earthquakes that 
# UNR identified to be large enough and close enough to produce a step 
# in a given GNSS station's data 
find_EQ_times = function(station_name) {
	steps = read.table('steps.txt')
	eq_times = rep(0, length = sum(steps[,1] == station_name))
	ind = which(steps[,1] == station_name)
	eq_times = steps[ind, 2]
	eq_times = unique(decimal_date(as.Date(eq_times, "%y%b%d")))
	
	filename = gsub("X", station_name, 'eq_times_X.csv')
	write.csv(eq_times, filename, row.names = FALSE)
	
	return(eq_times)
}

# --- data formatting --- #
# format_GPS_data will take the raw GNSS data, time shift the times and EQ times 
# based on the first time in the times column, cutoff the times and EQ times based on 
# the specified end time, then apply a running median filter to remove extreme outliers

format_GPS_data = function(data, eq_times, end_time) {
  #only want times & displacements before the designated end time
  eq_times = eq_times[eq_times < end_time]
  ind = which(data$times < end_time) 
  
  #apply a time shift so that the times and displacements start at zero
  start_time = data$times[1]
  eq_times = eq_times - start_time
  times = data$times[ind] - start_time
  
  #runmed used to remove extreme outliers due to instrument error
  E = runmed(data$east[ind], 15)
  N = runmed(data$north[ind], 15)
  
  #formats into a list saving displacements, times, eq times, and start time.
  GPS_data = list(disp = cbind(E, N), 
                  times = times, 
                  eq_times = eq_times, 
                  start_time = start_time)
  return(GPS_data)
}

# --- build G matrix --- #
# build_G_matrix is the overall function that calls
# build_BM_matrix, build_co_matrix, build_log_matrix, and build_D_matrix
# to build the larger G matrix for the inversion

build_BM_matrix = function(times) {
  # builds a lower triangular matrix of ones
  # each row represents a day of data
  # each column is the contribution of brownian motion to that day's motion
  Gbm = matrix(0, length(times), length(times))
  Gbm[lower.tri(Gbm, diag = TRUE)] = 1 
  return(Gbm)
}

build_co_matrix = function (times, eq_times) {
  # builds a matrix of ones, then replaces values in the ith row 
  # before the ith eq time with zeros
  Gco = matrix(1, nrow = length(times), ncol = length(eq_times))
  for(i in 1:length(eq_times)) {
    ind = times < eq_times[i]
    Gco[ind, i] = 0
  }
  return(Gco)
} 

build_log_matrix = function(times, eq_times, tR) {
  # builds a matrix of zeros with length times and ncol for number of eqs
  # if times are >= ith eq time, the Glog value is saved in the ith column
  Glog = matrix(0, nrow = length(times), ncol = length(eq_times))
  for(i in 1:length(eq_times)) {
    ind = times >= eq_times[i]
    Glog[ind, i] = log(1 + (times[ind] - eq_times[i])/tR)
  }
  return(Glog)
}

build_D_matrix = function(times, eq_times) {
  # The D matrix allows for smoothing of the Brownian motion component.
  n_pcol = length(eq_times)*2 + 2
  Deye = diag(length(times))
  Dzero = matrix(0, nrow = length(times), ncol = n_pcol)
  D = cbind(Deye, Dzero)
  return(D)
}

build_G_matrix = function(times, eq_times, tR, alpha) {
  # Gbm is the matrix corresponding to Brownian motion
  # The next two columns are linear model, one column of times, one column of ones
  # Then there are columns for the number of EQs for both Gco and Glog
  
  Gbm = build_BM_matrix(times)
  Glin = cbind(times, rep(1, length(times)))
  Gco = build_co_matrix(times, eq_times)
  Glog = build_log_matrix(times, eq_times, tR)
  D = build_D_matrix(times, eq_times)
  G = cbind(Gbm,Glin, Gco, Glog)
  G = rbind(G, alpha*D) #alpha is smoothing parameter for BM
  return(G)
}


# --- compute the model and fit --- #
# Adds on zeros to the data vector to make it the correct length
# Solves for the model vector (mhat)

compute_mhat = function(G, data) {
  data = c(data, rep(0, length(data)))
  mhat = solve(crossprod(G), crossprod(G, data))
  return(mhat)
}


compute_fits = function(G, mhat, times, eq_times) {
  L = length(times)
  G = G[1:L,] #removes the smoothing section of G
  
  #solves first for the full modeled data vector
  dhat = G %*% mhat
  
  #then solves for Brownian motion component
  brown = G[,1:L] %*% mhat[1:L]
  
  #then solves for the linear trend only
  trend = G[,(L+1):(L+2)] %*% mhat[(L+1):(L+2)]
  
  #then solves for the logarithmic component of the fit
  n = length(eq_times)
  log = matrix(nrow = L, ncol = n )
  m = length(mhat) - n
  for(i in 1:n) {
    log[,i] = G[,i + m] * mhat[i + m]
  }
  log = apply(log, 1, sum)
  
  fits = data.frame(dhat, trend, log, brown)
  return(fits)
  }


fit_GPS_data = function(station_name, end_time, tR, alpha) {
  #step 1: extract data from online
  raw_data = extract_data(station_name)
  
  #step 2: extract earthquake times
  orig_eq_times = find_EQ_times(station_name)
  
  #step 3: format, time shift, and filter the data
  data = format_GPS_data(raw_data, orig_eq_times, end_time)
  
  #step 4: build a G matrix with Gbm, Glin, Gco, and Glog components
  G = build_G_matrix(data$times, data$eq_times, tR, alpha)
  
  #step 5: solve for mhat
  mhat = apply(data$disp, 2, compute_mhat, G = G)
  
  #step 6: solve for the modeled data vector, linear trend, and log fit
  fits = apply(mhat, 2, compute_fits, G = G, times = data$times, eq_times = data$eq_times)
  
  #step 7: save the data and the fits for plotting and analysis purposes
  final = data.frame(times = data$times + data$start_time, data = data$disp, east = fits$E, north = fits$N)
  
  filename = gsub("X", station_name, 'fits_X_Y_Z.csv')
  filename = gsub("Y", tR, filename)
  filename = gsub("Z", alpha, filename)
  write.csv(final, filename, row.names = FALSE)
  
  mfilename = gsub("X", station_name, 'mhat_X_Y_Z.csv')
  mfilename = gsub("Y", tR, mfilename)
  mfilename = gsub("Z", alpha, mfilename)
  write.csv(mhat, mfilename, row.names = FALSE)
  return(head(final))
}

