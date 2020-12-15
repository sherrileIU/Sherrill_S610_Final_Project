# STAT-S 610
# FINAL PROJECT
# 2020-12-08

# TO RUN: testthat::test_dir('.')

# --- setup --- #

# load package functions and functions we want to test
library(readr)
library(lubridate)
library(testthat)

import::here(extract_data, find_EQ_times, fit_GPS_data, format_GPS_data, 
					build_G_matrix, compute_mhat, compute_fits, build_co_matrix, 
					build_log_matrix, build_BM_matrix, build_D_matrix,
					.from = 'sherrill_functions.R')

context("Check GPS data fitting functions")

################# data setup #####################################
# variables needed for running the functions
end_time = 2004
tR = 0.1
alpha = 20


# comparison dataset info
check_eq_times = c(1999.789)
check_data = read.csv('check_data.csv')

#compute these so they can be reused in checks below
check_GPS_data = format_GPS_data(check_data, check_eq_times, end_time)
check_G_matrix = build_G_matrix(check_GPS_data$times, check_GPS_data$eq_times, tR, alpha)


################ tests ##########################################

# --- test that data are downloaded and saved correctly --- #
test_that("GPS data are downloaded correctly", {
	# download data for test station SHOS
	test_data = extract_data('SHOS')
	
	# check that returned data has 3 columns
	expect_equal(dim(test_data), dim(check_data))
	
	# load the saved data file and check that its dimensions are the same as the returned data
	loaded_data = read.csv('data_SHOS.csv')
	expect_equal(dim(loaded_data), dim(test_data))
	
})



# --- test that number of earthquakes and earthquake times are correct --- #
test_that("earthquake times are correct", {
	# download data for test station SHOS
	test_eqs = find_EQ_times('SHOS')
	test_eqs = test_eqs[test_eqs <= end_time]
	
	# check that there are the same number of earthquakes
	expect_equal(length(test_eqs), length(check_eq_times))
	
	# the earthquake times are equal to the nearest tenth of a year
	expect_equal(round(test_eqs, digits = 3), round(check_eq_times, digits = 3))
		
})



# --- check that time shift and formatting are correct --- #
test_that("GPS data are formatted correctly", {
  # make example GPS data
  GPS_data = check_GPS_data
  start_time = GPS_data$start_time
  
  # check that the maximum time is not greater than the difference 
  # between end_time and start_time
  expect_true(max(GPS_data$times) <= end_time - start_time)
  
  # check that the time shifted times are >= 0
  expect_true(all(GPS_data$times >= 0))
  
  # check that formatted eq_times = eq_times - start_time
  expect_true(all(GPS_data$eq_times == check_eq_times - start_time))
  
  # check that formated times and displacements are the same length
  expect_equal(dim(GPS_data$disp)[1], length(GPS_data$times))
})


# --- checks that the G matrix is built correctly --- #
test_that("G matrix is built correctly", {
  # there are two earthquakes for the test data so there will be two columns corresponding 
  # to the linear trend, two columns for coseismic offset, and two columns for logarithmic postseismic decay
  
  GPS_data = check_GPS_data
  G = check_G_matrix
  L = length(GPS_data$times)
  
  # check dimensions of G are L + L rows (the D (smoothing) matrix adds on L rows) 
  # and L + 4 columns (4 added for the 2 linear, 1 coseismic, and 1 logarithmic columns)
  expect_equal(dim(G), c(L + L, L + 4))
 
   # check that the first column of the linear G matrix is equal to times
  expect_equal(G[1:L ,L+1], GPS_data$times)

  # check that the second column of the linear G matrix is all 1s
  expect_true(all(G[1:L ,L+2] == 1))
  
  # check that the coseismic columns are all 1 or 0
  expect_true(all(G[1:L, (L + 3)] == 1 | G[1:L, (L+ 3)] == 0))
  
  # check that all the values in the log columns are >= 0
  expect_true(all(G[1:L, (L + 4)] >= 0))
  
  # check that D has only 1s on the diagonal times alpha
  expect_equal(sum(G[(L + 1):(L + L), 1:L]), length(GPS_data$times)*alpha)
})



# --- checks that the final mhat and fit files are correct --- #
test_that("final products are correct", {
  # run full script to produce final files 
 fit_GPS_data('SHOS', end_time, tR, alpha)
 
 #load saved files
 fits = read.csv('fits_SHOS_0.1_20.csv')
 mhat = read.csv('mhat_SHOS_0.1_20.csv')
 
 # should have a row for every time and 
 # 11 columns (times as well as east and north displacements, dhat, trend, log, and brownian motion) 
 expect_equal(dim(fits), c(length(check_GPS_data$times), 11))
 
 # should have a row for every time corresponding to the Brownian motion,
 # two rows for the linear trend, one for the coseismic displacment, and one for logarithmic component
 # should have two columns, one for the east and one for the north
 expect_equal(dim(mhat), c((length(check_GPS_data$times) + 4), 2))
 })

