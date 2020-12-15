## final project
## plotting script
#library(ggplot2)

plot_GNSS_data = function(station_name, tR, alpha) {
	filename = gsub("X", station_name, 'fits_X_Y_Z.csv')
	filename = gsub("Y", tR, filename)
	filename = gsub("Z", alpha, filename)
	fit = read.csv(filename)

	eq_filename = gsub("X", station_name, 'eq_times_X.csv')
	eq_times = as.double(unlist(read.csv(eq_filename)))

	par(mfrow = c(1,2))
	plot(fit$times, fit$data.E - fit$east.trend, col = "blue", xlab = "time", ylab = "displacement")
	points(fit$times, fit$east.dhat - fit$east.trend, col = "red")
	lines(fit$times, fit$east.brown, col = "green")
	lines(fit$times, fit$east.log, col = "black")
	for(i in 1:length(eq_times)){
		lines(rep(eq_times[i], 2), c(-1,1), col = "black")
	}		
	
	plot(fit$times, fit$data.N - fit$north.trend, col = "blue", xlab = "time", ylab = "displacement")
	points(fit$times, fit$north.dhat - fit$north.trend, col = "red")
	lines(fit$times, fit$north.brown, col = "green")
	lines(fit$times, fit$north.log, col = "black")
	for(i in 1:length(eq_times)){
		lines(rep(eq_times[i], 2), c(-1, 1), col = "black")
	}
}	
	
	
	
