# tecan functions
library(ggplot2)
library(dplyr)

read_tecan <- function(file_path, measure.name = "OD"){
  # read_tecan takes in a path to a .csv file from the Tecan and reads it in
  # and converts it to long format, for easier use with ggplot
  # by default it names the measurement variable "OD"
  OD = read.table(file_path, sep = ",", as.is = TRUE, row.names = 1)
  OD = t(OD)
  OD = data.frame(OD)
  names(OD)[1:3] = c("cycle","time_s","temp_C")
  library(reshape2)
  OD = reshape2::melt(OD, id.vars = c("cycle","time_s","temp_C"),
                      variable.name = "well",value.name = measure.name)
  good_rows = !is.na(OD[,measure.name])
  OD = OD[good_rows,]
  return(OD)
}

apply_well_map <- function(results, map_path, new_name = "treatment"){
  # apply_well_map takes in a data.frame called results,
  # which should have a column called well, and also 
  # a path to a .csv containing a map of the treatments
  # e.g. [x   ecoli  ecoli  sal  x    methly ...]
  #      [x   sal    ecoli  x    x    methly ...]
  #      [... ...    ...    ...  ...  ...    ...]
  # it then matches the location in the map with the well name
  # and applies the treatment name in the right spot.
  # by default, the new column name containing the treatments is called
  # "treatment"
  treatments = read.csv(map_path,row.names = 1)
  # now the map is just a matrix with strings, with each different string
  # being a different treatment
  library(stringr)
  results[,new_name] = ""
  # basically we're just gonna travel through the map and assign that name
  # to a new column called treatment if the location matches the well name
  for (j in 1:nrow(treatments)){
    for (k in 1:ncol(treatments)){
      thisTreat = str_trim(treatments[j, k]) # this grabs the treatment at that loc and strips whitespace
      letter = LETTERS[j] # LETTERS is a constant vector of the alphabet letters, good for grabbing specific ones
      wellID = paste(letter,k,sep = "") # by putting the letter next to the col number, we get the well name
      results[OD$well == wellID,new_name] = thisTreat
    }
  }
  return(results)
}



get_growth_rate_per_hr = function(results, well_ID, response_var = "OD",
                                  window_size = 11, R2_cutoff = 0.9){
  library(dplyr)
  # get the data and log-transform the response variable
  curr_dat = results %>% filter(well == well_ID)
  x_dat = curr_dat$time_s / 3600
  y_dat = curr_dat[,response_var] + runif(length(x_dat), min = 1e-8, max = 2e-8) # to prevent multiple -Infs
  y_dat = log(y_dat - min(y_dat))
  bad_pos = which(y_dat == -Inf)
  y_dat[bad_pos] = mean(c(y_dat[bad_pos-1], y_dat[bad_pos+1]))
  
  slopes = rep(0, length(y_dat) - window_size)
  r2 = rep(0, length(y_dat) - window_size)
  for (k in 1:length(slopes)){
    x = x_dat[k:(k+window_size)]
    y = y_dat[k:(k+window_size)]
    m1 <- lm(y ~ x, na.action = na.omit)
    slopes[k] = coef(m1)[2]
    r2[k] = summary(m1)$r.squared
  }
  acceptable_slopes = which(r2 >= R2_cutoff)
  slopes = slopes[acceptable_slopes]
  growth_rate = max(slopes) # in per seconds
  cycle_of_max_slope = acceptable_slopes[which(slopes == growth_rate)]
  return(c(growth_rate = growth_rate,
           cycle = cycle_of_max_slope))
}

get_all_growth_rates_per_hr = function(results,  response_var, window_size, R2_cutoff,
                                       gr_name = "gr", cycle_name = "gr_cycle"){
  for (curr_well in unique(results$well)){
    gr_and_time = get_growth_rate_per_hr(results, curr_well, response_var, window_size, R2_cutoff)
    results[results$well == curr_well, gr_name] = gr_and_time[1]
    results[results$well == curr_well, cycle_name] = gr_and_time[2]
  }
  return(results)
}

generate_baranyi_prediction = function(results, r_name = "baranyi_r",
                                       lag_name = "baranyi_lag",
                                       ymax_name = "baranyi_ymax",
                                       y0_name = "barani_y0",
                                       new_name = "baranyi_pred"){
  for (curr_well in unique(results$well)){
    curr_dat = results %>%
      filter(well == curr_well)
    r = curr_dat[1, r_name]
    lag = curr_dat[1, lag_name]
    ymax = curr_dat[1, ymax_name]
    y0 = curr_dat[1, y0_name]
    
    x_dat = results$time_s[results$well == curr_well] / 3600
    y_pred = exp(baranyi(x_dat, r, lag, ymax, y0))
    results[results$well == curr_well, new_name] = y_pred
  }
  return(results)
}

get_all_baranyi_rates = function(results,  response_var,
                                 r_name = "baranyi_r", lag_name = "baranyi_lag",
                                 ymax_name = "baranyi_ymax", y0_name = "baranyi_y0",
                                 tries = 100){
  for (well_ID in unique(results$well)){
    gr_and_time = get_baranyi_rate(results, well_ID, response_var, tries)
    results[results$well == well_ID, r_name] = gr_and_time[1]
    results[results$well == well_ID, lag_name] = gr_and_time[2]
    results[results$well == well_ID, ymax_name] = gr_and_time[3]
    results[results$well == well_ID, y0_name] = gr_and_time[4]
  }
  return(results)
}


get_baranyi_rate = function(results, well_ID, response_var = "OD", tries = 100){
  curr_dat = results %>% filter(well == well_ID)
  curr_dat <- arrange(curr_dat, time_s)
  x = curr_dat$time_s / 3600
  y = log(curr_dat[,response_var])
  return(fit_baranyi(x, y, tries))
}



fit_baranyi = function(x, y, tries = 100){
  m1 <- NA
  
  #figure out start lag guess
  lag_guess = guess_half_max(x,y)
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ baranyi(x, r, lag, ymax, y0),
                        start = list(r = runif(1, 0.1, 0.4), 
                                     lag = lag_guess + sample(-20:20, 1), 
                                     ymax = max(y), 
                                     y0 = min(y))),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1)){
    return(c(NA, NA))
  }
  #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
  #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
  r = coef(m1)[1]
  lag = coef(m1)[2]
  ymax = coef(m1)[3]
  y0 = coef(m1)[4]
  return(c(r, lag, ymax, y0))
}

guess_half_max = function(x, y){
  # this function looks at a logistically-increasing time series
  # and guesses when the growth rate is at a maximum (which is also when
  # the trend is at half-max)
  
  # put in order
  y = y[sort(x, index.return = TRUE)$ix]
  x = sort(x)
  
  #find approximate time of max growth rate, using diffs, on smoothed y
  y2 = diff(y)
  y2 = zoo::rollmean(y2, 11, fill = NA, align = "center")
  half_max_idx = which(y2 == max(y2, na.rm = TRUE))[1]
  half_max = x[half_max_idx]
  return(half_max)
}

baranyi <- function(t, r, lag, logymax, logy0){
  At = t + (1 / r) * log(exp(-r * t) + exp(-r * lag) - exp(-r * (t + lag)))
  logy = logy0 + r * At - log(1 + (exp(r * At) - 1) / exp(logymax - logy0))
  return(logy)
}


plot_well <- function(results, well_ID, response_var = "OD"){
  curr_dat <- results %>% dplyr::filter(well == well_ID)
  curr_dat <- arrange(curr_dat, time_s)
  plot(curr_dat$time_s/3600, log(curr_dat[,response_var]), type = "l")
}