aphids2022A <- aggregate(aphids2022$pea_aphids, list(aphids2022$date, aphids2022$site), sum, na.rm = TRUE) # Aggregate transects

# Formalize dates and calculate julians 

aphids2022A$Date <- as.Date(aphids2022A$Date, format = "%m/%d/%Y")
aphids2022A$julian <- as.numeric(format(aphids2022A$Date, "%j"))


# Retrieve temperatures at sites

require(daymetr)

temp2022 <- download_daymet_batch(
  file_location = "locations2022.csv",
  start = 2022,
  end = 2022,
  internal = TRUE,
  force = FALSE,
  silent = FALSE,
  path = tempdir(),
  simplify = FALSE
)


# Calculate cumulative DDs for each site in farenheith

upperT <- c_to_f(30)
lowerT <- c_to_f(5.5)

tmax2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[7]))))
tmin2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[8]))))

seq()

hells_gate <- cumsum(calc_dd_vec(tmax = tmax2022[1:365], tmin = tmax2022[1:365], 
                                 lower_threshold = lowerT, upper_threshold = upperT, 
                                 cutoff = "horizontal"))

NJ_landing <- cumsum(calc_dd_vec(tmax = tmax2022[366:730], tmin = tmax2022[366:730], 
                                 lower_threshold = lowerT, upper_threshold = upperT, 
                                 cutoff = "horizontal"))


NJ_canyon <- cumsum(calc_dd_vec(tmax = tmax2022[731:1095], tmin = tmax2022[732:1095], 
                                 lower_threshold = lowerT, upper_threshold = upperT, 
                                 cutoff = "horizontal"))


alpowa_creek <- cumsum(calc_dd_vec(tmax = tmax2022[1096:1460], tmin = tmax2022[1098:1463], 
                                lower_threshold = lowerT, upper_threshold = upperT, 
                                cutoff = "horizontal"))


wawawai_park <- cumsum(calc_dd_vec(tmax = tmax2022[1461:1825], tmin = tmax2022[1464:1829], 
                                   lower_threshold = lowerT, upper_threshold = upperT, 
                                   cutoff = "horizontal"))


blyton_landing <- cumsum(calc_dd_vec(tmax = tmax2022[1826:2190], tmin = tmax2022[1830:2195], 
                                   lower_threshold = lowerT, upper_threshold = upperT, 
                                   cutoff = "horizontal"))


rose_creek <- cumsum(calc_dd_vec(tmax = tmax2022[2191:2555], tmin = tmax2022[2196:2561], 
                                     lower_threshold = lowerT, upper_threshold = upperT, 
                                     cutoff = "horizontal"))

sunshine_road <- cumsum(calc_dd_vec(tmax = tmax2022[2556:2920], tmin = tmax2022[2562:2927], 
                                 lower_threshold = lowerT, upper_threshold = upperT, 
                                 cutoff = "horizontal"))

grimes_road <- cumsum(calc_dd_vec(tmax = tmax2022[2921:3293], tmin = tmax2022[2928:3293], 
                                  lower_threshold = lowerT, upper_threshold = upperT, 
                                  cutoff = "horizontal"))

boyer_park <- cumsum(calc_dd_vec(tmax = tmax2022[3285:3659], tmin = tmax2022[3294:3659], 
                                  lower_threshold = lowerT, upper_threshold = upperT, 
                                  cutoff = "horizontal"))


evans_pond <- cumsum(calc_dd_vec(tmax = tmax2022[3649:4025], tmin = tmax2022[3660:4025], 
                                 lower_threshold = lowerT, upper_threshold = upperT, 
                                 cutoff = "horizontal"))
