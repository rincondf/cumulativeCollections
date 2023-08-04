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

upperT <- c_to_f(31)
lowerT <- c_to_f(10)

tmax2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[7]))))
tmin2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[8]))))

ind <- seq(1, length(tmax2022), 365)


DDs2022 <- rep(NA, length(tmax2022))

for(i in 1: length(locations2022$Site)) {
  DDs2022[ind[i]: (ind[i] + 364)] <- cumsum(calc_dd_vec(tmax = tmax2022[ind[i]: (ind[i] + 364)], tmin = tmax2022[ind[i]: (ind[i] + 364)], 
                                                        lower_threshold = lowerT, upper_threshold = upperT, 
                                                        cutoff = "horizontal"))
}

DDs2022 <- data.frame(site = rep(locations2022$Site, each = 365), DDs = DDs2022, julian  = rep(seq(1, 365), length(locations2022$Site)))

# assign DDs to each site and julian date in the datset

aphids2022A$DDs <- rep(NA, length(aphids2022A[, 1]))

for(i in 1: length(locations2022$Site)){
  aphids2022A$DDs[which(aphids2022A$Site == locations2022$Site[i])] <- DDs2022$DDs[which((DDs2022$site == locations2022$Site[i]) & 
                                                                                  (DDs2022$julian %in% 
                                                                                     aphids2022A$julian[which(aphids2022A$Site == locations2022$Site[i])]))]
}


aphids2022A$propT <- aphids2022A$Aphids / sum(aphids2022A$Aphids)
plot(aphids2022A$DDs[order(aphids2022A$DDs)], cumsum(aphids2022A$propT[order(aphids2022A$DDs)]))


aphids2022A$prop <- rep(NA, length(aphids2022A[, 1]))
for(i in 1: length(locations2022$Site)){
  aphids2022A$prop[which(aphids2022A$Site == locations2022$Site[i])] <- cumsum(aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])]) / sum(aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])])
}

plot(aphids2022A$DDs, aphids2022A$prop)

# 

