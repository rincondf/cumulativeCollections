aphids2022A <- aggregate(aphids2022$pea_aphids, list(aphids2022$date, aphids2022$site), sum, na.rm = TRUE)

aphids2022A$Date <- as.Date(aphids2022A$Date, format = "%m/%d/%Y")
aphids2022A$julian <- as.numeric(format(aphids2022A$Date, "%j"))






require(daymetr)

A3_bTs <- download_daymet(
  site = "Daymet",
  lat = 49.600334,
  lon = -119.681456,
  start = 1995,
  end = 2015,
  path = tempdir(),
  internal = TRUE,
  silent = FALSE,
  force = FALSE,
  simplify = FALSE
)