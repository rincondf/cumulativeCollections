# 1. Aggregate transects

aphids2022A <- aggregate(aphids2022$pea_aphids, list(aphids2022$date, aphids2022$site), sum, na.rm = TRUE) 

# 2. Formalize dates and calculate julians 

aphids2022A$Date <- as.Date(aphids2022A$Date, format = "%m/%d/%Y")
aphids2022A$julian <- as.numeric(format(aphids2022A$Date, "%j"))

# 3. Retrieve temperatures at sites

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

# 4. Calculate cumulative DDs for each site in Fahrenheit

# 4.1. Max and min thresholds in Fahrenheit

upperT <- c_to_f(30)
lowerT <- c_to_f(5.5)

# 4.2. Extract max and min temperatures of all sites and convert them to Fahrenheit

tmax2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[7]))))
tmin2022 <- as.numeric(unlist(lapply(temp2022, function(x) c_to_f(x$data[8]))))

# 4.3. Calculation of cumulative degree days using retrieved max and min temps and the DAS function for DDs calculation

ind <- seq(1, length(tmax2022), 365) #required indexation
DDs2022 <- rep(NA, length(tmax2022))

for(i in 1: length(locations2022$Site)) {
  DDs2022[ind[i]: (ind[i] + 364)] <- cumsum(calc_dd_vec(tmax = tmax2022[ind[i]: (ind[i] + 364)], tmin = tmax2022[ind[i]: (ind[i] + 364)], 
                                                        lower_threshold = lowerT, upper_threshold = upperT, 
                                                        cutoff = "horizontal"))
}

# Data frame with cumulative degree days per site

DDs2022 <- data.frame(site = rep(locations2022$Site, each = 365), DDs = DDs2022, julian  = rep(seq(1, 365), length(locations2022$Site)))

# 5. Assign DDs to each site and julian date in the dataset

aphids2022A$DDs <- rep(NA, length(aphids2022A[, 1]))

for(i in 1: length(locations2022$Site)){
  aphids2022A$DDs[which(aphids2022A$Site == locations2022$Site[i])] <- DDs2022$DDs[which((DDs2022$site == locations2022$Site[i]) & 
                                                                                  (DDs2022$julian %in% 
                                                                                     aphids2022A$julian[which(aphids2022A$Site == locations2022$Site[i])]))]
}

# 6. Proportion of aphids captured of the total 

aphids2022A$propT <- aphids2022A$Aphids / sum(aphids2022A$Aphids)

# 7. Cumulative proportion of aphids per site

aphids2022A$prop <- rep(NA, length(aphids2022A[, 1]))
for(i in 1: length(locations2022$Site)){
  aphids2022A$prop[which(aphids2022A$Site == locations2022$Site[i])] <- cumsum(aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])]) / sum(aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])])
}

# 8. Proportion of captured aphids per site

aphids2022A$prop1 <- rep(NA, length(aphids2022A[, 1]))
for(i in 1: length(locations2022$Site)){
  aphids2022A$prop1[which(aphids2022A$Site == locations2022$Site[i])] <- (aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])]) / sum(aphids2022A$Aphids[which(aphids2022A$Site == locations2022$Site[i])])
}


# 9. Parameter estimation by maximum likelihood

require(bbmle)
require(ExtDist)
require(SuppDists)

# This is a function that wraps mle2 (and optim) to find parameters of a Johnson SB for a given set of degree days (DDs) 
# and the corresponding proportion of captured insects (prs). Only works with constrained Broyden–Fletcher–Goldfarb–Shanno algorithm, and
# Nelder-Mead. The first is more robust, but convergence is difficult. The second allows an easier convergence.

estimat <- function(DDs, prs, method){
  LL1 <- function(gamma, delta, a, b) {
    -sum(pr * log(dJohnsonSB_ab(x = x, gamma = gamma, delta = delta, a = a, b = b)))
  }
  
  if(method == "L-BFGS-B"){
    MLL <- mle2(LL1, start = list(gamma = -0.5, delta = 2, a = min(DDs) - 1, b = max(DDs) + 1), 
                data = list(x = DDs, pr = prs),
                lower = list(gamma = -Inf, delta = 0, a = -Inf, b = max(DDs)), 
                upper = list(gamma = Inf, delta = Inf, a = min(DDs), b = Inf), method = "L-BFGS-B")
  }
  
  if(method == "Nelder-Mead"){
    MLL <- mle2(LL1, start = list(gamma = -0.5, delta = 2, a = min(DDs) - 1, b = max(DDs) + 1), 
                data = list(x = DDs, pr = prs),
                method = "Nelder-Mead")
  }
  
  MLL
}

# Model using proportions of totals from each site
MLL1 <- estimat(DDs = aphids2022A$DDs, prs = aphids2022A$prop1, method = "L-BFGS-B")

# Model using proportions of totals
MLL2 <- estimat(DDs = aphids2022A$DDs, prs = aphids2022A$propT, method = "Nelder-Mead")


summary(MLL1)
summary(MLL2)

funres <- function(DDs){
  xi = coef(MLL1)[3]
  lambda = coef(MLL1)[4] - coef(MLL1)[3]
  pJohnsonSB(DDs, gamma = coef(MLL1)[1],
             delta = coef(MLL1)[2],
             xi = xi,
             lambda = lambda)
}

funres2 <- function(DDs){
  xi = coef(MLL2)[3]
  lambda = coef(MLL2)[4] - coef(MLL2)[3]
  pJohnsonSB(DDs, gamma = coef(MLL2)[1],
             delta = coef(MLL2)[2],
             xi = xi,
             lambda = lambda)
}


par(mar = c(5, 5, 2, 2) + 0.1)
plot(aphids2022A$DDs, aphids2022A$prop, xlab = "Cumulative Degree Days", 
     ylab = "Proportion captured", cex.lab = 2, cex.axis = 2, lwd = 2)
curve(funres, from = 100, to = 5000, lwd = 2, add = TRUE, col = "blue")


par(mar = c(5, 5, 2, 2) + 0.1)
plot(aphids2022A$DDs[order(aphids2022A$DDs)], cumsum(aphids2022A$propT[order(aphids2022A$DDs)]), xlab = "Cumulative Degree Days", 
     ylab = "Proportion captured", cex.lab = 2, cex.axis = 2, lwd = 2)
curve(funres2, from = 100, to = 5000, lwd = 2, add = TRUE, col = "blue")