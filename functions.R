#' Calculate degree days using the single sine method
#'
#' Calculate degree days using the single sine method. This is based on the
#' \code{dd} function from the \href{https://github.com/wsu-das/das.service.models/blob/dev/app/Library/DegreedaysCalculator.php#L21-L104.}{DAS DegreedaysCalculator on GitHub}.
#' Also some parts of the code adapted from \href{https://rdrr.io/github/trenchproject/TrenchR/src/R/DDFunctions.R}{the TrenchR package}.
#' Per the GitHub documentation for DegreedaysCalculator.php,
#' "Calculates the degreeday value using the single-sine method
#' Code originally extracted from DAS v6.0 calculations."
#'
#' @param tmax A numeric maximum daily temperature value, Fahrenheit.
#' @param tmin A numeric minimum daily temperature value, Fahrenheit.
#' @param lower_threshold A numeric lower temperature threshold, in Fahrenheit, for the
#' developmental model being used.
#' @param upper_threshold A numeric upper temperature threshold, in Fahrenheit, for the
#' developmental model being used.
#' @param cuttoff A character string indicating the type of upper threshold cutoff
#' to be used in calculating degree days. Can be either horizontal or vertical.
#'
#' @return Returns a numeric value for the number of degree days for a day with
#' the temperature parameters and cutoff specified.
#'
#' @seealso The degree day explainer \href{http://ipm.ucanr.edu/WEATHER/ddconcepts.html}{from UC IPM.}\cr
#' The \href{https://rdrr.io/github/trenchproject/TrenchR/src/R/DDFunctions.R}{source code} from the TrenchR package.
#' \cr
#' Zalom, F. G., P. B. Goodell, L. T. Wilson, W. W. Barnett, and W. J. Bentley.
#' 1983. Degree-Days: The Calculation and Use of degree days in Pest Management.
#' University of California Division of Agriculture and Natural Resources Leaflet 21373.
#' \href{https://www.researchgate.net/publication/262337197_Degree-Days_the_Calculation_and_Use_of_Heat_Units_in_Pest_Management_University_of_California_Division_of_Agriculture_and_Natural_Resources_Leaflet_21373}{Link.}


calc_dd <- function(tmax, tmin, lower_threshold, upper_threshold, cutoff){
  
  # No degree day accumulation:
  if(tmin > tmax | tmax <= lower_threshold) {
    return(0)
  }
  
  sum_heat <- tmax + tmin
  diff_heat <- tmax - tmin
  
  # When vertical cutoff for upper threshold is used...
  if(cutoff == "vertical") {
    
    # No degree day accumulation: Options 1 & 2 from T5 in Zalom et al
    if(tmin >= upper_threshold | tmax <= lower_threshold) {
      return(0)
    }
    
    # Full temp range is within thresholds: Option 3 in T5 of Zalom et al
    if(lower_threshold <= tmin & upper_threshold >= tmax) {
      return((sum_heat / 2) - lower_threshold)
    }
    
    # Looks like "a" corresponds roughly to a portion of theta-1 in the Zalom diagrams
    a <- 2 * lower_threshold - sum_heat
    
    if(abs(diff_heat) > abs(a)) {
      b <- atan(a / sqrt(diff_heat * diff_heat - a * a))
    } else {
      b <- 0
    }
    
    # Looks like "c" corresponds roughly to a portion of theta-2 in the Zalom diagrams
    c <- 2 * upper_threshold - sum_heat
    
    if(abs(diff_heat) > abs(c)) {
      d <- atan(c / sqrt(diff_heat * diff_heat - c * c))
    } else {
      d <- 0
    }
    
    # Upper threshold limiting heat accumulation: Option 5 in T5 of Zalom et al
    if(lower_threshold <= tmin) {
      return((-diff_heat * cos(d) - a * (d + (0.5 * pi))) / (2 * pi))
    }
    
    # Both thresholds limiting heat accumulation: Option 6 in T5 of Zalom et al
    if(upper_threshold < tmax) {
      return((-diff_heat * (cos(d) - cos(b)) - a * (d - b)) / (2 * pi))
    }
    
    # Default to lower temp threshold limiting heat accumulation: Option 4
    # of T5 in Zalom et al
    return((diff_heat * cos(b) - a * ((0.5 * pi) - b)) / (2 * pi))
    
    
    # When a non-vertical cutoff method is used...
  } else {
    
    a <- sum_heat / 2 - lower_threshold
    
    # Full temp range is above both thresholds: Options 1 in T5 of Zalom et al
    if (tmin >= upper_threshold && tmax > upper_threshold) {
      return(upper_threshold - lower_threshold)}
    
    # Full temp range is within thresholds: Option 3 in T5 of Zalom et al
    if(lower_threshold <= tmin & upper_threshold >= tmax) {
      return(a)
    }
    
    # Part of theta-1
    b <- 2 * lower_threshold - sum_heat
    
    
    if(abs(diff_heat) > abs(b)) {
      c <- atan(b / sqrt(diff_heat * diff_heat - b * b))
      d <- (diff_heat * cos(c) - b * ((0.5 * pi) - c)) / (2 * pi)
      # Below both thresholds: Option 2 of T5 in Zalom et al
    } else {
      d <- 0
    }
    
    # Lower temp threshold limiting heat accumulation: Option 4 of T5 in Zalom
    # et al (or Option 2)
    if(upper_threshold >= tmax) {
      return(d)
    }
    
    # Part of theta-2
    e <- 2 * upper_threshold - sum_heat
    
    f <- atan(e / sqrt(diff_heat * diff_heat - e * e))
    
    g <- (diff_heat * cos(f) - e * ((0.5 * pi) - f)) / (2 * pi)
    
    
    # Both thresholds limiting heat accumulation: Option 6 in T5 of Zalom et al
    if(lower_threshold > tmin) {
      return((d - g))
    }
    
    # Default to upper threshold limiting heat accumulation: Option 5 in T5 of
    # Zalom et al
    return((a - g))
  }
  
}

# This does the same thing es the function above, but accepts vectors for tmax and tmin

calc_dd_vec <- function(tmax, tmin, lower_threshold, upper_threshold, cutoff) {
  ddss <- rep(NA, length(tmax))
  for(i in 1: length(tmax)) {
    ddss[i] <- calc_dd(tmax[i], tmin[i], lower_threshold, upper_threshold, cutoff)
  }
  ddss
}

# handy functions to convert Celsius to Fahrenheit and viceversa

c_to_f <- function(x) {
  (x * 9/5) + 32
}

f_to_c <- function(x) {
  (x - 32) / 1.8
}
