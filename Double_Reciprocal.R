## ------------------------------------------------------------------------------- ##
#           lb(data, unit_S, unit_v, weighting = F, plot = T)                       #
## ------------------------------------------------------------------------------- ##
#' Lineweaver-Burk Transformation
#' @description Obtain Km and Vm using double reciprocal transformation
#' @usage lb(data, unit_S = 'mM', unit_v = 'au', weighting = FALSE, plot = TRUE)
#' @param data a dataframe where the first column is the independent variable, [S], and the reamining columns (as many as experiment replicates) correspond to the dependent variable, v.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param weighting logical. When TRUE the weight 1/v^4 is employed.
#' @param plot logical. If TRUE the data and fitted line are plotted.
#' @detail
#' @return A double reciprocal plot and the Km and Vm computed using averaged 1/v (when more than one replicate is provided). In addition, this function returns a list of five elements. The first and second ones are vectors with the Km and Vm, respectively, computed individually for each replicate. The third one provides the R-squared values of the fits. The fourth element of the list gives the fitted Km and Vm. The last element of the list is a datafreme with the values of the transformed variables.
#' @author Juan Carlos Aledo
#' @examples double.rec(rMM(replicates = 3, plot = F, error = 'f', sd=0.05*Vm)[,c(1,3:5)])
#' @references J. Am. Chem. Soc.1934, 56, 3,658-666 (doi.org/10.1021/ja01318a036)
#' @seealso hw(), eh(), ecb()
#' @export

lb <- function(data, unit_S = 'mM', unit_v = 'au', weighting = FALSE, plot = TRUE){

  ## ---------------------------- Removing incomplete data ----------------------- ##
  data <- data[complete.cases(data), ]
  colnames(data)[1] <- 'S'

  ## -------------------------- Taking double reciprocal ----------------------------- ##
  tdata <- apply(data, MARGIN = 2, function(x)  1/x)
  tdata <- as.data.frame(round(tdata, 4))
  names(tdata)[1] <- 'inv_S'


  ## ----------------------- Replicates model fitting ----------------------------- ##
  Kms <- Vms <- R2 <- c()
  for (j in 2:ncol(tdata)){
    model <- lm(tdata[,j] ~ tdata$inv_S )
    Vm <- round(1/model$coefficients[1], 2)
    Km <- round(Vm * model$coefficients[2], 2)
    Vms <- c(Vms, Vm)
    Kms <- c(Kms, Km)
    R2 <- c(R2,  summary(model)$r.squared)
  }

  ## -------- Computing mean and sd if required of the transformed data ---------- ##
  if (ncol(tdata) > 2){
    tdata$inv_v <- round(apply(tdata[,-1], MARGIN = 1, mean), 4)
    tdata$sd <- round(apply(tdata[,-c(1, ncol(tdata))], MARGIN = 1, sd), 4)
  } else {
    names(tdata)[2] <- 'inv_v'
  }
  ## -------------------------- Mean model fitting -------------------------------- ##
  if (weighting){
    model <- lm(tdata$inv_v ~ tdata$inv_S, weights = 1/(tdata$inv_v)^4)
  } else {
    model <- lm(tdata$inv_v ~ tdata$inv_S)
  }

  Vm <- round(1/model$coefficients[1], 2)
  Km <- round(Vm * model$coefficients[2], 2)

  ## --------------------- Plotting the transformed data ------------------------- ##
  if (plot == TRUE){
    parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
    plot(tdata$inv_S, tdata$inv_v, ty = 'p',
         ylim = c(0, max(tdata[,-1]) + 0.1*max(tdata[,-1])),
         xlab = paste('1/[S] (1/', unit_S, ')', sep = ""),
         ylab = paste('1/v (1/', unit_v, ')', sep = ""), main = parameters)
    abline(model)

    if (ncol(tdata) > 2){
      arrows(tdata$inv_S, tdata$inv_v - tdata$sd,
             tdata$inv_S, tdata$inv_v + tdata$sd,
             length=0.05, angle=90, code=3)
    }
  }


  fitted_parameters <- c(Km, Vm)
  names(fitted_parameters) <- c('Km', 'Vm')
  output <- list(unname(Kms), unname(Vms), R2, fitted_parameters, tdata)
  names(output) <- c('Kms', 'Vms', 'R2s', 'fitted_parameters', 'inverse_data')
  return(output)
}
