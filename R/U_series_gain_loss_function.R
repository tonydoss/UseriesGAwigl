globalVariables(c("R234U_238U", "R230Th_238U", "Depth"))

#' U_series_gain_loss_function
#'
#' @param input_data Read a CSV file with the required data. For instance, if the CSV file is called "data.csv" and is in your working directory, have input_data = read.csv("data.csv")
#' @param nbit Number of iterations
#' @param logT_min Log10 of the minimum value allowed for the weathering age of the shallowest sample
#' @param logT_max Log10 of the maximum value allowed for the weathering age of the shallowest sample
#' @param Tmult_min Minimum value allowed for Tmult
#' @param logk238_min Log10 of the minimum value allowed for k238
#' @param logk238_max Log10 of the maximum value allowed for k238
#' @param k48_min Minimum value allowed for k234/k238
#' @param k48_max Maximum value allowed for k234/k238
#' @param logf238_k238_min Log10 of the minimum value allowed for f238/k238
#' @param logf238_k238_max Log10 of the maximum value allowed for f238/k238
#' @param f48_min Minimum value allowed for f234/f238
#' @param f48_max Maximum value allowed for f234/f238
#'
#' @return
#' @export
#'
#' @examples U_series_gain_loss_function (read.csv(system.file("extdata", "data.csv",
#' package = "UseriesGAwigl")),
#' logT_min = 1, logT_max = 7, Tmult_min = 0.5,
#' logk238_min = -7, logk238_max = -5,
#' k48_min = 0.1, k48_max = 10,
#' logf238_k238_min = -20, logf238_k238_max = 3,
#' f48_min = 0.1, f48_max = 10)
U_series_gain_loss_function <- function(
  input_data,
  nbit = 1000,
  logT_min = 1,
  logT_max = 7,
  Tmult_min = 0.5,
  logk238_min = -7,
  logk238_max = -5,
  k48_min = 1,
  k48_max = 2,
  logf238_k238_min = -8,
  logf238_k238_max = +0.5,
  f48_min = 0.1,
  f48_max = 10
){

  # # Example
  # U_series_gain_loss_function (read.csv("data.csv"),
  # logT_min = 1, logT_max = 7, Tmult_min = 0.5,
  # logk238_min = -7, logk238_max = -5, k48_min = 0.1, k48_max = 10,
  # logf238_k238_min = -20, logf238_k238_max = 3, f48_min = 0.1, f48_max = 10)

  # Import data -------------------------------------------------------------

  df <- input_data # this calls in our data we've stored in a csv format.
  df <- df[order(df$Depth..m.),] # sort data frame by increasing depth
  nb_sample <- nrow(df)-1 # nb of samples to solve - exclude deepest sample in solving, since used as initial conditions (so weathering age should be 0 a, which the model might have a hard time achieving)

  # Set variables for genetic algorithm -------------------------------------

  # nbrun <- 100 # the number of consecutive generations without any improvement in the best fitness value before the GA is stopped
  # elitism <- 10 # the number of best fitness individuals to survive at each generation. By default the top 5% individuals will survive at each iteration
  # maxfit <- -0.01
  # selection <- "gareal_lrSelection"
  # popsize <- 100

  # Set inital range of values -----------------------------------------------------

  Tmult_max <- 1

  # range of values allowed for (234U/238U) at t = 0
  r48i_min <- df$R234U_238U[which(df$Depth..m.==max(df$Depth..m.))]-2*df$R48_2SE[which(df$Depth..m.==max(df$Depth..m.))]
  r48i_max <- df$R234U_238U[which(df$Depth..m.==max(df$Depth..m.))]+2*df$R48_2SE[which(df$Depth..m.==max(df$Depth..m.))]
  # range of values allowed for (230Th/238U) at t = 0
  r08i_min <- df$R230Th_238U[which(df$Depth..m.==max(df$Depth..m.))]-2*df$R08_2SE[which(df$Depth..m.==max(df$Depth..m.))]
  r08i_max <- df$R230Th_238U[which(df$Depth..m.==max(df$Depth..m.))]+2*df$R08_2SE[which(df$Depth..m.==max(df$Depth..m.))]

  min_val <- c(logk238_min, k48_min, logf238_k238_min,
               f48_min, r48i_min, r08i_min, logT_min, rep(Tmult_min,
                                                          nb_sample - 1))

  max_val <- c(logk238_max, k48_max, logf238_k238_max,
               f48_max, r48i_max, r08i_max, logT_max, rep(Tmult_max,
                                                          nb_sample - 1))

  # Constants ---------------------------------------------------------------


  # values for 230Th and 232Th are set
  k230 <- 1e-10
  k232 <- 1e-10
  f230 <- 1e-10
  k232 <- 1e-10

  # decay constants
  l238 <- 0.1551e-9
  l234 <- 2.826e-6
  l232 <- 4.948e-11
  l230 <- 9.158e-6

  # U concentration at t=0 (in ppm); value doesn't matter
  U_init <- 2

  # Initialise vectors ------------------------------------------------------


  # initialise vectors and matrices used in loops
  R48_target <- vector(mode="numeric", length=nb_sample)
  R08_target <- vector(mode="numeric", length=nb_sample)
  R48_err <- vector(mode="numeric", length=nb_sample)
  R08_err <- vector(mode="numeric", length=nb_sample)
  depth <- vector(mode="numeric", length=nb_sample)
  T <- vector(mode="numeric", length=nb_sample)
  err_T <- vector(mode="numeric", length=nb_sample)
  R48_calc <- vector(mode="numeric", length=nb_sample)
  R08_calc <- vector(mode="numeric", length=nb_sample)
  f <- vector(mode="numeric", length=1)


  # set target ratios -------------------------------------------------------


  for (j in 1:nb_sample){
    R48_target[j] <- df$R234U_238U[j]
    R08_target[j] <- df$R230Th_238U[j]
    depth[j] <- as.matrix(df$Depth..m.[j])
    R48_err[j] <- df$R48_2SE[j]
    R08_err[j] <- df$R08_2SE[j]
  }

  target_ratios <- rbind(depth, R48_target, R48_err, R08_target, R08_err)
  target_ratios <- as.data.frame(t(target_ratios))


  # function to minimise ----------------------------------------------------


  fn <- function(x) {
    k238 <- 10^as.numeric(x[1])
    k234 <- k238*as.numeric(x[2])
    f238 <- k238*10^as.numeric(x[3])
    f234 <- f238*as.numeric(x[4])
    r48i <- as.numeric(x[5])
    r08i <- as.numeric(x[6])

    for (j in 1:nb_sample){
      if (j == 1){
        T[j] <-  10^as.numeric(x[7])}
      else {
        T[j] <- T[j-1]*as.numeric(x[6+j])
      }
    }

    r04i <- r08i/r48i
    U8i <- U_init*1e-6/238*6.02e23
    U4i <- r48i*l238/l234*U8i
    Th0i <- r04i*l234/l230*U4i

    F238 <- f238*l238*U8i
    F234 <- f234*l234*U4i
    F230 <- f230*l230*Th0i

    a4 <- l234 + k234
    a8 <- l238 + k238
    a0 <- l230 + k230

    A <- l230*l234/((a0 - a8)*(a4 - a8)) - l230*l234*F238/(l238*U8i*a8*(a4 - a8)*(a0 - a8))
    B <- l230/(a0 - a4)*r48i - l230*l234/((a0 - a4)*(a4 - a8)) - l230*F234/(l238*U8i*a4*(a0 - a4)) - l230*l234*F238/(l238*U8i*a4*(a0 - a4)*(a4 - a8))
    C <- r08i - l230/(a0 - a4)*r48i + l230*l234/((a0 - a4)*(a4 - a8)) - l230*l234/((a0 - a8)*(a4 - a8)) - F230/(l238*U8i*a0) - l230*F234/(l238*U8i*a0*a4) + l230*F234/(l238*U8i*a4*(a0 - a4)) - l230*l234*F238/(l238*U8i*a0*a4*a8) + l230*l234*F238/(l238*U8i*a4*a8*(a0 - a4)) - l230*l234*F238/(l238*U8i*a8*(a0 - a4)*(a4 - a8)) + l230*l234*F238/(l238*U8i*a8*(a0 - a8)*(a4 - a8))

    for (j in 1:nb_sample){
      R48_calc[j] <- ((l234/(a4 - a8) - l234*F238/(l238*U8i*a8*(a4 - a8)))*exp(-a8*T[j]) + (r48i - l234/(a4 - a8) - l234*F238/(l238*U8i*a4*(a4 - a8)) - F234/(l238*U8i*a4))*exp(-a4*T[j]) + F234/(l238*U8i*a4) + l234*F238/(l238*U8i*a4*a8))/((1 - F238/(l238*U8i*a8))*exp(-a8*T[j]) + F238/(l238*U8i*a8))
      R08_calc[j] <- (A*exp(-a8*T[j]) + B*exp(-a4*T[j]) + C*exp(-a0*T[j]) + l230*l234*F238/(l238*U8i*a0*a4*a8) + l230*F234/(l238*U8i*a0*a4) + F230/(l238*U8i*a0))/((1 - F238/(l238*U8i*a8))*exp(-a8*T[j]) + F238/(l238*U8i*a8))
    }

    f <- sum((R48_calc - R48_target)^2 + (R08_calc - R08_target)^2)
  }


  # solve -------------------------------------------------------------------
  # if(!require(GA)){
  #   install.packages("GA")
  #   require(GA)}

  GA <- GA::ga(type = "real-valued", fitness = function(x) -fn(x),
           lower = min_val, upper = max_val,
           # elitism = elitism, run = nbrun,
           # maxFitness = maxfit,
           # popSize = popsize,
           # selection = selection,
           maxiter = nbit, seed = 123, monitor = F)

  solution <- GA@solution


  # produce ratios ----------------------------------------------------------

  k238 <- 10^as.numeric(solution[1])
  k234 <- k238*as.numeric(solution[2])
  f238 <- k238*10^as.numeric(solution[3])
  f234 <- f238*as.numeric(solution[4])
  r48i <- as.numeric(solution[5])
  r08i <- as.numeric(solution[6])

  for (j in 1:nb_sample){
    if (j == 1){
      T[j] <-  10^as.numeric(solution[7])}
    else {
      T[j] <- T[j-1]*as.numeric(solution[6+j])
    }
  }

  r04i <- r08i/r48i
  U8i <- U_init*1e-6/238*6.02e23
  U4i <- r48i*l238/l234*U8i
  Th0i <- r04i*l234/l230*U4i

  F238 <- f238*l238*U8i
  F234 <- f234*l234*U4i
  F230 <- f230*l230*Th0i

  a4 <- l234 + k234
  a8 <- l238 + k238
  a0 <- l230 + k230
  A <- l230*l234/((a0 - a8)*(a4 - a8)) - l230*l234*F238/(l238*U8i*a8*(a4 - a8)*(a0 - a8))
  B <- l230/(a0 - a4)*r48i - l230*l234/((a0 - a4)*(a4 - a8)) - l230*F234/(l238*U8i*a4*(a0 - a4)) - l230*l234*F238/(l238*U8i*a4*(a0 - a4)*(a4 - a8))
  C <- r08i - l230/(a0 - a4)*r48i + l230*l234/((a0 - a4)*(a4 - a8)) - l230*l234/((a0 - a8)*(a4 - a8)) - F230/(l238*U8i*a0) - l230*F234/(l238*U8i*a0*a4) + l230*F234/(l238*U8i*a4*(a0 - a4)) - l230*l234*F238/(l238*U8i*a0*a4*a8) + l230*l234*F238/(l238*U8i*a4*a8*(a0 - a4)) - l230*l234*F238/(l238*U8i*a8*(a0 - a4)*(a4 - a8)) + l230*l234*F238/(l238*U8i*a8*(a0 - a8)*(a4 - a8))

  R48_calc <- vector(mode="numeric", length=nb_sample)
  R08_calc <- vector(mode="numeric", length=nb_sample)

  for (j in 1:nb_sample){
    R48_calc[j] <- ((l234/(a4 - a8) - l234*F238/(l238*U8i*a8*(a4 - a8)))*exp(-a8*T[j]) + (r48i - l234/(a4 - a8) - l234*F238/(l238*U8i*a4*(a4 - a8)) - F234/(l238*U8i*a4))*exp(-a4*T[j]) + F234/(l238*U8i*a4) + l234*F238/(l238*U8i*a4*a8))/((1 - F238/(l238*U8i*a8))*exp(-a8*T[j]) + F238/(l238*U8i*a8))
    R08_calc[j] <- (A*exp(-a8*T[j]) + B*exp(-a4*T[j]) + C*exp(-a0*T[j]) + l230*l234*F238/(l238*U8i*a0*a4*a8) + l230*F234/(l238*U8i*a0*a4) + F230/(l238*U8i*a0))/((1 - F238/(l238*U8i*a8))*exp(-a8*T[j]) + F238/(l238*U8i*a8))
  }

  data <- as.data.frame(cbind(R48_calc, R08_calc, T))
  data <- cbind(data, df$Depth..m.[1:nb_sample])
  colnames(data) <- c("R234U_238U", "R230Th_238U", "T", "Depth")


  # graphs ------------------------------------------------------------------

  theme_plots <- ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = (15), hjust = 0.5),
                       legend.title=ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_text(size=18),
                       axis.title.x = ggplot2::element_text(size=18),
                       legend.text=ggplot2::element_text(size=12)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())


  # (234U/238U) vs depth
  p48 <- ggplot2::ggplot(target_ratios, ggplot2::aes(R48_target,depth)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmax = R48_target+R48_err, xmin = R48_target-R48_err, height=0.02)) + # plot error bars
    ggplot2::geom_point(size = 5, colour = 'blue') +
    ggplot2:: geom_point(data=data, ggplot2::aes(R234U_238U,Depth), size = 5, colour = 'red') + ggplot2::scale_y_reverse() +
    theme_plots +  ggplot2::xlab(expression("("^234*"U/"^238*"U)")) +
    ggplot2::ylab('Depth (m)')

  # (230Th/238U) vs depth
  p08 <- ggplot2::ggplot(target_ratios, ggplot2::aes(R08_target,depth)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmax = R08_target+R08_err, xmin = R08_target-R08_err, height=0.02)) + # plot error bars
    ggplot2::geom_point(size = 5, colour = 'blue') +
    ggplot2::geom_point(data=data, ggplot2::aes(R230Th_238U,Depth), size = 5, colour = 'red') + ggplot2::scale_y_reverse() +
    theme_plots +  ggplot2::xlab(expression("("^230*"Th/"^238*"U)")) +
    ggplot2::ylab('Depth (m)')

  # Weathering age vs depth
  pt <- ggplot2::ggplot(data, ggplot2::aes(T/1000,Depth)) +
    ggplot2::geom_point(size = 5, colour = 'red') + ggplot2::scale_y_reverse() + theme_plots +
    ggplot2::xlab('Weathering age (ka)') +
    ggplot2::ylab('Depth (m)')

  # if(!(require(cowplot))){
  #   install.packages("cowplot")
  #   library(cowplot)
  # }

  print_plot <- cowplot::plot_grid(p48, p08, pt, ncol = 3)

  # finish ------------------------------------------------------------------
  # Regolith production rate (mm/kyr)
  rate <- (df$Depth..m.[nb_sample+1]-df$Depth..m.[1])*1e6/T[1]
  print(print_plot)
  solution_df <- as.data.frame(t(c(k238, k234/k238, f238/k238, f234/f238, r48i, r08i)))
  solution_list <- list(regolith_production_rate = rate,
                        parameters = solution_df,
                        calculated_ratios_and_ages = data,
                        graphs = print_plot, nb_Sample = nb_sample,
                        min_val = min_val, max_val = max_val)

  return(solution_list)
}

