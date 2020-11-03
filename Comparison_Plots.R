
# Initialization ----------------------------------------------------------

if (1) {
    setwd('~/Documents/MFE/FE5222Advanced Derivatives Pricing/')
    # 
    library(plyr)
    library(tidyverse)
    library(ggplot2)
    # 
    # load('')
}

# American Put by LSMC ----------------------------------------------------

AmericanPut_LSMC <- function(
    Spot = 100, Strike = 100, sigma = 0.2, r = 0.06, dr = 0, mT = 1,
    n = 1000, m = 252) {
    # 
    #' Pricing a vanilla American put using LSMC method, with GBM assumption.
    #'
    #' @param Spot Spot price of the underlying asset
    #' @param sigma Volatility of the underlying asset
    #' @param Strike Strike of the option
    #' @param r interest rate
    #' @param dr dividend yield of the underlying asset
    #' @param mT Maturity, in year
    #' @param n Number of paths simulated 
    #' @param m Number of time steps in the simulation
    #'
    set.seed(2020)
    seq_seed <- seq(n)
    # GBM matrix (n rows, m cols), simulated with reproducibility
    GBM <- llply(
        .data = seq_seed, .fun = function(iter_seed){
            set.seed(iter_seed)
            seq_t <- seq(from = 0, to = mT, length.out = m)
            seq_GBM <- (r - dr - .5 * sigma * sigma) * (mT / m) + 
                sigma * sqrt(mT/m) * rnorm(m, mean = 0, sd = 1)
            seq_GBM <- Spot * exp(cumsum(seq_GBM))
            return(seq_GBM)
        }
    )
    GBM <- Reduce(f = 'rbind', GBM)
    # ITM indicator matrix (n rows, m cols), where the underlying price S < Strike
    ITM <- ifelse(GBM < Strike, GBM, NA) 
    # Payoff matrix, at the current time step
    Payoff <- matrix(pmax(0, Strike - GBM), nrow = n, ncol = m)
    # Xsh matrix, underlying prices at the current time step
    # X2sh matrix, Xsh squared
    Xsh <- ITM[, -m]
    X2sh <- Xsh * Xsh
    # Y1 matrix, 1-step discounted payoff
    Y1 <- Payoff * exp(-1 * r * (mT/m))
    # Y2 matrix, PV of the cashflow if we hold on, iterated backwards for regression
    Y2 <- cbind((matrix(NA, nrow = n, ncol = m - 1)), Y1[, m])
    # Value matrix, value of the option, to be updated
    Value <- matrix(NA, nrow = n, ncol = m - 1)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    try(
        for (iter_i in ((m-1):1)) {
            # Regression
            iter_fit <- lm(Y2[, iter_i + 1] ~ Xsh[, iter_i] + X2sh[, iter_i])
            # 
            Value[, iter_i] <- 
                (iter_fit$coefficients)[1] + 
                (iter_fit$coefficients)[2] * Xsh[, iter_i] + 
                (iter_fit$coefficients)[3] * X2sh[, iter_i]
            Value[, iter_i] <- ifelse(
                test = is.na(Value[, iter_i]), 
                yes = 0, 
                no = Value[, iter_i]
            )
            Y2[, iter_i] <- ifelse(
                test = Payoff[, iter_i] > Value[, iter_i], 
                yes = Y1[, iter_i], 
                no = Y2[,iter_i + 1] * exp(-1 * r * (mT/m))
            )
        }, 
        silent = TRUE
    )
    Value <- cbind(ifelse(is.na(Value), 0, Value), 0)
    Value <- ifelse(Value > Payoff, 0, Payoff)
    # Pick the first positive value, and discount to time 0
    Value <-
        t(apply(Value, MARGIN = 1, FUN = function(vec){
            if(sum(vec>0)) {
                temp_index <- min(which(vec>0))
                return(
                    c(rep(0, temp_index-1), 
                      vec[temp_index] * exp(-1 * mT/m * r * temp_index), 
                      rep(0, length(vec) - temp_index))
                )
            } else {return(vec)}
        }))
    # Averaging among simulated paths
    Value <- mean(rowSums(Value))
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    return(Value)
}





# European Put by BSM -----------------------------------------------------

BSM <- function(Spot = 100, Strike = 100, sigma = .2, r = 0.06, mT = 1) {
    #' BSM for European Put, to directly compute the value at m - 1 step
    
    d1 <- (log(Spot / Strike) + (r + 0.5 * sigma * sigma) * mT) / (sigma * sqrt(mT))
    d2 <- d1 - sigma * sqrt(mT)
    Value <- pnorm(- d2) * Strike * exp(- r * mT) - pnorm(- d1) * Spot
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    return(Value)
}


# American Put by BBS -----------------------------------------------------

AmericanPut_BBS <- function(Spot = 100, Strike = 100, sigma = .2, r = .06, mT = 1, 
                            m = 100) {
    #' Binomial Black Scholes model for American Puts
    #' @param Spot Spot price of the underlying asset
    #' @param Strike Strike of the option
    #' @param sigma Volatility of the underlying asset
    #' @param r interest rate
    #' @param mT Maturity, in year
    #' @param m  Number of time steps in the simulation
    #' 
    u <- exp(sigma * sqrt(mT/m))
    d <- 1/u
    p <- (exp(r * mT/m) - d)/(u - d)
    # 
    GBM <- matrix(0, nrow = m + 1, ncol = m + 1)
    GBM[1,1] <- Spot
    # GBM matrix (n rows, m cols), simulated with reproducibility
    for (i in 2:(m+1)) {
        GBM[i,1] <- GBM[i-1,1] * u
        for (j in 2:i) {
            GBM[i,j] <- GBM[i-1,j-1] * d
        }
    }
    # Option value at final & m-1 node   
    Value <- matrix(0, nrow = m + 1, ncol = m + 1)
    Value[m + 1,] <- pmax(0, Strike - GBM[m + 1,])
    Value[m,1:m] <- BSM(Spot = GBM[m,1:m], Strike = Strike, sigma = sigma, r = r, mT = mT/m)
    # Backward induction
    for (i in (m-1):1) {
        for (j in 1:i) {
            Value[i,j] <- max(0, 
                              Strike - GBM[i, j],
                              (p * Value[i+1, j] + (1-p) * Value[i+1, j+1]) / exp(r * mT/m)
            )
        }
    }
    # Pick the Value at time 0
    Value <- Value[1,1]
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    return(Value)
}


# American Put by BBSR ----------------------------------------------------

AmericanPut_BBSR <- function(Spot = 100, Strike = 100, sigma = .2, r = .06, mT = 1, 
                             m = 100) {
    #' Binomial Black Scholes model with Richardson extrapolation, 
    #' for American Puts, as the benchmark in this study
    #' @param Spot Spot price of the underlying asset
    #' @param Strike Strike of the option
    #' @param sigma Volatility of the underlying asset
    #' @param r interest rate
    #' @param mT Maturity, in year
    #' @param m  Number of time steps in the simulation
    #' 
    Value <- 2 * 
        AmericanPut_BBS(Spot = Spot, Strike = Strike, sigma = sigma, r = r, mT = mT, 
                        m = 2 * m) - 
        AmericanPut_BBS(Spot = Spot, Strike = Strike, sigma = sigma, r = r, mT = mT, 
                        m = m)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    return(Value)
}



# Comparison and Testing --------------------------------------------------

ContourPlots <- function(
    temp_df,
    # iter_X, iter_Y, iter_Z, 
    fill_label = 'ValueDiff\n(LSMC - BBSR, in $)',
    xlab = 'X', ylab = 'Y', fig_title = '', fig_subtitle = '',
    scale_range = NULL,
    flag_X = F,
    flag_plot = F, flag_save_plot = T, save_folder = ''
) {
    #' Title
    #' @param iter_X                        x seq
    #' @param iter_Y                        y seq
    #' @param iter_Z                        z matrix
    #' @param fill_label                    legend title
    # Need to plot
    require(ggplot2)
    lut_16 <- matrix(
        data = c(0,0,0,    1,1,171,    1,1,224,    0,110,255,    1,171,254,
                 1,224,254,    1,254,1,    190,255,0,    255,255,0,    255,224,0,
                 255,141,0,    250,94,0,    245,0,0,    245,0,172,    222,180,222),
        ncol = 3, byrow = T)
    lut_16 <- rgb(lut_16, maxColorValue = 255)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # temp_df <- expand.grid(X = iter_X, Y = iter_Y)    
    # temp_df$Z <- c(iter_Z)
    seq_X <- sort(unique(temp_df$X))
    seq_Y <- sort(unique(temp_df$Y))
    temp_df$X <- temp_df$X/max(seq_X)
    temp_df$Y <- temp_df$Y/max(seq_Y)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    temp_p <- ggplot(data = temp_df, aes(x = X, y = Y, z = Z)) + 
        geom_raster(data = temp_df, aes(fill = Z)) +
        coord_fixed() +
        ggplot2::labs(title = fig_title, subtitle = fig_subtitle, x = xlab, y = ylab, fill = fill_label) +
        # ggplot2::labs(title = 'fig_title', subtitle = 'fig_subtitle', x = 'xlab', y = 'ylab', fill = 'fill_label') +
        ggplot2::theme_minimal()
    
    # ggplot2::scale_y_discrete(labels = unique(seq_Y)/max(seq_Y), labels = unique(seq_Y))
    # ggplot2::theme_bw()
    if (flag_X) {
        temp_seq_X <- seq_X[seq(from = 1, to = length(seq_X), by = 2)]
        temp_p <- temp_p + 
            ggplot2::scale_x_continuous(breaks = temp_seq_X/max(temp_seq_X), labels = temp_seq_X) + 
            ggplot2::scale_y_continuous(breaks = seq_Y/max(seq_Y), labels = seq_Y)
    } else {
        temp_p <- temp_p + 
            ggplot2::scale_x_continuous(breaks = seq_X/max(seq_X), labels = seq_X) + 
            ggplot2::scale_y_continuous(breaks = seq_Y/max(seq_Y), labels = seq_Y)
    }
    
    if (is.null(scale_range)) {
        temp_p <- temp_p +
            ggplot2::scale_fill_gradientn(
                colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                ), na.value = 'white'
            )
    }else if (length(scale_range) == 2) {
        temp_p <- temp_p +
            ggplot2::scale_fill_gradientn(
                limits = c(scale_range[1], scale_range[2]),
                colours = lut_16,
                # colours = colorRamps::matlab.like(300)
                guide = ggplot2::guide_colourbar(
                    raster = T, frame.colour = "black", frame.linewidth = 1
                ), na.value = 'white'
            )
    }
    
    if (flag_save_plot) {
        cat('Saving', paste0(save_folder, fig_title, '.png'), '...\n')
        ggplot2::ggsave(filename = paste0(save_folder, fig_title, '.png'), plot = temp_p)
        # , width = 7, height = 6.5)
    }
    if (flag_plot) {print(temp_p)}
    return(temp_p)
}



# Grid Comparison for pricing parameters ----------------------------------
if (1) {
    temp_seq <- paste0('seq_', c('Spot', 'Strike', 'sigma', 'r', 'mT'))
    temp_seq %>% print()
    combn(length(temp_seq), 2)
    # pricing parameters
    seq_Spot <- 1:20 * 10
    seq_Strike <- 1:20 * 10
    seq_sigma <- 1:20 / 40
    seq_r <- 1:20 / 200
    seq_mT <- 1:20 / 10
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 1 Spot & Strike
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Spot, seq_Strike), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Spot = temp_grid$X[iter_i], 
                             Strike = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Spot = temp_grid$X[iter_i], 
                             Strike = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Spot', ylab = 'Strike', 
        fig_title = paste0('Spot', ' & ', 'Strike')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}


if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 2 Spot & sigma
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Spot, seq_sigma), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Spot = temp_grid$X[iter_i], 
                             sigma = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Spot = temp_grid$X[iter_i], 
                             sigma = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Spot', ylab = 'sigma', 
        fig_title = paste0('Spot', ' & ', 'sigma')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 3 Spot & r
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Spot, seq_r), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Spot = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Spot = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Spot', ylab = 'r', 
        fig_title = paste0('Spot', ' & ', 'r')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 4 Spot & mT
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Spot, seq_mT), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Spot = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Spot = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Spot', ylab = 'mT', 
        fig_title = paste0('Spot', ' & ', 'mT') 
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 5 Strike & sigma
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Strike, seq_sigma), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Strike = temp_grid$X[iter_i], 
                             sigma = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Strike = temp_grid$X[iter_i], 
                             sigma = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Strike', ylab = 'sigma', 
        fig_title = paste0('Strike', ' & ', 'sigma')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 6 Strike & r
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Strike, seq_r), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Strike = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Strike = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Strike', ylab = 'r', 
        fig_title = paste0('Strike', ' & ', 'r')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 7 Strike & mT
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_Strike, seq_mT), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(Strike = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(Strike = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Strike', ylab = 'mT', 
        fig_title = paste0('Strike', ' & ', 'mT')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 8 sigma & r
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_sigma, seq_r), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(sigma = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(sigma = temp_grid$X[iter_i], 
                             r = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_X = T, 
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'sigma', ylab = 'r', 
        fig_title = paste0('sigma', ' & ', 'r')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 9 sigma & mT
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_sigma, seq_mT), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(sigma = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(sigma = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_X = T, 
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'sigma', ylab = 'mT', 
        fig_title = paste0('sigma', ' & ', 'mT')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 10 r & mT
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    temp_grid <- cbind(expand.grid(seq_r, seq_mT), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid) <- c('X', 'Y', 'Z')
    for (iter_i in 1:nrow(temp_grid)) {
        temp_grid$Z[iter_i] <- 
            # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
            AmericanPut_LSMC(r = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i]) - 
            AmericanPut_BBSR(r = temp_grid$X[iter_i], 
                             mT = temp_grid$Y[iter_i])
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    }
    ContourPlots(
        temp_grid, 
        fig_subtitle = 'Difference between LSMC & BBSR for American Puts',
        scale_range = c(-0.5, 6), flag_X = T, 
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'r', ylab = 'mT', 
        fig_title = paste0('r', ' & ', 'mT')
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    tictoc::toc()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
}


rm(iter_i, temp_seq)

# Grid Comparison for numerical parameters --------------------------------

if (1) {
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # 11 n & m, LSMC only
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    tictoc::tic()
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    # numerical parameters
    seq_n <- seq(from = 10, to = 1000, by = 30)
    seq_m <- seq(from = 3, to = 365, by = 10)
    temp_grid_nm <- cbind(expand.grid(seq_n, seq_m), NA)
    # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    colnames(temp_grid_nm) <- c('X', 'Y', 'Z')
    
    list_nm <- llply(
        .data = 1:nrow(temp_grid_nm), .progress = 'time', 
        .fun = function(iter_index) {
            print(iter_index)
            iter_result <- AmericanPut_LSMC(n = temp_grid_nm[iter_index, 1], 
                                            m = temp_grid_nm[iter_index, 2])
            names(iter_result) <- iter_index
            return(iter_result)
        }
    )
    
    temp_grid_nm$Z <- unlist(list_nm)
    
    # 
    AmericanPut_LSMC()
    AmericanPut_BBSR()
    
    temp_p <- ContourPlots(
        temp_grid_nm, 
        # fig_subtitle = 'LSMC for American Puts',
        scale_range = NULL, flag_X = T, 
        flag_plot = F, flag_save_plot = T, 
        save_folder = '~/Documents/MFE/FE5222Advanced Derivatives Pricing/',
        fill_label = 'LSMC Values, in $',
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
        xlab = 'Number of simulations (n)', ylab = 'Number of time steps (m)', 
        fig_title = 'LSMC for American Puts, n & m', 
        fig_subtitle = 'BBSR = 5.80'
        # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ 
    )
    ggsave(
        filename = 'LSMC for American Puts.png', plot = temp_p, device = 'png'
        )
    
    if (flag_save_plot) {
        cat('Saving', paste0(save_folder, fig_title, '.png'), '...\n')
        ggplot2::ggsave(filename = paste0(save_folder, fig_title, '.png'), plot = temp_p)
        # , width = 7, height = 6.5)
    }
    
    # fig_subtitle = 'Spot = 100, Strike = 100,\nsigma = 0.2, r = 0.06, mT = 1'
    # 
    if (1) {
        png(filename = 'LSMC n & m.png', width = 8, height = 6, units = 'in', res = 150)
        par(oma = rep(0,4))
        filled.contour(
            x = seq_n, y = seq_m, 
            z = matrix(temp_grid_nm$Z, nrow = length(seq_n), ncol = length(seq_m), byrow = T),
            # asp = 1, 
            color.palette = cm.colors
        )
        title(main = 'LSMC for American Puts, n & m', 
              xlab = 'Number of simulations (n)', ylab = 'Number of time steps (m)')
        dev.off()
    }
    # 
    tictoc::toc()
}

# m, BBSR only
if (1) {
    seq_m_BBSR <- round(seq(from = 3, to = 365, length.out = 100))
    seq_V_BBSR <- llply(
        .data = seq_m_BBSR, 
        .fun = function(iter_m){return(AmericanPut_BBSR(m = iter_m))}
    )
    png('BBSR m.png', width = 8, height = 6, units = 'in', res = 150)
    par(oma = rep(0,4))
    plot(seq_m_BBSR, seq_V_BBSR, type = 'b', 
         main = 'American put valued by BBSR, 
             using different numbers of time steps (m)',
         xlab = 'Number of time steps (m)', ylab = 'American Put Value')
    dev.off()
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# temp_seq <- paste0('seq_', c('Spot', 'Strike', 'sigma', 'r', 'mT'))
# combn(length(temp_seq), 2)
# temp_mat <- combn(length(temp_seq), 2)
# for (iter_i in 1:ncol(temp_mat)) {
#     temp_grid <- expand.grid(
#         eval(as.symbol(temp_seq[temp_mat[1,iter_i]])), 
#         eval(as.symbol(temp_seq[temp_mat[2,iter_i]]))
#     )
#     colnames(temp_grid) <- c(temp_seq[temp_mat[1,iter_i]], temp_seq[temp_mat[2,iter_i]])
#     for (iter_j in 1:nrow(temp_grid)) {
#         temp_grid$ValueDiff <- - AmericanPut_BBSR(Spot = )
#     }
# }
# AmericanPut_BBSR(Spot = )


