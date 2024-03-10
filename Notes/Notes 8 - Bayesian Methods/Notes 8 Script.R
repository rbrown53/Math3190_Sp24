##################
# Notes 8 Script #
##################

# First Bayesian example
# ======================
x <- c(1, 5, 3, 2, 4)
mu <- rnorm(500000, mean(x), sqrt(2.5 / 5))
quantile(mu, c(0.025, 0.975))


sigma_sq0 <- 1
mean_mu <- 1 / (1 / sigma_sq0 + 5 / 2.5) * (6 / sigma_sq0 + 5 * mean(x) / 2.5)
sigma_mu <- sqrt(1 / (1 / sigma_sq0 + 5 / 2.5))
mu <- rnorm(500000, mean_mu, sigma_mu)
quantile(mu, c(0.025, 0.975))

sigma_sq0 <- 10
mean_mu <- 1 / (1 / sigma_sq0 + 5 / 2.5) * (6 / sigma_sq0 + 5 * mean(x) / 2.5)
sigma_mu <- sqrt(1 / (1 / sigma_sq0 + 5 / 2.5))
mu <- rnorm(500000, mean_mu, sigma_mu)
quantile(mu, c(0.025, 0.975))

sigma_sq0 <- 1000
mean_mu <- 1 / (1 / sigma_sq0 + 5 / 2.5) * (6 / sigma_sq0 + 5 * mean(x) / 2.5)
sigma_mu <- sqrt(1 / (1 / sigma_sq0 + 5 / 2.5))
mu <- rnorm(500000, mean_mu, sigma_mu)
quantile(mu, c(0.025, 0.975))

# Second Bayesian Example
# =======================
# x data comes from Gam(1,1/40)
x <- c(18.9, 146.5, 12.4, 12.1, 80.4, 37.1, 37.0, 127.5, 4.2, 12.8, 25.9, 9.8)
shape <- 1 / 40 + 12 * 1
rate <- 1 + sum(x)
beta <- rgamma(100000, shape, rate)
quantile(beta, c(0.025, 0.975)) # 95% interval for beta
quantile(1 / beta, c(0.025, 0.975)) # 95% interval for the mean

hist(beta, breaks = 30, mgp = c(2.7, 1, 0), 
     cex.lab = 1.5, main = expression(paste(bold("Histogram of "), beta)), 
     xlab = expression(beta), cex.main = 1.5, freq = F)
hist(1 / beta, breaks = 30, mgp = c(2.7, 1, 0), cex.lab = 1.5, 
     main = expression(bold("Histogram of the Means")), xlab = "Mean", 
     cex.main = 1.5, freq = F)

t.test(x)$conf.int

set.seed(2024)
means <- c()
for(i in 1:10000) {
  index <- sample(1:length(x), length(x), replace = T)
  means[i] <- mean(x[index])
}
quantile(means, c(0.025, 0.975))


# Gibbs Sampling Example
# ======================
# x is generated from a normal(20,25) distribution
x <- c(
  23.3, 14.5, 19.0, 20.2, 5.4, 23.8, 21.3, 12.8, 17.6, 19.8,
  11.0, 21.5, 15.7, 22.1, 21.0, 13.7, 14.9, 13.8, 17.1, 11.3
)
nsamps <- 10000
mu <- rep(0, nsamps) # Initialize the mu vector
mu[1] <- 25 # Starting value of mu
tau <- rep(0, nsamps) # Initialize the tau vector
tau[1] <- 10 # Starting value of tau
for (i in 1:(nsamps - 1)) {
  mu[i + 1] <- rnorm(1, 1 / (1 / 100 + 20 * tau[i]) * (25 / 100 + 339.8 * tau[i]), 1 / (1 / 100 + 20 * tau[i]))
  tau[i + 1] <- rgamma(1, 15, 1 / 2 + 1 / 2 * sum((x - mu[i + 1])^2))
}
# Plots for all three parameters: mu, tau, and sigma^2
par(mfrow = c(2, 1))
plot(mu, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(mu), main = expression(paste("Sample Path of ", mu)), cex.main = 1.5)
hist(mu, breaks = 20, mgp = c(2.7, 1, 0), cex.lab = 1.5, xlab = expression(mu), main = expression(paste(bold("Histogram of "), mu)), freq = F, cex.main = 1.5)
par(mfrow = c(3, 1))
plot(tau, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(tau), main = expression(paste("Sample Path of ", tau)), cex.main = 1.5)
plot(tau[2:nsamps], type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(tau), xlab = "Index (First Value Removed)", main = expression(paste("Sample Path of ", tau, " (First Value Removed)")), cex.main = 1.5)
hist(tau[2:nsamps], xlab = expression(paste(tau, " (First Value Removed)")), mgp = c(2.7, 1, 0), cex.lab = 1.5, freq = F, main = expression(paste(bold("Histogram of "), tau)), cex.main = 1.5)
par(mfrow = c(2, 1))
plot(1 / tau[2:nsamps], type = "l", mgp = c(2.3, 1, 0), cex.lab = 1.5, ylab = expression(paste(sigma^2, " = ", 1 / tau)), xlab = "Index (First Value Removed)", main = expression(paste("Sample Path of ", sigma^2, " = ", 1 / tau)), cex.main = 1.5)
hist(1 / tau[2:nsamps], breaks = 20, xlab = expression(paste(sigma^2, " (First Value Removed)")), mgp = c(2.7, 1, 0), cex.lab = 1.5, freq = F, main = expression(paste(bold("Histogram of "), sigma^2)), cex.main = 1.5)
quantile(mu, c(0.025, 0.975))
quantile(tau, c(0.025, 0.975))
quantile(1 / tau, c(0.025, 0.975))


# =================== #
# Metropolis Sampling #
# =================== #
grocery <- c(
  1, 3, 4, 5, 5, 8, 11, 12, 15, 15, 16, 19, 21, 25, 26, 27, 30, 35, 35, 46, 50,
  55, 57, 72, 78, 85, 93, 137, 158, 212, 269
)
hist(grocery, xlab = "Grocery Prices", freq = F, breaks = 10, mgp = c(2.3, 1, 0), cex.lab = 1.5, cex.main = 1.5, main = "Histogram of Grocery Prices")
library(MASS)
gam_params <- fitdistr(grocery, "gamma")$estimate # Fits a gamma distribution to our data
lines(seq(0, 300, length = 1000), dgamma(seq(0, 300, length = 1000), gam_params[1], gam_params[2]), lwd = 3, col = "red")
legend(90, 0.01, expression(paste("Gamma Fit with ", alpha, " = 0.857 and ", beta, " = 0.016")), col = "red", lwd = 3)

# Use likelihood of gamma(77^2/sigma^2, 77/sigma)
# Prior on sigma: exponential(1/5)

sig <- 70
f1 <- (77 / sig^2)^(31 * 77^2 / sig^2) / gamma(77^2 / sig^2)^31 * prod(grocery^(77^2 / sig^2 - 1)) * exp(-sig / 50 - 77 / sig^2 * sum(grocery))
sig <- 64
f2 <- (77 / sig^2)^(31 * 77^2 / sig^2) / gamma(77^2 / sig^2)^31 * prod(grocery^(77^2 / sig^2 - 1)) * exp(-sig / 50 - 77 / sig^2 * sum(grocery))

sig <- 70
lf1 <- 31 * 77^2 / sig^2 * log(77 / sig^2) - 31 * log(gamma(77^2 / sig^2)) + (77^2 / sig^2 - 1) * sum(log(grocery)) - sig / 50 - 77 / sig^2 * sum(grocery)
sig <- 64
lf2 <- 31 * 77^2 / sig^2 * log(77 / sig^2) - 31 * log(gamma(77^2 / sig^2)) + (77^2 / sig^2 - 1) * sum(log(grocery)) - sig / 50 - 77 / sig^2 * sum(grocery)


# Metropolis with grocery data with known mu
# ==========================================
nsamps <- 10000
sigma <- rep(0, nsamps)
sigma[1] <- 64
mu <- 77
n <- length(grocery) # 31
for (i in 1:(nsamps - 1)) {
  sigma_star <- rnorm(1, sigma[i], 10)
  logf1 <- n * mu^2 / sigma_star^2 * log(mu / sigma_star^2) -
    n * log(gamma(mu^2 / sigma_star^2)) +
    (mu^2 / sigma_star^2 - 1) * sum(log(grocery)) -
    sigma_star / 50 - mu / sigma_star^2 * sum(grocery)
  logf2 <- n * mu^2 / sigma[i]^2 * log(mu / sigma[i]^2) -
    n * log(gamma(mu^2 / sigma[i]^2)) +
    (mu^2 / sigma[i]^2 - 1) * sum(log(grocery)) -
    sigma[i] / 50 - mu / sigma[i]^2 * sum(grocery)
  if (log(runif(1)) < (logf1 - logf2)) {
    sigma[i + 1] <- sigma_star
  } else {
    sigma[i + 1] <- sigma[i]
  }
}
par(mfrow = c(2, 1))
plot(sigma, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(sigma), xlab = "Index (Sample Number)", main = expression(paste("Sample Path of ", sigma)), cex.main = 1.5)
hist(sigma, breaks = 30, xlab = expression(sigma), mgp = c(2.7, 1, 0), cex.lab = 1.5, freq = F, main = expression(paste(bold("Histogram of "), sigma)), cex.main = 1.5)
quantile(sigma, c(0.025, 0.975))
acf(sigma, main = expression(paste("Autocorrelation of ", sigma, " Chain")), mgp = c(2.7, 1, 0), cex.lab = 1.3, lwd = 3)
1 + 2 * sum(abs(acf(sigma, lag.max = 100, plot = F)$acf))
quantile(sigma, c(0.025, 0.975))



# Change sigma_0 from 5 to 20 to 50 to 0.5
# ========================================
nsamps <- 10000
sigma <- rep(0, nsamps)
sigma[1] <- 64
mu <- 77
n <- length(grocery) # 31
par(mfrow = c(2, 2))
for (proposal_sig in c(5, 20, 50, 0.5)) {
  for (i in 1:(nsamps - 1)) {
    sigma_star <- rnorm(1, sigma[i], proposal_sig)
    logf1 <- n * mu^2 / sigma_star^2 * log(mu / sigma_star^2) -
      n * log(gamma(mu^2 / sigma_star^2)) +
      (mu^2 / sigma_star^2 - 1) * sum(log(grocery)) -
      sigma_star / 50 - mu / sigma_star^2 * sum(grocery)
    logf2 <- n * mu^2 / sigma[i]^2 * log(mu / sigma[i]^2) -
      n * log(gamma(mu^2 / sigma[i]^2)) +
      (mu^2 / sigma[i]^2 - 1) * sum(log(grocery)) -
      sigma[i] / 50 - mu / sigma[i]^2 * sum(grocery)
    if (log(runif(1)) < (logf1 - logf2)) {
      sigma[i + 1] <- sigma_star
    } else {
      sigma[i + 1] <- sigma[i]
    }
  }
  plot(sigma, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(sigma), xlab = "Index (Sample Number)", main = paste("Proposal SD = ", proposal_sig), cex.main = 1.5)
  print(1 + 2 * sum(abs(acf(sigma, lag.max = 100, plot = F)$acf)))
}


# Change mu from 77 to 100 to 50 to 30
# ====================================
nsamps <- 10000
sigma <- rep(0, nsamps)
sigma[1] <- 64
n <- length(grocery) # 31
par(mfrow = c(2, 2))
for (mu in c(77, 100, 50, 30)) {
  for (i in 1:(nsamps - 1)) {
    sigma_star <- rnorm(1, sigma[i], 10)
    logf1 <- n * mu^2 / sigma_star^2 * log(mu / sigma_star^2) -
      n * log(gamma(mu^2 / sigma_star^2)) +
      (mu^2 / sigma_star^2 - 1) * sum(log(grocery)) -
      sigma_star / 50 - mu / sigma_star^2 * sum(grocery)
    logf2 <- n * mu^2 / sigma[i]^2 * log(mu / sigma[i]^2) -
      n * log(gamma(mu^2 / sigma[i]^2)) +
      (mu^2 / sigma[i]^2 - 1) * sum(log(grocery)) -
      sigma[i] / 50 - mu / sigma[i]^2 * sum(grocery)
    if (log(runif(1)) < (logf1 - logf2)) {
      sigma[i + 1] <- sigma_star
    } else {
      sigma[i + 1] <- sigma[i]
    }
  }
  plot(sigma, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(sigma), xlab = "Index (Sample Number)", main = paste("True Mean = ", mu), cex.main = 1.5)
  print(quantile(sigma, c(0.025, 0.975)))
  print(mean(sigma))
  print(1 + 2 * sum(abs(acf(sigma, lag.max = 100, plot = F)$acf)))
}


# Metropolis with grocery data with unknown mu (two parameter sampling)
# =====================================================================
# sigma ~ exponential(1/50)
# mu ~ exponential(1/70)
nsamps <- 10000
sigma <- rep(0, nsamps)
sigma[1] <- 64
mu <- rep(0, nsamps)
mu[1] <- 77
n <- length(grocery) # 31
for (i in 1:(nsamps - 1)) {
  mu_star <- rnorm(1, mu[i], 8)
  sigma_star <- rnorm(1, sigma[i], 10)
  logf1 <- n * mu_star^2 / sigma_star^2 * log(mu_star / sigma_star^2) -
    n * log(gamma(mu_star^2 / sigma_star^2)) +
    (mu_star^2 / sigma_star^2 - 1) * sum(log(grocery)) -
    mu_star / 70 - sigma_star / 50 - mu_star / sigma_star^2 * sum(grocery)
  logf2 <- n * mu[i]^2 / sigma[i]^2 * log(mu[i] / sigma[i]^2) -
    n * log(gamma(mu[i]^2 / sigma[i]^2)) +
    (mu[i]^2 / sigma[i]^2 - 1) * sum(log(grocery)) -
    mu[i] / 70 - sigma[i] / 50 - mu[i] / sigma[i]^2 * sum(grocery)
  if (log(runif(1)) < (logf1 - logf2)) {
    mu[i + 1] <- mu_star
    sigma[i + 1] <- sigma_star
  } else {
    mu[i + 1] <- mu[i]
    sigma[i + 1] <- sigma[i]
  }
}
quantile(mu, c(0.025, 0.975))
mean(mu)
quantile(sigma, c(0.025, 0.975))
mean(sigma)
par(mfrow = c(2, 2))
plot(mu, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(mu), xlab = "Index (Sample Number)", main = expression(paste("Sample Path of ", mu)), cex.main = 1.5)
plot(sigma, type = "l", mgp = c(2.7, 1, 0), cex.lab = 1.5, ylab = expression(sigma), xlab = "Index (Sample Number)", main = expression(paste("Sample Path of ", sigma)), cex.main = 1.5)
hist(mu, breaks = 30, xlab = expression(mu), mgp = c(2.7, 1, 0), cex.lab = 1.5, freq = F, main = expression(paste(bold("Histogram of "), mu)), cex.main = 1.5)
hist(sigma, breaks = 30, xlab = expression(sigma), mgp = c(2.7, 1, 0), cex.lab = 1.5, freq = F, main = expression(paste(bold("Histogram of "), sigma)), cex.main = 1.5)
1 + 2 * sum(abs(acf(mu, lag.max = 100, plot = F)$acf)) # tau_int for mu
1 + 2 * sum(abs(acf(sigma, lag.max = 100, plot = F)$acf)) # tau_int for sigma


# Compare to bootstrapping
# ========================
library(boot)
b2 <- boot(grocery, function(x, i) mean(x[i]), 10000)
boot.ci(b2)
b3 <- boot(grocery, function(x, i) sd(x[i]), 10000)
boot.ci(b3)
