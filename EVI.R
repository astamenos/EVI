# Importing libraries
library(tidyverse)
library(rjags)
library(splines)

# Loading the data
df <- read.csv('EVI_Data.csv')
df['date'] <- as.Date(with(df,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")
df['t'] <- as.integer(df$date - min(df$date))
n <- length(df$t)
t <- df$t
Y <- df$EVI

################################### Model 1 ####################################

# Model definition
data1 <- list(n=n, t=t, Y=Y, pi=pi)
M1_str <- textConnection("model{
  # Likelihoods
  for(i in 1:n) {
    Y[i] ~ dbeta(r*mu[i], r*(1-mu[i]))
    logit(mu[i]) <- alpha + beta*sin((2*pi/365)*(t[i]-delta))
    like[i] <- dbeta(Y[i], r*mu[i], r*(1-mu[i]))
  }
  
  # Priors
  alpha ~ dnorm(0, 0.1)
  beta ~ dnorm(0, 0.1)
  delta ~ dunif(0, 365)
  r ~ dgamma(0.1, 0.1)
}")

# Model compilation
init1 <- list(alpha=0,beta=1.8, delta=280)
M1 <- jags.model(M1_str, data = data1, n.chains = 2, inits = init1, quiet = TRUE)
update(M1, 20000, progress.bar = 'none')
samples1 <- coda.samples(M1, variable.names = c('mu'),
                         n.iter = 20000, progress.bar = 'none')


# DIC calculation
DIC <- dic.samples(M1, n.iter = 20000, progress.bar = 'none')

# WAIC calculation
samples1_like <- coda.samples(M1, n.iter = 20000, 
                              variable.names = c('like'), progress.bar = 'none')
like1 <- rbind(samples1_like[[1]], samples1_like[[2]])
f_bar1 <- colMeans(like1)
p_W1 <- sum(apply(log(like1), 2, var))
WAIC1 <- -2*sum(log(f_bar1)) + 2*p_W1

# Summaries and convergence diagnostics
param_samples <- coda.samples(M1, variable.names = c('alpha', 'beta', 'delta', 'r'), 
                              n.iter = 20000, progress.bar = 'none')
gelman.diag(param_samples)
geweke.diag(param_samples)
plot(param_samples)
plot(samples1)
summary(samples1)

# Uncertainty Quantification of Mu by Day
sum1 <- summary(samples1[[1]])
quants <- sum1$quantiles
mu_t_025 <- quants[,1]
mu_t_median <- quants[,3]
mu_t_975 <- quants[,5]
df$date <- as.Date(df$date)
check <- data.frame(t=t[1:n], mu_t_025, mu_t_median, mu_t_975, Y=Y[1:n])
p <- ggplot(data = check) + 
  geom_line(aes(x=t, y=mu_t_median), size = 0.75) + geom_ribbon(data = check, aes(x=t, ymin = mu_t_025,ymax = mu_t_975), fill = 'steelblue', alpha = 0.5) + 
  geom_point(aes(x=t, y=Y), color = 'black', alpha = 0.5) + labs(x='Days since 06/10/1984', y='EVI')
p
################################### Model 2 ####################################

# Model definition
J <- 50
B <- bs(t, J) 
data2 <- list(n=n, J=J, B=B, Y=Y)
M2_str <- textConnection("model{
  # Likelihoods
  for(i in 1:n) {
    Y[i] ~ dbeta(r*mu[i], r*(1-mu[i]))
    logit(mu[i]) = alpha + inprod(B[i,], beta[])
    like[i] <- dbeta(Y[i], r*mu[i], r*(1-mu[i]))
  }
  
  # Priors
  for(j in 1:J) {
    beta[j] ~ dnorm(0, tau_e*tau_b)
  }
  tau_e~ dgamma(0.1, 0.1)
  tau_b ~ dgamma(0.1, 0.1)
  r ~ dgamma(0.1, 0.1)
  alpha ~ dnorm(0, 0.01)
}")

# Model compilation
M2 <- jags.model(M2_str, data = data2, n.chains = 2, quiet = TRUE)
update(M2, 5000, progress.bar = 'none')
samples2 <- coda.samples(M2, variable.names = c('mu'), 
                         n.iter = 5000, progress.bar = 'none')

# Summaries and diagnostics
plot(samples2)
summary(samples2)
param_samples <- coda.samples(M2, variable.names = c('beta'), 
                              n.iter = 10000, progress.bar = 'none')
gelman.diag(param_samples)
geweke.diag(param_samples)
plot(param_samples)

# DIC calculation
DIC <- dic.samples(M2, n.iter = 5000, progress.bar = 'none')

# WAIC calculation
samples2_like <- coda.samples(M2, n.iter = 5000, 
                              variable.names = c('like'), progress.bar = 'none')
like2 <- rbind(samples2_like[[1]], samples2_like[[2]])
f_bar2 <- colMeans(like2)
p_W2 <- sum(apply(log(like2), 2, var))
WAIC2 <- -2*sum(log(f_bar2)) + 2*p_W2

# Uncertainty Quantification of Mu by Day
sum2 <- summary(samples2[[1]])
quants <- sum2$quantiles
mu_t_025 <- quants[,1]
mu_t_median <- quants[,3]
mu_t_975 <- quants[,5]
df$date <- as.Date(df$date)
check <- data.frame(t=t[1:n], mu_t_025, mu_t_median, mu_t_975, Y=Y[1:n])
p <- ggplot(data = check) + geom_point(aes(x=t, y=Y), color = 'black', alpha = 0.5) + labs(x='Days since 06/10/1984', y='EVI')
p <- p + geom_line(aes(x=t, y=mu_t_median), size = 0.75) + geom_ribbon(data = check, aes(x=t, ymin = mu_t_025,ymax = mu_t_975), fill = 'steelblue', alpha = 0.5)
p
################################################################################
# Preparing the data
df <- df %>% filter(DOY < 218)
X <- df[,c(1,4,5)]
DOY <- X$DOY
EVI <- X$EVI
year <- X$Year-min(X$Year) + 1
n <- length(EVI)
m <- length(unique(year))

data3 <- list(year=year, DOY=DOY, EVI=EVI, n=n, m=m)

M3_str <- textConnection("model{
    # Likelihoods
    for(i in 1:n) {
      EVI[i] ~ dnorm(mu[i], tau[year[i]])
      mu[i] = 0.2 + beta[year[i]]*exp(-(1/gamma[year[i]])*(DOY[i]-delta[year[i]])^2)
    }
    
    for(j in 1:m) {
      delta[j] ~ dunif(sqrt(-gamma[j]*log(0.3/beta[j])), 365 + sqrt(-gamma[j]*log(0.3/beta[j])))
      beta[j] ~ dunif(0.3, 0.8)
      gamma[j] ~ dgamma(0.1, 0.1)
      tau[j] ~ dgamma(0.1, 0.1)
      gut[j] <- delta[j]-sqrt(-gamma[j]*log(0.3/beta[j]))
    }
    
  }")

M3_str <- textConnection("model{
    # Likelihoods
    for(i in 1:n) {
      EVI[i] ~ dnorm(mu[i], tau[year[i]])
      mu[i] = alpha[year[i]]/(1 + exp(-gamma[year[i]]*(DOY[i]-delta[year[i]])))
    }
    
    for(j in 1:m) {
      alpha[j] ~ dunif(0.50, 1)
      gamma[j] ~ dbeta(1, 10)
      delta[j] ~ dunif(0 + (1/gamma[j])*log(2*alpha[j]-1), 365 + (1/gamma[j])*log(2*alpha[j]-1))
      tau[j] ~ dgamma(0.1, 0.1)
      gut[j] <- delta[j]-(1/gamma[j])*log(2*alpha[j]-1)
    }
    
  }")



# Model compilation
#init3 <- list(gamma = 5000)
M3 <- jags.model(M3_str, data = data3, n.chains = 1, quiet = TRUE)
update(M3, 70000, progress.bar = 'none')
samples <- coda.samples(M3, 
                        variable.names = c('gut'), 
                        n.iter = 10000, progress.bar = 'none')  
gut_sum <- as.data.frame(summary(samples)$quantiles)
gut_sum <- cbind(gut_sum, unique(X$Year))
colnames(gut_sum) <- c('q025', 'q25', 'q50', 'q75', 'q975', 'year')

samples2 <- coda.samples(M3, 
                        variable.names = c('mu'), 
                        n.iter = 10000, progress.bar = 'none')
# Uncertainty Quantification of Mu by Day
sum2 <- summary(samples2[[1]])
quants <- sum2$quantiles
mu_t_025 <- quants[,1]
mu_t_median <- quants[,3]
mu_t_975 <- quants[,5]
df$date <- as.Date(df$date)
check <- data.frame(t=t[1:n], mu_t_025, mu_t_median, mu_t_975, Y=Y[1:n])
p <- ggplot(data = check) + geom_point(aes(x=t, y=Y), color = 'black', alpha = 0.5) + labs(x='Days since 06/10/1984', y='EVI')
p <- p + geom_line(aes(x=t, y=mu_t_median), size = 0.75) + geom_ribbon(data = check, aes(x=t, ymin = mu_t_025,ymax = mu_t_975), fill = 'steelblue', alpha = 0.5)
p

gut_plot <- ggplot(gut_sum, aes(year)) + 
  geom_line(aes(y=q50, color = 'Posterior Median'), size = 1.25) + 
  geom_ribbon(data = gut_sum, aes(x=year, ymin = q025, ymax = q975, fill = '95% CI'), alpha = 0.5) + 
  geom_ribbon(data = gut_sum, aes(x=year, ymin = q25,ymax = q75, fill = '50% CI'), alpha = 0.3) +
  scale_fill_manual(name = '', values = c('50% CI' = 'black', '95% CI' = 'steelblue')) +
  scale_colour_manual(name = '', values = c('Posterior Median' = 'black'))
gut_plot <- gut_plot + labs(x = 'Year', y = 'GUT')
gut_plot
################################ GUT Analysis ##################################

# Preparing the data
X <- df[,c(1,4,5)]
temp <- split(X, f=X$Year)
data2 <- list(temp=temp)

# Model definition
output <- lapply(temp, function(x) {
  DOY <- x$DOY
  EVI <- x$EVI
  n <- length(EVI)
  data3 <- list(DOY=DOY, EVI=EVI, n=n)
  
  M3_str <- textConnection("model{
    # Likelihoods
    for(i in 1:n) {
      EVI[i] ~ dnorm(mu[i], tau)
      mu[i] = exp(-(1/gamma)*(DOY[i]-delta)^2)
    }
    
    # Priors
    gamma ~ dgamma(0.1, 0.1)
    delta ~ dunif(0, 365)
    tau ~ dgamma(0.1, 0.1)

  }")
  
  # Model compilation
  init3 <- list(gamma = 5000)
  M3 <- jags.model(M3_str, data = data3, n.chains = 1, inits = init3, quiet = TRUE)
  update(M3, 50000, progress.bar = 'none')
  samples <- coda.samples(M3, 
                          variable.names = c('gamma', 'delta'), 
                          n.iter = 50000, progress.bar = 'none')  
  return(samples)
})

# Calculating GUT
gut <- lapply(output, function(x) {
  t <- x[[1]][,1]-sqrt(x[[1]][,2]*log(2))
})

# Comparison of GUT 95% credible intervals for each year
gut_sum <- lapply(gut, function(x) round(summary(x)$quantiles, 4))
gut_sum <- data.frame(matrix(unlist(gut_sum), nrow=length(gut_sum), byrow=TRUE))
colnames(gut_sum) <- c('q025', 'q25', 'q50', 'q75', 'q975')
gut_sum <- cbind(year = unique(df$Year), gut_sum)
gut_plot <- ggplot(gut_sum, aes(year)) + 
  geom_line(aes(y=q50, color = 'Posterior Median'), size = 1.25) + 
  geom_ribbon(data = gut_sum, aes(x=year, ymin = q025, ymax = q975, fill = '95% CI'), alpha = 0.5) + 
  geom_ribbon(data = gut_sum, aes(x=year, ymin = q25,ymax = q75, fill = '50% CI'), alpha = 0.3) +
  scale_fill_manual(name = '', values = c('50% CI' = 'black', '95% CI' = 'steelblue')) +
  scale_colour_manual(name = '', values = c('Posterior Median' = 'black'))
gut_plot <- gut_plot + labs(x = 'Year', y = 'GUT')
gut_plot

# Parsing model results
gut_sum[gut_sum$q975==max(gut_sum$q975),]
probs <- lapply(gut, function(x) mean(x>=355.45))
probs <- data.frame(matrix(unlist(probs), nrow=length(probs), byrow=TRUE))
colnames(probs) <- 'tail_prob'
probs['year'] <- seq(1984, 2019)
probs['ratio'] <- probs$tail_prob/min(probs$tail_prob)

# Time-trend visualization
gut_comp <- ggplot(probs, aes(year)) + 
  geom_line(aes(y = ratio, color = "Ratio of P(GUT > 355.45) 
  to Min P(GUT > 355.45)"), size = 1.25) + 
  geom_hline(yintercept = 1, size = 1) +
  scale_colour_manual(name = '', values = c("Ratio of P(GUT > 355.45) 
  to Min P(GUT > 355.45)" = 'steelblue')) +
  ylim(0.95, 1.25)
gut_comp
################################################################################
