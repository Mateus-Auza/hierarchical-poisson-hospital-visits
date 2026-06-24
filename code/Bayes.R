#######################################################################
###############Pre processing the data #################################
#######################################################################

#Libraries
library(coda)
library(knitr)
library(ggplot2)
library(rjags)
library(tidyr)
library(dplyr)

set.seed(219)
# Load the dataset
data = read.table("HospitalVisits.txt", header = TRUE)

# Prepare data
n = nrow(data)
group = data$hospital
g = as.numeric(as.factor(group))
G = length(unique(g)) #Or max(g)-same thing
y = data$y
age = data$age
chronic = data$chronic
X = cbind(1, age, chronic)  # Design matrix

# Initial values
beta = c(0, 0, 0)
v = rep(1, G)
alpha = 1

# Prior parameters
a.0 = 0.01
b.0 = 0.01

#################################################
#################Question 3######################
#################################################

#Variances of the proposal
sd.beta=0.5
sd.alpha=0.5

# Monte Carlo Markov Chain setup
n_iter = 5000
samples = matrix(0, nrow=n_iter, ncol=length(beta) + G + 1)
colnames(samples) = c(paste0("beta", 0:2), paste0("v", 1:G), "alpha")

# Log-posterior for beta
log.post.beta = function(beta, v, alpha) {
  eta = X %*% beta
  mu = v[g] * exp(eta)
  lik = sum(y * eta - mu)
  lp.b = sum(-beta^2 / 200)
  return(lik + lp.b)
}

# Log-posterior for alpha
log.post.alpha = function(alpha, v) {
  if (alpha <= 0) return(-Inf)
  G = length(v)
  lp.a = (a.0 - 1) * log(alpha) - b.0 * alpha
  lp.a = lp.a + sum((alpha - 1) * log(v) - alpha * v)+ G*(alpha*log(alpha) - lgamma(alpha))
  return(lp.a)
}

# MCMC sampler
for (iter in 1:n_iter) {
  
  # 1. Update beta via Metropolis
  for (j in 1:3) {
    beta.prop = beta
    beta.prop[j] = rnorm(1, beta[j], sd = sd.beta)
    log.r = log.post.beta(beta.prop, v, alpha) - log.post.beta(beta, v, alpha)
    if (log(runif(1)) < log.r) beta = beta.prop
  }
  
  # 2. Update v via Gibbs
  for (h in 1:G) {
    idx = which(g == h)
    eta_h = X[idx, ] %*% beta
    shape = alpha + sum(y[idx])
    rate = alpha + sum(exp(eta_h))
    v[h] = rgamma(1, shape, rate)
  }
  
  # 3. Update alpha via Metropolis-Hastings on log-scale!
  log.alpha.prop = rnorm(1, log(alpha), sd = sd.alpha)
  alpha.prop = exp(log.alpha.prop)
  
  log.r = log.post.alpha(alpha.prop, v) - log.post.alpha(alpha, v) - log.alpha.prop + log(alpha)
  if (log(runif(1)) < log.r) alpha = alpha.prop
  
  # Save samples
  samples[iter, ] = c(beta, v, alpha)
}

# Basic diagnostics
par(mfrow=c(2,2))
plot(samples[, "beta0"], type="l", main="Traceplot beta0",ylab="beta0")
plot(samples[, "beta1"], type="l", main="Traceplot beta1",ylab="beta1")
plot(samples[, "beta2"], type="l", main="Traceplot beta2",ylab="beta2")
plot(samples[, "alpha"], type="l", main="Traceplot alpha",ylab="alpha")


cat("Posterior means:\n")
print(colMeans(samples))

########Computing the Geweke statistic for the chains alpha and beta

convergence_mcmc = as.mcmc(samples[, c("beta0", "beta1", "beta2","alpha")])
geweke.diag(convergence_mcmc)

##Plots of the convergence
geweke.plot(convergence_mcmc)

###Autocorrolation plots
autocorr.plot(convergence_mcmc)

# -------------------------------
# EXERCISE 4: Posterior Analysis
# -------------------------------



# (a) Numerical summaries and credible intervals for beta and alpha
burn.in = 1000  # Adjust if needed

# Extract posterior samples after burn-in
post.samples = samples[(burn.in + 1):n_iter, ]

# Separate beta and alpha
post.beta = post.samples[, c("beta0", "beta1", "beta2")]
post.alpha = post.samples[, "alpha"]

# Summary function
summary.stats = function(samples) {
  c(
    Mean    = round(mean(samples), 4),
    SD      = round(sd(samples), 4),
    Median  = round(median(samples), 4),
    Min     = round(min(samples), 4),
    Max     = round(max(samples), 4),
    `2.5%`  = round(quantile(samples, 0.025), 4),
    `97.5%` = round(quantile(samples, 0.975), 4)
  )
}

# Combine summaries
Matrix = rbind(
  summary.stats(post.beta[, 1]),
  summary.stats(post.beta[, 2]),
  summary.stats(post.beta[, 3]),
  summary.stats(post.alpha)
)

rownames(Matrix) = c("Beta0", "Beta1", "Beta2", "Alpha")

# Print nicely
kable(Matrix, format = "markdown", caption = "Posterior Summaries for β and α")

# -------------------------------
# (b) Graphical summaries for v_g
# -------------------------------

# Extract v_g samples (columns 4 to 4 + G - 1)
post.v = post.samples[, paste0("v", 1:G)]

# Traceplots for v_g (30 hospitals)
par(mfrow = c(5, 6), mar = c(2, 2, 2, 1))
for (g in 1:G) {
  plot(post.v[, g], type = "l", main = paste0("v", g), xaxt = 'n', yaxt = 'n')
}

# (b) Evaluate the posterior predictive performance of the fitted model

# Sample S posterior draws for prediction
S = 1000
sample.indices = sample(1:nrow(post.beta), S)

# Matrix to store synthetic outcomes
y.simul = matrix(0, nrow = S, ncol = length(y))

# Generate posterior predictive draws
for (s in 1:S) {
  beta.s = post.beta[sample.indices[s], ]
  v.s = post.v[sample.indices[s], ][g]  # g = hospital ID index from earlier
  eta.s = beta.s[1] + beta.s[2] * age + beta.s[3] * chronic
  y.simul[s, ] = rpois(length(y), lambda = v.s * exp(eta.s))
}

########################################
######################Visualisation of the posterior density compared with the true model##################################################
########################################

# Compute average prediction per observation
y.pred.mean = colMeans(y.simul)

# Create a combined data frame
dens.df = bind_rows(
  data.frame(value = y, type = "Observed y"),
  data.frame(value = y.pred.mean, type = "Prediction mean y")
)

# Plot with ggplot2
ggplot(dens.df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Density: Observed data vs Predictions",
    x = "y",
    y = "Density",
    fill = "Legend"
  ) +
  scale_fill_manual(values = c("Observed y" = "red", "Prediction mean y" = "blue"))


### Numerical comparison
# Compare mean
cat("Observed mean of y:", mean(y), "\n")
cat("Posterior predictive mean of y:", mean(apply(y.simul, 1, mean)), "\n")

# Compare variance
cat("Observed SD of y:", sd(y), "\n")
cat("Average posterior predictive SD:", mean(apply(y.simul, 1, sd)), "\n")

#From these numerical and graphical summaries we can see that the model tends to locate quite accurately the mean of the observed distribution but underestimates the variability of the observed y.


#####################################################
################Computing the DIC####################
#####################################################
#Now we compare the hierachical model (with a hyperparameter v.g) and a non hierarchical model

  
###DIC for the hierarchical model ###
deviance_samples = numeric(nrow(post.beta))
for (s in 1:nrow(post.beta)) {
  beta_s = post.beta[s, ]
  v_s = post.v[s, ]
  mu_s = exp(X %*% beta_s)*v_s[group]
  deviance_samples[s] <- -2 * sum(dpois(y, mu_s, log = TRUE))
}

# DIC
D_bar = mean(deviance_samples)
beta_mean = colMeans(post.beta)
v_mean = colMeans(post.v)
mu_mean = exp(X %*% beta_mean)*v_mean[group]
D_hat = -2 * sum(dpois(y, mu_mean, log = TRUE))
pD = D_bar - D_hat
DIC = D_bar + pD

# Modèle Poisson sans effets aléatoires
beta_mean_simple = colMeans(post.beta)  # Moyenne a posteriori des beta
mu_mean_simple = exp(X %*% beta_mean_simple)
deviance_simple = -2 * sum(dpois(y, mu_mean_simple, log = TRUE))
pD_simple = 3  # Nombre de paramètres (beta0, beta1, beta2)
DIC_simple = deviance_simple + pD_simple

#--------------------------------------------
# Comparison of the two models 
#--------------------------------------------
# Results for the hierarchical model
cat("DIC (hiérarchique) =", round(DIC, 1), "\npD =", round(pD, 1))
# Results for the non hierarchical model
cat("DIC (simple) =", round(DIC_simple, 1), "\npD =", pD_simple)


##############################Simple model
# Sample S posterior draws for prediction
S = 1000
sample.indices = sample(1:nrow(post.beta), S)

# Matrix to store synthetic outcomes
y.simul.simple = matrix(0, nrow = S, ncol = length(y))

# Generate posterior predictive draws
for (s in 1:S) {
  beta.s = post.beta[sample.indices[s], ]
  eta.s = beta.s[1] + beta.s[2] * age + beta.s[3] * chronic
  y.simul.simple[s, ] = rpois(length(y), lambda =exp(eta.s))
}

# Compute average prediction per observation
y.pred.simple.mean = colMeans(y.simul.simple)

# Create a combined data frame
dens.df = bind_rows(
  data.frame(value = y, type = "Observed y"),
  data.frame(value = y.pred.mean, type = "Prediction mean y"),
  data.frame(value= y.pred.simple.mean, type="Prediction simple mean y")
)

# Plot with ggplot2
ggplot(dens.df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Density: Observed data vs Predictions",
    x = "y",
    y = "Density",
    fill = "Legend"
  ) +
  scale_fill_manual(values = c("Observed y" = "red", "Prediction mean y" = "blue", "Prediction simple mean y"="green"))



###############################################
## Exercice 5 #################################
###############################################

#(a) Do the same as exercice 3 but using JAGS (Just Another Gibs Sampler)

#Putting the model into a string model
model.string = "
model {
  for (i in 1:N) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- log(v[hospital[i]]) + beta0 + beta1 * age[i] + beta2 * chronic[i]
  }

  for (g in 1:G) {
    v[g] ~ dgamma(alpha, alpha)
  }

  # Priors
  beta0 ~ dnorm(0, 10)
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)
  alpha ~ dgamma(a0, b0)
}
"

#Using JAGS
data.jags = list(
  y = y,
  age = age,
  chronic = chronic,
  hospital = data$hospital,
  N = nrow(data),
  G = G,
  a0 = a.0,
  b0 = b.0
)
initial.values = function() {
  list(
    beta0 = 0, beta1 = 0, beta2 = 0,
    alpha = 1,
    v = rep(1, data.jags$G)
  )
}
parameters=c("beta0","beta1","beta2","alpha","v")

#Model that stores everything-Interface between R and JAGS 

model = jags.model(textConnection(model.string),
                   data = data.jags,
                   inits = initial.values,
                   n.chains = 1,
                   n.adapt = 1000)

#We update the model such that the model forgets the iterations before the burn-in
update(model,burn.in)

#We take the samples
samples.jags=coda.samples(model, variable.names = parameters,n_iter)

#Storing this in a matrix

samples.mat=as.matrix(samples.jags)

post.beta.jags = samples.mat[, c("beta0", "beta1", "beta2")]
post.alpha.jags = samples.mat[, "alpha"]
post.v.jags = samples.mat[, grep("^v\\[", colnames(samples.mat))]

#(b) Compare the posterior summaries obtained in 3 #####################

###Post.beta vs Post.beta.jags #################

comparison.beta = data.frame(
  Parameter = c("beta0", "beta1", "beta2"),
  Mean.R = colMeans(post.beta),
  SD.R   = apply(post.beta, 2, sd),
  Mean.JAGS   = colMeans(post.beta.jags),
  SD.JAGS     = apply(post.beta.jags, 2, sd)
)
kable(comparison.beta, format = "markdown", caption = "Comparison between the beta parameters in R and JAGS")

###Alpha ########################################

comparison.alpha=data.frame(
  Parameter = "alpha",
  Mean.R = mean(post.alpha),
  SD.R   = sd(post.alpha),
  Mean.JAGS   = mean(post.alpha.jags),
  SD.JAGS     = sd(post.alpha.jags)
)
kable(comparison.alpha, format = "markdown", caption = "Comparison between the alpha in R and JAGS")


### Comparison v ################################

comparison.v= data.frame(
  Mean.R = colMeans(post.v),
  SD.R   = apply(post.v, 2, sd),
  Mean.JAGS   = colMeans(post.v.jags),
  SD.JAGS     = apply(post.v.jags, 2, sd)
)
kable(comparison.v, format = "markdown", caption = "Comparison between the v parameters in R and JAGS")

############################
#####Visual representation #
############################

###Beta
dif.beta = rbind(
  data.frame(value = as.vector(post.beta), 
             parameter = rep(c("beta0", "beta1", "beta2"), each = nrow(post.beta)),
             model = "Basic R"),
  data.frame(value = as.vector(post.beta.jags), 
             parameter = rep(c("beta0", "beta1", "beta2"), each = nrow(post.beta.jags)),
             model = "JAGS")
)

ggplot(dif.beta, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior distributions of beta parameters",
       x = "Value", y = "Density")

###Alpha ##############################################3

dif.alpha = rbind(
  data.frame(value = as.vector(post.alpha), 
             parameter = "alpha",
             model = "Basic R"),
  data.frame(value = as.vector(post.alpha.jags), 
             parameter = "alpha",
             model = "JAGS")
)


ggplot(dif.alpha, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Posterior distributions of alpha",
       x = "Value", y = "Density")


###v ###########################################


# Add group index
comparison.v = comparison.v %>%
  mutate(Group = 1:n())


v.small = bind_rows(
  comparison.v %>%
    transmute(
      Group,
      Model = "R",
      Mean = Mean.R,
      Lower = Mean.R - 1.96 * SD.R,
      Upper = Mean.R + 1.96 * SD.R
    ),
  comparison.v %>%
    transmute(
      Group,
      Model = "JAGS",
      Mean = Mean.JAGS,
      Lower = Mean.JAGS - 1.96 * SD.JAGS,
      Upper = Mean.JAGS + 1.96 * SD.JAGS
    )
)

ggplot(v.small, aes(x = Group, y = Mean, color = Model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  theme_minimal() +
  labs(title = "Posterior Means and 95% CIs for v by Model",
       x = "Group (hospital)",
       y = "Mean ± 95% CI")

###########################################################
################ Exercice 6- Robustness to the prior choice
##########################################################

#Replacing the prior

model.string2 = "
model {
  for (i in 1:N) {
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) <- log(v[hospital[i]]) + beta0 + beta1 * age[i] + beta2 * chronic[i]
  }

  for (g in 1:G) {
    log.v[g] ~ dnorm(0, tau)
    v[g] <- exp(log.v[g])
  }

  # Priors
  beta0 ~ dnorm(0, 10)
  beta1 ~ dnorm(0, 10)
  beta2 ~ dnorm(0, 10)
  tau ~ dgamma(a0, b0)
}
"
initial.values2 = function() {
  list(beta0 = 0,beta1 = 0,beta2 = 0,tau = 1,log.v = rep(0, G))}

parameters2=c("beta0","beta1","beta2","tau","v")

#Model that stores everything-Interface between R and JAGS 

model2 = jags.model(textConnection(model.string2),
                    data = data.jags,
                    inits = initial.values2,
                    n.chains = 1,
                    n.adapt = 1000)

update(model2,burn.in)

#We take the samples
samples.jags2=coda.samples(model2, variable.names = parameters2,n_iter)

#Storing this in a matrix
samples.mat2=as.matrix(samples.jags2)

post.beta.jags2 = samples.mat2[, c("beta0", "beta1", "beta2")]
post.tau.jags = samples.mat2[, "tau"]
samples.mat2
post.v.jags2 = samples.mat2[, grep("^v\\[", colnames(samples.mat2))]

##Differences between the two JAGS models

###Post.beta.jags vs Post.beta.jags2

comparison.beta.jags = data.frame(
  Parameter = c("beta0", "beta1", "beta2"),
  Mean.JAGS = colMeans(post.beta.jags),
  SD.JAGS   = apply(post.beta.jags, 2, sd),
  Mean.JAGS2   = colMeans(post.beta.jags2),
  SD.JAGS2     = apply(post.beta.jags2, 2, sd)
)
kable(comparison.beta.jags, format = "markdown", caption = "Comparison between the v parameters in R and JAGS")


### Comparison v

comparison.v.jags= data.frame(
  Group = 1:ncol(post.v.jags),
  Mean.JAGS = colMeans(post.v.jags),
  SD.JAGS   = apply(post.v.jags, 2, sd),
  Mean.JAGS2   = colMeans(post.v.jags2),
  SD.JAGS2     = apply(post.v.jags2, 2, sd)
)
kable(comparison.v.jags, format = "markdown", caption = "Comparison between the v parameters in R and JAGS")


###Visualisation of the differences in the parameters

##Beta
dif.beta.jags = rbind(
  data.frame(value = as.vector(post.beta.jags), 
             parameter = rep(c("beta0", "beta1", "beta2"), each = nrow(post.beta.jags)),
             model = "JAGS"),
  data.frame(value = as.vector(post.beta.jags2), 
             parameter = rep(c("beta0", "beta1", "beta2"), each = nrow(post.beta.jags2)),
             model = "JAGS2")
)

ggplot(dif.beta.jags, aes(x = value, fill = model)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Posterior distributions of beta parameters",
       x = "Value", y = "Density")


##Vg

#We add the confidence intervals to our data.framme of differences
comparison.v.jags$Lower.JAGS = comparison.v.jags$Mean.JAGS - 1.96 * comparison.v.jags$SD.JAGS
comparison.v.jags$Upper.JAGS = comparison.v.jags$Mean.JAGS + 1.96 * comparison.v.jags$SD.JAGS
comparison.v.jags$Lower.JAGS2 = comparison.v.jags$Mean.JAGS2 - 1.96 * comparison.v.jags$SD.JAGS2
comparison.v.jags$Upper.JAGS2 = comparison.v.jags$Mean.JAGS2 + 1.96 * comparison.v.jags$SD.JAGS2


v.long = bind_rows(
  comparison.v.jags %>%
    mutate(Group = 1:n()) %>%
    transmute(
      Group,
      Model = "JAGS",
      Mean = Mean.JAGS,
      Lower = Mean.JAGS - 1.96 * SD.JAGS,
      Upper = Mean.JAGS + 1.96 * SD.JAGS
    ),
  comparison.v.jags %>%
    mutate(Group = 1:n()) %>%
    transmute(
      Group,
      Model = "JAGS2",
      Mean = Mean.JAGS2,
      Lower = Mean.JAGS2 - 1.96 * SD.JAGS2,
      Upper = Mean.JAGS2 + 1.96 * SD.JAGS2
    )
)

ggplot(v.long, aes(x = Group, y = Mean, color = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  theme_minimal() +
  labs(title = "Posterior Means and 95% CIs for v.g by Model",
       x = "Group (hospital)",
       y = "Mean ± 95% CI")










