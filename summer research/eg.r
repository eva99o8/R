N0 <- 1000
true.betas <- c(-1, 1, .5, -.5) 
P <- 4
X0 <- matrix(NA, nrow = N0, ncol = P)
for(i in 1:P) X0[, i] <- rnorm(N0)
sy <- 2
y0 <- rnorm(N0, mean = X0%*%true.betas, sd = sy)

summary(lm(y0 ~ -1 + X0))

library(rstan)
rstan_options(auto_write = TRUE)

compiled.model.prior <- stan_model("simple_linear_regression_NIG_prior.stan")

as <- .5
bs <- 2
vb <- 1.5
lm.data <- list(
  N_0 = N0,
  X_0 = X0,
  y_0 = y0,
  mu_beta = rep(0, P),
  lambda_0 = solve(vb * diag(P)),
  alpha0 = as,
  beta0 = bs,
  a_0 = .65
)

prior.lm <- sampling(compiled.model.prior, data = lm.data)
prior.lm
pairs(prior.lm, pars = "beta")

############
### Posterior
### will draw from the a similar process
N <- 100
X <- matrix(NA, nrow = N, ncol = P)
new.betas <- true.betas + rnorm(P, 0, sd = .2)
for(i in 1:P) X[, i] <- rnorm(N)
y <-  rnorm(N, mean = X%*%new.betas, sd = sy)

compiled.model.posterior <- stan_model("simple_linear_regression_NIG_posterior.stan")

lm.data.forposterior <- list(
  N_0 = lm.data$N0,
  X_0 = lm.data$X0,
  y_0 = lm.data$y0,
  mu_beta = lm.data$mu_beta,
  lambda_0 = lm.data$lambda_0,
  alpha0 = lm.data$alpha0,
  beta0 = lm.data$beta0,
  X = X,
  N = N,
  y = y,
  a_0 = lm.data$a_0
)
posterior.lm <- sampling(compiled.model.posterior, data = lm.data.forposterior)
posterior.lm
pairs(posterior.lm, pars = c("beta"))
