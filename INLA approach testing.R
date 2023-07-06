library(mvtnorm)
library(INLA)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
num_core = detectCores()-1
source("C:\\Users\\yzwan\\Desktop\\Research problems\\Spatial temporal outbreak detection\\Spatial temporal Bayesian hierarchical model\\Iterative Fitting(Current)\\Loading_real_data.R")
source("C:\\Users\\yzwan\\Desktop\\Research problems\\Spatial temporal outbreak detection\\Spatial temporal Bayesian hierarchical model\\Iterative Fitting(Current)\\user_defined_funcs.R")
source("C:\\Users\\yzwan\\Desktop\\Research problems\\Spatial temporal outbreak detection\\Spatial temporal Bayesian hierarchical model\\Iterative Fitting(Current)\\Y_sim_generate.R")
source("C:\\Users\\yzwan\\Desktop\\Research problems\\Spatial temporal outbreak detection\\Spatial temporal Bayesian hierarchical model\\Iterative Fitting(Current)\\simulation toy code\\user_funcs_INLA.R")

logoffset = log(Xpop)
Xtime = log(1:T+0.01) # b_time
levels(DOW) = c("Weekend",rep("Weekday",4),rep("Weekend",2))
Xdow = model.matrix(~DOW) # b_dow


Xden = log(pop_dc_MA$pop /pop_dc_MA$`dt$area`)

### All the parameters
rho_temp =.5; tau_temp = 1/0.1
tau_spat = 1/.1; rho_spat =.5
b_dow = c(-11.5,.2,.2,.2,.2,0,0)
b_dow0 = -11.5
b_dow1 = .2
b_time = 0.1
b_den = 0.2
delta = 0
# cps = rpois(n_region,lambda = (T*.7))
cps = c(123, 141, 120, 145, 142, 144, 152, 157, 134, 139, 144, 128, 137) # cps fixed for simulation


## spatial effect
D_w = diag(rowSums(W))
Sigma0 = 1/tau_spat*solve(D_w - rho_spat*W)
err_spat = rmvnorm(1,mean = rep(0,n_region),sigma = Sigma0)


## temporal effect
err_temp = arima.sim(model = list(ar=rho_temp),sd=tau_temp**-.5,n = T,n.start = 100)

## Vectorized quanties
logoffset_vec= rep(logoffset,T)
Xtime_vec = rep(Xtime,each = n_region)
Xden_vec = rep(Xden,T)
DOW_vec = rep(DOW,each =n_region)
err_temp_vec = rep(err_temp,each = n_region)
err_spat_vec =rep(err_spat,T)
lambda_vec =exp(logoffset_vec+b_time*Xtime_vec+
                  b_den*Xden_vec+
                  b_dow0+b_dow1*(as.numeric(DOW_vec)-1)+
                  err_temp_vec+err_spat_vec)
Y0 = rpois(T*n_region,lambda = lambda_vec)


data_tbl_all = tibble(
  Y = c(Y),
  day = rep(1:T,each=n_region),
  county = rep(1:n_region,T),
  Xtime = Xtime_vec,
  Xden=Xden_vec,
  logoffset =logoffset_vec,
  XDOW = DOW_vec
)

data_tbl=data_tbl_all

data_tbl_all%>%
  ggplot()+
  geom_line(aes(x=day,y=Y))+
  facet_wrap(~county)


inla_obj = inla(Y~Xtime+Xden+XDOW+f(day,model='ar1')+f(county,model=CAR.model),
                family = "poisson",data = data_tbl,E=exp(logoffset),
                control.compute = list(config=TRUE)
                      )

summary(inla_obj)

inla_obj$summary.random$day%>%
  mutate(true = err_temp)%>%
  ggplot()+
  geom_line(aes(x=ID,y=mean))+
  geom_ribbon(aes(x=ID,ymin=`0.025quant`,ymax=`0.975quant`),fill='red',alpha=0.5)+
  geom_line(aes(x=ID,y=true),linetype=2,color='blue')

rbind(c(rho_temp,tau_temp),c(inla_obj$summary.hyperpar[2,'mean'],
                             inla_obj$summary.hyperpar[1,"mean"]))

rbind(c(b_dow0,b_time,b_den,b_dow1),(inla_obj$summary.fixed[,"mean"]))
rbind(c(tau_spat,rho_spat),c(exp(inla_obj$summary.hyperpar[3,'mean']),
                             1/(1+exp(inla_obj$summary.hyperpar[4,'mean']))))
err_spat
mean(inla_obj$summary.random$county[,'mean'])

err_spat_sample=inla.posterior.sample(1000,inla_obj,selection = list(county=1:13))

cl = makeCluster(num_core)
registerDoParallel(cl)
err_spat_collect=function(i,samples){
  return(samples[[i]]$latent[,1])
}
err_spat_sampled_vec = foreach(i = 1:1000,.combine = c) %dopar%
  err_spat_collect(i,samples = err_spat_sample)
stopCluster(cl)

err_spat_tbl=tibble(
  err_spat_sampled = err_spat_sampled_vec,
  county = rep(1:n_region,1000),
  err_spat_true = rep(err_spat,1000)
)
err_spat_true_tbl = tibble(
  err_spat_true = err_spat,
  county = 1:n_region
)
ggplot()+
geom_histogram(data=err_spat_tbl,aes(x=err_spat_sampled,group=county),fill="blue")+
geom_vline(data=err_spat_true_tbl, aes(xintercept = err_spat,group=county),color='red')+
facet_wrap(~county)




  