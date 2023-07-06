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



inla_obj = inla(Y~Xtime+Xden+XDOW+f(day,model='ar1',
                                    hyper = list(prec=list(param=c(1,0.0001))))+
                  f(county,model=CAR.model),
                family = "poisson",data = data_tbl,E=exp(logoffset),
                control.compute = list(config=TRUE),
                control.predictor = list(compute=TRUE),
                control.fixed = list(mean=list(Xtime=0,Xden=0,XDOW=0),
                                     prec=(list(Xtime=1,Xden=1,XDOW=1)))
                  
)

summary(inla_obj)


fitted_val_tbl=(inla_obj$summary.fitted.values)
fitted_val_tbl%>%
  mutate(lam_fitted = exp(rep((logoffset),T))*mean,
         lam_fitted_lower = exp(rep((logoffset),T))*`0.025quant`,
         lam_fitted_upper = exp(rep((logoffset),T))*`0.975quant`,
         Y_fitted_upper =qpois(.95,lam_fitted_upper) ,
         Y=c(Y),
         county= rep(1:n_region,T),
         day = rep(1:T,each = n_region))%>%
  ggplot()+
  # geom_line(aes(x=day,y=lam_fitted,group=county),linetype=2,color='green')+
  geom_ribbon(aes(x=day,ymin=0,ymax=Y_fitted_upper,group=county),
              fill="green",alpha=.9)+
  geom_line(aes(x=day,y=Y,group=county),color="red",alpha=.5)+
  facet_wrap(~county)
  
  
  
