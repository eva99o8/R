# A single simulation of spatio-temporal data

# gen_Y = function(){
#   print("sim")
# }
# 

## Xpop
## T
## n_region
## W
## Xtime
## Xdow
## Xden
## ran_temp
## ran_spat
library(mvtnorm)
library(rstan)
library(dplyr)
options(mc.cores = parallel::detectCores())
source("Loading_real_data.R")
source(("user_defined_funcs.R"))
source(("Y_sim_generate.R"))
###
logoffset = log(Xpop)
Xtime = log(1:T+0.01) # b_time
Xdow = model.matrix(~DOW) # b_dow
Xden = log(pop_dc_MA$pop /pop_dc_MA$`dt$area`)

### All the parameters
rho_temp=0.7; tau_temp = 1/0.3 
tau_spat = 1/0.6; rho_spat =.7 
b_dow = c(-13,-0.1,-0.2,-0.2,-0.1,-.4,-.4) 
b_time = 0.4 
b_den = 0.3
delta = 1
# cps = rpois(n_region,lambda = (T*.7))
cps = c(123, 141, 120, 145, 142, 144, 152, 157, 134, 139, 144, 128, 137) # cps fixed for simulation

set.seed(0614)
list_sim_Y_cp = gen_Y_sim(Y=Y,W=W,Xpop=Xpop,Xdow=Xdow,Xden=Xden,rho_temp =rho_temp, tau_temp =tau_temp,
                          tau_spat = tau_spat, rho_spat =rho_spat,
                          b_dow = b_dow,
                          b_time = b_time,
                          b_den = b_den,delta= delta,cps=cps)
gen_diseaseplot(list_sim_Y_cp$Y_sim,region_names = rownames(W))
# gen_diseaseplot(Y,region_names = rownames(W))
set.seed(0613)
list_sim_Y_cp = gen_Y_sim(Y=Y,W=W,Xpop=Xpop,Xdow=Xdow,Xden=Xden,rho_temp =rho_temp, tau_temp =tau_temp,
                          tau_spat = tau_spat, rho_spat =rho_spat,
                          b_dow = b_dow,
                          b_time = b_time,
                          b_den = b_den,delta= delta,cps=cps)
gen_diseaseplot(list_sim_Y_cp$Y_sim,region_names = rownames(W))



stan_code_path = "stan code"
spt_model_simple = stan_model(file =paste(stan_code_path,'\\CAR.stan',sep=''))
# saveRDS(spt_model_simple,file =paste(stan_code_path,'\\CAR_stan.rds',sep='') )
# spt_model_simple = readRDS(file =paste(stan_code_path,'\\CAR_stan.rds',sep=''))


time_window =30
sim_num=48
sim_save_path = paste("simulated rds\\sim_",
                      sim_num,sep='')
stan_samples_path = paste(sim_save_path,'\\stan_samples',sep='')
dir.create(sim_save_path)
dir.create(stan_samples_path)
alpha_thres = 0.05

Y_original = list_sim_Y_cp$Y_sim ## Copy of the true count
Y_adapt = list_sim_Y_cp$Y_sim ## The true count might be replaced as the loop goes on
day_last = ncol(Y_original)-time_window
# day_last =100 # For testing
T_window = time_window
Xtime_tw = log(1:(T_window+1)+0.01) # tw:time window


Y_pred_collection = Y_original[,1:T_window]
Compute_time =c()
# EY_pred_collection = array(NA,dim = c(n_iter*2,n_region,day_last))
day_b=1
day_e = day_b + time_window-1
Xdow_temp = Xdow[day_b:(day_e+1),]
Y_temp=cbind(Y_adapt[,day_b:day_e])
data0 = list(
  T = T_window,
  n = n_region,
  Y = Y_temp,
  # T0 = T_window,## For STAR model
  # n0 = n_region,## For STAR model
  # n = (T_window+1)*n_region, ## For STAT model
  # Y =c(Y_temp), ## For STAR model
  Xdow = Xdow_temp,
  Xtime = Xtime_tw,
  Xden=Xden,
  logoffset=logoffset,
  W =W,
  W_n = sum(W)/2
  # W=W_extended,## For STAR model
  # W_n = sum(W_extended)/2 ## For STAR model
)
time_temp0 = Sys.time()
samples_temp = sampling(spt_model_simple,data=data0,chains=4, iter =2000,control = list(adapt_delta =0.9,max_treedepth = 15,metric='diag_e'))
time_diff_temp =as.numeric(Sys.time()-time_temp0)
Compute_time = c(Compute_time,time_diff_temp)
cat(paste("We are making prediction on day ",day_b+T_window," now!\n",sep=''))
check_hmc_diagnostics(samples_temp)
cat(paste("The smallest eff size of the current model is ",min(summary(samples_temp)$summary[,'n_eff']),"\n",sep=''))
filename_temp = paste(sim_save_path,"\\stan_samples\\samples_day_",day_b,'.rds',sep='')
saveRDS(samples_temp,file=filename_temp)
pred_temp = extract(samples_temp,pars = c('EY_pred'))
## Collecting predictions
EY_pred_temp = exp(apply(pred_temp$EY_pred,2,mean))
Y_pred_collection = cbind(Y_pred_collection,EY_pred_temp)

### update the fitting for the next round
Y_next_true = Y_original[,day_e+1]
crit_bool_temp = ppois(Y_next_true,lambda = EY_pred_temp,lower.tail = FALSE)>alpha_thres
Y_update_vec = Y_next_true*crit_bool_temp + round(EY_pred_temp)*(1-crit_bool_temp) ## Adjusting for large p-value
Y_adapt[,day_e+1] = Y_update_vec

data_involved_collection = list(Y_raw = Y_original,Y_adapt = Y_adapt,
                                Y_pred = Y_pred_collection,
                                Compute_time = Compute_time
)
saveRDS(data_involved_collection,file=paste(sim_save_path,'\\data_involved_collection.rds',sep=''))





