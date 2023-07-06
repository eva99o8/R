
gen_Y_sim =function(Y=Y,W=W,Xpop=Xpop,Xdow=Xdow,Xden=Xden,rho_temp =rho_temp, tau_temp =tau_temp,
                    tau_spat = tau_spat, rho_spat =rho_spat,
                    b_dow = b_dow,
                    b_time = b_time,
                    b_den = b_den,delta= delta,cps=cps)
  {
  logoffset = log(Xpop)
  T = ncol(Y)
  n_region  = ncol(W)
  Xtime = log(1:T+0.01) # b_time


  cps_mat = c()
  for(cps_i in seq_along(cps)){
    cps_temp = c(rep(0,cps[cps_i]-1),rep(1,T-cps[cps_i]+1))*delta
    cps_mat = rbind(cps_mat,cps_temp)
  }
  
  ran_temp_vec = arima.sim(model = list(ar=c(rho_temp)),sd = 1/tau_temp,n=T,n.start = 100)
  ran_temp = matrix(rep(ran_temp_vec,n_region),byrow = TRUE,nrow = n_region)
  # 
  # ran_temp =matrix(0,ncol=T,nrow=n_region)
  # for(n_region_i in seq_along(1:n_region)){
  #   for(T_i in seq_along(1:T)){
  #     if(T_i==1){
  #       ran_temp[n_region_i,T_i]= rho_temp*ran_temp[n_region_i,1] + rnorm(1,0,1/tau_temp)
  #     }
  #     else{
  #       ran_temp[n_region_i,T_i]= rho_temp*ran_temp[n_region_i,T_i-1] + rnorm(1,0,1/tau_temp)}
  # 
  #   }
  # }
  
  
  D_w = diag(rowSums(W))
  Sigma0 = 1/tau_spat^2*solve(D_w - rho_spat*W)
  ran_spat_vec = rmvnorm(1,mean = rep(0,n_region),sigma = Sigma0)
  ran_spat = matrix(rep(ran_spat_vec,T),ncol = T)
  
  fixed_dow_mat = matrix(rep(c(Xdow%*%b_dow),n_region),ncol=T,byrow = TRUE)
  fixed_den_mat = matrix(rep(Xden*b_den,T),ncol=T)
  fixed_time_mat = matrix(rep(Xtime*b_time,n_region),ncol=T,byrow = TRUE)
  logoffset_mat = matrix(rep(logoffset,T),ncol=T)
  
  EY =logoffset_mat+fixed_time_mat+fixed_den_mat+fixed_dow_mat +ran_temp+ran_spat+cps_mat
  Y_sim = matrix(rpois(n=length(c(EY)),lambda = c(exp(EY))),ncol=T)
  # gen_diseaseplot(Y_sim,region_names = rownames(W)) # plot the simulated disease
  return(list(Y_sim=Y_sim))
}
