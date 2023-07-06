## Loading raw data
path_data= "D:\\WPI\\22summer\\Iterative Fitting(Current)\\data"
pop_dc_MA=readRDS('data\\pop_dc_MA_temp_area.rds')
W=readRDS('data\\ma_mat.rds')
# Y_raw=as.matrix(pop_dc_MA[,3:ncol(pop_dc_MA)])
Y_raw=as.matrix(pop_dc_MA[,4:ncol(pop_dc_MA)]) ## With new variable population density
Y_date=seq(as.Date('2020-01-22'),length.out = ncol(Y_raw),by='days')
colnames(Y_raw)=format(Y_date,'%m/%d/%y')
Y_raw[,c('11/26/20','12/25/20',"01/01/21")]=round(Y_raw[,c('11/27/20','12/26/20',"01/02/21")]/2)
Y_raw[,c('11/27/20','12/26/20',"01/02/21")]=round(Y_raw[,c('11/27/20','12/26/20',"01/02/21")]/2)

newdate_s=which(colnames(Y_raw)=='06/15/20')
# newdate_e=which(colnames(Y_raw)==colnames(Y_raw)[ncol(Y_raw)]) # default
newdate_e=which(colnames(Y_raw)=="12/31/20")
(newdate_e-newdate_s)
Y=Y_raw[,newdate_s:newdate_e]
gen_diseaseplot(Y,region_names = rownames(W))

## time
T = ncol(Y);
## Number of counties
n_region = nrow(Y);n_region
## population
Xpop = (pop_dc_MA$pop)


# Y_inci=c()
# denom_incident = 10^6
# for(n_region_i in seq_len(n_region)){
#   Y_inci_row_temp = round(Y[n_region_i,]/Xpop[n_region_i]*denom_incident)
#   Y_inci = rbind(Y_inci,Y_inci_row_temp)
# }

# rownames(Y_inci)=pop_dc_MA$subregion
DOW = factor(weekdays(as.Date(colnames(Y),'%m/%d/%y')),levels = c("Monday","Tuesday",   "Wednesday", "Thursday",  "Friday",    "Saturday",  "Sunday"))
