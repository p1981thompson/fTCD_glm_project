#=========================================================================================================#
# Simulation code - gamma HRF (fastGLS)
#=========================================================================================================#

# 14-JAN-2020 

#Testing out a different method of simulating the time series data.

#=========================================================================================================#
# https://www.flutterbys.com.au/stats/tut/tut8.3a.html - reference for simulation code.
#=========================================================================================================#

#Packages Required

library(tidyverse)
library(nlme)
library(cladoRcpp) # turbo charges the 'convolve' function used to convolve the stimulus and HRF (super slow previously). Based on the function 'fmri.stimulus' from the fmri package.
library(Rcpp)
library(RcppArmadillo)

set.seed(54321)

N = 100

N = 100

.gammaHRF <- function(t, par = NULL) {
  th <- 0.242 * par[1]
  1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
}


obs_per_sec <- 2

stim1_s <- rep(c(rep(0,5*obs_per_sec),rep(c(1,0),times=c(15*obs_per_sec,35*obs_per_sec))),15) #insert from other programs simulted stim and convolve with gamma. 
stim2_s <- rep(c(rep(0,20*obs_per_sec),rep(c(1,0),times=c(5*obs_per_sec,30*obs_per_sec))),15) #same as above

par <- floor((30/28)*4)

y <- .gammaHRF(0:30, par) 

stim1_c <- rcpp_convolve(a=stim1_s, b=rev(y))
stim2_c <- rcpp_convolve(a=stim2_s, b=rev(y))

t1 <- seq(0.5,length(stim1_c)/obs_per_sec,by=0.5)

stim1 <- c(stim1_c,stim1_c)
stim2 <- c(stim2_c,stim2_c)
t <- c(t1,t1)

N <- length(stim1_c)

signal <- rep(c(0,1),each=N)
stim1_signal <- signal*stim1

b0 <- 100
b1 <- 2.26
b2 <- -2.18
b3 <- 0.00002
b4 <- -0.00000004
b5 <- -0.000000005
b6 <- -0.5
b7 <- -0.5

sigma <- 2

rho <- 0.7

mod1<-mod2<-mod3<-se1<-se2<-se3<-data.frame(matrix(NA,nrow=100,ncol=8))
names(mod1)<-names(mod2)<-names(mod3)<-names(se1)<-names(se2)<-names(se3)<-c('(Intercept)', 'stim1', 'stim2','t', 'I(t^2)', 'I(t^3)', 'signal1', 'stim1_signal')


for(i in 1:100)
{
  eps1<-eps2<-vector(mode='numeric',length=length(stim1_c))
  ## define a constructor for a first-order
  ## correlation structure
  ar1 <- corAR1(form = ~t|signal, value = rho)
  ## initialize this constructor against our data
  AR1 <- Initialize(ar1, data = data.frame(t,signal))
  ## generate a correlation matrix
  V <- corMatrix(AR1)
  ## Cholesky factorization of V
  
  V1<-make.positive.definite(V[[1]])
  V2<-make.positive.definite(V[[2]])
  
  Cv1 <- base::chol(V1)
  Cv2 <- base::chol(V2)
  
  
  ## simulate AR1 errors
  eps[1:length(stim1_c)] <- t(Cv1) %*% rnorm(length(stim1_c), 0, sigma)  # cov(e) = V * sig^2
  eps[(1+length(stim1_c)):(length(stim1_c)+length(stim1_c))] <- t(Cv2) %*% rnorm(length(stim1_c), 0, sigma) 
  ## generate response
  
  y <- b0 + b1*stim1 + b2*stim2 + b3*t + b4*(t^2) + b5*(t^3) + b6*signal + b7*stim1_signal + eps
  sim_data = data.frame(y = y, t=t, signal = signal,stim1=stim1,stim2=stim2,stim1_signal=stim1_signal)
  sim_data$signal<-as.factor(sim_data$signal)
  #ggplot(sim_data,aes(x=t,y=y))+geom_point(aes(colour=signal))+theme_bw()+theme(legend.position = 'top') + xlab('Time')+ylab('Blood flow velocity')
  
  #=========================================================================================================#
  
  #---------------------------------------------------------------------------------------------------------------#
  # Generalised linear model for each individual.
  myfit<-glm(y ~ stim1+stim2+t+I(t^2)+I(t^3)+signal+stim1_signal,data=sim_data,family=gaussian())
  
  
  #---------------------------------------------------------------------------------------------------------------#
  # 28-10-2019 added gls addition via:https://stackoverflow.com/questions/52741610/matrix-inversion-and-transpose-r-vs-c
  
  # Create Armadillo function
  Rcpp::cppFunction(depends = "RcppArmadillo", code = '
Rcpp::List fastgls(const arma::mat &X, const arma::colvec &y, const arma::mat &Sigma) {
   arma::mat Sigi = inv(Sigma);
   arma::mat xtxi = inv(trans(X) * Sigi * X);
   arma::mat B = xtxi * trans(X) * Sigi * y;
   return Rcpp::List::create( Rcpp::Named("Estimates") = B);
}')
  
  
  myx <- model.matrix(myfit)
  Sigma <- diag(dim(myx)[1])
  Sigma1 <- (cor(myfit$residuals[-1],myfit$residuals[-length(myfit$residuals)]))^abs(row(Sigma)-col(Sigma))
  #Sigma1<-cor(myfit$residuals)
  
  beta<-fastgls(X=myx,y=sim_data$y,Sigma=Sigma1)
  
  mod1[i,] <- beta$Estimates
  
  #=========================================================================================================#
  
  mymod2<- nlme::gls(y ~ stim1 + stim2 + t + I(t^2) + I(t^3) + signal + stim1_signal,data=sim_data, correlation=corAR1(form=~t|signal))
  mod2[i,] <- mymod2$coefficients
  se2[i,] <- sqrt(diag(vcov(mymod2)))
  
  #=========================================================================================================#
  mymod3 <- lm(y ~ stim1 + stim2 + t + I(t^2) + I(t^3) + signal + stim1_signal,data=sim_data)
  mod3[i,] <- mymod3$coefficients
  se3[i,] <- sqrt(diag(vcov(mymod3)))
  
  #=========================================================================================================#
  
}


names(mod1)<-c("(Intercept)", "stim1", "stim2", "t", "I(t^2)", "I(t^3)", "signal1", "stim1_signal")

write.csv(mod1,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/FASTgls_gamma_results1_inc_corErr0.7_newmethod.csv")

#=========================================================================================================#
#=========================================================================================================#

long_sim_data_mod1<-gather(mod1,key='variable','estimate')
CI_mod1 <- long_sim_data_mod1 %>% dplyr::group_by(variable) %>% dplyr::summarise(CI_low=quantile(estimate,probs = c(0.025)),mean=mean(estimate),CI_upp=quantile(estimate,probs = c(0.975)))

write.csv(CI_mod1,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/FASTgls_gamma_results1_CIs_inc_corErr0.7_newMethod.csv")

vline.dat <- data.frame(variable=unique(long_sim_data_mod1$variable), vl=c(b0,b1,b2,b3,b4,b5,b6,b7))
ggplot(long_sim_data_mod1,aes(x=estimate))+geom_histogram()+geom_vline(aes(xintercept=vl), data=vline.dat,colour='red') +facet_wrap(.~variable,scales = 'free')+theme_bw()
ggsave('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/FASTgls_gamma_results1_estimate_dist_inc_corErr0.7_newMethod.png')
#=========================================================================================================#
#=========================================================================================================#

long_sim_data_mod2<-gather(mod2,key='variable','estimate')
CI_mod2 <- long_sim_data_mod2 %>% dplyr::group_by(variable) %>% dplyr::summarise(CI_low=quantile(estimate,probs = c(0.025)),mean=mean(estimate),CI_upp=quantile(estimate,probs = c(0.975)))

#write.csv(CI_mod2,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/gls_gamma_results1_CIs_inc_corErr0.7_newMethod.csv")

vline.dat2 <- data.frame(variable=unique(long_sim_data_mod2$variable), vl=c(b0,b1,b2,b3,b4,b5,b6,b7))
ggplot(long_sim_data_mod2,aes(x=estimate))+geom_histogram()+geom_vline(aes(xintercept=vl), data=vline.dat2,colour='red') +facet_wrap(.~variable,scales = 'free')+theme_bw()
#ggsave('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/gls_gamma_results1_estimate_dist_inc_corErr0.7_newMethod.png')



write.csv(mod3,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/lm_gamma_results1_inc_corErr0.7_newmethod.csv")
write.csv(se3,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/lm_gamma_results1_inc_corErr0.7_SE_newMethod.csv")

long_sim_data_mod3<-gather(mod3,key='variable','estimate')
CI_mod3 <- long_sim_data_mod3 %>% dplyr::group_by(variable) %>% dplyr::summarise(CI_low=quantile(estimate,probs = c(0.025)),mean=mean(estimate),CI_upp=quantile(estimate,probs = c(0.975)))

write.csv(CI_mod3,"/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/lm_gamma_results1_CIs_inc_corErr0.7_newMethod.csv")


vline.dat <- data.frame(variable=unique(long_sim_data_mod3$variable), vl=c(b0,b1,b2,b3,b4,b5,b6,b7))
ggplot(long_sim_data_mod3,aes(x=estimate))+geom_histogram()+geom_vline(aes(xintercept=vl), data=vline.dat,colour='red') +facet_wrap(.~variable,scales = 'free')+theme_bw()
ggsave('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/lm_gamma_results1_estimate_dist_inc_corErr0.7_newMethod.png')

#=========================================================================================================#
#=========================================================================================================#

res_all<-do.call("rbind", list(CI_mod1,CI_mod2,CI_mod3))

res_all$type<-rep(c('FASTgls','gls','ols'),each=8)


vline.dat2 <- data.frame(variable=unique(long_sim_data_mod2$variable), vl=c(b0,b1,b2,b3,b4,b5,b6,b7))


ggplot(data=res_all,
       aes(x = type,y = mean, ymin = CI_low, ymax = CI_upp))+
  geom_pointrange(aes(col=type))+
  geom_hline(aes(yintercept=vl), data=vline.dat2,linetype=2)+
  xlab('Type (FASTgls, GLS or OLS)')+ ylab("mean estimate")+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_upp,col=type),width=0.5,cex=1)+ 
  facet_wrap(~variable,strip.position="left",nrow=9,scales = "free") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()

ggsave('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/simulations_newMethod/FASTgls_gamma_results1_estimate_boot CI_inc_corErr0.7_newMethod.png')

#ggplot(sim_data,aes(x=t,y=y))+geom_point(colour='grey')+theme_bw()+theme(legend.position = 'top') + geom_line(aes(y=predict(mod1),x=t,colour=signal))+ geom_line(aes(y=predict(mod2),x=t,colour=signal),linetype='dashed')

#=========================================================================================================#
#=========================================================================================================#

