#--------------------------------------------------------------------#
#Processing script for fMRI and fTCD glm data.
#--------------------------------------------------------------------#

#15-08-2019

library(tidyverse)

get_data=0

if(get_data==1)
{
user_renviron = path.expand(file.path("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project", ".Renviron"))
if(!file.exists(user_renviron)) # check to see if the file already exists
  file.create(user_renviron)
file.edit(user_renviron) #paste in the PAT token created earlier into file like this: OSF_PAT='insert token token'

osf_auth(token=Sys.getenv("OSF_PAT"))

lb_project <- osf_retrieve_node("gw4en")

osf_retrieve_file("https://osf.io/gw4en/") %>%
  osf_download('Raw Data/fMRI_leftdata.csv')  

osf_retrieve_file("https://osf.io/gw4en/") %>%
  osf_download('Raw Data/fMRI_rightdata.csv')  
}

#--------------------------------------------------------------------#

key_data<-read.csv('Key_fMRI_fTCD.csv')
key_data2<-key_data[is.na(key_data$fMRI_row)==FALSE,]

which(is.na(key_data2$fTCD_row)==TRUE)

zoe_testL<-read.csv('fMRI_leftdata.csv')
zoe_testR<-read.csv('fMRI_rightdata.csv')


ts_summary_data<-read.csv("ts_summary_data.csv")


for(i in 1:dim(my_results)[1])
{
my_results$Lparam1_perc <- my_results$Lparam1*(1/ts_summary_data$bmeanL)*100
my_results$Rparam1_perc <- my_results$Rparam1*(1/ts_summary_data$bmeanR)*100
}


#--------------------------------------------------------------------#

zoe_testL_gg<-as.data.frame(t(zoe_testL))
colnames(zoe_testL_gg)<-key_data2$fTCD_ID
zoe_testL_gg<-zoe_testL_gg[,-which(is.na(key_data2$fTCD_row)==TRUE)]

zoe_testR_gg<-as.data.frame(t(zoe_testR))
colnames(zoe_testR_gg)<-key_data2$fTCD_ID
zoe_testR_gg<-zoe_testR_gg[,-which(is.na(key_data2$fTCD_row)==TRUE)]

IDs<-colnames(zoe_testL_gg)

medianL<-apply(zoe_testL_gg,2,median)
medianR<-apply(zoe_testR_gg,2,median)

dop_est<-as.data.frame(matrix(NA,nrow =length(medianL),ncol=2))
names(dop_est)<-c("Left_dop_param","right_dop_param")

#--------------------------------------------------------------------#

library(ggpubr)

pdf("dist_plot_params.pdf")
for(k in seq_along(IDs))
{
p1<-ggplot(data=zoe_testL_gg)+geom_histogram(aes_string(x=as.name(IDs[k])),colour="grey",fill="grey",alpha=0.5) +geom_vline(aes_string(xintercept = my_results$Lparam1_perc[my_results$ID==IDs[k]]),col='red',size=1)+theme_bw()+geom_vline(aes_string(xintercept = medianL[which(IDs==IDs[k])]),col='blue',size=1)
 
p2<-ggplot(data=zoe_testR_gg)+geom_histogram(aes_string(x=as.name(IDs[k])),colour="grey",fill="grey",alpha=0.5) +geom_vline(aes_string(xintercept = my_results$Rparam1_perc[my_results$ID==IDs[k]]),col='red',size=1)+theme_bw()+geom_vline(aes_string(xintercept = medianR[which(IDs==IDs[k])]),col='blue',size=1)
  
p3<-ggarrange(p1,p1,nrow = 2,labels=c('left','right'))
print(p3)

dop_est[k,1]<-my_results$Lparam1_perc[my_results$ID==IDs[k]]
dop_est[k,2]<-my_results$Rparam1_perc[my_results$ID==IDs[k]]
}
dev.off()

#--------------------------------------------------------------------#

cor(dop_est[,1],medianL)

cor(dop_est[,2],medianR)

plot(dop_est[,1],medianL)
plot(dop_est[,2],medianR)

cor(medianL-medianR,dop_est[,1]-dop_est[,2])

