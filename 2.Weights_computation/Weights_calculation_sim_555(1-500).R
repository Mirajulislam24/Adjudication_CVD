library(SoftBart);library(dplyr)
load("data_train_feature_event.RData")
df_train_event=data_train_feature_event
load("data_test_feature_event.RData")
df_test_event=data_test_feature_event
load("data_train_feature_censored.RData")
df_train_censored=data_train_feature_censored
load("data_test_feature_censored.RData")
df_test_censored=data_test_feature_censored
load("Simdata_test_withTrue_EstCVD_3R.RData")
data_sim_test=data_sim_test

data_sim_test_withADJ<-list()
for (n_sim in 1:500){
  s=5
#AMONG UNADJ==EVENT
df_A1<-df_train_event[[n_sim]];df_A2=df_test_event[[n_sim]]
data_A1<-df_A1[,c("id","CVD_DTH_ADJ_sim","Age_death","base_age","BMI",
                  "SEX","RACE","V1","V2","V3",
                  "LH","AH")]#,"RXHYP","SMOKER1")]#only necessary variables
data_A1$H<-ifelse(data_A1$LH == 0& data_A1$AH == 0, 1, 0) 

library(dplyr)
data_A2<-list()
for (i in 1:500){
  data_A2[[i]]=data_A1[seq(i, nrow(data_A1),100), ] 
}

data_A3<-df_A2[,c("id","CVD_DTH_ADJ_sim","Age_death","base_age","BMI",
                  "SEX","RACE","V1","V2","V3",
                  "LH","AH")]#,"RXHYP","SMOKER1")]#only necessary variables
data_A3$H<-ifelse(data_A3$LH == 0& data_A3$AH == 0, 1, 0) 

data_A4<-list()
for (i in 1:500){
  data_A4[[i]]=data_A3[seq(i, nrow(data_A3),100), ] 
}

test_data<-list()
P1<-matrix(0,ncol=s,nrow=sum(data_A4[[1]]$CVD_DTH_ADJ_sim==1))
P0<-matrix(0,ncol=s,nrow=sum(data_A4[[1]]$CVD_DTH_ADJ_sim==0))
varimp<-matrix(0,nrow=ncol(data_A4[[1]][,-c(1,2)]),ncol=s)
P<-list()
for (i in 1:s){
  data1=data_A2[[i]] #train data
  data2=data_A4[[i]]#test data
  trn_data <- data.frame(X = data1[,c(3:13)], Y = factor(data1[,2], levels = c(0,1))) 
  tst_data <- data.frame(X = data2[,c(3:13)],Y=factor(data2[,2], levels = c(0,1)))
  fitted_probit <- softbart_probit(Y ~ ., data = trn_data, test_data = tst_data, verbose = TRUE,num_tree=30,
                                   opts = Opts(num_burn =500, num_save =s))
  
  test_data[[i]]=tst_data
  test_data[[i]]$event=tst_data$Y
  P[[i]]=fitted_probit$p_test
  test_data[[i]]$P=fitted_probit$p_test_mean
  P1[,i]<-test_data[[i]]$P[test_data[[i]]$event==1]
  P0[,i]<-test_data[[i]]$P[test_data[[i]]$event==0]
  varimp[,i]<-posterior_probs(fitted_probit)$varimp
}
P1=cbind(t(P[[1]]),t(P[[2]]),t(P[[3]]),t(P[[4]]),t(P[[5]]))#N row and s*s column
N=length(data_A4[[1]]$CVD_DTH_ADJ_sim)#total individual in testdata

V<-list()
for (i in 1:(s*s)){
  V[[i]]=matrix(0,nrow=N,ncol=s)
  for(k in 1:s){
    U<-runif(N)
    V[[i]][,k]<-ifelse(U<P1[,i],1,0)
  } }

ADJ_event<-matrix(unlist(V),ncol=s*s*s,nrow=N)




#AMONG UNADJ==CENSORED

df_B1<-df_train_censored[[n_sim]];df_B2=df_test_censored[[n_sim]]
data_B1<-df_B1[,c("id","CVD_DTH_ADJ_sim","Age_death","base_age","BMI",
                  "SEX","RACE","V1","V2","V3",
                  "LH","AH")]#,"RXHYP","SMOKER1")]#only necessary variables
data_B1$H<-ifelse(data_B1$LH == 0& data_B1$AH == 0, 1, 0) 

library(dplyr)
data_B2<-list()
for (i in 1:500){
  data_B2[[i]]=data_B1[seq(i, nrow(data_B1),100), ] 
}

data_B3<-df_B2[,c("id","CVD_DTH_ADJ_sim","Age_death","base_age","BMI",
                  "SEX","RACE","V1","V2","V3",
                  "LH","AH")]#,"RXHYP","SMOKER1")]#only necessary variables
data_B3$H<-ifelse(data_B3$LH == 0& data_B3$AH == 0, 1, 0) 

library(dplyr)
data_B4<-list()
for (i in 1:500){
  data_B4[[i]]=data_B3[seq(i, nrow(data_B3),100), ] 
}


s=5
test_data<-list()
P<-list()
for (i in 1:s){
  data1=data_B2[[i]] #train data
  data2=data_B4[[i]]#test data
  trn_data <- data.frame(X = data1[,c(3:13)], Y = factor(data1[,2], levels = c(0,1))) 
  tst_data <- data.frame(X = data2[,c(3:13)],Y=factor(data2[,2], levels = c(0,1)))
  fitted_probit <- softbart_probit(Y ~ ., data = trn_data, test_data = tst_data, verbose = TRUE,num_tree=30,
                                   opts = Opts(num_burn =500, num_save =s))
  
  test_data[[i]]=tst_data
  test_data[[i]]$event=tst_data$Y
  P[[i]]=fitted_probit$p_test
}

P1=cbind(t(P[[1]]),t(P[[2]]),t(P[[3]]),t(P[[4]]),t(P[[5]]))#N row and s*s column
N=length(data_B4[[1]]$CVD_DTH_ADJ_sim)#total individual in testdata

V<-list()
for (i in 1:(s*s)){
  V[[i]]=matrix(0,nrow=N,ncol=s)
  for(k in 1:s){
    U<-runif(N)
    V[[i]][,k]<-ifelse(U<P1[,i],1,0)
  } }

ADJ_censored<-matrix(unlist(V),ncol=s*s*s,nrow=N)

ID=c(data_A4[[1]]$id,data_B4[[1]]$id)
Ind_Adj<-rbind(ADJ_event,ADJ_censored)
data_adj_id<-data.frame(id=ID,CVD_ADJ_sim_est1=Ind_Adj)
data_sim_test_withADJ[[n_sim]]=merge(data_sim_test[[n_sim]],data_adj_id,by="id",all.x=TRUE)
}
save(data_sim_test_withADJ,file="data_sim_test_withADJ_555.RData")
