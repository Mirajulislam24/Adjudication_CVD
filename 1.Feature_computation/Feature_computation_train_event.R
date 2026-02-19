
library(dplyr);library(JM);library(splines);library(rstan)

load("Simdata_train_withTrue_EstCVD_3R.RData")
data_train_feature_event=list()
for (i in 1:500){
  data_train<-data_sim_train[[i]]
  data_train<-subset(data_train,CVD_DTH_UNADJ_sim==1)
  data_train$new_id <- match(data_train$id, unique(data_train$id))
  data_train.id<-data_train[!duplicated(data_train$new_id), ]# to get baseline info
  table(data_train.id$CVD_DTH_ADJ_sim)
  
  data=data_train
  data.id=data_train.id
  T=15 #maximum time
  legP2<-function(t){
    #T<-max(t)
    cbind((2*t/T-1),(-0.5+1.5*(2*t/T-1)^2))
  }
  
  lmeObject1<- lme(y1 ~legP2(time), data = data,
                   random=list(new_id=pdDiag(form=~legP2(time))), control = lmeControl(opt = "optim"))
  
  lmeObject2<- lme(y2 ~legP2(time), data = data,
                   random=list(new_id=pdDiag(form=~legP2(time))), control = lmeControl(opt = "optim"))
  
  lmeObject3<- lme(y3 ~legP2(time), data = data,
                   random=list(new_id=pdDiag(form=~legP2(time))), control = lmeControl(opt = "optim"))
  
  ### Set
  timeVar <- "time"
  lag <- 0
  survMod <- "spline-PH"
  Time <- data.id$Time###the survival time
  
  # for the continuous longitudinal outcome create the design matrices
  id <- data$new_id
  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  
  # 1st longitudinal outcome
  formYx1 <- formula(lmeObject1)
  TermsX1 <- lmeObject1$terms
  mfX1 <- model.frame(TermsX1, data = data)
  X1 <- model.matrix(formYx1, mfX1)
  
  formYz1 <- formula(lmeObject1$modelStruct$reStruct[[1]])
  mfZ1 <- model.frame(terms(formYz1), data = data)
  TermsZ1 <- attr(mfZ1, "terms")
  Z1 <- model.matrix(formYz1, mfZ1)
  
  
  data.id[[timeVar]] <- pmax(Time - 0, 0)
  
  mfX.id1 <- model.frame(TermsX1, data = data.id) 
  mfZ.id1 <- model.frame(TermsZ1, data = data.id)  
  Xtime1 <- model.matrix(formYx1, mfX.id1)
  Ztime1 <- model.matrix(formYz1, mfZ.id1)
  
  ###All design matrices are same
  X=X1;Z=Z1;Xtime=Xtime1;Ztime=Ztime1
  ncZ=ncol(Z)
  nb <- 3*ncZ
  b <- cbind(data.matrix(ranef(lmeObject1)),data.matrix(ranef(lmeObject2)),
             data.matrix(ranef(lmeObject3)))
  
  nY <- nrow(b)

  
  #################################
  Data <- list(N = nY,N1=nrow(X),P_x=ncol(X),
               P_z=ncol(Z),  offset = offset,
               y= cbind(data$y1,data$y2,data$y3), 
               X = X, Z = Z, Xtime = Xtime,
               nb=nb,
               n_long=3) #number of risk factors
  scode <- "functions{
vector evaluate_eta(matrix X_f, matrix Z_f, vector beta_f,
matrix bMat_f, int N, int[] offset_f) {
int N_f = rows(X_f); // num rows in design matrix
vector[N_f] eta;
for (i in 1:N)
  for (j in offset_f[i]:(offset_f[i + 1] - 1))
eta[j] = X_f[j,] * beta_f+ bMat_f[i,]* Z_f[j,]';

return eta;
}
}
data{
  int<lower=1> N; // number of individuals
  int<lower=1> N1; // number of total observations
  int<lower=1> P_x;// ncols of original fixed design matrix
  int<lower=1> P_z;// ncols of original random design matrix
  int offset[N+1];
  int<lower=1> n_long;
  matrix[N1,n_long] y; // long outcomes
  matrix[N1, P_x] X;
  matrix[N, P_x] Xtime;
  matrix[N1, P_z] Z;
  int<lower=1>  nb;// ncols of random effects
}
parameters {
  matrix[P_x,n_long] betas;
  vector <lower=0> [n_long] sigma;

  // group level params (random effects params)   
  vector<lower=0>[nb] sigma_b;
  matrix[nb,N] z_b;   // unscaled group level params 
  cholesky_factor_corr[nb] L_corr;// cholesky factor of corr matrix
}

transformed parameters {
  matrix[N, nb] b;                                // actual random effects
  matrix[nb, nb] L_b;
  L_b = diag_pre_multiply(sigma_b, L_corr);
  b[,1:nb] = (L_b * z_b)';                    // note: transpose for row vector
}
model{
// declare linear predictors
vector[N1] y1_eta;
vector[N1] y2_eta;
vector[N1] y3_eta;
// evaluate linear predictor for each long. submodel
y1_eta = evaluate_eta(X, Z, betas[,1], b[,1:3],N,offset);
y2_eta = evaluate_eta(X, Z, betas[,2], b[,4:6],N,offset);
y3_eta = evaluate_eta(X, Z, betas[,3], b[,7:9],N,offset);
// long models
y[,1]~normal(y1_eta, sqrt(sigma[1]));
y[,2]~normal(y2_eta, sqrt(sigma[2]));
y[,3]~normal(y3_eta, sqrt(sigma[3]));

to_vector(z_b) ~ normal(0, 1);                 // prior for latent variables
  sigma_b ~ normal(0, 1);                      // weakly informative prior
  L_corr ~ lkj_corr_cholesky(1);                 // prior on correlations
}
generated quantities {
    vector[N] V1;
    vector[N] V2;
    vector[N] V3;
  V1 = rows_dot_product(b[,1:3],Xtime[,1:3])+Xtime[,1:3]*betas[,1];
  V2 = rows_dot_product(b[,4:6],Xtime[,1:3])+Xtime[,1:3]*betas[,2];
  V3 = rows_dot_product(b[,7:9],Xtime[,1:3])+Xtime[,1:3]*betas[,3];
  
}
"
nec_par=c("V1","V2","V3")

#setwd("D:/A_UF STUDY_PHD/1.DISSERTATION/Adjudication/Bi-level using normal mix")


library(rstan)
options(mc.cores=parallel::detectCores())
system.time(
  fit2 <- stan(model_code = scode, 
               data = Data,
               pars= nec_par,
               warmup=400,
               iter =500, 
               chains = 1, 
               verbose = FALSE,
               thin=1,
               cores = 2,
               init_r=1,
               seed=7865
  )
)

fit_summary <- summary(fit2)
#traceplot(fit2,pars="betas2")
round((fit_summary$summary[1:50,c(1:4,8:10)]),2)

s=100
n=length(Data$Xtime[,1])
VS1<-rstan::extract(fit2,  pars=c("V1"),permuted = FALSE, inc_warmup = FALSE,include=TRUE)
V1<-as.vector(apply(VS1,2,c))
VS2<-rstan::extract(fit2,  pars=c("V2"),permuted = FALSE, inc_warmup = FALSE,include=TRUE)
V2<-as.vector(apply(VS2,2,c))
VS3<-rstan::extract(fit2,  pars=c("V3"),permuted = FALSE, inc_warmup = FALSE,include=TRUE)
V3<-as.vector(apply(VS3,2,c))

data_pst=data.id[rep(seq_len(nrow(data.id)), each = s), ]

data_train_feature_event[[i]]<-cbind(data_pst,V1=V1,V2=V2,V3=V3)
}
save(data_train_feature_event,file="data_train_feature_event.RData")
