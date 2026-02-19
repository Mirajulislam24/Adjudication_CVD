
idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
s=500#Number of posterior sample data
parameters <-1:s
myparam <- parameters[idx]

Sample_our<-list()

library(dplyr);library(JM);library(splines)
load("/home/mirajulislam/ADJUDICATION_ARIC/Simulation/Final_new/1.Sim_500/Sample_1000/Our_approach_sim/2.Weights_computation//data_sim_test_withADJ_555.RData")

data<-data_sim_test_withADJ[[myparam]]
data$new_id <- match(data$id, unique(data$id))
data.id<-data[!duplicated(data$new_id), ]# to get baseline info
table(data.id$CVD_DTH_ADJ_sim)
T=15 #maximum time
legP2<-function(t){
  #T<-max(t)
  cbind((2*t/T-1),(-0.5+1.5*(2*t/T-1)^2))
}

ctrl <- lmeControl(opt = "optim", msMaxIter = 200, msVerbose = TRUE)
lmeObject1<- lme(y1 ~legP2(time), data = data,
                 random=list(new_id=pdDiag(form=~legP2(time))),control=ctrl)

lmeObject2<- lme(y2 ~legP2(time), data = data,
                 random=list(new_id=pdDiag(form=~legP2(time))),control=ctrl)

lmeObject3<- lme(y3 ~legP2(time), data = data,
                 random=list(new_id=pdDiag(form=~legP2(time))),control=ctrl)

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


#################################
# survival submodel
# design matrices for the survival submodel
#data.id$RACE<-data.id$group
WD <- model.matrix(~-1+base_age+BMI+LH+AH+SEX+RACE, data=data.id)
nT <- length(Time)
zeros <- numeric(nT)

x <- list(X1 = X1, Z1 = Z1, WD =WD)

###################################
# for the longitudinal outcomes - design matrices for the 15-point Gauss-Kronrod quadrature rule approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)
P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))


# 1st longitudinal outcome
mfX1 <- model.frame(TermsX1, data = data.id2)  
mfZ1 <- model.frame(TermsZ1, data = data.id2)    
Xs1 <- model.matrix(formYx1, mfX1)
Zs1 <- model.matrix(formYz1, mfZ1)


#################################
# set MCMC details
con <- list( K = 100,C = 5000, knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4)

# design matrices for the baseline hazard (a B-splines baseline hazard function is asssumed)
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- if (con$ObsTimes.knots) {
    Time
  } else {Time[event == 1]    }
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr

W2D <- splineDesign(rr, Time, ord = con$ordSpline)
# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation

W2sD <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(W2D = W2D, W2sD = W2sD))

#################################
ncX <- ncol(X1)
ncZ <- ncol(Z1)
ncWD <- ncol(x$WD)
ncW2D <- ncol(x$W2D)
C <- con$C
nb <- 3*ncZ
b <- cbind(data.matrix(ranef(lmeObject1)),data.matrix(ranef(lmeObject2)),
           data.matrix(ranef(lmeObject3)))

nY <- nrow(b)

#################################
# priors/hyperpriors
mu01=mu02=mu03<- rep(0, 3)




###All design matrices are same
X=X1;Z=Z1;Xtime=Xtime1;Ztime=Ztime1;Xs=Xs1;Zs=Zs1
#################################
#################################
# 2 stage approach - calculation of the mean and sd of the longitudinal outcomes
# 1st longitudinal outcome
muV1 <- mean(Xtime%*%fixef(lmeObject1) + rowSums(Ztime*ranef(lmeObject1)))
stdV1 <- sd(Xtime%*%fixef(lmeObject1) + rowSums(Ztime*ranef(lmeObject1)))

# 2nd longitudinal outcome
muV2 <- mean(Xtime%*%fixef(lmeObject2) + rowSums(Ztime*ranef(lmeObject2)))
stdV2 <- sd(Xtime%*%fixef(lmeObject2) + rowSums(Ztime*ranef(lmeObject2)))

# 3rd longitudinal outcome
muV3 <- mean(Xtime%*%fixef(lmeObject3) + rowSums(Ztime*ranef(lmeObject3)))
stdV3 <- sd(Xtime%*%fixef(lmeObject3) + rowSums(Ztime*ranef(lmeObject3)))

#################################
#################################
##Stan Code
library("rstan")
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
vector evaluate_eta1(matrix X_f, matrix Z_f, vector beta_f,
matrix bMat_f) {
int N_f = rows(X_f); // num rows in design matrix
vector[N_f] eta;
for (i in 1:N_f)
eta[i] = dot_product(X_f[i,],beta_f)+ dot_product(bMat_f[i,], Z_f[i,]);
return eta;
}
}
data{
  int<lower=1> N; // number of individuals
  int<lower=1> N1; // number of total observations
  int<lower=1> N_s; // number of total observations of Gauss kronrod
  int<lower=1> P_x;// ncols of original fixed design matrix
  int<lower=1> P_z;// ncols of original random design matrix
  int<lower=1> K;// 15point Gauss-quardrature
  int offset[N+1];
  int<lower=0> P_w;// ncols of baseline covariates CoxModel
  int<lower=0> P_w2;// ncols of baseline hazard design matrix
  int<lower=0> P_w2s;// ncols of baseline hazard design matrix Gauss Quadrature
  int<lower=1> n_long;
  int<lower=1>  nb;// ncols of random effects
  real<lower=0> P[N];
  vector[K] wk;
  matrix[N1,n_long] y; // long outcomes
  matrix[N1, P_x] X;
  matrix[N, P_x] Xtime;
  matrix[N_s, P_x] Xs;
  matrix[N1, P_z] Z;
  matrix[N, P_z] Ztime;
  matrix[N_s, P_z] Zs;
  vector[N] eventD;
  matrix[N, P_w] WD;
  matrix[N, P_w2] W2D;
  matrix[N_s, P_w2s] W2sD;
  vector[n_long] Mu;
  vector<lower=0>[n_long] SD;
  vector[n_long*P_x] Means;
  real <lower=0> prior_scale_long;
}
parameters {
  matrix[P_x,n_long] z_betas;// primitive coefs in long. submodels
  vector <lower=0> [n_long] z_sigma;
  // group level params (random effects params)   
  vector<lower=0>[nb] b_sd1;
  matrix[nb,N] z_b_mat1;   // unscaled group level params 
  cholesky_factor_corr[nb] b_cholesky1;// cholesky factor of corr matrix
  vector[P_w] gammasD;// primitive coefs in surv. submodel
  vector[P_w2s] Bs_gammasD;
  //VSpriors parameters
 vector[n_long] alpha;
}

transformed parameters{
matrix[P_x,n_long] betas;
  vector <lower=0> [n_long] sigma;
  matrix[N, nb] b;
  vector[N] etaBaselineD;
  vector[N] log_h0_TD;
  matrix[N,n_long] f_T;
  matrix[N,n_long]  f_T_derivY;
  matrix[N,n_long]  fA_T_derivY;
  vector[N] log_hazardD;
  matrix[K,N] log_h0_sD;
  matrix[K,N] f_s1;
  matrix[K,N] f_s2;
  matrix[K,N] f_s3;
  matrix[K,N] SurvLongD;
  vector[N] log_survivalD;
  b[,1:nb] = (diag_pre_multiply(b_sd1, b_cholesky1) * z_b_mat1)';
  betas[,1]=prior_scale_long*z_betas[,1]+Means[1:3];
  betas[,2]=prior_scale_long*z_betas[,2]+Means[4:6];
  betas[,3]=prior_scale_long*z_betas[,3]+Means[7:9];
  sigma=prior_scale_long*z_sigma;
//baseline linear predictor and features calculation for event model 
    etaBaselineD = WD*gammasD;
    log_h0_TD = W2D*Bs_gammasD;
    f_T[,1] = evaluate_eta1(Xtime,Ztime,betas[,1],b[,1:P_z]);
    f_T[,2] = evaluate_eta1(Xtime,Ztime,betas[,2],b[,(P_z+1):(2*P_z)]);
    f_T[,3] = evaluate_eta1(Xtime,Ztime,betas[,3],b[,(2*P_z+1):(3*P_z)]);
                      
    
     log_hazardD= log_h0_TD + etaBaselineD+ 
     alpha[1] * ((f_T[,1] - Mu[1])/SD[1]) +alpha[2]* ((f_T[,2] - Mu[2])/SD[2])+ alpha[3]* ((f_T[,3] - Mu[3])/SD[3]);
    
    for (i in 1:N) {
      log_h0_sD[,i] = W2sD[K * (i - 1) + 1:K * (i - 1) + 15,]*Bs_gammasD;
      f_s1[,i] = Xs[K * (i - 1) + 1:K * (i - 1) + 15,]*betas[,1]+
                 Zs[K * (i - 1) + 1:K * (i - 1) + 15,]*b[i, 1:3]';
      f_s2[,i] = Xs[K * (i - 1) + 1:K * (i - 1) + 15,]*betas[,2]+
                 Zs[K * (i - 1) +1:K * (i - 1) + 15,]*b[i, 4:6]';
      f_s3[,i] = Xs[K * (i - 1) + 1:K * (i - 1) + 15,]*betas[,3]+
                 Zs[K * (i - 1) +1:K * (i - 1) + 15,]*b[i, 7:9]';
     
      SurvLongD[,i] =  wk.* exp(log_h0_sD[, i] + 
                      alpha[1] * ((f_s1[, i] - Mu[1])/SD[1])+alpha[2]* ((f_s2[, i] - Mu[2])/SD[2])+
                      alpha[3]* ((f_s3[, i] - Mu[3])/SD[3])); 
    log_survivalD[i] = -(exp(etaBaselineD[i]) * P[i] *sum(SurvLongD[,i]));
    }
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
//----- Log-priors in long_submodel
z_betas[,1] ~ normal(0,1);
z_betas[,2] ~ normal(0,1);
z_betas[,3] ~ normal(0,1);
z_sigma~normal(0,1);
 b_sd1~normal(0, 1);
 to_vector(z_b_mat1) ~ normal(0,1);
 target += lkj_corr_cholesky_lpdf(b_cholesky1 | 1);
//Event submodel
target += dot_product(eventD,log_hazardD)+sum(log_survivalD);

//----- Log-priors in surv_submodel
// long. submodels
// baseline survival parameters
gammasD~normal(0,10);
Bs_gammasD~normal(0,10);
//VS priors
alpha~normal(0,10);
}

generated quantities {
vector[n_long] Talpha;
  for (i in 1:n_long){
   Talpha[i]=alpha[i]/(SD[i]);
   }
   
}"
nec_par=c("Talpha","gammasD")


Data <- list(N = nY,N1=nrow(X),N_s=nrow(Xs),P_x=ncol(X),
             P_z=ncol(Z), K = K, offset = offset,
             y= cbind(data$y1,data$y2,data$y3), 
             X = X, Xtime = Xtime,
             Xs = Xs,
             Z = Z, Ztime = Ztime,
             Zs = Zs,
             WD = WD,
             P_w = ncol(WD),
             
             W2D = W2D, P_w2 = ncol(W2D),
             W2sD = W2sD, P_w2s = ncol(W2sD), 
             
             P = P,
             wk = wk,
             Mu =c(muV1,muV2,muV3),
             SD=c(stdV1,stdV2,stdV3),
             nb=nb,
             n_long=3, #number of risk factors
             Means=c(mean(data$y1),rep(0,2),mean(data$y2),rep(0,2),mean(data$y3),rep(0,2)),
             prior_scale_long=10)

npost=300;nwarm=700;niter=nwarm+npost;nrefresh=500
mat_samp=matrix(0,nrow=125*npost,ncol=10)

for (i in 1:125){
  #i=73
  Data$eventD<-data.id[,(15+i)]
  
  #options(mc.cores=parallel::detectCores())
  system.time(
    fit <- stan(model_code = scode, 
                data = Data,
                pars= nec_par,
                warmup=nwarm,
                iter =niter, 
                thin=1,
                chains = 1, 
                verbose = FALSE,
                cores = 2,
                init_r=1,
                seed=786,
                refresh=nrefresh
    )
  )
  
  fit_summary <- summary(fit)
  (S=round(fit_summary$summary,3))
  if(abs(S[2,1])>5){
    fit <- stan(model_code = scode, 
                data = Data,
                pars= nec_par,
                warmup=nwarm,
                iter =niter, 
                thin=1,
                chains = 1, 
                verbose = FALSE,
                cores = 2,
                init_r=1,
                seed=786,
                refresh=nrefresh,
                control = list(max_treedepth=12)#adapt_delta = 0.99)
    )
  }
  #fit_summary1 <- summary(fit)
  #(round(fit_summary1$summary,3))
  samples<-rstan:::extract(fit)
  mat_samp[((i-1)*npost+1):(npost*i),]<-round(cbind(myparam,samples$Talpha,samples$gammasD),4)
}
Sample_our[[1]]<-mat_samp
filename1 = 'Sample_ADJ_EST.txt'
write.table(Sample_our,file=filename1, append=TRUE)