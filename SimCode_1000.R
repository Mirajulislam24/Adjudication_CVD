x=c("JMbayes","splines","parallel","runjags","coda","dplyr")
#install.packages(x)
lapply(x, library, character.only = TRUE)
load("ARICdata_1000.RData")
data.id <- data[!duplicated(data$ID),]
SEX=data.id$SEX1;RACE=data.id$RACE1;LH=data.id$LH
AH=data.id$AH;BMI=data.id$BMI
CVD_DTH_UNADJ=data.id$CVD_DTH_ICD
CVD_DTH_ADJ=data.id$ADJUDICATED_CVDDTH
STRK_UNADJ_num=data.id$STRK_ICD_num
MI_UNADJ_num=data.id$MI_ICD_num
HF_UNADJ_num=data.id$HF_ICD_num
##Mixed model params
load("LM_Est_ARIC_1200_3R.RData")
betas1=data_long$betas1;betas2=data_long$betas2
betas3=data_long$betas3
sigma1=data_long$sigma1;sigma2=data_long$sigma2
sigma3=data_long$sigma3-300
D=data_long$D
alpha1=0.003;alpha2=-0.05;alpha3=0
#alpha1=0;alpha2=-.14;alpha3=0.03
gamma0=0.002#baseline_age
gamma1=0.020;gamma2=0.035;gamma3=-0.489;gamma4=0.474;gamma5=-0.470

param=c(1001,1004,1010,1011,1012,1013,1015,1019,1022,1025,1030,1031,1034,1036,1039,1040,1041,1043,1049,1051,
       1053,1054,1055,1057,1058,1059,1060,1061,1063,1066,1067,1070,1071,1074,1075,1076,1077,1079,1083,1084,
       1085,1086,1087,1088,1095,1097,1100,1102,1104,1105,1107,1116,1118,1120,1123,1126,1127,1128,1131,1132, 
       1135,1136,1138,1141,1142,1143,1145,1147,1149,1150,1152,1155,1159,1160,1161,1163,1164,1165,1166,1168,
       1169,1170,1171,1174,1178,1179,1181,1182,1183,1184,1188,1189,1191,1192,1194,1198,1199,1201,1203,1204, 
       1205, 1206, 1207, 1208, 1210, 1211, 1214, 1215, 1216, 1218, 1221, 1222, 1223, 1224, 1227, 1228, 1229, 1231, 1234, 1235,
       1237, 1238, 1241, 1243, 1248, 1249, 1250, 1251, 1252, 1254, 1258, 1260, 1262, 1263, 1265, 1269, 1271, 1272, 1274, 1276, 
       1277, 1279, 1280, 1281, 1284, 1285, 1286, 1291, 1293, 1294, 1296, 1297, 1303, 1305, 1307, 1308, 1309, 1313, 1319, 1321, 
       1326, 1328, 1329, 1331, 1334, 1335, 1338, 1340, 1341, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1352, 1354, 1357, 1358, 
       1361, 1362, 1363, 1365, 1367, 1371, 1372, 1374, 1375, 1377, 1378, 1379, 1380, 1382, 1383, 1389, 1390, 1391, 1395, 1397, 
       1399, 1400, 1403, 1405, 1413, 1417, 1419, 1420, 1421, 1427, 1428, 1429, 1430, 1431, 1432, 1434, 1439, 1440, 1441, 1442, 
       1444, 1449, 1450, 1454, 1456, 1457, 1459, 1460, 1461, 1462, 1463, 1468, 1469, 1471, 1474, 1475, 1476, 1477, 1478, 1479, 
       1480, 1482, 1484, 1485, 1486, 1487, 1488, 1491, 1494, 1495, 1497, 1499, 1501, 1504, 1508, 1512, 1514, 1518, 1522, 1526, 
       1529, 1530, 1531, 1534, 1535, 1538, 1540, 1541, 1542, 1543, 1545, 1547, 1548, 1551, 1553, 1554, 1557, 1558, 1560, 1561, 
       1562, 1563, 1564, 1566, 1569, 1570, 1571, 1572, 1573, 1576, 1577, 1580, 1585, 1587, 1588, 1589, 1591, 1592, 1595, 1597, 
       1598, 1599, 1600, 1601, 1605, 1606, 1607, 1608, 1611, 1612, 1613, 1615, 1616, 1625, 1626, 1627, 1629, 1632, 1634, 1635, 
       1637, 1638, 1642, 1643, 1645, 1646, 1647, 1649, 1654, 1656, 1658, 1659, 1660, 1661, 1665, 1666, 1667, 1668, 1670, 1671, 
       1673, 1675, 1677, 1678, 1680, 1681, 1683, 1684, 1687, 1689, 1690, 1694, 1696, 1697, 1699, 1700, 1701, 1705, 1706, 1707, 
       1708, 1713, 1714, 1715, 1716, 1722, 1725, 1726, 1727, 1728, 1731, 1732, 1736, 1737, 1742, 1743, 1747, 1748, 1752, 1754, 
       1757, 1760, 1761, 1764, 1768, 1770, 1772, 1774, 1775, 1776, 1777, 1778, 1779, 1782, 1784, 1786, 1789, 1790, 1792, 1794, 
       1795, 1798, 1802, 1803, 1804, 1806, 1809, 1810, 1811, 1812, 1814, 1815, 1816, 1823, 1826, 1828, 1831, 1833, 1835, 1837, 
       1838, 1840, 1846, 1850, 1851, 1854, 1856, 1858, 1863, 1864, 1865, 1866, 1867, 1868, 1869, 1870, 1874, 1876, 1877, 1878, 
       1879, 1881, 1885, 1886, 1887, 1889, 1890, 1891, 1895, 1897, 1898, 1900, 1902, 1903, 1905, 1918, 1922, 1924, 1926, 1927, 
       1929, 1930, 1931, 1932, 1933, 1934, 1937, 1938, 1939, 1942, 1944, 1946, 1947, 1951, 1952, 1953, 1955, 1956, 1957, 1958, 
       1959, 1963, 1964, 1966, 1969, 1970, 1973, 1974, 1976, 1977, 1978, 1981, 1984, 1986, 1987, 1988, 1995, 1996, 1997, 1998)

data_sim_test<-list()
data_sim_train<-list()

for (i in 1:length(param)){
  seed=param[i]
set.seed(seed)
#alpha1=0.003;alpha2=-0.05;alpha3=0
n=1000
  K <- 5 #10  # number of planned repeated measurements per subject, per outcome
  t.max <- 15 # maximum follow-up time
  times <- list()#c(replicate(n, c(0,3,6,9,12)))
  time0<-sample(0:19,n,replace=T)#baseline age
  for (j in 1:n){
    times[[j]]<-seq(0,12,by=3)
  }
  times<-unlist(times)
  t_m=t.max
  legP2<-function(t){
    #T<-max(t)
    cbind((2*t/t_m-1),(-0.5+1.5*(2*t/t_m-1)^2))
  }
  
  #Av.BP
  betas1 <- betas1
  sigma.y1 <- sqrt(sigma1) # measurement error standard deviation
  #Glucose
  betas2 <- betas2
  sigma.y2 <- sqrt(sigma2) # measurement error standard deviation
  #TOTCHL
  betas3 <- betas3
  sigma.y3 <- sqrt(sigma3) # measurement error standard deviation
  
  D<-data_long$D
  # design matrices for the longitudinal measurement model
  
  DF <- data.frame(base_age=rep(time0,each=K),time = times,CVD_DTH_UNADJ=rep(CVD_DTH_UNADJ,each=K),CVD_DTH_ADJ=rep(CVD_DTH_ADJ,each=K),
                   BMI=rep(BMI,each=K),LH=rep(LH,each=K),AH=rep(AH,each=K),SEX=rep(SEX,each=K),RACE=rep(RACE,each=K),
                   STRK_UNADJ_num=rep(STRK_UNADJ_num,each=K),MI_UNADJ_num=rep(MI_UNADJ_num,each=K),
                   HF_UNADJ_num=rep(HF_UNADJ_num,each=K))
  
  
  X=Z <- model.matrix(~legP2(times), data = DF)
  
  
  # design matrix for the survival model
  W <- cbind(time0,BMI,LH,AH,SEX,RACE)
  
  ################################################
  
  # simulate random effects
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  
  
  # simulate longitudinal responses
  id <- rep(1:n, each = K)
  eta.y1 <- as.vector(X %*% betas1 + rowSums(Z * b[id,1:3]))
  y1 <- rnorm(length(id), eta.y1, sigma.y1)
  
  eta.y2 <- as.vector(X %*% betas2 + rowSums(Z * b[id,4:6]))
  y2 <- rnorm(length(id), eta.y2, sigma.y2)
  
  eta.y3 <- as.vector(X %*% betas3 + rowSums(Z * b[id,7:9]))
  y3 <- rnorm(length(id), eta.y3, sigma.y3)
  
  
  # parameters for the survival model
  gammas <- c("baseline_age"=gamma0, "BMI"=gamma1,"LH"=gamma2,"AH"=gamma3, "SEX"=gamma4,"RACE" = gamma5)
  shape=2;scale=18
  # simulate CVD_ADJ times
  eta.t <- as.vector(W %*% gammas)
  
  invS <- function (t, u, j) {
    h <- function (s) {
      P2 <-legP2(s)
      XX=ZZ<- cbind(1,P2[,1],P2[,2])
      
      fV1 <- as.vector(XX %*% betas1 + rowSums(ZZ * b[rep(j, nrow(ZZ)),1:3]))
      fV2 <- as.vector(XX %*% betas2 + rowSums(ZZ * b[rep(j, nrow(ZZ)),4:6]))
      fV3 <- as.vector(XX %*% betas3 + rowSums(ZZ * b[rep(j, nrow(ZZ)),7:9]))
      
      exp(log(shape) + (shape - 1) * log(s)-shape*log(scale)+ eta.t[j] +
          fV1 * alpha1+fV2 * alpha2+fV3 * alpha3)
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  
  u <- runif(n)
  trueTimes<- numeric(n)
  for (j in 1:n) {
    Up <- 5000
    tries <- 5
    Root <- try(uniroot(invS, interval = c(0.00001, Up), u = u[j], j = j)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 5000
      Root <- try(uniroot(invS, interval = c(0.0001, Up), u = u[j], j = j)$root, TRUE)
    }
    trueTimes[j] <- if (!inherits(Root, "try-error")) Root else NA
  }
  #(Root <- try(uniroot(invS, interval = c(0, 30), u = u[1], j = 1)$root, TRUE))
  
  na.ind <- !is.na(trueTimes)
  trueTimes <- trueTimes[na.ind]
  summary(trueTimes);length(trueTimes)
  ###Features
  
  W <- W[na.ind, , drop = FALSE]
  long.na.ind <- rep(na.ind, each = K)
  y1 <- y1[long.na.ind]
  y2 <- y2[long.na.ind]
  y3 <- y3[long.na.ind]
  X <- X[long.na.ind, , drop = FALSE]
  Z <- Z[long.na.ind, , drop = FALSE]
  DF <- DF[long.na.ind, ]
  n <- length(trueTimes)
  
  #mean.Cens<-8
  Ctimes <- t.max#runif(n, 0, 2 * mean.Cens)
  Time <- pmin(trueTimes, Ctimes)
  CVD_DTH_ADJ_sim <- as.numeric(trueTimes <= Ctimes) # CVD_ADJ indicator
  table(CVD_DTH_ADJ_sim)
  
################################################
  
  # keep the nonmissing cases, i.e., drop the longitudinal measurements
  # that were taken after the observed CVD_ADJ time for each subject.
  ind <- times[long.na.ind] <= rep(Time, each = K)
  y1 <- y1[ind]
  length(y1)
  X <- X[ind, , drop = FALSE]
  Z <- Z[ind, , drop = FALSE]
  y2 <- y2[ind]
  y3 <- y3[ind]
  id <- id[long.na.ind][ind]
  id <- match(id, unique(id))
  
  
  dat <- DF[ind,]
  dat$id <- id
  dat$y1 <- y1
  dat$y2 <- y2
  dat$y3 <- y3
  dat$Time <- Time[id]
  dat$CVD_DTH_ADJ_sim <- CVD_DTH_ADJ_sim[id]
  
  dat$Age_death<-dat$Time
  dat.id <- dat %>%group_by(id) %>%filter(time == min(time))
  dat1=subset(dat,time>=6&CVD_DTH_ADJ_sim==1)
  table(dat.id$CVD_DTH_ADJ_sim);table(dat$CVD_DTH_ADJ_sim);length(unique(dat1$id))
  
  
  
  
  
  t_m=15
  Xtime=Ztime=model.matrix(~legP2(Time),data=dat.id)
  dat.id$V1<- Xtime%*% betas1 + rowSums(Ztime,b[id,1:3])
  dat.id$V2<- Xtime%*% betas2 + rowSums(Ztime,b[id,4:6])
  dat.id$V3<- Xtime%*% betas3 + rowSums(Ztime,b[id,7:9])
  
#table(dat.id$CVD_DTH_ADJ_sim);table(dat.id$CVD_DTH_UNADJ)
#length(dat.id$time);length(dat$time)

data.id1=subset(dat.id,CVD_DTH_ADJ==1)#TRUE CVD
#table(data.id1$CVD_DTH_UNADJ)
data.id0=subset(dat.id,CVD_DTH_ADJ==0)#TRUE NOT
#table(data.id0$CVD_DTH_UNADJ)
library(rstanarm)

# Fit the model
fit1 <- stan_glm(CVD_DTH_UNADJ ~ Age_death+BMI+LH+AH+SEX+RACE+V1+V2+V3+HF_UNADJ_num+MI_UNADJ_num+STRK_UNADJ_num, 
                 data = data.id1, family = binomial(), chains = 1, iter = 1000)
# Summary
#summary(fit1)
# Extract posterior linear predictions (on logit scale)
eta_draws <- posterior_linpred(fit1, transform = FALSE)
# Convert to probabilities using logistic transformation
prob_draws <- plogis(eta_draws)  # matrix: [iterations x observations]
# For each observation, average across posterior draws
probs1 <- colMeans(prob_draws) 
set.seed(seed)
U<-runif(length(probs1),0,1)
data.id1$CVD_DTH_UNADJ_sim=as.numeric(probs1>U)
#table(data.id1$CVD_DTH_UNADJ);table(data.id1$CVD_DTH_UNADJ_sim)

# Fit the model
fit2 <- stan_glm(CVD_DTH_UNADJ ~ BMI+LH+AH+SEX+RACE+V1+V2+V3+HF_UNADJ_num+MI_UNADJ_num+STRK_UNADJ_num, 
                 data = data.id0, family = binomial(), chains = 1, iter = 1000)
# Summary
#summary(fit2)
# Extract posterior linear predictions (on logit scale)
eta_draws <- posterior_linpred(fit2, transform = FALSE)
# Convert to probabilities using logistic transformation
prob_draws <- plogis(eta_draws)  # matrix: [iterations x observations]
# For each observation, average across posterior draws
probs0 <- colMeans(prob_draws) 
set.seed(seed)
U<-runif(length(probs0),0,1)
data.id0$CVD_DTH_UNADJ_sim=as.numeric(probs0>U)
#table(data.id0$CVD_DTH_UNADJ);table(data.id0$CVD_DTH_UNADJ_sim)

data.id_EST=rbind(data.id1,data.id0)
#table(data.id_EST$CVD_DTH_UNADJ_sim)
data.id_EST1=data.id_EST[,c("id","CVD_DTH_UNADJ_sim")]


simdata<-dplyr:::left_join(dat,data.id_EST1,by="id")
simdata=simdata[,c("id","base_age","time","Time","BMI","LH","AH","SEX","RACE","y1","y2","y3",
                   "CVD_DTH_ADJ_sim","CVD_DTH_UNADJ_sim","Age_death")]
simdata.id <- simdata[!duplicated(simdata$id),]

N_A=length(simdata.id$id)
ID_train=simdata.id$id[1:(0.75*N_A)];ID_test=simdata.id$id[(round(0.75*N_A)+1):N_A]
#length(ID_train);length(ID_test)

data_sim_train[[i]]<-subset(simdata,id%in% ID_train)#
data_sim_test[[i]]<-subset(simdata,id%in% ID_test)#
}

table(data_sim_test[[3]]$CVD_DTH_UNADJ_sim)


