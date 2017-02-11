library(rstan); library(lavaan)
library(MCMCvis)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rstanmulticore)

xx = list()
for(i in 1:100){
  xx[i] = paste("x", i, sep = "")
}
uu = list()
for(i in 1:length(xx)){
  uu[i] = paste("0*",xx[i],sep="")
}


reg.list = paste(uu,collapse="+")

sim.list = list()


sim.list[[1]] = "f1 =~ y1 + 0.8*y2 + 1.2*y3 + 0.8*y4 + 0.5*y5 + 1.5*y6
"

sim.list[[2]] = "
f1 ~ 1*x1001 + 1*x1002 + 1*x1003 + 1*x1004 + 1*x1005
f1 ~ .2*x1006 + .2*x1007 + .2*x1008 + .2*x1009 + .2*x1010
"
sim.list[[3]] = paste(paste("f1"," ~ "), reg.list)

pop.mod = " "
for(i in 1:length(sim.list)){
  pop.mod = paste(pop.mod,sim.list[[i]],sep="\n")
}

reg.list2 = paste(xx,collapse="+")
run.list = list()
run.list[[1]] = "f1 =~ y1 + y2 + y3 + y4 + y5 + y6
"
run.list[[2]] = paste(paste("f1"," ~ "), reg.list2)
run.list[[3]] = "
f1 ~ x1001 + x1002 + x1003 + x1004 + x1005 + x1006 + x1007 + x1008 + x1009 + x1010
"
run.mod = " "
for(k in 1:length(run.list)){
  run.mod = paste(run.mod,run.list[[k]],sep="\n")
}

dat <- simulateData(pop.mod,sample.nobs=200,fixed.x=TRUE)

head(dat)


X <- dat[,1:6]
cov <- dat[,7:116]


X = as.matrix(X)
cov= as.matrix(cov)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  cov = cov,
  npen=110)


mod.stan <-"
data{
int npen; // number penalized parameters
int N; // sample size
matrix[N,6] X; // data matrix of order [N,P]
matrix[N,npen] cov; // data matrix of order [N,P]
}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
vector<lower=0>[6] sigma;
vector[5] lam;
vector[6] alpha;
vector[npen] beta;
real<lower=0> psi;
real<lower=0> lambda;

vector<lower=0>[npen] tau;
vector<lower=0>[npen] ptau;
}

transformed parameters{

vector[6] mu[N];
vector[N] mu2;
vector<lower=0>[npen] pen;
real<lower=0> ppsi;

ppsi = pow(lambda,-1);
for(j in 1:npen){
  pen[j] = ppsi*ptau[j];
}



for(i in 1:N){
mu[i,1] = alpha[1] + 1*FS[i];
mu[i,2] = alpha[2] + lam[1]*FS[i];
mu[i,3] = alpha[3] + lam[2]*FS[i];
mu[i,4] = alpha[4] + lam[3]*FS[i];
mu[i,5] = alpha[5] + lam[4]*FS[i];
mu[i,6] = alpha[6] + lam[5]*FS[i];


  mu2[i] =  beta'*cov[i,]';


}
}

model{
sigma ~ gamma(2,2);
alpha ~ normal(0,1);
psi ~ gamma(2,2);

lambda ~ gamma(1,.05);
tau ~ gamma(1,.1);
ptau ~ gamma(1,tau/2);

for(j in 1:npen){
beta[j] ~ normal(0,pen[j]);
}


/////

for(i in 1:N){  
  for(j in 1:6){
    X[i,j] ~ normal(mu[i,j],pow(sigma[j],0.5));
  }
FS[i] ~ normal(mu2[i],psi);
}
}
"


init.list <- function(){
  list(sigma=rep(1,6),
  psi=1,
  lam=rep(1,5),
  beta=rep(0,110),
  alpha = rep(0,6))
}

# https://www.ariddell.org/horseshoe-prior-with-stan.html

fa.model=stan(model_code=mod.stan,
              data = dat,chains=1,init=init.list,
              pars=c("sigma","lam","beta","psi","alpha"))

print(fa.model)
MCMCplot(fa.model,params="beta")
