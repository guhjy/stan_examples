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
  npen=110,
  nu=1)


mod.stan <-"
functions {
// square root of a vector (elementwise)
vector sqrt_vec(vector x) {
vector[dims(x)[1]] res;
for (m in 1:dims(x)[1]){
res[m] <- sqrt(x[m]);
}
return res;
}
}

data{
int npen; // number penalized parameters
int N; // sample size
matrix[N,6] X; // data matrix of order [N,P]
matrix[N,npen] cov; // data matrix of order [N,P]
real<lower=1> nu; // degrees of freedom for the half t-priors

}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
vector<lower=0>[6] sigma;
vector[5] lam;
vector[6] alpha;
real<lower=0> psi;

// auxiliary variables for the variance parameters
vector[npen] z;
real<lower=0> r1_global;
real<lower=0> r2_global;
vector<lower=0>[npen] r1_local;
vector<lower=0>[npen] r2_local;

}



transformed parameters{

vector[6] mu[N];
vector[N] mu2;
real<lower=0> tau;
vector<lower=0>[npen] lambda;
vector[npen] beta;

tau = r1_global * sqrt(r2_global);
lambda = r1_local .* sqrt_vec(r2_local);
beta = z .* lambda*tau;



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


// https://arxiv.org/pdf/1508.02502.pdf
// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
z ~ normal(0, 1);
r1_local~ normal(0.0, 1.0);
r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
// half cauchy for tau
r1_global ~ normal(0.0, 1.0);
r2_global ~ inv_gamma(0.5, 0.5);



for(i in 1:N){  
for(j in 1:6){
X[i,j] ~ normal(mu[i,j],pow(sigma[j],0.5));
}
FS[i] ~ normal(mu2[i],psi);
}
}
"


fa.model=stan(model_code=mod.stan,
              data = dat,chains=1,
              pars=c("sigma","lam","beta","psi","alpha"))

print(fa.model)
MCMCplot(fa.model,params="beta")
