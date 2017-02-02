library(rstan); library(lavaan)

HS <- HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),]


mod <- "
f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f1 ~ sex + grade + ageyr
"
out <- sem(mod, HS)
summary(out)

X <- HS[,7:15]
cov <- HS[,c(2,3,6)]


X = as.matrix(X)
cov= as.matrix(cov)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  cov = cov)


mod.stan <-"
data{
int N; // sample size
vector[9] X[N]; // data matrix of order [N,P]
vector[3] cov[N]; // data matrix of order [N,P]
}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
real SD;
vector[9] sigma;
vector[8] lam;
vector[3] beta;
}

transformed parameters{

vector[9] mu[N];
vector[N] mu2;


for(i in 1:N){
mu[i,1] = 1*FS[i];
mu[i,2] = lam[1]*FS[i];
mu[i,3] = lam[2]*FS[i];
mu[i,4] = lam[3]*FS[i];
mu[i,5] = lam[4]*FS[i];
mu[i,6] = lam[5]*FS[i];
mu[i,7] = lam[6]*FS[i];
mu[i,8] = lam[7]*FS[i];
mu[i,9] = lam[8]*FS[i];

mu2[i] = beta[1]*cov[i,1] + beta[2]*cov[i,2] + beta[3]*cov[i,3];

}
}

model{
SD ~ gamma(2,2);
sigma ~ gamma(2,2);
beta ~ normal(0,1);

for(i in 1:N){  
  for(j in 1:9){
    X[i] ~ normal(mu[i],pow(sigma[j],0.5));
  }
FS[i] ~ normal(mu2[i],SD);
}
}
"


fa.model=stan(model_code=mod.stan,
              data = dat,chains=1,
              pars=c("sigma","lam","beta","SD"))