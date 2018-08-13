library(rstan); library(lavaan)

HS <- HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),]


mod <- "
f1 =~ NA*x1 + x2 + x3 + 1*x4 + x5 + x6 + x7 + x8 + x9
f1 ~ sex + grade + ageyr
#f1~~1*f1
"
out <- sem(mod, HS,meanstructure=T)
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
matrix[N,9] X; // data matrix of order [N,P]
matrix[N,3] cov; // data matrix of order [N,P]
}

parameters{
vector[N] FS; // factor scores, matrix of order [N,D]
vector<lower=0>[9] sigma;
vector[8] lam;
vector[9] alpha;
vector[3] beta;
real<lower=0> psi;
}

transformed parameters{

vector[9] mu[N];
vector[N] mu2;


for(i in 1:N){
mu[i,1] = alpha[1] + lam[1]*FS[i];
mu[i,2] = alpha[2] + lam[2]*FS[i];
mu[i,3] = alpha[3] + lam[3]*FS[i];
mu[i,4] = alpha[4] + 1*FS[i];
mu[i,5] = alpha[5] + lam[4]*FS[i];
mu[i,6] = alpha[6] + lam[5]*FS[i];
mu[i,7] = alpha[7] + lam[6]*FS[i];
mu[i,8] = alpha[8] + lam[7]*FS[i];
mu[i,9] = alpha[9] + lam[8]*FS[i];

mu2[i] =  beta' * cov[i,]';

}
}

model{
sigma ~ gamma(2,2);
beta ~ normal(0,10);
alpha ~ normal(2,10);
psi ~ gamma(2,2);
lam ~ normal(.5,10);

for(i in 1:N){  
  for(j in 1:9){
    X[i,j] ~ normal(mu[i,j],pow(sigma[j],0.5));
  }
FS[i] ~ normal(mu2[i],psi);
}
}
"


fa.model=stan(model_code=mod.stan,
              data = dat,chains=3,cores=3,
              pars=c("sigma","lam","beta","psi","alpha"))

print(fa.model)
