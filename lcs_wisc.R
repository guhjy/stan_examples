library(blavaan)

wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


X <- wisc[,c("V1","V2","V4","V6")]
data = list()
#data$alpha <- c(1,1,2,2)
X = as.matrix(X)
N <- nrow(X)

dat <- list(
  N = N,
  X = X,
  t = 4)








stanmodel <- "
data {
  int N; 
  int t; 
  vector[t] X[N]; 
} 
parameters {
  vector[2] beta;  
  //real pi;
  real<lower=0> sigma;
  cov_matrix[2] phi;
} 
transformed parameters {
  vector[2] s[N];
  vector[t] y[N];
  vector[t-1] d[N];

for (i in 1:N){
  y[i,1] = s[i,2];

  for (tt in 2:t){
    d[i,tt-1] = s[i,1];
    y[i,tt] = d[i,tt-1]+y[i,tt-1];
  }
}
  
} 

model {

for (i in 1:N){
  s[i,1:2] ~ multi_normal_prec(beta, phi);
  X[i,1] ~ normal(y[i,1],pow(sigma, -0.5));
  

  for (tt in 2:t){
    X[i,tt] ~ normal(y[i,tt], pow(sigma, -0.5));
  }
}

sigma ~ gamma(.001,.001);
//pi~normal(0,2);
beta[1]~normal(-2,10);
beta[2]~normal(20,.10);
phi~wishart(2,diag_matrix(rep_vector(1.0,2)));
}
"


library(rstan)

fit <- stan(model_code = stanmodel, model_name = "LCS", 
          data = dat,chains=1,
        pars=c("beta","sigma"))
