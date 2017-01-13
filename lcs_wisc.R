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
  "N" = N,
  "X" = X,
  "T" = 4)








stanmodel <- "
data {
int N; 
int T; 
matrix[N,T] X; 
} 
parameters {
vector[2] beta;  
real pi;
real<lower=0> sigma; 
cov_matrix[2] phi;
} 
//transformed parameters {
//cov_matrix[2] phi;
//phi = diag_matrix(rep_vector(1.0,2)); 
//} 

model {
matrix[N,2] s;
matrix[N,T] y;
matrix[N,T-1] d;

for (i in 1:N){
s[i,1:2] ~ multi_normal(beta, phi);
X[i,1]~ normal(y[i,1],pow(sigma, -0.5));
y[i,1] = s[i,2];

for (t in 2:T){
X[i,t] ~ normal(y[i,t], pow(sigma, -0.5));
d[i,t-1] = pi*y[i,t-1]+s[i,1];
y[i,t] = d[i,t-1]+y[i,t-1];
}
}

sigma ~ gamma(1,10);
pi~normal(0,10);
beta[1]~normal(-2,10);
beta[2]~normal(20,10);
phi~wishart(10,diag_matrix(rep_vector(1.0,2)));
}
"


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(model_code = stanmodel, 
        model_name = "LCS", 
          data = dat,chains=1)
