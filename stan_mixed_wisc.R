wisc <- read.table("C:/Users/RJacobucci/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
wisc <- read.table("C:/Users/jacobucc/Documents/GitHub/EDM_Labs/2015/wisc4vpe.dat")
names(wisc)<- c("V1","V2","V4","V6","P1","P2","P4", "P6", "Moeducat")


library(rstan)

wisc.verb <- wisc[,c(1:4,9)]

# create subset for plotting
ntot <- nrow(wisc.verb)    # total number of observations
wisc.verb.sel <- wisc.verb[sample(ntot, 30), ]

wisc.long <- reshape(wisc.verb, varying = c("V1", "V2", "V4", "V6"), v.names = "verbal",
                     times = c(1, 2, 4, 6), direction = "long")

wisc.long.sel <- reshape(wisc.verb.sel, varying = c("V1", "V2", "V4", "V6"),
                         v.names = "verbal", times = c(1, 2, 4, 6),
                         direction = "long")
head(wisc.long,3)
names(wisc.long)[2] <- "grade"
names(wisc.long.sel)[2] <- "grade"

head(wisc.long)




library(nlme)
mix1 <- lme(fixed = verbal ~ grade, random = ~ grade | id, data = wisc.long, method="ML" )
summary(mix1) # get same estimates as in LGM, notice SD not VAR









dat <- list(
  verbal = wisc.long[,"verbal"],
  grade = wisc.long[,"grade"],
  N = nrow(wisc.long),
  R = diag(2),
  m = c(0,0))

stanmodel <- "
data {
  int N; 
  vector[2] m;
  matrix[2,2] R;
  real verbal[N]; 
  real grade[N];
} 
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  corr_matrix[2] psi;
}
model{
  real mu;
  for (i in 1:N){
    mu = beta[1] + beta[2]*grade[i];
    verbal ~ normal(mu,sigma);
  }
  beta ~ multi_normal(m,psi);
  psi ~ lkj_corr(2);
}
"

mixed.model=stan(model_code=stanmodel,
               data = dat,chains=1,
               pars=c("beta","sigma","psi"))
print(mixed.model)




