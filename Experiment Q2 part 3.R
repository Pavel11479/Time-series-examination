data <- read.csv("Downloads/infl.csv")
inflFR<-data$inflFR
year<-data$year
T1 <- length(inflFR)
infl <- inflFR[2:T1]
laginfl <- inflFR[1:T1-1]
year <- year[2:T1]
T <- length(year)
infl <- infl-mean(infl)
laginfl <- laginfl-mean(laginfl)
summary(lm(infl$inflFR ~ infl$year))

kalman <- function(mu0, Sigma0, H, Q, R, y) {
  mu.p <- rep(0, T)
  mu.f <- rep(0, T)
  Sigma.p <- rep(0,T)
  Sigma.f <- rep(0, T)
  S <- rep(0,T)
  yhat <- rep(0,T)
  
  #First Iteration
  mu.p[1] = mu0
  Sigma.p[1] = Sigma0 + Q
  S[1] = H[1]*Sigma.p[1]*H[1]+R
  mu.f[1] = mu.p[1] + Sigma.p[1]*H[1]*(1/S[1])*(y[1]-H[1]*mu.p[1])
  Sigma.f[1] = (1 - Sigma.p[1]*H[1]*(1/S[1])*H[1])*Sigma.p[1]
  yhat[1] = H[1]*mu.p[1]
  
  #Following Iterations:
  for (t in (2:T)) {
    mu.p[t] = mu.f[t-1]
    Sigma.p[t] = Sigma.f[t-1] + Q
    S[t] = H[t]*Sigma.p[t]*H[t]+R
    mu.f[t] = mu.p[t] + (Sigma.p[t]*H[t]*(1/S[t]))*(y[t]-H[t]*mu.p[t])
    Sigma.f[t] = (1 - Sigma.p[t]*H[t]*(1/S[t])*H[t])*Sigma.p[t]
    yhat[t] = H[t]*mu.p[t]
  }
  
  llk <- sum(dnorm(y, mean = yhat, sd = sqrt(S), log = TRUE))
  return(list(mu.f = mu.f, Sigma.f = Sigma.f, yhat = yhat, S = S, llk = llk, 
              Sigma.p = Sigma.p , mu.p = mu.p))
}

result <- kalman(0, 1, laginfl, 0.01, 4, infl)
result$llk
result$mu.f

llkfin <- rep(0, 1000000)
t=1
for (i in seq(0.0000001, 1, 0.001)) {
  for (j in seq(0.0000001, 10, 0.01)){
    llkfin[t] <- kalman(0, 1, laginfl, i, j, infl)$llk
    t=t+1
  }
}
max(llkfin)
match(max(llkfin), llkfin)
kalman(0, 1, laginfl, 0.01, 4.01, infl)$llk
