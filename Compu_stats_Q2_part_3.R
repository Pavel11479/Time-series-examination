data <- read.csv("Downloads/infl.csv")
inflFR<-data$inflFR
year<-data$year
T1 <- length(inflFR)
infl <- as.matrix(inflFR[2:T1])
laginfl <- as.matrix(inflFR[1:T1-1])
year <- as.matrix(year[2:T1])
T <- length(year)
infl <- infl-mean(infl)
laginfl <- laginfl-mean(laginfl)



y = t(infl)
lagy = t(laginfl)
yhat <- matrix(0, nrow = 1, ncol = T)
Shat <- matrix(0, nrow = 1, ncol = T)

#-------------------------------------------------------------------------------
#Kalman Filter

kalman = function(y, F, G, Q, H, R, mu0, Sigma0) {
  dy = nrow(y)
  T = ncol(y)
  dx = length(mu0)
  I = diag(dx)
  ## INITIALIZATION ##
  
  
  mu.p = matrix(0, nrow = dx, ncol = T)
  Sigma.p = array(0, c(dx, dx, T))
  mu.f = matrix(0, nrow = dx, ncol = T)
  Sigma.f = array(0, c(dx, dx, T))
  mu.s = matrix(0, nrow = dx, ncol = T)
  Sigma.s = array(0, c(dx, dx, T))
  
  
  ## FORWARD RECURSION ## Time 1
  mu.p[, 1] = F %*% mu0
  Sigma.p[, , 1] = F %*% Sigma0 %*% t(F) + G %*% Q %*% t(G)
  nu = y[, 1] - H[, 1] %*% mu.p[, 1]
  S = H[, 1] %*% Sigma.p[, , 1] %*% t(H[, 1]) + R
  K = Sigma.p[, , 1] %*% t(H[, 1]) %*% solve(S)
  mu.f[, 1] = mu.p[, 1] + K %*% nu
  Sigma.f[, , 1] = (I - K %*% H[, 1]) %*% Sigma.p[, , 1]
  yhat[, 1] = H[, 1] %*% mu.p[, 1]
  Shat[, 1] = H[, 1] %*% Sigma.p[, , 1] %*% t(H[, 1]) + R
  
  # Time 2:T
  for (t in (2:T)) {
    # Prediction
    mu.p[, t] = F %*% mu.f[, t - 1]
    Sigma.p[, , t] = F %*% Sigma.f[, , t - 1] %*% t(F) + G %*% Q %*% t(G)
    # Update
    nu = y[, t] - H[, t] %*% mu.p[, t]
    S = H[, t] %*% Sigma.p[, , t] %*% t(H[, t]) + R
    K = Sigma.p[, , t] %*% t(H[, t]) %*% solve(S)
    mu.f[, t] = mu.p[, t] + K %*% nu
    Sigma.f[, , t] = (I - K %*% H[, t]) %*% Sigma.p[, , t]
    yhat[, t] = H[, t] %*% mu.p[, t]
    Shat[, t] = H[, t] %*% Sigma.p[, , t] %*% t(H[, t]) + R
  }
  ## BACKWARD RECURSION ##
#  mu.s[, T] = mu.f[, T]
#  Sigma.s[, , T] = Sigma.f[, , T]
#  for (t in (T - 1):1) {
#    J = Sigma.f[, , t] %*% t(F) %*% solve(Sigma.p[, , t + 1])
#    mu.s[, t] = mu.f[, t] + J %*% (mu.s[, t + 1] - mu.p[, t + 1])
#    Sigma.s[, , t] = Sigma.f[, , t] + J %*% (Sigma.s[, , t + 1] - Sigma.p[,
#                                                                          , t + 1]) %*% t(J)
#  }
  llk <- sum(dnorm(y, mean = yhat, sd = sqrt(Shat), log = TRUE))
  return(list(mu.f = mu.f, Sigma.f = Sigma.f, mu.p = mu.p, Sigma.p = Sigma.p,
              mu.s = mu.s, Sigma.s = Sigma.s, yhat = yhat, Shat = Shat, 
              llk = llk))
}

#-------------------------------------------------------------------------------
#a)

Q = 0.01
R = 4
F = 1
G = 1 
H = lagy
Sigma0 = 1
mu0 = 0 

results.KF = kalman(y, F, G, Q, H, R, mu0, Sigma0)

mu.f = results.KF$mu.f
Sigma.f = results.KF$Sigma.f
se.f = sqrt(Sigma.f[,,])
mu.f0 = c(mu0,mu.f)
se.f0 = c(sqrt(Sigma0),se.f)
alpha = 0.05
cv95 = qnorm(1-alpha/2)
CIupper = mu.f0+cv95*se.f0
CIlower = mu.f0-cv95*se.f0
time0 = c(1956,year)





#-------------------------------------------------------------------------------
#Checking the model works as intended (part b)

results.KFnew = kalman(y, F, G, Q, H, R, mu0, Sigma0)
yhatideal <- results.KFnew$yhat
Shatideal <- results.KFnew$Shat
results.KFnew$llk
plot(t(y))
points(t(yhatideal), col = 'blue')
points(mu.f0, col = 'red')



y
yhatideal
mu.f0
#-------------------------------------------------------------------------------
#Numerical Search for optimal Parameters

llkfin <- rep(0, 15000)
t=1
for (sn in seq(0.000001, 3, 0.001)) {
  for (st in seq(0.000001, 5, 0.1)) {
    llkfin[t] <- kalman(y, F, G, sn, H, st, mu0, Sigma0)$llk
    t=t+1
  }
}
max(llkfin)

#-------------------------------------------------------------------------------
# Re-writing the algoritm with scalars instead of matrices

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

#-------------------------------------------------------------------------------
# Plot
yhatideal0 = c(NA,yhatideal)
plot(year,y,cex=1,col='darkgreen',pch=15,ylim=c(-13,15), 
     xlab = "Year", ylab = "De-meaned Inflation rate", cex.lab=1.3, 
     cex.axis = 1.3, cex.main = 1.3,
     main = "Filtering mean for Beta coefficient")
points(time0,mu.f0,cex=0.8,col='red',pch=16)
points(time0,CIupper,col='blue',type="l",lty=2,lwd=2)
points(time0,CIlower,col='blue',type="l",lty=2,lwd=2)
#points(time0,yhatideal0,cex=1,col='red',pch=16)
legend("topright", legend = c("Y", "Filtering mean", "95% CI"), 
       col = c("darkgreen", "red", "blue"), pch = c(15,16,NA),
       lty = c(NA, NA, 2), cex = 1
)
#legend("topright", legend = c("Y", "Yhat"), 
#      col = c("darkgreen", "red"), pch = c(15,16), cex = 1
#)

max(mu.f0)
points(time0, yhatideal0, cex=0.8, col='orange', type = "l")

