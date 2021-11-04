
#Packages installation 
install.packages("sandwich")
library("sandwich")

#--------------------------------
#Data preparation

dvis <- read.csv("Downloads/dvis (2).csv")
dvisreg <- dvis
dvisreg$female <- as.factor(dvis$female)
dvisreg$age <- as.numeric(dvis$age)
dvisreg$hhninc <- as.numeric(dvis$hhninc)
dvisreg$hhkids <- as.factor(dvis$hhkids)
dvisreg$educyrs <- as.numeric(dvis$educyrs)
dvisreg$docvis <- as.numeric(dvis$docvis)
dvisreg$addins <- as.factor(dvis$addins)
str(dvisreg)

#---------------------------------

#a)

model_glm <- glm(docvis ~ female + hhninc + hhkids + age + educyrs + addins, 
    family=poisson(link=log), data=dvisreg)
summary(model_glm)
std_se <- summary(model_glm)$coefficients[,2]
hc_se <- sqrt(diag(vcovHC(model_glm)))
std_se ; hc_se
tsglm <- summary(model_glm)$coefficients[,1]/summary(model_glm)$coefficients[,2]
tssandwich <- summary(model_glm)$coefficients[,1]/hc_se


tsglm
2*pnorm(-abs(tsglm))
tssandwich
2*pnorm(-abs(tssandwich))
#-------------------------------------------------------------------------------

#b)
n = length(dvisreg$docvis)
B = 10000
intercept.boot.paired <- rep(0,B)
beta.female.boot.paired <- rep(0,B)
beta.age.boot.paired <- rep(0,B)
beta.hhninc.boot.paired <- rep(0,B)
beta.hhkids.boot.paired <- rep(0,B)
beta.educyrs.boot.paired <- rep(0,B)
beta.addins.boot.paired <- rep(0,B)
for (j in 1:B) {
  ind <- sample(1:n, n, replace=TRUE)
  female.boot.paired <- dvisreg$female[ind]
  age.boot.paired <- dvisreg$age[ind]
  hhninc.boot.paired <- dvisreg$hhninc[ind]
  hhkids.boot.paired <- dvisreg$hhkids[ind]
  educyrs.boot.paired <- dvisreg$educyrs[ind]
  addins.boot.paired <- dvisreg$addins[ind]
  docvis.boot.paired <- dvisreg$docvis[ind]
  fit.boot.paired <- glm(docvis.boot.paired ~ 
                           female.boot.paired + age.boot.paired + 
                           hhninc.boot.paired + hhkids.boot.paired + 
                           educyrs.boot.paired + addins.boot.paired, 
                         family=poisson(link=log))
  intercept.boot.paired[j]  <- fit.boot.paired$coefficients[1]
  beta.female.boot.paired[j]  <- fit.boot.paired$coefficients[2]
  beta.age.boot.paired[j]  <- fit.boot.paired$coefficients[3]
  beta.hhninc.boot.paired[j]  <- fit.boot.paired$coefficients[4]
  beta.hhkids.boot.paired[j]  <- fit.boot.paired$coefficients[5]
  beta.educyrs.boot.paired[j]  <- fit.boot.paired$coefficients[6]
  beta.addins.boot.paired[j]  <- fit.boot.paired$coefficients[7]
}

intercept.boot.paired.sq <- rep(0, B)
  for (i in 1:B) {
    intercept.boot.paired.sq[i] <- (intercept.boot.paired[i] - 
                                      mean(intercept.boot.paired))^2
  }
stderr.intercept.boot.paired <- sqrt(mean(intercept.boot.paired.sq))
stderr.intercept.boot.paired

beta.female.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.female.boot.paired.sq[i] <- (beta.female.boot.paired[i] - 
                                    mean(beta.female.boot.paired))^2
}
stderr.beta.female.boot.paired <- sqrt(mean(beta.female.boot.paired.sq))
stderr.beta.female.boot.paired

beta.age.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.age.boot.paired.sq[i] <- (beta.age.boot.paired[i] - 
                                      mean(beta.age.boot.paired))^2
}
stderr.beta.age.boot.paired <- sqrt(mean(beta.age.boot.paired.sq))
stderr.beta.age.boot.paired

beta.hhninc.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.hhninc.boot.paired.sq[i] <- (beta.hhninc.boot.paired[i] - 
                                      mean(beta.hhninc.boot.paired))^2
}
stderr.beta.hhninc.boot.paired <- sqrt(mean(beta.hhninc.boot.paired.sq))
stderr.beta.hhninc.boot.paired

beta.hhkids.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.hhkids.boot.paired.sq[i] <- (beta.hhkids.boot.paired[i] - 
                                      mean(beta.hhkids.boot.paired))^2
}
stderr.beta.hhkids.boot.paired <- sqrt(mean(beta.hhkids.boot.paired.sq))
stderr.beta.hhkids.boot.paired

beta.educyrs.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.educyrs.boot.paired.sq[i] <- (beta.educyrs.boot.paired[i] - 
                                      mean(beta.educyrs.boot.paired))^2
}
stderr.beta.educyrs.boot.paired <- sqrt(mean(beta.educyrs.boot.paired.sq))
stderr.beta.educyrs.boot.paired

beta.addins.boot.paired.sq <- rep(0, B)
for (i in 1:B) {
  beta.addins.boot.paired.sq[i] <- (beta.addins.boot.paired[i] - 
                                       mean(beta.addins.boot.paired))^2
}
stderr.beta.addins.boot.paired <- sqrt(mean(beta.addins.boot.paired.sq))
stderr.beta.addins.boot.paired

stderr.boot.paired <- c(stderr.intercept.boot.paired, 
                        stderr.beta.female.boot.paired, 
                        stderr.beta.hhninc.boot.paired, 
                        stderr.beta.hhkids.boot.paired, 
                        stderr.beta.age.boot.paired, 
                        stderr.beta.educyrs.boot.paired, 
                        stderr.beta.addins.boot.paired)
stderr.boot.paired
std_se ; hc_se
tsboot <- summary(model_glm)$coefficients[,1]/stderr.boot.paired
tssandwich; tsboot
2*pnorm(-abs(tssan))
#-------------------------------------------------------------------------------

#c)
coeff_m_glm <- matrix(c(model_glm$coefficients[4], model_glm$coefficients[6]), 
                      nrow = 2, ncol =1)
cov_m_glm <- matrix(c(vcov(model_glm)[4,4], vcov(model_glm)[6,4], 
                      vcov(model_glm)[4,6], vcov(model_glm)[6,6]), 
                    nrow = 2, ncol =2)
Wald_glm <- t(coeff_m_glm)%*%solve(cov_m_glm)%*%coeff_m_glm


cov_m_hc <- matrix(c(vcovHC(model_glm)[4,4], vcovHC(model_glm)[6,4], 
                      vcovHC(model_glm)[4,6], vcovHC(model_glm)[6,6]), 
                    nrow = 2, ncol =2)
Wald_hc <- t(coeff_m_glm)%*%solve(cov_m_hc)%*%coeff_m_glm

B <- 10000 
sqsum.boot.paired <- rep(0,B)
for (i in 1:B) {
  sqsum.boot.paired[i] <- ((beta.hhkids.boot.paired[i] - 
                             mean(beta.hhkids.boot.paired))*
    (beta.educyrs.boot.paired[i] - mean(beta.educyrs.boot.paired)))
}
cov.boot.paired <- sum(sqsum.boot.paired)/9999
cov.boot.paired
cov_m_boot <- matrix(c((stderr.beta.hhkids.boot.paired^2), cov.boot.paired, 
                       cov.boot.paired, (stderr.beta.educyrs.boot.paired^2)), 
                   nrow = 2, ncol =2)
Wald_boot <- t(coeff_m_glm)%*%solve(cov_m_boot)%*%coeff_m_glm
Wald_glm; Wald_hc; Wald_boot
qchisq(.95, df=2)
pchisq(Wald_glm, df=2, lower.tail=FALSE)
pchisq(Wald_hc, df=2, lower.tail=FALSE)
pchisq(Wald_boot, df=2, lower.tail=FALSE)

#-------------------------------------------------------------------------------

#d) 
B <- 10000 
teststat.boot.paired <- rep(0,B)
for (i in 1:B) {
  teststat.boot.paired[i] <- t(matrix(c(beta.hhkids.boot.paired[i], 
                                        beta.educyrs.boot.paired[i]), 
                                      nrow = 2, ncol = 1) - coeff_m_glm)%*%
    solve(cov_m_hc)%*%
    (matrix(c(beta.hhkids.boot.paired[i], 
             beta.educyrs.boot.paired[i]), 
           nrow = 2, ncol = 1)- coeff_m_glm)
}
p.hat = mean(Wald_hc[1,1] <= teststat.boot.paired)
p.hat
pchisq(Wald_hc, df=2, lower.tail=FALSE)
Wald_hc

