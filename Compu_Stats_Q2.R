#1) 
library("ggplot2")
library("dplyr")
install.packages("Hmisc")
library("Hmisc")

infl <- read.csv("Downloads/infl.csv")
#Plot of inflation
ggplot(infl, aes(x=year, y=inflFR)) +
  geom_line() +
  geom_hline(yintercept = mean(infl$inflFR), color="blue")+
  labs(title = "Inflation rate", x = "Year", y = "Inflation Rate")+
#  ggtitle("Inflation rate")+
#  xlab("Year")+
#  ylab("Inflation Rate")+
  theme_classic()+
  theme(text=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))
#No Trend - as expected

hist(infl$inflFR, main = "Inflation rate histogram", xlab = "Inflation rate", 
     cex.lab=1.3, cex.axis = 1.3, cex.main = 1.3)
qqnorm(infl$inflFR, main = "Normal QQ plot of inflation", cex.lab=1.3, 
       cex.axis = 1.3, cex.main = 1.3)
qqline(infl$inflFR)
mean(infl$inflFR)
sd(infl$inflFR)
acf(infl$inflFR, main = "Autocorrelation plot of inflation", cex.lab=1.3, 
    cex.axis = 1.3, cex.main = 1.3)
pacf(infl$inflFR, main = "Partial autocorrelation plot of inflation", 
     cex.lab=1.3, 
     cex.axis = 1.3, cex.main = 1.3)

#Checking For a unit root
infl$inflFR.ch <- infl$inflFR - infl$inflFR.L
infl$inflFR.ch.L <- Lag(infl$inflFR.ch, +1)
infl$inflFR.ch.L2 <- Lag(infl$inflFR.ch, +2)
infl$inflFR.ch.L3 <- Lag(infl$inflFR.ch, +3)
infl$inflFR.ch.L4 <- Lag(infl$inflFR.ch, +4)
infl

urtest <- lm(inflFR.ch ~ inflFR.L + inflFR.ch.L + inflFR.ch.L2 + 
               inflFR.ch.L3 + inflFR.ch.L4, data = infl)
summary(urtest)

#-------------------------------------------------------------------------------
#2)

infl$inflFR.L <- Lag(infl$inflFR, +1)
infl
AR1 <- lm(infl$inflFR ~ infl$inflFR.L)
summary(AR1)
plot(year, resid(AR1), type = "l", main = "Residual plot of AR(1)", 
     cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3, xlab = "Year", ylab = "Residual")
acf(resid(AR1), main = "Autocorrelation plot of residuals", 
    cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)
pacf(resid(AR1), main = "Partial autocorrelation plot of residuals", 
    cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)

qqnorm(resid(AR1), main = "Normal QQ plot of residuals", cex.lab=1.3, cex.axis = 1.3, 
       cex.main = 1.3)
qqline(resid(AR1))


#Now look at autocorrelation in the residuals (can I do it, given the series
#is non stationary?)
#Breusch-Godfrey test
resid <- resid(AR1)
resid.L <- Lag(resid, +1)
sercorr <- lm(resid ~ resid.L)
summary(sercorr)
1-pchisq(0.007527*63, 1)
0.007527*63
pchisq(-500, 1, lower.tail = FALSE)











