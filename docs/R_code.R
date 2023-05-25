
# R code for "Estimating and presenting flexible interaction effects on the additive and multiplicative scale from a Cox model"
# May 2023

library(survival)
library(rms)
library(readxl)
library(interactionRCS)
library(flexsurv)
library(eventglm)
library(riskRegression)
library(pec)
library(ggplot2)

risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = time)$surv)
}


# fit Cox model with a single binary predictor (curr_smoke)
fit1<-coxph(Surv(days2miistr, miistrfu) ~ curr_smoke,
            data=data, ties="breslow",y=TRUE,x=TRUE)
summary(fit1)


# Direct estimation of absolute risk difference (riskRegression package)
riskreg <- ate(fit1, data = data, treatment = "curr_smoke", times = 1096.75)
as.numeric(summary(riskreg)$diffRisk$estimate)

# Predict survival probabilities
data$surv<- 1- predictSurvProb(fit1,newdata=data,time=1096.75)

# Average risk over levels of curr_smoke
risk0<-mean(data$surv[data$curr_smoke==0])
risk1<-mean(data$surv[data$curr_smoke==1])
arr<-risk1-risk0
arr


### 1) model with interaction and no splines

fit_int1 <- coxph(Surv(days2miistr, miistrfu) ~ curr_smoke*age,
                  data=data, ties="breslow",y=TRUE,x=TRUE)

newdat = expand.grid(age = seq(min(data$age), max(data$age), by = .1),
                     curr_smoke = unique(data$curr_smoke) )

# predict on the grid
newdat$surv<- 1- predictSurvProb(fit_int1,newdata=newdat,time=1096.75)

newdat<-merge(newdat[newdat$curr_smoke==1,],newdat[newdat$curr_smoke==0,],by="age")
newdat$ard<-newdat$surv.x-newdat$surv.y

# Plot ARD over leves of age

ggplot() +
  geom_line(data = newdat, aes(x = age, y = ard)) +
  ylab("Absolute Risk Difference at 3 years") +
  xlab("age")

## 2) add splines
fit_int2 <- coxph(Surv(days2miistr, miistrfu) ~ curr_smoke*rcs(age,3),
                  data=data, ties="breslow",y=TRUE,x=TRUE)

newdat = expand.grid(age = seq(50, 70, by = .1),
                     curr_smoke = unique(data$curr_smoke) )

# predict on the grid
newdat$surv<- 1- predictSurvProb(fit_int2,newdata=newdat,time=1096.75)

newdat<-merge(newdat[newdat$curr_smoke==1,],newdat[newdat$curr_smoke==0,],by="age")
newdat$ard<-newdat$surv.x-newdat$surv.y


ggplot() +
  geom_line(data = newdat, aes(x = age, y = ard)) +
  ylab("Absolute Risk Difference at 3 years") +
  xlab("age")

# 3) add CIs with boostrap

library(boot)
grid.ad <- expand.grid(age = seq(50, 70, by = .1),
                       curr_smoke = unique(data$curr_smoke) )

boot.ad <- function(data, i, grid) {
  dat.boot <- data[i, ]
  fit_int2 <- coxph(Surv(days2miistr, miistrfu) ~ curr_smoke*rcs(age,3),
                    data=dat.boot, ties="breslow",y=TRUE,x=TRUE)
  
  grid.ad$surv <- risk(fit_int2,newdata=grid.ad,time=1096.75)
  grid.ad2 <- merge(grid.ad[grid.ad$curr_smoke==1,],grid.ad[grid.ad$curr_smoke==0,],by="age")
  grid.ad2$surv.x-grid.ad2$surv.y
}

boot.out <- boot(data, boot.ad, R = 1000, grid = grid.ad)
newdat<-cbind(newdat,
              t(apply(boot.out$t, 2, quantile, probs = c(0.025, 0.975)))) # 95% quantile-based CIs


# plot
ggplot() +
  geom_line(data = newdat, aes(x = age, y = ard)) +
  geom_line(data = newdat, aes(x = age, y = `2.5%`),linetype="dashed") +
  geom_line(data = newdat, aes(x = age, y = `97.5%`),linetype="dashed")+
  ylab("Absolute Risk Difference at 3 years") +
  xlab("age")

# 4) add CIs based on SE

risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, se.fit = T, conf.int = F), times = time)$surv)
}
se = function(model, newdata, time) {
  as.numeric(summary(survfit(model, newdata = newdata, se.fit = T, conf.int = F), times = time)$std.err)
}


fit_int2 <- coxph(Surv(days2miistr, miistrfu) ~ curr_smoke*rcs(age,3),
                  data=data, ties="breslow",y=TRUE,x=TRUE)

newdat = expand.grid(age = seq(50, 70, by = .1),
                     curr_smoke = unique(data$curr_smoke) )

# predict on the grid
newdat$surv<- risk(fit_int2,newdat,1096.75)
newdat$se<-se(fit_int2,newdat,1096.75)


newdat<-merge(newdat[newdat$curr_smoke==1,],newdat[newdat$curr_smoke==0,],by="age")
newdat$ard<-newdat$surv.x-newdat$surv.y
newdat$sediff<-sqrt(newdat$se.x^2+newdat$se.y^2)
newdat$lb<-newdat$ard-1.96*newdat$sediff
newdat$ub<-newdat$ard+1.96*newdat$sediff

# plot
ggplot() +
  geom_line(data = newdat, aes(x = age, y = ard)) +
  geom_line(data = newdat, aes(x = age, y = lb),linetype="dashed") +
  geom_line(data = newdat, aes(x = age, y = ub),linetype="dashed")+
  ylab("Absolute Risk Difference at 3 years") +
  xlab("age")


