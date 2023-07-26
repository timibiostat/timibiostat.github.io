#############################################################################
## AUTHOR:  Andrea Bellavia
## DATE:    7/25/2023
## ----------------------------------------------------------------
## PURPOSE: Data simulations and figures for JSM 2023 poster:  FLEXIBLE SPLINES MODELING OF ABSOLUTE RISK IN SURVIVAL ANALYSIS
## ----------------------------------------------------------------
############################################################################## */

# packages
pacman::p_load(dplyr,pec,survival, survminer, ggplot2,rms,simsurv,eventglm)


set.seed(12345)
cov <- data.frame(id = 1:10000,
                  x1 = rnorm(10000, 0, 1),
                  x2 = rnorm(10000, 0, 1))

cov$x2sq<-cov$x2^2


# Simulate the event times
# X1 with log-linear effect, x2 with quadratic

dat <- simsurv(lambdas = 0.1,
               gammas = 1,
               betas = c(x1=0.2,x2sq = 0.15,x2=0.2),
               x = cov,
               maxt = 5)


# Merge the simulated event times onto covariate data frame
dat <- merge(cov, dat)


dat$age<-dat$x1*5+55
dat$bmi<-dat$x2*2+32
#write_xlsx(dat,"dat_10000_new_v2.xlsx")

hist(dat$age)
hist(dat$bmi)

risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata), times = time)$surv)
}

low = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata), times = time)$lower)
}

up = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata), times = time)$upper)
}

## Example 1: age
# Cox model

fit.coxph <- coxph(Surv(eventtime, status) ~ rcs(age,3) , ties="breslow", data = dat,x=TRUE,y=TRUE)
summary(fit.coxph)

fit.coxph2 <- cph(Surv(eventtime, status) ~ rcs(age,3) , data = dat,surv=T)


newdat = expand.grid(age = seq(45, 65, by = .1))

## 1) function based on survfit
newdat$surv<- risk(fit.coxph,newdata=newdat,time=3)
newdat$upper<- up(fit.coxph,newdata=newdat,time=3)
newdat$lower<- low(fit.coxph,newdata=newdat,time=3)

## 2) rms prediction
newdat$surv3<-1-survest(fit.coxph2, newdat,
                        times=3, conf.int=.95)$surv
newdat$lower3<-1-survest(fit.coxph2, newdat,
                         times=3, conf.int=.95)$upper
newdat$upper3<-1-survest(fit.coxph2, newdat,
                         times=3, conf.int=.95)$lower



# GLM (pseudo-values/ eventGLM)
fit.glm <- cumincglm(Surv(eventtime, status) ~ rcs(age,3), time = 3, data = dat)
summary(fit.glm)



critval <- 1.96 ## approx 95% CI
preds<-predict.glm(fit.glm, newdat=newdat, type="link", se.fit = TRUE)

newdat$upr <- preds$fit + (critval * preds$se.fit)
newdat$lwr <- preds$fit - (critval * preds$se.fit)
newdat$surv_glm <- preds$fit

# Figure for poster


cols    <- c( "c2" = "red", "c1" = "blue")
ggplot() +
  geom_line(data = newdat, aes(x = age, y = surv3,col="c1"),size=1.05) +
  geom_line(data = newdat, aes(x = age, y = surv_glm,col="c2"),size=1.05) +
  geom_ribbon(data = newdat, aes(x=age, ymin = lower3, ymax = upper3), fill = "dodgerblue",alpha=0.4)+
  geom_ribbon(data = newdat, aes(x=age, ymin = lwr, ymax = upr), fill = "orange",alpha=0.2)+
  scale_x_continuous("Age (years)")+scale_y_continuous("Absolute Risk of Overall Mortality at 3 years")+
  scale_color_manual(name = "",
                     breaks = c("c1", "c2"),
                     values = cols,
                     labels = c("AR from Cox Model", "AR from GLM (pseudo-values)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 15),
        axis.title = element_text(size = 15))



## Example 2: BMI
# Cox model

fit.coxph <- coxph(Surv(eventtime, status) ~ rcs(bmi,3) , ties="breslow", data = dat,x=TRUE,y=TRUE)
summary(fit.coxph)

fit.coxph2 <- cph(Surv(eventtime, status) ~ rcs(bmi,3) , data = dat,surv=T)


newdat = expand.grid(bmi = seq(28, 38, by = .1))

## 1) function based on survfit

newdat$surv<- risk(fit.coxph,newdata=newdat,time=3)
newdat$upper<- up(fit.coxph,newdata=newdat,time=3)
newdat$lower<- low(fit.coxph,newdata=newdat,time=3)

## 2) rms prediction
newdat$surv3<-1-survest(fit.coxph2, newdat,
                        times=3, conf.int=.95)$surv
newdat$lower3<-1-survest(fit.coxph2, newdat,
                         times=3, conf.int=.95)$upper
newdat$upper3<-1-survest(fit.coxph2, newdat,
                         times=3, conf.int=.95)$lower



# GLM (pseudo-values/ eventGLM)

fit.glm <- cumincglm(Surv(eventtime, status) ~ rcs(bmi,3), time = 3, data = dat)
summary(fit.glm)



critval <- 1.96 ## approx 95% CI
preds<-predict.glm(fit.glm, newdat=newdat, type="link", se.fit = TRUE)

newdat$upr <- preds$fit + (critval * preds$se.fit)
newdat$lwr <- preds$fit - (critval * preds$se.fit)
newdat$surv_glm <- preds$fit

# figure for poster
cols    <- c( "c2" = "red", "c1" = "blue")

ggplot() +
  geom_line(data = newdat, aes(x = bmi, y = surv3,col="c1"),size=1.2) +
  geom_line(data = newdat, aes(x = bmi, y = surv_glm,col="c2"),size=1.2) +
  geom_ribbon(data = newdat, aes(x=bmi, ymin = lower3, ymax = upper3), fill = "dodgerblue",alpha=0.4)+
  geom_ribbon(data = newdat, aes(x=bmi, ymin = lwr, ymax = upr), fill = "orange",alpha=0.2)+
  scale_x_continuous("BMI (kg/m2)")+scale_y_continuous("Absolute Risk of Overall Mortality at 3 years")+
  scale_color_manual(name = "",
                     breaks = c("c1", "c2"),
                     values = cols,
                     labels = c("AR from Cox Model", "AR from GLM (pseudo-values)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 15),
        axis.title = element_text(size = 15))
