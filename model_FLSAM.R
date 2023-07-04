# model_SAM.R - DESC
# sol.27.4_benchmark/model_SAM.R

# Copyright Iago MOSQUEIRA (WMR), 2023
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("model")


library(FLSAM)
library(icesAdvice)

load('data/data.rda')

# SETUP control

control <- FLSAM.control(stock, indices)

# ADD survey correlation

# BTS
control@cor.obs[2, 1:9] <- c(rep(101, 2), rep(102, 7))
# SNS
control@cor.obs[3, 1:5] <- c(rep(201,2), rep(202,3))
control@cor.obs.Flag[2:3] <- as.factor("AR")

# F random walks
control@f.vars["catch unique",] <- c(1,2,rep(3,6),rep(4,2))

# F correlation structure
control@cor.F <- 2

# Observation variances
control@obs.vars["catch unique",] <- c(1,2,3,3,rep(4,6))
control@obs.vars["BTS", ac(1:10)] <- c(101, rep(102,4),rep(103,2),rep(104,3))
control@obs.vars["SNS", ac(1:6)] <- c(rep(201,3),rep(202,3))

# Catchabilities
control@catchabilities["BTS",ac(1:10)] <- c(1,1,2,3,4,5,5,6,6,6)
control@catchabilities["SNS",ac(1:6)] <- c(1,2,3,4,4,4) + 101

control <- update(control)

# RUN

fit <- FLSAM(stock, indices, control) 

run <- fit + stock

# RETRO

control@residuals <- FALSE
retro <- FLStocks(lapply(retro(stock, indices, control, retro=5), '+', stock))

monhsrho <- data.frame(qname=c("F", "SSB", "REC"),
  mr=format(c(mohn(mohnMatrix(retro, metric=fbar)),
    mohn(mohnMatrix(retro, metric=ssb)),
    mohn(mohnMatrix(retro, metric=rec))), digits=3),
  y=c(0.05, 5000, 25000))

# RESIDUALS

resids <- residuals(fit)
resids$std.res[which(is.na(resids$std.res))] <- 0

save(fit, run, resids, retro, monhsrho, file="model/flsam.rda", compress="xz")

# SAVE control to stockassessment

load("data/sam.rda")

conf <- ctrl2conf(control, data)
par <- stockassessment::defpar(data, conf)

save(data, conf, par, file="data/sam.rda", compress="xz")

# RUN using stockassessment

samfit <- stockassessment::sam.fit(data, conf, par)
samresid <- residuals(samfit)

save(data, conf, par, samfit, samresid, file="model/sam.rda", compress="xz")




# stock
plot(fit, futureYrs=F)

# retro w/error & Mohn's rho

dat <- cbind(rbindlist(setNames(lapply(c("ssb", "rec", "fbar"), function(x)
  do.call(x, list(fit))), c("SSB", "Rec", "F")), idcol="qname"), run="base")
colnames(dat)[3] <- "data"

rmets <- rbindlist(lapply(retro, function(x) data.table(as.data.frame(metrics(x,
  list(SSB=ssb, F=fbar, Rec=rec)), drop=TRUE))), idcol="run")
err <- rbindlist(list(rmets, dat), fill=TRUE)

monhsrho <- data.frame(qname=c("F", "SSB", "Rec"),
  mr=format(c(mohn(mohnMatrix(retro, metric=fbar)),
    mohn(mohnMatrix(retro, metric=ssb)),
    mohn(mohnMatrix(retro, metric=rec))), digits=3),
  y=c(0.05, 5000, 25000))

ggplot(err[year > 2000], aes(x=year)) +
  geom_line(aes(y=data, colour=run, group=run)) + ylim(c(0,NA)) +
  geom_ribbon(data=err[year > 2000 & run == "base", ],
    aes(ymin=lbnd, ymax=ubnd), alpha=0.3, colour="grey") +
  facet_grid(qname~., scales="free_y") +
  theme(legend.position="none") +
  geom_text(data=monhsrho, aes(x=2016, y=y, label=mr, fill=NULL, colour=NULL))

# residuals
plot(rsd)

# fits

# selectivities
ggplot(catch.sel(res), aes(x=age, y=data)) + geom_line() +
  facet_wrap(~year)

catch <- catchabilities(run)

print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
         scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
         type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
         subset=fleet %in% names(tune),
         main="Survey catchability parameters",ylab="Catchability",xlab="Age"))

# Fs

dev.off()
