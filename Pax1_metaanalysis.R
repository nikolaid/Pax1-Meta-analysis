library(metafor)

Pax1 <- data.frame (study = c(1, 2, 3, 4, 5, 6,7), author = c("Huang", "Kan", "Lai", "Lai", "Lai", "Lin","Wang"), year = c(2010, 2014, 2008, 2010, 2014, 2011,2014), tpos = c(22, 32, 143, 70, 48, 15,25), tneg = c(24, 11, 9, 20, 44, 27,17), cpos = c(0, 6, 0, 0, 19, 13,0), cneg = c(17, 47, 41, 52, 180, 107,39))
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=Pax1)
dat

res <- rma(yi, vi, data=dat, method="FE")
res
predict(res, transf=exp, digits=3)
res <- rma(yi, vi, data=dat, method="DL")
res
predict(res, transf=exp, digits=3)

args(rma)
result.md <-rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, method="DL", measure="OR",data=Pax1)
names(result.md)

contributions <- 1/result.md$vi/sum(1/result.md$vi) * 100
cbind(contributions)

contributions <- 1/result.md$vi/sum(1/result.md$vi) * 100
par(mar = c(5, 10, 5, 5))
barplot(contributions, names = Pax1$study,xlim = c(0, 50), las = 2, horiz = T,col = "royalblue")

forest(result.md)
args(forest.rma)

leave1out(result.md)
cases <- leave1out(result.md)
which(cases$I2 == min(cases$I2))
sum(cases$I2 < 30)

cbind(exp(cases$estimate), cases$pval < 0.05)

result.mdr <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg,data = Pax1, mods =Pax1[,"year"],measure = "OR", method = "ML")
summary(result.mdr)

exp(c(result.md$b, result.mdmr$b[1]))
c(result.md$I2, result.mdmr$I2)
(result.md$I2 - result.mdmr$I2)/result.md$I2 * 100

funnel(result.md)
result.rd <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg,data = Pax1, measure = "OR",method = "ML")
trimfill(result.rd)
funnel(trimfill(result.rd))
funnel(result.rd)

value <- fsn(y = result.md$yi, v = result.md$vi)
value

regtest(res)
ranktest(res)


library (mada)

SROC <- data.frame(study = c(1, 2, 3, 4, 5, 6,7), author = c("Huang", "Kan", "Lai", "Lai", "Lai", "Lin","Wang"), year = c(2010, 2014, 2008, 2010, 2014, 2011,2014), TP = c(22,32,143,70,48,15,25), FN = c(28, 11,9,20,44,27,17), FP = c(0,6,0,0,19,13,0), TN = c(17,47,41,52,180,107,39))
SROC

forest(madauni(SROC))
forest(SROC.d, type = "spec", xlab = "Specificity")
forest(SROC.d, type = "sens", xlab = "Sensitivity")

madad(SROC, TP, FN, FP, TN, level = 0.95, correction = 0.5, correction.control = "all", method = "wilson", yates = TRUE,suppress = TRUE)

SROC.d<- madad(SROC)
print(SROC.d, digits=2)

madauni(SROC, type = "DOR", method = "DSL", suppress = TRUE)

SROC.uni <- madauni(SROC)
summary(SROC.uni)

SROC.uni_low <- madauni(SROC, correction = 0.1)
SROC.uni_single <- madauni(SROC,correction.control = "single")
confint(SROC.uni_single)

ROCellipse(SROC)

ROCellipse(SROC)
mslSROC(SROC, add = TRUE)
msl_SROC <- mslSROC(SROC, col = 3, lwd = 3, lty = 3)
msl_SROC$A2

phm(data = SROC, subset=NULL,TP="TP", FN="FN", FP="FP", TN="TN", correction = 0.5, correction.control = "all", hetero = TRUE, estimator = "APMLE", l = 100)

(fit <- phm(SROC))
summary(fit)
plot(fit)

print(SROC) 
summary(SROC, level = 0.95)

sroc(fit, fpr = 1:99/100)

sroc.SROC <- sroc(fit)
plot(sroc.SROC, type = "l")

reitsma(SROC, subset=NULL, formula = NULL,TP="TP", FN="FN", FP="FP", TN="TN", alphasens = 1, alphafpr = 1, correction = 0.5, correction.control = "all", method = "reml", control = list())
reitsma(SROC, subset=NULL, formula = NULL,TP="TP", FN="FN", FP="FP", TN="TN", alphasens = 1, alphafpr = 1, correction = 0.5, correction.control = "all", method = "ml", control = list())


plot(fit)
cr.SROC <- ROCellipse(fit)
sroc.SROC <- sroc(fit)
plot(cr.SROC$ROCellipse, type = "l", xlim = c(0,1), ylim = c(0,1))
points(cr.SROC$fprsens)
lines(sroc.SROC)

ROCellipse(SROC)
rsSROC(Dementia, add = TRUE)
rsSROC(SROC, add = TRUE)
rs_SROC <- rsSROC(SROC, col = 3, lwd = 3, lty = 3,plotstudies = TRUE)
rs_SROC$lambda
rs_SROC$aa
rs_SROC$bb

mcmc_sum <- SummaryPts(fit, n.iter = 10^3)
mcmc_sum
