library(survival)
library(timeROC)


# we evaluate our proximity as a prognostic biomarker.
ROC.proximity<-timeROC(T=data_for_roc$TimeSurv,
                  delta=data_for_roc$vitalStatus,
                  marker=as.numeric(data_for_roc$proximity_score),
                  cause=1,
                  weighting="marginal",
                  times=quantile(data_for_roc$TimeSurv,probs=seq(0.2,0.8,0.2)),
                  iid=TRUE)
# we evaluate ELN as a prognostic biomarker.
ROC.ELN<-timeROC(T=data_for_roc$TimeSurv,
                  delta=data_for_roc$vitalStatus,
                  marker=as.numeric(data_for_roc$ELN_score),
                  cause=1,
                  weighting="marginal",
                  times=quantile(data_for_roc$TimeSurv,probs=seq(0.2,0.8,0.2)),
                  iid=TRUE)

print(ROC.proximity)
print(ROC.ELN)

compare(ROC.proximity,ROC.ELN) #compute p-values of comparison tests

pdf(file = './plots/AUC_curve.pdf')
# plot AUC curve for albumin only with pointwise confidence intervals
# and simultaneous confidence bands
plotAUCcurve(ROC.proximity,conf.int=TRUE,conf.band=TRUE)
# plot AUC curve for albumin and bilirunbin with pointwise confidence intervals
plotAUCcurve(ROC.proximity,conf.int=TRUE,col="red")
plotAUCcurve(ROC.ELN,conf.int=TRUE,col="blue",add=TRUE)
legend("bottomright",c("Proxmity","ELN 2017"),col=c("red","blue"),lty=1,lwd=2)
dev.off()
