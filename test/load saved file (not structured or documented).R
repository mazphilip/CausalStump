#load saves

# full and limited #####
library(xtable)

load("/Users/philip/Downloads/H11_100_27_august.RData")

#for(i in 1:nr_models){
#  results[c(1,2,12,13) + (i-1)*nr_metrix,]=log(results[c(1,2,12,13) + (i-1)*nr_metrix,])
#}

load("/Users/philip/Downloads/H11_100_limited_27_august.RData")

#results[(1:nr_metrix) + (2-1)*nr_metrix,] = results_CF
#results[(1:nr_metrix) + (3-1)*nr_metrix,] = results_VTRF
#results[(1:nr_metrix) + (4-1)*nr_metrix,] = results_CFRF
#

mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train",
                    "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test"),
                  c("BART","CF","VT-RF","CF-RF","GP","BCF","BART.PS","GP.PS","TP2.PS","TP2000.PS"))


nr_metrix = length(mydimnames[[1]])
nr_models = length(mydimnames[[2]])

out = round(cbind(matrix(apply(results,1,mean),nr_metrix,nr_models),sqrt(matrix(apply(results,1,var),nr_metrix,nr_models))),2)#,matrix(apply(results,1,min),11,7),matrix(apply(results,1,max),11,7)),2)
out2 = matrix(NaN,nr_metrix*2,nr_models)
colnames(out2) <- mydimnames[[2]]
for(i in  1:2){
  for(j in 1:nr_models)
    #    all "means" for method "j" , <-
    out2[(1:nr_metrix -1) * 2 + i, j] <- c(out[,j + nr_models*(i-1) ])
}
colnames(out2) <- mydimnames[[2]]
rownames(out2) <- rep(mydimnames[[1]],each=2)

out2
out.BART.PS = out2[,"BART.PS"]
out.GP = out2[,"GP"]
xtable(out2)


mycol80 <- rgb(191, 191, 191, max=255, alpha = (100-80)*255/100)
mycol60 <- rgb(191, 191, 191, max=255, alpha = (100-60)*255/100)
mycol30 <- rgb(191, 191, 191, max=255, alpha = (100-30)*255/100)
mycol10 <- rgb(191, 191, 191, max=255, alpha = (100-10)*255/100)
mycol5 <- rgb(191, 191, 191, max=255, alpha = (100-5)*255/100)

#PEHE.all train BART
hist(results[1,])
i=5 #2: CF, 5: GP
hist(results[1+nr_metrix*(i-1),])

#PEHE.all test BART
hist(results[12,])
i=5 #2: CF, 5: GP
hist(results[12+nr_metrix*(i-1),])

#ATE bias test BART
hist(results[16+nr_metrix*(i-1),])
i=5 #2: CF, 5: GP, 6: BCF


#load("/Users/philip/Downloads/H11_100_27_august.RData")
mymax=0
mymin=0
## ATE bias #####
for(i in 1:nr_models){
  tmp =  max(results[16+nr_metrix*(i-1),])
  if (tmp>mymax){ mymax = tmp }
  tmp =  min(results[16+nr_metrix*(i-1),])
  if (tmp<mymin){ mymin = tmp }
}
mymax=ceiling(mymax)
mymin=floor(mymin)
myseq = seq(mymin,mymax,by=0.2)

window=2
plot_w = 3; plot_h = 1.3

for(i in 1:nr_models){

  {
    hist.test  = hist(results[16+nr_metrix*(i-1),],breaks = myseq)
    hist.train = hist(results[5+nr_metrix*(i-1),],breaks = myseq)#,breaks = seq(-6,6,by=0.2))
  }

  tikz(file = sprintf("IHDP_ATEhisto_full_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.train, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.test, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[16+nr_metrix*(i-1),]),lty=2)
  abline(v=mean(results[5+nr_metrix*(i-1),]),lty=3)

  #legend("topleft",c("train","test","mean"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA),col=c(NA,NA,1),border=c(NA,mycol30,NA),density=c(1000,20,NA),lty=c(NA,NA,2))
  dev.off()
}

mymax=0
mymin=0
## ATT bias #####
for(i in 1:nr_models){
  tmp =  max(results[17+nr_metrix*(i-1),])
  if (tmp>mymax){ mymax = tmp }
  tmp =  min(results[17+nr_metrix*(i-1),])
  if (tmp<mymin){ mymin = tmp }
}
mymax=ceiling(mymax)
mymin=floor(mymin)
myseq = seq(mymin,mymax,by=0.2)

window=2
plot_w = 3; plot_h = 1.3 #1.5

for(i in 1:nr_models){
  {
    hist.test  = hist(results[17+nr_metrix*(i-1),],breaks = myseq)
    hist.train = hist(results[6+nr_metrix*(i-1),],breaks = myseq)#,breaks = seq(-6,6,by=0.2))
  }

  tikz(file = sprintf("IHDP_ATThisto_full_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.train, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.test, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[17+nr_metrix*(i-1),]),lty=2)
  #legend("topleft",c("train","test","mean"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA),col=c(NA,NA,1),border=c(NA,mycol30,NA),density=c(1000,20,NA),lty=c(NA,NA,2))
  dev.off()
}

#### limited ####
load("/Users/philip/Downloads/H11_100_limited_27_august.RData")
#load("/Users/philip/Downloads/H11_100_limited_29_august.RData")
mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train",
                    "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test"),
                  c("BART","CF","VT-RF","CF-RF","GP","BCF","BART.PS","GP.PS","TP2.PS","TP2000.PS"))

limited_results = results

nr_metrix = length(mydimnames[[1]])
nr_models = length(mydimnames[[2]])
mymax=0
mymin=0
## ATE bias #####
for(i in 1:nr_models){
  tmp =  max(results[16+nr_metrix*(i-1),])
  if (tmp>mymax){ mymax = tmp }
  tmp =  min(results[16+nr_metrix*(i-1),])
  if (tmp<mymin){ mymin = tmp }
}
for(i in 1:nr_models){
  tmp =  max(results[5+nr_metrix*(i-1),])
  if (tmp>mymax){ mymax = tmp }
  tmp =  min(results[5+nr_metrix*(i-1),])
  if (tmp<mymin){ mymin = tmp }
}
mymax=ceiling(mymax)
mymin=floor(mymin)
myseq = seq(mymin,mymax,by=0.2)

window=2
plot_w = 3; plot_h = 1.3

for(i in 1:nr_models){

    hist.test  = hist(results[16+nr_metrix*(i-1),],breaks = myseq)
    hist.train = hist(results[5+nr_metrix*(i-1),],breaks = myseq)#,breaks = seq(-6,6,by=0.2))


  tikz(file = sprintf("IHDP_ATEhisto_limited_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.train, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="ATE error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.test, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[16+nr_metrix*(i-1),]),lty=2)
  #legend("topleft",c("train","test","mean"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA),col=c(NA,NA,1),border=c(NA,mycol30,NA),density=c(1000,20,NA),lty=c(NA,NA,2))
  dev.off()
}

mymax=0
mymin=0
## ATT bias #####
for(i in 1:nr_models){
  tmp =  max(results[17+nr_metrix*(i-1),])
  if (tmp>mymax){ mymax = tmp }
  tmp =  min(results[17+nr_metrix*(i-1),])
  if (tmp<mymin){ mymin = tmp }
}
mymax=ceiling(mymax)
mymin=floor(mymin)
myseq = seq(mymin,mymax,by=0.2)

window=5
plot_w = 3; plot_h = 1.3

for(i in 1:nr_models){


    hist.test  = hist(results[17+nr_metrix*(i-1),],breaks = myseq)
    hist.train = hist(results[6 +nr_metrix*(i-1),],breaks = myseq)#,breaks = seq(-6,6,by=0.2))


  tikz(file = sprintf("IHDP_ATThisto_limited_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.train, col=mycol60, xlim=c(-window,window),ylim=c(0,100),main=NULL,xlab="ATE error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.test, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[17+nr_metrix*(i-1),]),lty=2)
  #legend("topleft",c("train","test","mean"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA),col=c(NA,NA,1),border=c(NA,mycol30,NA),density=c(1000,20,NA),lty=c(NA,NA,2))
  dev.off()
}

# weights #####

load("/Users/philip/Downloads/H11_100_limited_27_august.RData")
mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train",
                    "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test"),
                  c("BART","CF","VT-RF","CF-RF","GP","BCF","BART.PS","GP.PS","TP2.PS","TP2000.PS"))

nr_metrix = length(mydimnames[[1]])
nr_models = length(mydimnames[[2]])
i=5
mymax=ceiling(max(results[16+nr_metrix*(i-1),]))
mymin=floor(min(results[16+nr_metrix*(i-1),]))
mymax=max(c(mymax,abs(mymin))); mymin=-mymax
myseq = seq(mymin,mymax,by=0.2)

hist.gp.test  = hist(results[16+nr_metrix*(i-1),],breaks=myseq)
old_mean = mean(results[16+nr_metrix*(i-1),])

#now load robust simulations
#mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train","discard percent train",
#                    "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test","discard percent test"),
#                  c("BART w R1","BART w R2","BART w R3","GP w R1","GP w R2","GP weights","GP balancing","GP proper weights","GP proper balancing","BART w GP-R2"))
#nr_metrix = length(mydimnames[[1]])
#nr_models = length(mydimnames[[2]])
#
#load("/Users/philip/Downloads/H11_100_weights_add_29_august.RData")



#bart.add = results[216:(nr_metrix*nr_models),]

#load("/Users/philip/Downloads/H11_100_weights_29_august.RData")
#load("/Users/philip/Downloads/H11_33_weights_add_29_august.RData")
load("/Users/philip/Downloads/H11_100_weights_29_august.RData")
#load("/Users/philip/Downloads/H11_100_weights_5_sept.RData")

mydimnames = list(c("PEHE.all.train","PEHE.tt.train","RMSE.all.train","RMSE.tt.train","E.ate.train","E.att.train","SE.ate.train","SE.att.train","coverage.ite.train","coverage.ate.train","range.ate.train","discard percent train",
                    "PEHE.all.test","PEHE.tt.test","RMSE.all.test","RMSE.tt.test","E.ate.test","E.att.test","SE.ate.test","SE.att.test","coverage.ite.test","coverage.ate.test","range.ate.test","discard percent test"),
                  c("BART w R1","BART w R2","BART w R3","GP w R1","GP w R2","GP w R3","GP weights","GP bal","GP var1","GP ivar2","BART w GP-R2","BART GP-R3","GP BART R2"))

nr_metrix = length(mydimnames[[1]])
nr_models = length(mydimnames[[2]])


nr_metrix_old = length(out.BART.PS)/2
out.BART.PS_rep = rep(NA,2*nr_metrix)
out.BART.PS_rep[1:(nr_metrix-2)] = out.BART.PS[1:nr_metrix_old]
out.BART.PS_rep[(nr_metrix+1):(2*nr_metrix-2)] = out.BART.PS[(nr_metrix_old+1):(2*nr_metrix_old)]

out.GP_rep = rep(NA,2*nr_metrix)
out.GP_rep[1:(nr_metrix-2)] = out.GP[1:nr_metrix_old]
out.GP_rep[(nr_metrix+1):(2*nr_metrix-2)] = out.GP[(nr_metrix_old+1):(2*nr_metrix_old)]

out = round(cbind(matrix(apply(results,1,mean),nr_metrix,nr_models),sqrt(matrix(apply(results,1,var),nr_metrix,nr_models))),2)#,matrix(apply(results,1,min),11,7),matrix(apply(results,1,max),11,7)),2)
out2 = matrix(NaN,nr_metrix*2,nr_models)
colnames(out2) <- mydimnames[[2]]
for(i in  1:2){
  for(j in 1:nr_models)
    #    all "means" for method "j" , <-
    out2[(1:nr_metrix -1) * 2 + i, j] <- c(out[,j + nr_models*(i-1) ])
}
colnames(out2) <- mydimnames[[2]]
rownames(out2) <- rep(mydimnames[[1]],each=2)
out2 = cbind(out2,BART.PS =out.BART.PS_rep ,GP = out.GP_rep)
out2


#
xtable(out2[,c(14,1,2,3,11,12,15,4,5,6,13)])
#weighted methods
xtable(out2[,c(15,8,7,9,10)])

#limited_results
i=7 #BART with PS
bart_ate_limited_results = limited_results[c(5+nr_metrix_old*(i-1),16+nr_metrix_old*(i-1)),]
bart_att_limited_results = limited_results[c(6+nr_metrix_old*(i-1),17+nr_metrix_old*(i-1)),]
i=5 #BART with PS
gp_ate_limited_results = limited_results[c(5+nr_metrix_old*(i-1),16+nr_metrix_old*(i-1)),]
gp_att_limited_results = limited_results[c(6+nr_metrix_old*(i-1),17+nr_metrix_old*(i-1)),]

#ATE
myseq = seq(-10,10,by=0.2)
window=2
plot_w = 3; plot_h = 1.3

for(i in c(7,8,9,10)){

  hist.gp = hist(gp_ate_limited_results[1,],breaks = myseq)
  hist.weight = hist(results[5+nr_metrix*(i-1),],breaks = myseq)#5,6,17,18

  tikz(file = sprintf("IHDP_ATEhisto_weight_train_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.gp, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.weight, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[5+nr_metrix*(i-1),]),lty=3)
  abline(v=mean(gp_ate_limited_results[1,]),lty=2)
  #if(i==7){
  #  legend("topleft",c("rule","robust","mean","robust"),cex=0.8,
  #         fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
  #}
  dev.off()
}

for(i in c(7,8,9,10)){

  hist.gp = hist(gp_ate_limited_results[2,],breaks = myseq)
  hist.weight = hist(results[17+nr_metrix*(i-1),],breaks = myseq)#5,6,17,18

  tikz(file = sprintf("IHDP_ATEhisto_weight_test_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.gp, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.weight, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[17+nr_metrix*(i-1),]),lty=3)
  abline(v=mean(gp_ate_limited_results[2,]),lty=2)
  #legend("topleft",c("CS-GP","robust","mean","robust"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
  dev.off()
}

for(i in c(7,8,9,10)){

  hist.gp = hist(gp_att_limited_results[1,],breaks = myseq)
  hist.weight = hist(results[6+nr_metrix*(i-1),],breaks = myseq)#5,6,17,18

  tikz(file = sprintf("IHDP_ATThisto_weight_train_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.gp, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.weight, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[6+nr_metrix*(i-1),]),lty=3)
  abline(v=mean(gp_att_limited_results[1,]),lty=2)
  #legend("topleft",c("CS-GP","robust","mean","robust"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
  dev.off()
}

for(i in c(7,8,9,10)){

  hist.gp = hist(gp_att_limited_results[2,],breaks = myseq)
  hist.weight = hist(results[18+nr_metrix*(i-1),],breaks = myseq)#5,6,17,18

  tikz(file = sprintf("IHDP_ATThisto_weight_test_%d.tex",i), width = plot_w, height = plot_h)
  par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
  plot( hist.gp, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
  plot( hist.weight, col=mycol30, density=20, border=TRUE, add=T)
  #abline(v=0,lty=3,col="gray")
  abline(v=mean(results[18+nr_metrix*(i-1),]),lty=3)
  abline(v=mean(gp_att_limited_results[2,]),lty=2)
  #legend("topleft",c("CS-GP","robust","mean","robust"),cex=0.8,
  #       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
  dev.off()
}

# discarding ####

myseq = seq(-10,10,by=0.2)

# plots for discarding ###
window=2
plot_w = 3; plot_h = 1.3

####### ATE BART w G5% train
i=11
hist.orig = hist(bart_ate_limited_results[1,],breaks = myseq)
hist.new = hist(results[5+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_train_bart.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[5+nr_metrix*(i-1),]),lty=3)
abline(v=mean(bart_ate_limited_results[1,]),lty=2)
legend("topleft",c("no rule","rule","no rule","rule"),cex=0.8,
       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATE BART w G5% test
i=11
hist.orig = hist(bart_ate_limited_results[2,],breaks = myseq)
hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_test_bart.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[17+nr_metrix*(i-1),]),lty=3)
abline(v=mean(bart_ate_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATE GP w G5% train
i=5
hist.orig = hist(gp_ate_limited_results[1,],breaks = myseq)
hist.new = hist(results[5+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_train_gp.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[5+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_ate_limited_results[1,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATE GP w G5% test
i=5
hist.orig = hist(bart_ate_limited_results[2,],breaks = myseq)
hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_test_gp.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[17+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_ate_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATE GP w B10% train
i=13
hist.orig = hist(gp_ate_limited_results[1,],breaks = myseq)
hist.new = hist(results[5+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_train_gp_b10.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[5+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_ate_limited_results[1,]),lty=2)
#legend("topleft",c("orig.","B10pc","orig.","B10pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATE GP w B10% test
i=13
hist.orig = hist(bart_ate_limited_results[2,],breaks = myseq)
hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATEhisto_test_gp_b10.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[17+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_ate_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","B10pc","orig.","B10pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT BART w G5% train
i=11
hist.orig = hist(bart_att_limited_results[1,],breaks = myseq)
hist.new = hist(results[6+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_train_bart.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[6+nr_metrix*(i-1),]),lty=3)
abline(v=mean(bart_att_limited_results[1,]),lty=2)
legend("topleft",c("no rule","rule","no rule","rule"),cex=0.8,
       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT BART w G5% test
i=11
hist.orig = hist(bart_att_limited_results[2,],breaks = myseq)
hist.new = hist(results[18+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_test_bart.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[18+nr_metrix*(i-1),]),lty=3)
abline(v=mean(bart_att_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT GP w G5% train
i=5
hist.orig = hist(gp_att_limited_results[1,],breaks = myseq)
hist.new = hist(results[6+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_train_gp.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[6+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_att_limited_results[1,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT GP w G5% test
i=5
hist.orig = hist(bart_att_limited_results[2,],breaks = myseq)
hist.new = hist(results[18+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_test_gp.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[18+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_att_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","G5pc","orig.","G5pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT GP w B10% train
i=13
hist.orig = hist(gp_att_limited_results[1,],breaks = myseq)
hist.new = hist(results[6+nr_metrix*(i-1),],breaks = myseq)
#hist.new = hist(results[17+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_train_gp_b10.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[6+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_att_limited_results[1,]),lty=2)
#legend("topleft",c("orig.","B10pc","orig.","B10pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

####### ATT GP w B10% test
i=13
hist.orig = hist(bart_att_limited_results[2,],breaks = myseq)
hist.new = hist(results[18+nr_metrix*(i-1),],breaks = myseq)

tikz(file = sprintf("IHDP_discard_ATThisto_test_gp_b10.tex",i), width = plot_w, height = plot_h)
par(mfrow=c(1,1),mar=c(2.5,2.5,0.1,0.1),mgp=c(1.3,0.3,0),cex.lab=0.8,cex.axis=0.8,tcl=-0.3)
plot( hist.orig, col=mycol60, xlim=c(-window,window),ylim=c(0,50),main=NULL,xlab="Error",border=FALSE,cex=0.8)  # first histogram
plot( hist.new, col=mycol30, density=20, border=TRUE, add=T)
#abline(v=0,lty=3,col="gray")
abline(v=mean(results[18+nr_metrix*(i-1),]),lty=3)
abline(v=mean(gp_att_limited_results[2,]),lty=2)
#legend("topleft",c("orig.","B10pc","orig.","B10pc"),cex=0.8,
#       fill=c(mycol60,mycol30,NA,NA),col=c(NA,NA,1,1),border=c(NA,mycol30,NA,NA),density=c(1000,20,NA,NA),lty=c(NA,NA,2,3))
dev.off()

