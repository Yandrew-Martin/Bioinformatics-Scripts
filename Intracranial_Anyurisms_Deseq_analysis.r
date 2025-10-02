dat <- Table(gds5186)

#remove rows with all NA
dat <- dat[rowSums(is.na(dat))<23,]

#filter columns by factor
dat.rupt <- dat.table[,dat.factor=="ruptured intracranial aneurysm"]
dat.unrupt <- dat.table[,dat.factor=="unruptured intracranial aneurysm"]
dat.sta <- dat.table[,dat.factor=="superficial temporal artery"]
dat.ia <- dat.table[,dat.factor=="unruptured intracranial aneurysm" | dat.factor=="ruptured intracranial aneurysm"]

#find average gene expression for all genes in each factor
dat.rupt.m <- apply(dat.rupt,1,mean,na.rm=T)
dat.unrupt.m <- apply(dat.unrupt,1,mean,na.rm=T)
dat.sta.m <- apply(dat.sta,1,mean,na.rm=T)
dat.ia.m <- apply(dat.ia,1,mean,na.rm=T)

dat.cor <- cor(dat.table)

#plot correlation
par(oma=c(1,7,0,1), mar=c(9,5,2,1)) 
cx <- rev(colorpanel(25,"yellow","black","blue")) 
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10) 
image(dat.cor,main="Correlation plot Normal/Tumor data",axes=F,col=cx) 
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dat.ann,cex.axis=0.9,las=2) 
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dat.ann,cex.axis=0.9,las=2) 
image(as.matrix(leg),col=cx,axes=F) 
tmp <- round(leg,2) 
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)

#transpose and plot cluster dendrogram
dat.t <- t(dat.table)
dat.dist <- dist(dat.t,method="euclidean") 
dat.clust <- hclust(dat.dist,method="single") 
plot(dat.clust,labels=names(dat.t),cex=0.75, xlab="Distance" )

#CVM plot
dat.mean <- apply(dat.table,2,mean) 
dat.sd <- sqrt(apply(dat.table,2,var)) 
dat.cv <- dat.sd/dat.mean 
plot(dat.mean,dat.cv,main="Intracranial Aneurysms dataset CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n") 
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21) 
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

#Average correlation plot
dat.avg <- apply(dat.cor,1,mean) 
par(oma=c(3,0.1,0.1,0.1)) 
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Tumor/Normal samples",axes=F) 
points(dat.avg,bg="red",col=1,pch=21,cex=1.25) 
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat.table)[[2]],las=2,cex.lab=0.4,cex.axis=0.6) 
axis(2) 
abline(v=seq(0.5,62.5,1),col="grey")

#remove suspected outliers from all three plots
dat.out <- within(dat.table, rm("GSM1306900", "GSM1306894", "GSM1306889"))
dat.unrupt <- within(dat.unrupt, rm("GSM1306894"))
dat.rupt <- within(dat.rupt, rm("GSM1306889"))
dat.sta <- within(dat.sta, rm("GSM1306900"))
dat.ia <- within(dat.ia, rm("GSM1306894", "GSM1306889"))

#calculate rowmeans again
dat.rupt.m <- apply(dat.rupt,1,mean,na.rm=T)
dat.unrupt.m <- apply(dat.unrupt,1,mean,na.rm=T)
dat.sta.m <- apply(dat.sta,1,mean,na.rm=T)
dat.ia.m <- apply(dat.ia,1,mean,na.rm=T)

#student’s t-test on all genes between control and IA groups and unruptured and ruptured groups
pv.sta.ia <- apply(dat.out,1,t.test.all.genes,s1=colnames(dat.sta),s2=colnames(dat.ia))
pv.unrupt.rupt <- apply(dat.out,1,t.test.all.genes,s1=colnames(dat.unrupt),s2=colnames(dat.rupt))

#p-value distribution histogram
hist(pv.sta.ia,col="lightblue",xlab="p-values",main="P-value dist’n between\nControl and IA groups",cex.main=0.9)

hist(pv.unrupt.rupt,col="lightblue",xlab="p-values",main="P-value dist’n between\nUnruptured and Ruptured groups",cex.main=0.9)

#bonferoni thresholds
bonf <- 0.5/length(pv.sta.ia)

#fold change between control and IA; UIA vs RIA
fold.sta.ia <- dat.sta.m - dat.ia.m
fold.unrupt.rupt <- dat.unrupt.m - dat.rupt.m

#filter significant fold and p-values (bonferoni correction)
fold.sta.ia.sig <- names(fold.sta.ia[abs(fold.sta.ia)>1])
pv.sta.ia.sig <- names(pv.sta.ia[pv.sta.ia<bonf])
fold.unrupt.rupt.sig <- names(fold.unrupt.rupt[abs(fold.unrupt.rupt)>1])

#use p<0.05: not as many significant p-values as the other comparison
pv.unrupt.rupt.sig <- names(pv.unrupt.rupt[pv.unrupt.rupt<0.05])

#volcano plot of each comparison
plot(range(pv.trans.sta.ia),range(fold.sta.ia),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\nControl and IA group differences') 
points(pv.trans.sta.ia,fold.sta.ia,col='black',pch=21,bg=1) 
points(pv.trans.sta.ia[(pv.trans.sta.ia> -log10(bonf.sta.ia)&fold.sta.ia>1)],fold.sta.ia[(pv.trans.sta.ia> -log10(bonf.sta.ia)&fold.sta.ia>1)],col=1,bg=2,pch=21) 
points(pv.trans.sta.ia[(pv.trans.sta.ia> -log10(bonf.sta.ia)&fold.sta.ia<(-1))],fold.sta.ia[(pv.trans.sta.ia> -log10(bonf.sta.ia)&fold.sta.ia<(-1))],col=1,bg=3,pch=21) 
abline(v= -log10(bonf.sta.ia)) 
abline(h= (-1)) 
abline(h=1)

#UIA vs RIA
plot(range(pv.trans.unrupt.rupt),range(fold.unrupt.rupt),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\nUnruptured and Ruptured group differences') 
points(pv.trans.unrupt.rupt,fold.unrupt.rupt,col='black',pch=21,bg=1) 
points(pv.trans.unrupt.rupt[(pv.trans.unrupt.rupt > -log10(0.05) & fold.unrupt.rupt > 1)],fold.unrupt.rupt[(pv.trans.unrupt.rupt > -log10(0.05)&fold.unrupt.rupt>1)],col=1,bg=2,pch=21)  
points(pv.trans.unrupt.rupt[(pv.trans.unrupt.rupt > log10(.05) & fold.unrupt.rupt < (-1))],fold.unrupt.rupt[(pv.trans.unrupt.rupt > -log10(.05) & fold.unrupt.rupt < (-1))],col=1,bg=3,pch=21) 
abline(v= -log10(0.05)) 
abline(h= (-1)) 
abline(h=1)

 #subset data tables
dat.sig.sta.ia <- dat.out[(pv.sta.ia>bonf&abs(fold.sta.ia)>1),]
dat.sig.unrupt.rupt <- dat.ia[(pv.unrupt.rupt>0.05&abs(fold.unrupt.rupt)>1),]

#heatmap of significantly expressed genes
hm.rg <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")
heatmap(as.matrix(dat.sig.sta.ia),col=hm.rg,xlab="Samples",ylab="Genes",main="Heatmap of significantly regulated genes\nSTA vs IA data
",cexCol = 0.5)

heatmap(as.matrix(dat.sig.unrupt.rupt),col=hm.rg,xlab="Samples",ylab="Genes",main="Heatmap of significantly regulated genes\nUIA vs RIA data", cexCol = 0.5)

dat.pca.sta.ia <- prcomp(t(dat.sig.sta.ia),cor=F)
dat.pca.unrupt.rupt <- prcomp(t(dat.sig.unrupt.rupt),cor=F)
dat.loadings <- dat.pca.sta.ia$x[,1:2]
par(xpd=NA) 
plot(range(dat.loadings[,1]), range(dat.loadings[,2]), type="n", xlab='p1',ylab='p2',main='PCA plot of Sotiriou Data\np2 vs. p1') 
points(dat.loadings[,1][colnames(dat.sta)], dat.loadings[,2][colnames(dat.sta)], col=1,bg='red',pch=21,cex=1.5) 
points(dat.loadings[,1][colnames(dat.ia)], dat.loadings[,2][colnames(dat.ia)],col=1,bg='blue',pch=21,cex=1.5) 
legend("bottomright", legend=c("STA","IA"), pch=21, pt.bg=c("red","blue"), title="site")
 
dat.loadings2 <- dat.pca.unrupt.rupt$x[,1:2]
par(xpd=NA) 
plot(range(dat.loadings2[,1]), range(dat.loadings2[,2]), type="n", xlab='p1',ylab='p2',main='PCA plot of Sotiriou Data\np2 vs. p1') 
points(dat.loadings2[,1][colnames(dat.unrupt)], dat.loadings2[,2][colnames(dat.unrupt)], col=1,bg='red',pch=21,cex=1.5) 
points(dat.loadings2[,1][colnames(dat.rupt)], dat.loadings2[,2][colnames(dat.rupt)],col=1,bg='blue',pch=21,cex=1.5) 
legend("bottomright", legend=c("UIA","RIA"), pch=21, pt.bg=c("red","blue"), title="site")

#Scree plots
dat.pca.var1 <- round(dat.pca.sta.ia$sdev^2 / sum(dat.pca.sta.ia$sdev^2)*100,2) 
plot(c(1:length(dat.pca.var1)),dat.pca.var1,type="b",xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1.5) 
title("Scree plot showing % variability explained by each eigenvalue\nSTA/IA dataset")

dat.pca.var2 <- round(dat.pca.unrupt.rupt$sdev^2 / sum(dat.pca.unrupt.rupt$sdev^2)*100,2) 
plot(c(1:length(dat.pca.var2)),dat.pca.var2,type="b",xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1.5) 
title("Scree plot showing % variability explained by each eigenvalue\nUIA/RIA dataset")

#LDA prediction and plot
training <- cbind(dat.pca.sta.ia$x[,1:3], as.numeric(dat.out.factor))
dat.train <- data.frame(training, t(dat.out))
dat.lda <- lda(tissue ~ PC1 + PC2 + PC3, data = dat.train)
dat.pred <- predict(dat.lda, dat.train)

training2 <- cbind(dat.pca.unrupt.rupt$x[,1:3], as.numeric(dat.ia.factor))
colnames(training2)[4] <- "tissue"
dat.train2 <- data.frame(training2, t(dat.ia))
dat.lda2 <- lda(tissue ~ PC1 + PC2 + PC3, data = dat.train2)
dat.pred2 <- predict(dat.lda2, dat.train2)

plot(dat.pred$x,col=dat.out.factor,pch=as.numeric(dat.pred$class), ylab="Discriminant function",xlab="Score",main="Discriminant function for STA/IA dataset") 
legend("bottomright", legend=c("predRIA","predUIA","predSTA","RIA","UIA","STA"), pch=c(1,3,2,21,21,21), pt.bg=c(1,3,2))

plot(dat.pred2$x,col=dat.ia.factor,pch=as.numeric(dat.pred2$class), ylab="Discriminant function",xlab="Score",main="Discriminant function for UIA/RIA dataset")
legend("bottomright", legend=c("predRIA","predUIA","RIA","UIA"), pch=c(1,2,21,21), pt.bg=c(1,2))

#find top 5 upregulated and downregulated genes
fold <- fold.sta.ia[rownames(dat.sig.sta.ia)]
rank <- order(fold)
rank <- fold[rank]
fold2 <- fold.unrupt.rupt[rownames(dat.sig.unrupt.rupt)]
rank2 <- order(fold2)
rank2 <- fold2[rank2]
discrim1 <- cbind(2^(rank[c(1:5,2588:2593)]), dat[names(rank[c(1:5,2588:2593)]),2])
discrim2 <- cbind(2^(rank2[c(1:5,594:598)]), dat[names(rank2[c(1:5,594:598)]),2])