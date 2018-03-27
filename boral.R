#boral

library(boral)
library(Matrix)

##########
## (2) Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM like in equation 3 of main text, i.e. two latent variables, (random) site effect, and a quadratic effect of soil. 
## Also note that this can fit models to more data (e.g. try full dataset) but takes longer to fit it - for full spdier it will still be a few min, for full aravo dataset over half an hour.
##########   

## Covariates need to be stored as a matrix for boral. no need for intercept
covX <- hmscXe[,-1]

## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#start 1:26pm end 3:08pm
fit.lvmlo <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", row.eff = "random", save.model = TRUE, calc.ics = T, hypparams = c(20,20,20,20))#calc.ics = F, use ics (information criteria at your own risk)
summary(fit.lvm) # To look at estimated parameter values
fit.lvm$hpdintervals # 95% credible intervals for model parameters.


#check convergence
plot(get.mcmcsamples(fit.lvm)[,4])

mcmchi<-get.mcmcsamples(fit.lvm)
head(mcmchi)
dim(mcmchi)

gew.pvals <- 2*pnorm(abs(unlist(fit.lvm$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, th trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

fit.lvm$geweke.diag
#(4th species) ece270a38485359f4d1e2a33f964a983 for snowdepth did not converge, but looking at the trace plot, it seems fine

## compare model coefficients with GLMM:
cbind(fit.lvm$lv.coefs.mean,fit.lvm$X.coefs.mean)

## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lvmlo)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lvm)

## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc...
res.corslo <- get.residual.cor(fit.lvmlo)
corrplot(res.cors$correlation, diag = F, type = "lower", title = "Residual correlations from LVM", mar=c(3,0.5,2,1), tl.srt = 45)#method = "color",
corrplot(res.cors$sig.correlaton, type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45)

env.cors<-get.enviro.cor(fit.lvm)
corrplot(env.cors$cor, type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45)
corrplot(env.cors$sig.cor, type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45)

#Finally, one approach to quantify how much of the species co-occurrence is explained by covariates (that is, how well the predictor variables describe the species assemblage) is through differences in the trace of the estimated residual covariance matrix induced by the latent variables (Warton et al. 2015). From the above code, this can be obtained as rescors$trace. For the spider data set, when we compared a pure latent variable model (similar to the equation 1 but without site effects) to the correlated response model, the trace decreased from 17892 to 10792. This implies that environmental covariates accounted for approximately 40% of the covariation between species.
(107.92-178.92)/178.92

#Need to fit a model with only latent variables
res.cors$trace

str(res.cors)
dim(res.corslo$sig.correlaton)
res.corslo$sig.correlaton[1:4,1:4]
length(which(res.corslo$sig.correlaton!=0))-dim(res.corslo$sig.correlaton)[1]
colMat <- matrix(NA, nrow = nrow(res.cors$sig.correlaton), ncol = ncol(res.cors$sig.correlaton))
colMat[which(res.cors$sig.correlaton > 0.63, arr.ind = TRUE)] <- "red"
colMat[which(res.cors$sig.correlaton < -0.63, arr.ind = TRUE)] <- "blue"
chordDiagram(res.cors$sig.correlaton, symmetric = TRUE,annotationTrack = c("name", "grid"),grid.col = "grey",col=colMat)

colMat <- matrix(NA, nrow = nrow(res.cors$correlation), ncol = ncol(res.cors$correlation),dimnames = list(rownames(res.cors$correlation),rownames(res.cors$correlation)))
colMat[which(res.cors$correlation > 0.8, arr.ind = TRUE)] <- 1
colMat[which(res.cors$correlation < -0.8, arr.ind = TRUE)] <- -1
colMat[which(res.cors$correlation<.8&res.cors$correlation>(-.8))]<-0

which(res.cors$correlation>.8&res.cors$correlation<1,arr.ind = T)

#creating sparse matrix
colMat<-res.corslo$sig.correlaton
colMat[which(res.corslo$sig.correlaton>0)]<-1
colMat[which(res.corslo$sig.correlaton<0)]<- 0

colMat[1:20,1:20]

#colMat2<-as(colMat, "dgCMatrix")
graph1<-graph_from_adjacency_matrix(colMat, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
E(graph1)$weight
myedgelist<-data.frame(as_edgelist(graph1),weight=E(graph1)$weight) #just the edges

#graph2<-simplify(graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE))
graph2<-graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE)
E(graph2)

#(sum(colMat,na.rm=T)-322)/2 #(134 interactions)
#(24*24-24)/2 #every pairwise interaction (276)

graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)

plot(graph2,vertex.size=4,edge.curved=T,vertex.label=NA,edge.color=ifelse(myedgelist$weight==1,"red","blue"),vertex.color=colorgraph$color)

tax_table(datEukN3)[rownames(tax_table(datEukN3))=="6f914ead2160e51670d3dc70c25e107b"]
