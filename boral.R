#boral

library(boral)
library(Matrix)

##########
## Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM like in equation 3 of main text, i.e. two latent variables, (random) site effect, and a quadratic effect of soil. 
##########   

## Covariates need to be stored as a matrix for boral. no need for intercept
covX <- hmscXe[,-1]

## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#start 1:26pm, end 3:08pm for lo with 4 latent variables
#start 11:20pm, end 1:30am, for hi with 8 latent
#start 11:30am end 1:47pm, to hi with 8 or more occurrances, 4 latent variables, fit.lvmhiocc8
fit.lvmhiocc6 <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", row.eff = "random", save.model = TRUE, calc.ics = T, hypparams = c(20,20,20,20))#calc.ics = F, use ics (information criteria at your own risk)
save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis3.Rdata")  #alternate between 1 and 2

summary(fit.lvmhiocc8) # To look at estimated parameter values
fit.lvmhiocc8$hpdintervals # 95% credible intervals for model parameters.


#check information criteria
fit.lvm$ics
fit.lvmhi8$ics


#check convergence
plot(get.mcmcsamples(fit.lvmhiocc8)[,5])

mcmchi<-get.mcmcsamples(fit.lvm)
head(mcmchi)
dim(mcmchi)

gew.pvals <- 2*pnorm(abs(unlist(fit.lvmhiocc8$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, th trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

fit.lvmhiocc8$geweke.diag
#(4th species) ece270a38485359f4d1e2a33f964a983 for snowdepth did not converge, but looking at the trace plot, it seems fine

## compare model coefficients with GLMM:
cbind(fit.lvm$lv.coefs.mean,fit.lvm$X.coefs.mean)

## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lvmhiocc8)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lvm)


## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc...
#res.cors (hi with 4latent)
#res.corslo (lo with 4)
#res.corshi8 (hi with 8) - hardly any correlations significant or strong
#res.corshiocc8 (hi with 8 occurrences)

res.corshiocc8 <- get.residual.cor(fit.lvmhiocc8)
corrplot(res.corshiocc8$correlation[1:50,1:50], diag = F, type = "lower", title = "Residual correlations from LVM", mar=c(3,0.5,2,1), tl.srt = 45,method = "color")#
corrplot(res.cor90$sig.correlaton[1:50,1:50], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")

#use modified function to get 90% CI
res.corshiocc8.90<-get.residual.cor2(fit.lvmhiocc8)


### Environmental correlations
env.cors<-get.enviro.cor(fit.lvm)
corrplot(env.cors$cor, type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45)
corrplot(env.cors$sig.cor, type="lower", diag=F, title="Environmental correlations", mar=c(3,0.5,2,1), tl.srt=45)

#Finally, one approach to quantify how much of the species co-occurrence is explained by covariates (that is, how well the predictor variables describe the species assemblage) is through differences in the trace of the estimated residual covariance matrix induced by the latent variables (Warton et al. 2015). From the above code, this can be obtained as rescors$trace. For the spider data set, when we compared a pure latent variable model (similar to the equation 1 but without site effects) to the correlated response model, the trace decreased from 17892 to 10792. This implies that environmental covariates accounted for approximately 40% of the covariation between species.
(107.92-178.92)/178.92

#Need to fit a model with only latent variables
res.corshi8$trace


#Chord diagrams
str(res.cors)
dim(res.cors$sig.correlaton)
res.corslo$sig.correlaton[1:4,1:4]
(length(which(res.corslo$sig.correlaton!=0))-dim(res.corslo$sig.correlaton)[1])/2
colMat <- matrix(NA, nrow = nrow(res.cors$sig.correlaton), ncol = ncol(res.cors$sig.correlaton))
colMat[which(res.cors$sig.correlaton > 0.63, arr.ind = TRUE)] <- "red"
colMat[which(res.cors$sig.correlaton < -0.63, arr.ind = TRUE)] <- "blue"
chordDiagram(res.cors$sig.correlaton, symmetric = TRUE,annotationTrack = c("name", "grid"),grid.col = "grey",col=colMat)

colMat <- matrix(NA, nrow = nrow(res.corshiocc8$correlation), ncol = ncol(res.corshiocc8$correlation),dimnames = list(rownames(res.corshiocc8$correlation),rownames(res.corshiocc8$correlation)))
colMat[which(res.corshiocc8$correlation > 0.5, arr.ind = TRUE)] <- 1
colMat[which(res.corshiocc8$correlation < -0.5, arr.ind = TRUE)] <- -1
colMat[which(res.corshiocc8$correlation<.5&res.corshiocc8$correlation>(-.5))]<-0

length(which(res.cor90$correlation>.7&res.cor90$correlation<1,arr.ind = T))
length(which(res.corshi8$correlation>.7&res.corshi8$correlation<1,arr.ind = T))
hist(res.cors$correlation)
hist(res.corshi8$correlation)


#Plotting with igraph
#creating sparse matrix
colMat<-res.corshiocc8.90$sig.correlaton
colMat[which(res.corshiocc8.90$sig.correlaton>0)]<-1
colMat[which(res.corshiocc8.90$sig.correlaton<0)]<- -1

colMat<-res.corslo$sig.correlaton
colMat[which(res.corslo$sig.correlaton>0)]<-1
colMat[which(res.corslo$sig.correlaton<0)]<- -1

colMat<-res.cors$sig.correlaton
colMat[which(res.cors$sig.correlaton>0)]<-1
colMat[which(res.cors$sig.correlaton<0)]<- -1

colMat<-res.corshiocc8$sig.correlaton
colMat[which(res.corshiocc8$sig.correlaton>0)]<-1
colMat[which(res.corshiocc8$sig.correlaton<0)]<- -1


colMat[1:20,1:20]

#colMat2<-as(colMat, "dgCMatrix")
graph1<-graph_from_adjacency_matrix(colMat, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelist<-data.frame(as_edgelist(graph1),weight=E(graph1)$weight) #just the edges

#graph2<-simplify(graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE))
graph2<-graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE)
graph2

#(sum(colMat,na.rm=T)-322)/2 #(134 interactions)
#(24*24-24)/2 #every pairwise interaction (276)

graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist$weight==1,"red","blue"),vertex.color=colorgraph$color)

10 or higher 580/134
8 or higher 630/158



##Exploring nematode correlations in hi
#use colorgraph from below
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
colorgraph2
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%colorgraph2$oldotu]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%colorgraph2$oldotu]
temp<-c("a64f0e90f25d607e4c34fffa870f19e4","71f3f9709801ee03a20e33c46a8b797d","a7a7353d6231f298e75b69f4161926a1","ee45b87e6ac04257510a6d99644760cd","eafa1d1d561a47042f099e619cd25693","e4ae90d4c05ef00793fc6b93fb6a9af7","f7b389ab203fda4cc3cf1fffe23004a6","9bd3df7f714986d788a2cce65be56d10","7a97813269020725286a63989363ceef","d3898987a65cc50f9226529374935e74","854124b2d8e0badeb138eb9e0e1808c2")#lo 10
temp<-c("f4e7445ced25e7028c593551658c091a","f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8","479660228dd139e20ad8ff42698a123e","65cb8aa63bed1f44c9a576ae2e94e0ab")#hi 8
temp<-c("f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8")#hi 10
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$oldotu%in%temp),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"),]

#tax_table(datEukN3)[rownames(tax_table(datEukN3))=="6f914ead2160e51670d3dc70c25e107b"]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph2<-graph.edgelist(as.matrix(myedgelist2[,1:2]),directed=FALSE)
graph2
graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph4<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"red","blue"),vertex.color=colorgraph4$color)

#using occ8 with correlation >.5 or <-.5, I get 10 mesofauna taxa
#using occ10 significnat correlations, I get 5 mesofauna and 4 heterotrophic eukaryotes. 67/51
#using occ8 significnat correlations, I get 5 mesofauna and 7 heterotrophic eukaryotes. 53/45
#lo 10, significant, 1 mesofuana, 10 heterotrophic euks 128/81
#using occ8 with 90%CI, 7 mesofauna (didn't look for euks)




#90% CI
get.residual.cor2<-function (object, est = "median", prob = .9) 
{
  if (is.null(object$jags.model)) 
    stop("MCMC samples not found")
  fit.mcmc <- get.mcmcsamples(object)
  y <- object$y
  X <- object$X
  num.lv <- object$num.lv
  if (length(grep("lvs", colnames(fit.mcmc))) == 0) 
    stop("Cannot find MCMC samples corresponding to latent variables")
  n <- nrow(y)
  p <- ncol(y)
  sig_rescor_mat <- rescor_mat <- rescov_mat <- matrix(0, nrow = p, 
                                                       ncol = p)
  sig_respres_mat <- respres_mat <- matrix(0, nrow = p, ncol = p)
  if (is.null(colnames(y))) 
    colnames(y) <- 1:ncol(y)
  rownames(rescor_mat) <- colnames(rescor_mat) <- rownames(sig_rescor_mat) <- colnames(sig_rescor_mat) <- colnames(y)
  rownames(rescov_mat) <- colnames(rescov_mat) <- colnames(y)
  rownames(respres_mat) <- colnames(respres_mat) <- rownames(sig_respres_mat) <- colnames(sig_respres_mat) <- colnames(y)
  all_rescor_mat <- all.rescov_mat <- all.respres_mat <- array(0, 
                                                               dim = c(nrow(fit.mcmc), p, p))
  all_trace_rescor <- numeric(nrow(fit.mcmc))
  for (k0 in 1:nrow(fit.mcmc)) {
    lv.coefs <- matrix(fit.mcmc[k0, grep("lv.coefs", colnames(fit.mcmc))], 
                       nrow = p)
    lambdalambdaT <- tcrossprod(as.matrix(lv.coefs[, 2:(num.lv + 
                                                          1)]))
    all.rescov_mat[k0, , ] <- (lambdalambdaT)
    all_trace_rescor[k0] <- sum(diag(lambdalambdaT))
    all_rescor_mat[k0, , ] <- cov2cor(all.rescov_mat[k0, 
                                                     , ])
    all.respres_mat[k0, , ] <- ginv(all_rescor_mat[k0, , 
                                                   ])
  }
  for (j in 1:p) {
    for (j2 in 1:p) {
      if (est == "median") {
        rescov_mat[j, j2] <- median(all.rescov_mat[, 
                                                   j, j2])
        rescor_mat[j, j2] <- median(all_rescor_mat[, 
                                                   j, j2])
        respres_mat[j, j2] <- median(all.respres_mat[, 
                                                     j, j2])
      }
      if (est == "mean") {
        rescov_mat[j, j2] <- mean(all.rescov_mat[, j, 
                                                 j2])
        rescor_mat[j, j2] <- mean(all_rescor_mat[, j, 
                                                 j2])
        respres_mat[j, j2] <- mean(all.respres_mat[, 
                                                   j, j2])
      }
      sig_rescor_mat[j, j2] <- rescor_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all_rescor_mat[, 
                                                         j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_rescor_mat[j, j2] <- 0
      sig_respres_mat[j, j2] <- respres_mat[j, j2]
      get.hpd.cors <- HPDinterval(as.mcmc(all.respres_mat[, 
                                                          j, j2]), prob = .9)
      if (0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
        sig_respres_mat[j, j2] <- 0
    }
  }
  if (est == "median") 
    final_trace <- median(all_trace_rescor)
  if (est == "mean") 
    final_trace <- mean(all_trace_rescor)
  return(list(correlation = rescor_mat, sig.correlaton = sig_rescor_mat, 
              covariance = rescov_mat, precision = respres_mat, sig.precision = sig_respres_mat, 
              trace = final_trace))
}
