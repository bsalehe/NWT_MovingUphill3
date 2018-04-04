#boral

library(boral)
library(Matrix)

##########
## Fitting a Latent Variable Model using the boral package ##
## Use boral to fit a LVM with 4 latent variables, a random site effect, and the effect of fixed environmental variables
##########   

## Covariates need to be stored as a matrix for boral. no need for intercept
covX <- hmscXe[,-1]

## fit the LVM using boral and look at results.  Need version 0.7 or later, available on CRAN.
## NOTE: boral uses JAGS so you need to first install JAGS-3.0.0.exe or higher from http://sourceforge.net/projects/mcmc-jags/files/JAGS/

#Using the default mcmc parameters, the models take about 2 hrs to fit
#fit.lvmhiocc6 
#load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis3b.Rdata")#to load fit.lvmhiocc6
#fit.hilv4occ8exp5
#fit.lolv4occ8exp5
#fit.hilv4occ8exp4 (and b)
#fit.melv4occ8exp4

fit.lolv4occ8exp4 <- boral(y = hmscYe, X = covX, num.lv = 4, family = "negative.binomial", row.eff = "random", save.model = TRUE, calc.ics = T,mcmc.control = list(n.burnin = 1000, n.iteration = 4000, n.thin = 3, seed = 123), hypparams = c(20,20,20,20))#calc.ics = F, use ics (information criteria at your own risk)

#mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123)

res.corslolv4occ8exp4 <- get.residual.cor(fit.lolv4occ8exp4) #this will not fit with a cuttoff of occ6, R studio crashes on my laptop
save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis3.Rdata") 

summary(fit.lvmhiocc8) # To look at estimated parameter values
fit.lvmhiocc8$hpdintervals # 95% credible intervals for model parameters.


#check information criteria
fit.lvm$ics
fit.lvmhi8$ics


#check convergence
plot(get.mcmcsamples(fit.hilv4occ8exp4)[,5])

mcmchi<-get.mcmcsamples(fit.hilv4occ8exp4)
head(mcmchi)
dim(mcmchi)

gew.pvals <- 2*pnorm(abs(unlist(fit.hilv4occ8exp4$geweke.diag[[1]])), lower.tail = FALSE)
length(gew.pvals)
gew.pvals[which(gew.pvals<.05)] #technically these did not converge, however, the trace plots look fine to me
p.adjust(gew.pvals, method = "holm")

fit.hilv4occ8exp4$geweke.diag
#(4th species) ece270a38485359f4d1e2a33f964a983 for snowdepth did not converge, but looking at the trace plot, it seems fine

## compare model coefficients with GLMM:
cbind(fit.lvm$lv.coefs.mean,fit.lvm$X.coefs.mean)

## Dunn-Smyth residual plots  to check model assumption, outliers etc...
plot(fit.lolv4occ8exp5)

## Residual ordination plot of the sites (please see Figure 2b in the main text for relevant explanation)
## Please note that boral currently does not automatically produce biplots, although species loadings can be obtained from fit.lvm$lv.coefs.median[,2:3]
lvsplot(fit.lolv4occ8exp5)


## Use corrplot package to plot residual correlations between species, e.g. possibly due to species interaction etc...
#res.cors (hi with 4latent)
#res.corslo (lo with 4 latent)
#res.corshi8 (hi with 8 latent) - hardly any correlations significant or strong
#res.corshiocc8 (hi with 8 occurrences)
#res.corshiocc8.90 (hi with 8 occurrences and 90CI)
#res.corshiocc6 (hi 6 occ)
#res.corshilv4occ8exp5 (hi with 5 environmetnal variables)
#res.corslolv4occ8exp5

res.corshiocc8 <- get.residual.cor(fit.lvmhiocc8)
res.corshiocc6 <- get.residual.cor(fit.lvmhiocc6)
#save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill3_Workspace_Analysis4.Rdata")  #alternate between 1 and 2

corrplot(res.corshiocc8$sig.correlaton[1:100,1:100], diag = F, type = "lower", title = "Residual correlations from LVM", mar=c(3,0.5,2,1), tl.srt = 45,method = "color")#
corrplot(res.corshilv4occ8exp5$sig.correlaton[1:100,1:100], type="lower", diag=F, title="Residual correlations", mar=c(3,0.5,2,1), tl.srt=45,method = "color")

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
colMat<-res.corshilv4occ8exp5$sig.correlaton
colMat[which(res.corshilv4occ8exp5$sig.correlaton>0)]<-1
colMat[which(res.corshilv4occ8exp5$sig.correlaton<0)]<- -1

colMat<-res.corslolv4occ8exp5$sig.correlaton
colMat[which(res.corslolv4occ8exp5$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp5$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(res.corshilv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corshilv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corsmelv4occ8exp4$sig.correlaton
colMat[which(res.corsmelv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corsmelv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(res.corslolv4occ8exp4$sig.correlaton>0)]<-1
colMat[which(res.corslolv4occ8exp4$sig.correlaton<0)]<- -1

colMat<-res.corshilv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

colMat<-res.corslolv4occ8exp4$sig.correlaton
colMat[which(colMat>.6)]<-1
colMat[which(colMat<(-.6))]<- -1
colMat[which(colMat<.6&colMat>(-.6))]<-0

#hist(colMat)
#dim(colMat)

#colMat2<-as(colMat, "dgCMatrix")
graph1<-graph_from_adjacency_matrix(colMat, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelist<-data.frame(as_edgelist(graph1),weight=E(graph1)$weight) #just the edges

length(which(myedgelist$weight==1))
length(which(myedgelist$weight==-1))
length(which(myedgelist$weight==1))/(length(which(myedgelist$weight==1))+length(which(myedgelist$weight==-1)))

#graph2<-simplify(graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE))
graph2<-graph.edgelist(as.matrix(myedgelist[,1:2]),directed=FALSE)
graph2

#(sum(colMat,na.rm=T)-322)/2 #(134 interactions)
#(24*24-24)/2 #every pairwise interaction (276)

graph2$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph2))))
colorgraph<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,layout=l3

unique(tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%colorgraph$oldotu])
temp<-myedgelist[which(myedgelist$X1%in%c("8f0fb68673a8b3c26f3fd210e29b7916")|myedgelist$X2%in%c("8f0fb68673a8b3c26f3fd210e29b7916")),]
8f0fb68673a8b3c26f3fd210e29b7916 #clostridium bowmanii, byproducts of acetic and butyric acid
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

temp<-myedgelist[which(myedgelist$X1%in%c("6472eb8b1e09f892aca2f23182962903")|myedgelist$X2%in%c("6472eb8b1e09f892aca2f23182962903")),]
8f0fb68673a8b3c26f3fd210e29b7916  #bradyrhizobium fixes N
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%c(as.character(temp2$X1),as.character(temp2$X2))] #one chlorophyceae - possible syntrophic relationship
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]


"631c2b290edb6b71d9e6b5577a0613ae"# syntrophobacteriaceae
"8a8895d160b2fda424ad5446e9bf8f4e"#Desulfosporosinus sulfate reducing
temp<-myedgelist[which(myedgelist$X1%in%c("8a8895d160b2fda424ad5446e9bf8f4e")|myedgelist$X2%in%c("8a8895d160b2fda424ad5446e9bf8f4e")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

"b8f64099026ffffde7cd2c6aa4642240" #nitrosomonadaceae, ammonia oxidizers
temp<-myedgelist[which(myedgelist$X1%in%c("b8f64099026ffffde7cd2c6aa4642240")|myedgelist$X2%in%c("b8f64099026ffffde7cd2c6aa4642240")),]
temp2<-temp[which(temp$weight==1),]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]
tax_table(datITSS3)[rownames(tax_table(datITSS3))%in%c(as.character(temp2$X1),as.character(temp2$X2))]

myedgelist[myedgelist$X1%in%temp,]
myedgelist[myedgelist$X2%in%temp,]
tax_table(datBacS3)[rownames(tax_table(datBacS3))%in%c(temp)]
colorgraph[which(colorgraph$oldotu=="854124b2d8e0badeb138eb9e0e1808c2"),]

#Creating a subgraph
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Bacteria"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph3$oldotu]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-subgraph.edges(graph2, eids=which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu), delete.vertices = F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraph$color)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)

#investigating nematode relationships
#lo
temp<-"b85db42310af5ddb08354eef2427cc8e" #omnivore nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="b85db42310af5ddb08354eef2427cc8e"|myedgelist$X2=="b85db42310af5ddb08354eef2427cc8e"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#5 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#2 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#6
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#1
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#6
tax_table(datEukN)[rownames(tax_table(datEukN))%in%temp2$X2]#0

#hi
6f914ead2160e51670d3dc70c25e107b __Aphelenchida fungal feeder
e485c9c4bdda247bc813c36275d99dce __Plectidae bacterial feeder
temp<-"6f914ead2160e51670d3dc70c25e107b" #ff nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="6f914ead2160e51670d3dc70c25e107b"|myedgelist$X2=="6f914ead2160e51670d3dc70c25e107b"),]
#only positive rels
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#2 bacteria
temp<-"e485c9c4bdda247bc813c36275d99dce" #bf nematode in lo
temp1<-myedgelist[which(myedgelist$X1=="e485c9c4bdda247bc813c36275d99dce"|myedgelist$X2=="e485c9c4bdda247bc813c36275d99dce"),]
temp2<-temp1[temp1$weight==-1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#0 euks
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]#0
temp2<-temp1[temp1$weight==1,]
tax_table(datBacS)[rownames(tax_table(datBacS))%in%temp2$X2]#3 bacteria
tax_table(datEukS)[rownames(tax_table(datEukS))%in%temp2$X2]#
tax_table(datITSS)[rownames(tax_table(datITSS))%in%temp2$X2]
# 2 plants 





#trying to keep the layout the same and plot only partial graph
l3<-layout_in_circle(graph2)#layout_with_fr(graph2) #I don't want to run this accidentally
rownames(l3) <- V(graph2)$name
l3b<-layout.norm(l3,xmin=-1,xmax=1,ymin=-1,ymax=1)
#plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph3,layout=l3b,rescale=F,xlim=c(-1,1),ylim=c(-1,1))




##Exploring nematode correlations in hi
#use colorgraph from below
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$labels=="Eukaryota"),]
colorgraph2<-colorgraph[which(colorgraph$labels=="Mesofauna"),]
colorgraph2
tax_table(datEukN3)[rownames(tax_table(datEukN3))%in%colorgraph2$oldotu]
tax_table(datEukN)[rownames(tax_table(datEukN))%in%colorgraph2$oldotu]
tax_table(datEukS3)[rownames(tax_table(datEukS3))%in%colorgraph2$oldotu]
temp<-c("a64f0e90f25d607e4c34fffa870f19e4","71f3f9709801ee03a20e33c46a8b797d","a7a7353d6231f298e75b69f4161926a1","ee45b87e6ac04257510a6d99644760cd","eafa1d1d561a47042f099e619cd25693","e4ae90d4c05ef00793fc6b93fb6a9af7","f7b389ab203fda4cc3cf1fffe23004a6","9bd3df7f714986d788a2cce65be56d10","7a97813269020725286a63989363ceef","d3898987a65cc50f9226529374935e74","854124b2d8e0badeb138eb9e0e1808c2")#lo 10
temp<-c("f4e7445ced25e7028c593551658c091a","f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8","479660228dd139e20ad8ff42698a123e","65cb8aa63bed1f44c9a576ae2e94e0ab")#hi 8
temp<-c("f8d2629f619aea2ec1a4a2571c91f1af","01a1a6cfe5266fe585bed622ab7f64c4","124f8bbd2e10abcb816e84a974280321","3094f58213b725653eb04f35aea790f8")#hi 10
colorgraph3<-colorgraph[which(colorgraph$labels=="Mesofauna"|colorgraph$oldotu%in%temp),]
colorgraph3<-colorgraph[which(colorgraph$labels=="Fungi"),]

#tax_table(datEukN3)[rownames(tax_table(datEukN3))=="6f914ead2160e51670d3dc70c25e107b"]

myedgelist2<-myedgelist[which(myedgelist[,"X1"]%in%colorgraph3$oldotu|myedgelist[,"X2"]%in%colorgraph3$oldotu),]
graph3<-graph.edgelist(as.matrix(myedgelist2[,1:2]),directed=FALSE)
graph3
graph3$layout <- layout_in_circle
verticesgraph<-data.frame(oldotu=rownames(as.matrix(V(graph3))))
colorgraph4<-merge(verticesgraph,labelsall,"oldotu",all.y=F,all.x=F,sort=F)
plot(graph3,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelist2$weight==1,"red","blue"),vertex.color=colorgraph$color,layout=l3b)#,rescale=F,xlim=c(-1,1),ylim=c(-1,1)



#using hi occ8 with correlation >.5 or <-.5, I get 10 mesofauna taxa
#using hi occ10 significnat correlations, I get 5 mesofauna and 4 heterotrophic eukaryotes. 67/51
#using hi occ8 significnat correlations, I get 5 mesofauna and 7 heterotrophic eukaryotes. 53/45
#lo 10, significant, 1 mesofuana, 10 heterotrophic euks 128/81
#using hi occ8 with 90%CI, 7 mesofauna (didn't look for euks)









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
