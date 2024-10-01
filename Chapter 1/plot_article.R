###libraries####
library(vcfR)
library(dartR)
library(adegenet)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(ggrepel)


###data import###
setwd("~/Desktop/pinctadapt/article_1")
vcf <- "chip_pinctadapt_data_cleaned_outgroup.vcf"


obj.vcf <- read.vcfR(vcf)
obj.gl <- vcfR2genlight(obj.vcf)

ID <- ((dimnames(obj.vcf@gt))[[2]])[-1]
ID_df <- as.data.frame(ID)
pop <- strsplit(as.character(ID_df$ID), split = "_")
POP <- as.data.frame(t(as.data.frame(pop)))

POP$V1[POP$V1 == "WildJuvIndo"] <- "INDO"
POP$V1[POP$V1 == "MP"] <- "NHV-S"
POP$V1[POP$V1 == "MG"] <- "NHV-P"
POP$V1[POP$V1 == "MarSud"] <- "MARS"
POP$V1[POP$V1 == "WC1"] <- "PEI"
POP$V1[POP$V1 == "WC2"] <- "MAT"
POP$V1[POP$V1 == "WC3"] <- "RIK"

pop_list <- as.factor(POP$V1)
pop(obj.gl) <- pop_list
#obj.gl <- obj.gl[pop(obj.gl)!= "MARS",]
obj.glx <- gl.compliance.check(obj.gl)
obj.glx <- gl.drop.ind(obj.glx, "NHV-S_08")
obj.gi <- gl2gi(obj.glx)

table(pop(obj.gl))



###admixture####

df_ID <- read.csv("admixture_PF_indo_marq/chip_pinctadapt_data_cleaned_PF_indo_marq.fam", sep = "", header = F)
labels <- as.data.frame(cbind(paste(df_ID$V1,df_ID$V2, sep = "_"), df_ID$V1))

labels$V2[labels$V2 == "WildJuvIndo"] <- "INDO"
labels$V2[labels$V2 == "MP"] <- "NHV-S"
labels$V2[labels$V2 == "MG"] <- "NHV-P"
labels$V2[labels$V2 == "MarSud"] <- "MARS"
labels$V2[labels$V2 == "WC1"] <- "PEI"
labels$V2[labels$V2 == "WC2"] <- "MAT"
labels$V2[labels$V2 == "WC3"] <- "RIK"

# Name the columns
names(labels)<-c("ind","pop")

populations_order <- "INDO,AUQ,NHV-S,NHV-P,UAH,SCI,MOP,AHE,MAN,TKP,ARA,KAU,RAR,ANA,KAT,TAH,MOT,TEA,TAK,RVV,MARS,MOR,RIK,GMBW,PEI,MAT"

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(populations_order,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))
rep<-as.vector(table(labels$n))
rep_pop <-0
for (i in 1:length(levels(as.factor(labels$pop)))) {
  rep_pop <- c(rep_pop,
               rep((as.vector(strsplit(populations_order,",")))[[1]][i],
                   rep[i]))
}
rep_pop_K1 <- rep(rep_pop[-1],1)
rep_pop_K2 <- rep(rep_pop[-1],2)
rep_pop_K3 <- rep(rep_pop[-1],3)
rep_pop_K4 <- rep(rep_pop[-1],4)
rep_pop_K5 <- rep(rep_pop[-1],5)

# read in the different admixture output files
directory_prefix <- "admixture_PF_indo_marq//chip_pinctadapt_data_cleaned_PF_indo_marq"
minK=1
maxK=5
tbl<-lapply(minK:maxK, function(x) read.table(paste0(directory_prefix,".",x,".Q")))

#K1
K1_df <- as.data.frame(cbind(tbl[[1]][order(labels$n),],
                             Categorie=1:dim(tbl[[1]])[1]))
K1.df.long <- tidyr:: gather(K1_df, "Cluster", "Proportion", -Categorie)
K1.df.long.pop <- as.data.frame(cbind(K1.df.long,
                                      Population = factor(rep_pop_K1,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K1 <- ggplot(K1.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K1 <- K1 + geom_bar(stat='identity',width=1)
K1 <- K1 + scale_fill_manual(values = c("red","green","blue"),name="cluster") 
K1 <- K1 + facet_grid(~K1.df.long.pop$Population, scales = "free")
K1 <- K1 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K1 <- K1 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"))

#K2
K2_df <- as.data.frame(cbind(tbl[[2]][order(labels$n),],
                             Categorie=1:dim(tbl[[2]])[1]))
K2.df.long <- tidyr:: gather(K2_df, "Cluster", "Proportion", -Categorie)
K2.df.long.pop <- as.data.frame(cbind(K2.df.long,
                                      Population = factor(rep_pop_K2,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K2 <- ggplot(K2.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K2 <- K2 + geom_bar(stat='identity',width=1)
K2 <- K2 + scale_fill_manual(values = c("red","green","blue", "violet"),name="cluster") 
K2 <- K2 + facet_grid(~K2.df.long.pop$Population, scales = "free")
K2 <- K2 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K2 <- K2 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"))

#K3
K3_df <- as.data.frame(cbind(tbl[[3]][order(labels$n),],
                             Categorie=1:dim(tbl[[3]])[1]))
K3.df.long <- tidyr:: gather(K3_df, "Cluster", "Proportion", -Categorie)
K3.df.long.pop <- as.data.frame(cbind(K3.df.long,
                                      Population = factor(rep_pop_K3,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K3 <- ggplot(K3.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K3 <- K3 + geom_bar(stat='identity',width=1)
K3 <- K3 + scale_fill_manual(values = c("red","green","blue", "violet", "light blue"),name="cluster") 
K3 <- K3 + facet_grid(~K3.df.long.pop$Population, scales = "free")
K3 <- K3 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K3 <- K3 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"))

#K4
K4_df <- as.data.frame(cbind(tbl[[4]][order(labels$n),],
                             Categorie=1:dim(tbl[[4]])[1]))
K4.df.long <- tidyr:: gather(K4_df, "Cluster", "Proportion", -Categorie)
K4.df.long.pop <- as.data.frame(cbind(K4.df.long,
                                      Population = factor(rep_pop_K4,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K4 <- ggplot(K4.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K4 <- K4 + geom_bar(stat='identity',width=1)
K4 <- K4 + scale_fill_manual(values = c("red","green","blue", "violet", "light blue"),name="cluster") 
K4 <- K4 + facet_grid(~K4.df.long.pop$Population, scales = "free")
K4 <- K4 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K4 <- K4 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"))

#K5
K5_df <- as.data.frame(cbind(tbl[[5]][order(labels$n),],
                             Categorie=1:dim(tbl[[5]])[1]))
K5.df.long <- tidyr:: gather(K5_df, "Cluster", "Proportion", -Categorie)
K5.df.long.pop <- as.data.frame(cbind(K5.df.long,
                                      Population = factor(rep_pop_K5,
                                                          levels = (as.vector(strsplit(populations_order,",")))[[1]])))


K5 <- ggplot(K5.df.long.pop, aes(x=Categorie, y=Proportion, fill=Cluster))  + ylab("Assignment") +xlab("Population")
K5 <- K5 + geom_bar(stat='identity',width=1)
K5 <- K5 + scale_fill_manual(values = c("red","green","blue", "violet", "light blue"),name="cluster") 
K5 <- K5 + facet_grid(~K5.df.long.pop$Population, scales = "free")
K5 <- K5 + theme_gray(base_size = 12) + guides(fill=guide_legend(ncol=1))
K5 <- K5 + theme(axis.text.x = element_blank(),
                 axis.ticks = element_blank(),
                 panel.spacing.x=unit(0.1, "lines"))


admixture_1 <- ggarrange(K2,K3, ncol = 1, nrow = 2, labels = c("K=2",
                                                               "K=3"))
admixture_2 <- ggarrange(K4,K5, ncol = 1, nrow = 2, labels = c("K=4",
                                                               "K=5"))
admixture_1
admixture_2

####ACP####
obj.gl_corrected <- gl.drop.ind(obj.glx, "NHV-S_08")
obj.gi_cor <- gl2gi(obj.gl_corrected)
#prepare for the pca and to it
X <- scaleGen(obj.gi_cor, NA.method ="mean")
#pca1 <- dudi.pca(X)
pca1 <- dudi.pca(X, scannf = FALSE, nf = 2)

#plot the acp and the variables
#eigencalues
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

#Avec ggplot
population <- as.character(pop(obj.gi_cor))
tuamotus <- c("AHE", "ANA", "ARA", "KAT", "KAU", "MAN", "MOT", "RAR", "TAH", "TAK", "TEA", "TKP")
marquesas <- c("NHV-S", "NHV-P", "UAH", "UAP", "AUQ")
gambier <- c("GMBW", "MOR", "MARS", "PEI", "MAT", "RIK")
society <- c("SCI", "MOP")
australe <- c("RVV")
indo <- c("INDO")

for (i in tuamotus) {
  population[population == i] <- "Tuamotus"
}
for (i in marquesas) {
  population[population == i] <- "Marquesas"
}

for (i in gambier) {
  population[population == i] <- "Gambier"
}

for (i in society) {
  population[population == i] <- "Society"
}

for (i in australe) {
  population[population == i] <- "Austral"
}

for (i in indo) {
  population[population == i] <- "Indonesia"
}

new<-cbind(pca1$li,population)
str(new)
names(new)<-c("PC1","PC2","pop")

centroids <- aggregate(cbind(PC1, PC2)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")


coly <- funky(length(unique(test$pop)))
a <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = coly)+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (2.23%)", y="PC2 (0.89%)")
plot(a)
#####


obj.gl_restricted <- gl.drop.pop(obj.gl_corrected, c("NHV-S", "NHV-P", "UAH", "INDO"))
obj.glx <- gl.compliance.check(obj.gl_restricted)             
obj.gi_restricted <- gl2gi(obj.glx)


X2 <- scaleGen(obj.gi_restricted, NA.method ="mean")
#pca1 <- dudi.pca(X)
pca2 <- dudi.pca(X2, scannf = FALSE, nf = 2)

#plot the acp and the variables
#eigencalues
barplot(pca2$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))


population_restricted <- as.character(pop(obj.gi_restricted))

for (i in tuamotus) {
  population_restricted[population_restricted == i] <- "Tuamotus"
}
for (i in marquesas) {
  population_restricted[population_restricted == i] <- "Marquesas"
}

for (i in gambier) {
  population_restricted[population_restricted == i] <- "Gambier"
}

for (i in society) {
  population_restricted[population_restricted == i] <- "Society"
}

for (i in australe) {
  population_restricted[population_restricted == i] <- "Austral"
}

for (i in indo) {
  population_restricted[population_restricted == i] <- "Indonesia"
}

new<-cbind(pca2$li,population_restricted)
str(new)
names(new)<-c("PC1","PC2","pop")

centroids <- aggregate(cbind(PC1, PC2)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")


colx <- coly[c(-3,-4)]
b <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = colx)+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (0.30%)", y="PC2 (0.27%)")
plot(b)

ACP <- ggarrange(a,b , ncol = 2, nrow = 1, labels = c("A","B"))
ACP

#### FSt ####
archipelago <- ifelse(pop(obj.glx) == "AHE", "Tuamotus", pop(obj.glx))
archipelago <- ifelse(pop(obj.glx) == "ANA", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "ARA", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "GMBW", "Gambier", archipelago)
archipelago <- ifelse(pop(obj.glx) == "KAT", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "KAU", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MAN", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "NHV-P", "Marquesas", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MOP", "Society", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MOR", "Gambier", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MOT", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "NHV-S", "Marquesas", archipelago)
archipelago <- ifelse(pop(obj.glx) == "RAR", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "RVV", "Autrales", archipelago)
archipelago <- ifelse(pop(obj.glx) == "SCI", "Society", archipelago)
archipelago <- ifelse(pop(obj.glx) == "TAH", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "TAK", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "TEA", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "TKP", "Tuamotus", archipelago)
archipelago <- ifelse(pop(obj.glx) == "UAH", "Marquesas", archipelago)
archipelago <- ifelse(pop(obj.glx) == "PEI", "Gambier", archipelago)
archipelago <- ifelse(pop(obj.glx) == "RIK", "Gambier", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MAT", "Gambier", archipelago)
archipelago <- ifelse(pop(obj.glx) == "INDO", "Indonesia", archipelago)
archipelago <- ifelse(pop(obj.glx) == "MARS", "Gambier", archipelago)

obj.gl.arch <- obj.glx
pop(obj.gl.arch) <- archipelago
summary(pop(obj.gl.arch))

obj.glx.arch <- gl.compliance.check(obj.gl.arch)

data_fst_test_arch <- gl.fst.pop(obj.glx.arch, verbose = 0)
data_fst_test_arch[["Fsts"]]
data_fst_test_arch[["Pvalues"]]

obj.gl.tuam <- obj.glx[pop(obj.gl.arch)== "Tuamotus",]
obj.gl.gamb  <- obj.glx[pop(obj.gl.arch)== "Gambier",]
obj.gl.marq  <- obj.glx[pop(obj.gl.arch)== "Marquesas",]
obj.gl.soc  <- obj.glx[pop(obj.gl.arch)== "Society",]

obj.glx.tuam <- gl.compliance.check(obj.gl.tuam)
obj.glx.gamb <- gl.compliance.check(obj.gl.gamb)
obj.glx.marq <- gl.compliance.check(obj.gl.marq)
obj.glx.soc <- gl.compliance.check(obj.gl.soc)

obj.gi.tuam <- gl2gi(obj.gl.tuam)
obj.gi.gamb  <- gl2gi(obj.gl.gamb)
obj.gi.marq  <- gl2gi(obj.gl.marq)
obj.gi.soc <- gl2gi(obj.gl.soc)

X.tuam <- scaleGen(obj.gi.tuam, NA.method ="mean")
X.gamb <- scaleGen(obj.gi.gamb, NA.method ="mean")
X.marq <- scaleGen(obj.gi.marq, NA.method ="mean")
X.soc <- scaleGen(obj.gi.soc, NA.method ="mean")

pca.tuam <- dudi.pca(X.tuam, scannf = FALSE, nf = 2)
pca.gamb <- dudi.pca(X.gamb, scannf = FALSE, nf = 2)
pca.marq <- dudi.pca(X.marq, scannf = FALSE, nf = 2)
pca.soc <- dudi.pca(X.soc, scannf = FALSE, nf = 2)

coord.tuam<-cbind(pca.tuam$li,pop(obj.gi.tuam))
coord.gamb<-cbind(pca.gamb$li,pop(obj.gi.gamb))
coord.marq<-cbind(pca.marq$li,pop(obj.gi.marq))
coord.soc<-cbind(pca.soc$li,pop(obj.gi.soc))

names(coord.tuam)<-c("PC1","PC2","pop")
names(coord.gamb)<-c("PC1","PC2","pop")
names(coord.marq)<-c("PC1","PC2","pop")
names(coord.soc)<-c("PC1","PC2","pop")

centroids.tuam <- aggregate(cbind(PC1, PC2)~pop,coord.tuam,mean)
centroids.gamb <- aggregate(cbind(PC1, PC2)~pop,coord.gamb,mean)
centroids.marq <- aggregate(cbind(PC1, PC2)~pop,coord.marq,mean)
centroids.soc <- aggregate(cbind(PC1, PC2)~pop,coord.soc,mean)

pl.tuam <- merge(coord.tuam, centroids.tuam, by = "pop")
pl.gamb <- merge(coord.gamb, centroids.gamb, by = "pop")
pl.marq <- merge(coord.marq, centroids.marq, by = "pop")
pl.soc <- merge(coord.soc, centroids.soc, by = "pop")

colnames(pl.tuam) <- colnames(pl.gamb) <- colnames(pl.marq) <- colnames(pl.soc) <- c("pop", "PC1", "PC2", "PC1_ctr", "PC2_ctr")

a <- ggplot(data=pl.tuam,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids.tuam,size=3, alpha = 1) +
  geom_segment(data = pl.tuam, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  ggtitle("tuamotu")+
  scale_x_continuous(limits = c(-25,25))+
  scale_y_continuous(limits = c(-25,25))

b <- ggplot(data=pl.gamb,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids.gamb,size=3, alpha = 1) +
  geom_segment(data = pl.gamb, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  ggtitle("gambier")+
  geom_label_repel(data=centroids.gamb, label = centroids.gamb$pop, size = 5, )

c <- ggplot(data=pl.marq,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids.marq,size=3, alpha = 1) +
  geom_segment(data = pl.marq, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  ggtitle("marquesas")+
  geom_label_repel(data=centroids.marq, label = centroids.marq$pop, size = 5, )

d <- ggplot(data=pl.soc,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids.soc,size=3, alpha = 1) +
  geom_segment(data = pl.soc, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  ggtitle("society")+
  geom_label_repel(data=centroids.soc, label = centroids.soc$pop, size = 5, )

ggarrange(a,b,c,d, ncol = 2, nrow = 2)

#### relatedness ####

obj.glx.arch <- obj.glx
pop(obj.glx.arch) <- archipelago
obj.glx.arch <- gl.compliance.check(obj.glx.arch)

t1 <- gl.filter.callrate(obj.glx.arch)
res <- gl.grm(t1)
res2 <- gl.grm.network(res, t1, relatedness_factor = 0.004, 'mds')

#by archipelago
obj.glx.marq <- obj.glx[archipelago == "Marquesas",]
obj.glx.soc <- obj.glx[archipelago == "Society",]
obj.glx.tuam <- obj.glx[archipelago == "Tuamotus",]
obj.glx.gamb<- obj.glx[archipelago == "Gambier",]
obj.gl.wout.marq <- obj.glx[archipelago != "Marquesas",]

obj.glx.marq <- gl.compliance.check(obj.glx.marq)
obj.glx.soc <- gl.compliance.check(obj.glx.soc)
obj.glx.tuam <- gl.compliance.check(obj.glx.tuam)
obj.glx.gamb <- gl.compliance.check(obj.glx.gamb)
obj.gl.wout.marq <- gl.compliance.check(obj.gl.wout.marq)

t1.wm <- gl.filter.callrate(obj.gl.wout.marq)
res.wm <- gl.grm(t1.wm)
res2.wm <- gl.grm.network(res.wm, t1.wm, 'mds', relatedness_factor = 0)

t1.m <- gl.filter.callrate(obj.glx.marq)
res.m <- gl.grm(t1.m)
res2.m <- gl.grm.network(res.m, t1.m, 'mds', relatedness_factor = 0)

t1.s <- gl.filter.callrate(obj.glx.soc)
res.s <- gl.grm(t1.s)
res2.s <- gl.grm.network(res.s, t1.s, 'mds', relatedness_factor = 0)

t1.t <- gl.filter.callrate(obj.glx.tuam)
res.t <- gl.grm(t1.t)
res2.t <- gl.grm.network(res.t, t1.t, 'mds', relatedness_factor = 0.089)

t1.g <- gl.filter.callrate(obj.glx.gamb)
res.g <- gl.grm(t1.g)
res2.g <- gl.grm.network(res.g, t1.g, 'mds', relatedness_factor = 0)




#### DAPC outliers extraction ####
obj.glx.arch <- obj.glx
population <- as.character(pop(obj.glx))
tuamotus <- c("AHE", "ANA", "ARA", "KAT", "KAU", "MAN", "MOT", "RAR", "TAH", "TAK", "TEA", "TKP")
marquesas <- c("NHV-S", "NHV-P", "UAH", "UAP", "AUQ")
gambier <- c("GMBW", "MOR", "MARS", "PEI", "MAT", "RIK")
society <- c("SCI", "MOP")
australe <- c("RVV")
indo <- c("INDO")
obj.glx.arch <- obj.glx

for (i in tuamotus) {
  population[population == i] <- "Tuamotus"
}
for (i in marquesas) {
  population[population == i] <- "Marquesas"
}

for (i in gambier) {
  population[population == i] <- "Gambier"
}

for (i in society) {
  population[population == i] <- "Society"
}

for (i in australe) {
  population[population == i] <- "Austral"
}

for (i in indo) {
  population[population == i] <- "Indonesia"
}

pop(obj.glx.arch) <- population
obj.gl.wo.indo <- gl.drop.pop(obj.glx.arch, "Indonesia")
obj.glx.wo.indo <- gl.compliance.check(obj.gl.wo.indo)
obj.gi.wo.indo <- gl2gi(obj.glx.wo.indo)
dapc_marq <-dapc(obj.gi.wo.indo, pop = pop(obj.gi.wo.indo), n.da = 2, n.pca = 600)
axis1_marq_top <- head(order(dapc_marq$pca.loadings[,1]),200)
snp_axis1_marq_top <- attributes(dapc_marq$pca.loadings)$dimnames[[1]][axis1_marq_top]
snp_marq_top <- snp_axis1_marq_top
snp_marq_top <- strsplit(as.character(snp_marq_top), split = ".A")
snp_marq_top <- strsplit(as.character(snp_marq_top), split = ".T")
snp_marq_top <- unique(as.character(snp_marq_top))


obj.gl.pol.clus <- gl.drop.pop(obj.glx.arch, c("Indonesia","Marquesas", "Austral"))
obj.glx.pol.clus <- gl.compliance.check(obj.gl.pol.clus)
obj.gi.pol.clus <- gl2gi(obj.glx.pol.clus)

dapc_rest <- dapc(obj.gi.pol.clus, pop = pop(obj.gi.pol.clus), n.da = 2, n.pca = 600)
axis1_rest_top <- head(order(dapc_rest$pca.loadings[,1]),1500)
axis2_rest_top <- head(order(dapc_rest$pca.loadings[,2]),1500)
snp_axis1_rest_top <- attributes(dapc_rest$pca.loadings)$dimnames[[1]][axis1_rest_top]
snp_axis2_rest_top <- attributes(dapc_rest$pca.loadings)$dimnames[[1]][axis2_rest_top]
snp_rest_top <- c(snp_axis1_rest_top, snp_axis2_rest_top)
snp_rest_top <- strsplit(as.character(snp_rest_top), split = ".A")
snp_rest_top <- strsplit(as.character(snp_rest_top), split = ".T")
snp_rest_top <- unique(as.character(snp_rest_top))

snp_selected <- unique(c(snp_marq_top, snp_rest_top))

snp_selected.obj.gi <- obj.gi[loc= snp_selected, drop =TRUE]

test <- find.clusters(snp_selected.obj.gi, "AIC")

#control of the selection with an acp
pop(snp_selected.obj.gi) <- population
X_selection <- scaleGen(snp_selected.obj.gi, NA.method ="mean")
pca_selected <- dudi.pca(X_selection, scannf = FALSE, nf = 3)

new<-cbind(pca_selected$li,population)
str(new)
names(new)<-c("PC1","PC2", "PC3","pop")

centroids <- aggregate(cbind(PC1, PC2, PC3)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC3","PC1_ctr", "PC2_ctr", "PC3_ctr")


coly <- funky(length(unique(test$pop)))
a <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = coly)+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (4.70%)", y="PC2 (1.46%)") +
  ggtitle("all archipelagos")
plot(a)

b <- ggplot(data=test,aes(x=PC2,y=PC3,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = coly)+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC2, y = PC3, xend = PC2_ctr, yend = PC3_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC2 (1.46%)", y="PC3 (0.95%)")+
  ggtitle("")
plot(b)

obj.gl_rest <- obj.glx
pop(obj.gl_rest) <- population
obj.gl_rest <- gl.drop.pop(obj.gl_rest, c("Indonesia", "Marquesas"))
obj.gi_rest <- gl2gi(obj.gl_rest)
obj.gi_rest <- obj.gi_rest[loc= snp_selected, drop =TRUE]


X_selection_pol <- scaleGen(obj.gi_rest, NA.method ="mean")
pca_selected_pol <- dudi.pca(X_selection_pol, scannf = FALSE, nf = 3)

new<-cbind(pca_selected_pol$li,pop(obj.gi_rest))
str(new)
names(new)<-c("PC1","PC2", "PC3","pop")

centroids <- aggregate(cbind(PC1, PC2, PC3)~pop,new,mean)
test <- merge(new, centroids, by = "pop")
colnames(test) <- c("pop", "PC1", "PC2", "PC3","PC1_ctr", "PC2_ctr", "PC3_ctr")

coly <- funky(length(unique(test$pop)))
c <- ggplot(data=test,aes(x=PC1,y=PC2,col=pop)) + 
  geom_point(alpha=0.5,size=2) + 
  scale_color_manual(values = coly)+
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none") +
  stat_ellipse() +
  geom_point(data=centroids,size=3, alpha = 1) +
  geom_segment(data = test, aes(x = PC1, y = PC2, xend = PC1_ctr, yend = PC2_ctr), alpha = 0.3)+
  geom_label_repel(data=centroids, label = centroids$pop, size = 5, )+
  labs(x="PC1 (1.28%)", y="PC2 (1.07%)")+
  ggtitle("Only polynesian cluster")
plot(c)

ggarrange(a,b,c)
