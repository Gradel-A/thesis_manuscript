library(ggplot2)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

outliers.lim <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
}


load("/Users/antoinegradel/Downloads/fit1d_marq_tuam.RData")
load("/Users/antoinegradel/Downloads/fit1d_marq_soc.RData")
load("/Users/antoinegradel/Downloads/fit1d_marq_gamb.RData")
load("/Users/antoinegradel/Downloads/fit1d_pol.RData")
load("/Users/antoinegradel/Downloads/matrice_all.RData")
chr_pos_info <- colnames(matrice.marq_soc)
rm(matrice.marq_soc)

chromosome <- strsplit(as.character(chr_pos_info), split = "_")
chr_list <- "init"
start_pos <- seq(1,length(chromosome),10000)
for (i in start_pos[-length(start_pos)]) {
     chr_list <- c(chr_list, chromosome[[i]][1])  
}

index <- 1:length(fit1d.marq_gamb$points[,1])

data_ggplot.MG<- as.data.frame(cbind(index=index, pos = fit1d.marq_gamb$points[,1],chr=chr_list[-1]))
data_ggplot.MS<- as.data.frame(cbind(index=index, pos = fit1d.marq_soc$points[,1],chr=chr_list[-1]))
data_ggplot.MT<- as.data.frame(cbind(index=index, pos = fit1d.marq_tuam$points[,1],chr=chr_list[-1]))
data_ggplot.pol<- as.data.frame(cbind(index=index, pos = fit1d.pol$points[,1],chr=chr_list[-1]))


data_ggplot.MG$pos <- as.numeric(data_ggplot.MG$pos)
data_ggplot.MS$pos <- as.numeric(data_ggplot.MS$pos)
data_ggplot.MT$pos <-as.numeric(data_ggplot.MT$pos)
data_ggplot.pol$pos <-as.numeric(data_ggplot.pol$pos)

lims.MG <- outliers.lim(data_ggplot.MG$pos, 1.96)
lims.MS <- outliers.lim(data_ggplot.MS$pos, 1.96)
lims.MT <- outliers.lim(data_ggplot.MT$pos, 1.96)
lims.pol <- outliers.lim(data_ggplot.pol$pos, 1.96)

data_ggplot.MG$outliers <- rep(FALSE, nrow(data_ggplot.MG))
data_ggplot.MS$outliers <- rep(FALSE, nrow(data_ggplot.MS))
data_ggplot.MT$outliers <- rep(FALSE, nrow(data_ggplot.MT))
data_ggplot.pol$outliers <- rep(FALSE, nrow(data_ggplot.pol))

data_ggplot.MG$outliers[data_ggplot.MG$pos < lims.MG[1] | data_ggplot.MG$pos > lims.MG[2]] <- 0.05
data_ggplot.MS$outliers[data_ggplot.MS$pos < lims.MS[1] | data_ggplot.MS$pos > lims.MS[2]] <- 0.05
data_ggplot.MT$outliers[data_ggplot.MT$pos < lims.MT[1] | data_ggplot.MT$pos > lims.MT[2]] <- 0.05
data_ggplot.pol$outliers[data_ggplot.pol$pos < lims.pol[1] | data_ggplot.pol$pos > lims.pol[2]] <- 0.05

lims.MG <- outliers.lim(data_ggplot.MG$pos, 2.58)
lims.MS <- outliers.lim(data_ggplot.MS$pos, 2.58)
lims.MT <- outliers.lim(data_ggplot.MT$pos, 2.58)
lims.pol <- outliers.lim(data_ggplot.pol$pos, 2.58)

data_ggplot.MG$outliers[data_ggplot.MG$pos < lims.MG[1] | data_ggplot.MG$pos > lims.MG[2]] <- 0.01
data_ggplot.MS$outliers[data_ggplot.MS$pos < lims.MS[1] | data_ggplot.MS$pos > lims.MS[2]] <- 0.01
data_ggplot.MT$outliers[data_ggplot.MT$pos < lims.MT[1] | data_ggplot.MT$pos > lims.MT[2]] <- 0.01
data_ggplot.pol$outliers[data_ggplot.pol$pos < lims.pol[1] | data_ggplot.pol$pos > lims.pol[2]] <- 0.01


data_ggplot.MG$comp <- rep("Gambier", nrow(data_ggplot.MG))
data_ggplot.MS$comp <- rep("Society", nrow(data_ggplot.MS))
data_ggplot.MT$comp <- rep("Tuamotu", nrow(data_ggplot.MT))
data_ggplot.pol$comp <- rep("Polynesia", nrow(data_ggplot.pol))



rbind(data_ggplot.MG,data_ggplot.MS,data_ggplot.MT, data_ggplot.pol) -> data_ggplot

ggplot(data = data_ggplot.MG, aes(x=as.numeric(index), y=abs(pos)))+
  geom_point(aes(fill = as.factor(chr), color = ifelse(outliers, "grey", "red")),
size = 2, shape = 21, stroke = 1)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)))+
  scale_color_manual(values = c("red", "white"))

ggplot(data = data_ggplot.MS, aes(x=as.numeric(index), y=abs(pos)))+
  geom_point(aes(fill = as.factor(chr), color = ifelse(outliers, "grey", "red")),
             size = 2, shape = 21, stroke = 1)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)))+
  scale_color_manual(values = c("red", "white"))

ggplot(data = data_ggplot.MT, aes(x=as.numeric(index), y=abs(pos)))+
  geom_point(aes(fill = as.factor(chr), color = as.factor(outliers)),
             size = 2, shape = 21, stroke = 1)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)))+
  scale_color_manual(values = c("white","brown", "red"))


ggplot()+
  geom_point(data = data_ggplot.MT, aes(x=as.numeric(index), y=abs(pos), fill = as.factor(chr), color = as.factor(outliers)),
            size = 3, shape = 21, stroke = 1.5)+
  geom_point(data = data_ggplot.MG, aes(x=as.numeric(index), y=abs(pos),fill = as.factor(chr), color = as.factor(outliers)),
            size = 3, shape = 22, stroke = 1.5)+
  geom_point(data = data_ggplot.MS, aes(x=as.numeric(index), y=abs(pos),fill = as.factor(chr), color = as.factor(outliers)),
            size = 3, shape = 23, stroke = 1.5)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)), )+
  scale_color_manual(values = c("white","green", "red"), name = "Significativity", labels = c("","0.01", "0.05"))+
  guides(fill=FALSE) +
  scale_shape_manual(values = c(21,22,23), labels = c("MT" ,"MG", "MS"), name = "comparison with")+
  theme_bw()+
  xlab("SNPs window position")+
  ylab("MDS 1")+
  theme(axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_text(size=16, colour="black", vjust=-0.1),
        axis.title.y=element_text(size=16, colour="black"),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour="black"),
        legend.text = element_text(size=16, colour="black"),
        legend.background = element_rect(fill = NA , colour = NA),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "bottom") 

a <- ggplot()+
    geom_point(data = data_ggplot.MS, aes(x=as.numeric(index), y=abs(pos), fill = as.factor(chr), color = as.factor(outliers)),
            size = 3, shape = 21, stroke = 1.5)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)), )+
  scale_color_manual(values = c("white","green", "red"), name = "Significativity", labels = c("","0.01", "0.05"))+
  guides(fill=FALSE) +
  theme_bw()+
  xlab("SNPs window position")+
  ylab("MDS 1")+
  ggtitle("Marquesas-Society")+
  theme(axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, colour="black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour="black"),
        legend.text = element_text(size=16, colour="black"),
        legend.background = element_rect(fill = NA , colour = NA),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") 
b <- ggplot()+
    geom_point(data = data_ggplot.MT, aes(x=as.numeric(index), y=abs(pos), fill = as.factor(chr), color = as.factor(outliers)),
            size = 3, shape = 21, stroke = 1.5)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)), )+
  scale_color_manual(values = c("white","green", "red"), name = "Significativity", labels = c("","0.01", "0.05"))+
  guides(fill=FALSE) +
  theme_bw()+
  ggtitle("Marquesas-Tuamotu")+
  ylab("MDS 1")+
  theme(axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, colour="black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour="black"),
        legend.text = element_text(size=16, colour="black"),
        legend.background = element_rect(fill = NA , colour = NA),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") 

c <- ggplot()+
  geom_point(data = data_ggplot.MG, aes(x=as.numeric(index), y=abs(pos),fill = as.factor(chr), color = as.factor(outliers)),
             size = 3, shape = 21, stroke = 1.5)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)), )+
  scale_color_manual(values = c("white","green", "red"), name = "Significativity", labels = c("","0.01", "0.05"))+
  guides(fill=FALSE) +
  theme_bw()+
  xlab("SNPs window position")+
  ylab("MDS 1")+
  ggtitle("Marquesas-Gambier")+
  theme(axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, colour="black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour="black"),
        legend.text = element_text(size=16, colour="black"),
        legend.background = element_rect(fill = NA , colour = NA),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") 

d <- ggplot()+
  geom_point(data = data_ggplot.pol, aes(x=as.numeric(index), y=abs(pos),fill = as.factor(chr), color = as.factor(outliers)),
             size = 3, shape = 21, stroke = 1.5)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)), )+
  scale_color_manual(values = c("white","green", "red"), name = "Significativity", labels = c("","0.01", "0.05"))+
  guides(fill=FALSE) +
  theme_bw()+
  xlab("SNPs window position")+
  ylab("MDS 1")+
  ggtitle("Society-Gambier-Tuamotu")+
  theme(axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, colour="black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        plot.title = element_text(size=16, colour="black"),
        legend.title = element_text(size=16, colour="black"),
        legend.text = element_text(size=16, colour="black"),
        legend.background = element_rect(fill = NA , colour = NA),
        axis.ticks.x = element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") 
  
ggarrange(a,b,c,d, ncol = 1)

model <- aov(data = data_ggplot, abs(pos)~comp)
myhsd <- HSD.test(model, "comp", group = TRUE, alpha = 0.001)
myhsd$groups


# et si tu veux qu'un chromosome
liste_chr <- unique(datatoprocess.marq_soc$chr)
for (i in liste_chr){
a <- ggplot(data = datatoprocess.marq_soc[datatoprocess.marq_soc$chr == i,], aes(x=as.numeric(index), y=pos))+
  geom_point(aes(fill = as.factor(chr), color = ifelse(outliers, "grey", "red")),
             size = 2, shape = 21, stroke = 1)+
  scale_fill_manual(values = c(rep(c("black", "grey"),7)))+
  scale_color_manual(values = c("red", "white"))

plot(a)
}




