library(lostruct)
library(ggplot2)
library(dartR)


load("/home1/scratch/agradel/wall_genome_vcf/matrice_marq_gamb.RData") #attention il faut les transposer pour les rentrer dans lostruct

pcs <- eigen_windows(t(matrice.marq_gamb),k=2, win=10000)
pcdist <- pc_dist(pcs,npc=2)
fit1d.marq_gamb <-cmdscale(pcdist,eig=TRUE, k=1)

save(fit1d.marq_gamb, file= "wall_genome_vcf/fit1d_marq_gamb.RData")
