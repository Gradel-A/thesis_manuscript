library(vcfR)
library(dartR)

obj.vcf.test <- read.vcfR("wall_genome_vcf/allgenome_pop_filtered_10percent_NA.vcf.recode.vcf") 

ID <- ((dimnames(obj.vcf.test@gt))[[2]])[-1]
ID_df <- as.data.frame(ID)
pop <- strsplit(as.character(ID_df$ID), split = "_")
POP <- as.data.frame(t(as.data.frame(pop)))
pop_list <- as.factor(POP$V1)
table(pop_list)

obj.gl.test <- vcfR2genlight(obj.vcf.test)
rm(obj.vcf.test)
obj.gi.test <- gl2gi(obj.gl.test)
pop(obj.gi.test) <- pop_list
rm(obj.gl.test)

obj.gi.marq_soc <- obj.gi.test[pop=c("MARQ","SOC")]
obj.gi.marq_tuam <- obj.gi.test[pop=c("MARQ","TUAM")]
obj.gi.marq_gamb <- obj.gi.test[pop=c("MARQ","GAMB")]
obj.gi.pol <- obj.gi.test[pop=c("SOC","GAMB","TUAM")]

rm(obj.gi.test)

ligne.marq_soc <- seq(1,dim(obj.gi.marq_soc@tab)[2],2) 
ligne.marq_tuam <- seq(1,dim(obj.gi.marq_tuam@tab)[2],2) 
ligne.marq_gamb <- seq(1,dim(obj.gi.marq_gamb@tab)[2],2) 
ligne.pol <- seq(1,dim(obj.gi.pol@tab)[2],2) 

matrice.marq_soc <- obj.gi.marq_soc@tab[,ligne.marq_soc]
matrice.marq_tuam <- obj.gi.marq_tuam@tab[,ligne.marq_tuam]
matrice.marq_gamb <- obj.gi.marq_gamb@tab[,ligne.marq_gamb]
matrice.pol <- obj.gi.pol@tab[,ligne.pol]

save(matrice.marq_soc, file="/home1/scratch/agradel/wall_genome_vcf/matrice_marq_soc.RData")
save(matrice.marq_tuam, file="/home1/scratch/agradel/wall_genome_vcf/matrice_marq_tuam.RData")
save(matrice.marq_gamb, file="/home1/scratch/agradel/wall_genome_vcf/matrice_marq_gamb.RData")
save(matrice.pol, file="/home1/scratch/agradel/wall_genome_vcf/matrice_pol.RData")




