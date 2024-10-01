library(vcfR)
library(dartR)

obj10.vcf <- read.vcfR("./wall_genome_vcf/allgenome_pop_filtered_10percent_NA.vcf.recode.vcf")
save(obj10.vcf, file = "./wall_genome_vcf/obj10_vcf_10.RData")
load("./wall_genome_vcf/obj10_vcf.RData")


obj10.gl <- vcfR2genlight(obj10.vcf)
ID <- ((dimnames(obj10.vcf@gt))[[2]])[-1]
ID_df <- as.data.frame(ID)
pop <- strsplit(as.character(ID_df$ID), split = "_")
POP <- as.data.frame(t(as.data.frame(pop)))
pop_list <- as.factor(POP$V1)
pop(obj10.gl) <- pop_list
table(pop_list)

save(obj10.gl, file = "./wall_genome_vcf/obj10_gl_10.RData")
load("./wall_genome_vcf/obj10_gl.RData")
rm(obj10.vcf)

obj10.glx <- gl.compliance.check(obj10.gl)
rm(obj10.gl)

save(obj10.glx, file = "./wall_genome_vcf/obj10_glx_10.RData")
load("./wall_genome_vcf/obj10_glx_10.RData")
data_heteroZ_pop_10 <- gl.report.heterozygosity(obj10.glx)
data_heteroZ_ind_10 <- gl.report.heterozygosity(obj10.glx, method = 'ind')


save(data_heteroZ_pop_10, file = "./wall_genome_vcf/het_pop_10.RData")
save(data_heteroZ_ind_10, file = "./wall_genome_vcf/het_ind_10.RData")

data_fst_test_10 <- gl.fst.pop(obj10.glx, nboots = 0, percent = 95, verbose = 3)
data_fst_test_10[["Fsts"]]
data_fst_test_10[["Pvalues"]]

save(data_fst_test_10, file = "./wall_genome_vcf/fst_info_10.RData")
