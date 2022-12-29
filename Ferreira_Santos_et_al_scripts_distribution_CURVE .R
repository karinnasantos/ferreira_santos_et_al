#############################################################################
#############################################################################
# 0.0 IMPORTS
require(vegan)
require(BiodiversityR)

#############################################################################
#############################################################################
# 0.2 DATA COLLECTION 
setwd(choose.dir())
getwd()

# abundance rank of all sampling units
abundance <-read.table ("abundance.txt",h=T)

# species abundance rank in each topographic habitat
abundance_valley <-read.table ("valley.txt",h=T)
abundance_slope <-read.table ("slope.txt",h=T)
abundance_ridge <-read.table ("ridge.txt",h=T)

# abundance rank of the 50 most abundant sampling units
ab2 <-read.table ("sampling.2.txt",h=T)
ab3 <-read.table ("sampling.3.txt",h=T)
ab4 <-read.table ("sampling.4.txt",h=T)
ab5 <-read.table ("sampling.5.txt",h=T)
ab6 <-read.table ("sampling.6.txt",h=T)
ab7 <-read.table ("sampling.7.txt",h=T)
ab8 <-read.table ("sampling.8.txt",h=T)
ab9 <-read.table ("sampling.9.txt",h=T)
ab10 <-read.table ("sampling.10.txt",h=T)
ab11 <-read.table ("sampling.11.txt",h=T)
ab12 <-read.table ("sampling.12.txt",h=T)
ab13 <-read.table ("sampling.13.txt",h=T)
ab14 <-read.table ("sampling.14.txt",h=T)
ab15 <-read.table ("sampling.15.txt",h=T)
ab16 <-read.table ("sampling.16.txt",h=T)
ab17 <-read.table ("sampling.17.txt",h=T)
ab18 <-read.table ("sampling.18.txt",h=T)
ab19 <-read.table ("sampling.19.txt",h=T)
ab20 <-read.table ("sampling.20.txt",h=T)
ab21 <-read.table ("sampling.21.txt",h=T)
ab22 <-read.table ("sampling.22.txt",h=T)
ab23 <-read.table ("sampling.23.txt",h=T)
ab24 <-read.table ("sampling.24.txt",h=T)
ab25 <-read.table ("sampling.25.txt",h=T)
ab26 <-read.table ("sampling.26.txt",h=T)
ab27 <-read.table ("sampling.27.txt",h=T)
ab28 <-read.table ("sampling.28.txt",h=T)
ab29 <-read.table ("sampling.29.txt",h=T)
ab30 <-read.table ("sampling.30.txt",h=T)
ab31 <-read.table ("sampling.31.txt",h=T)
ab32 <-read.table ("sampling.32.txt",h=T)
ab33 <-read.table ("sampling.33.txt",h=T)
ab34 <-read.table ("sampling.34.txt",h=T)
ab35 <-read.table ("sampling.35.txt",h=T)
ab36 <-read.table ("sampling.36.txt",h=T)
ab37 <-read.table ("sampling.37.txt",h=T)
ab38 <-read.table ("sampling.38.txt",h=T)
ab39 <-read.table ("sampling.39.txt",h=T)
ab40 <-read.table ("sampling.40.txt",h=T)
ab41 <-read.table ("sampling.41.txt",h=T)
ab42 <-read.table ("sampling.42.txt",h=T)
ab43 <-read.table ("sampling.43.txt",h=T)
ab44 <-read.table ("sampling.44.txt",h=T)
ab45 <-read.table ("sampling.45.txt",h=T)
ab46 <-read.table ("sampling.46.txt",h=T)
ab47 <-read.table ("sampling.47.txt",h=T)
ab48 <-read.table ("sampling.48.txt",h=T)
ab49 <-read.table ("sampling.49.txt",h=T)
ab50 <-read.table ("sampling.50.txt",h=T)
ab51 <-read.table ("sampling.51.txt",h=T)

#############################################################################
#############################################################################
# 1.0 abundance rank of all sampling units

rank.abundance <-rankabundance(abundance)

#############################################################################
#############################################################################
# 2.0 species abundance rank in each topographic habitat

rank.abundance_valley <-rankabundance(abundance_valley)
rank.abundance_valley
########################################
rank.abundance_slope <-rankabundance(abundance_slope)
rank.abundance_slope
########################################
rank.abundance_ridge <-rankabundance(abundance_ridge)
rank.abundance_ridge

#############################################################################
#############################################################################
# 3.0 abundance rank of the 50 most abundant sampling units

rank.mf2 <-rankabundance(ab2)
rank.mf2 
########################################
rank.mf3 <-rankabundance(ab3)
rank.mf3 
########################################
rank.mf4 <-rankabundance(ab4)
rank.mf4 
########################################
rank.mf5 <-rankabundance(ab5)
rank.mf5 
########################################
rank.mf6 <-rankabundance(ab62)
rank.mf6 
########################################
rank.mf7 <-rankabundance(ab7)
rank.mf7 
########################################
rank.mf8 <-rankabundance(ab8)
rank.mf8 
########################################
rank.mf9 <-rankabundance(ab9)
rank.mf9 
########################################
rank.mf10 <-rankabundance(ab10)
rank.mf10 
########################################
rank.mf11 <-rankabundance(ab11)
rank.mf11 
########################################
rank.mf12 <-rankabundance(ab12)
rank.mf12 
########################################
rank.mf13 <-rankabundance(ab13)
rank.mf13 
########################################
rank.mf14 <-rankabundance(ab14)
rank.mf14 
########################################
rank.mf15 <-rankabundance(ab15)
rank.mf15 
########################################
rank.mf16 <-rankabundance(ab16)
rank.mf16 
########################################
rank.mf17 <-rankabundance(ab17)
rank.mf17 
########################################
rank.mf18 <-rankabundance(ab18)
rank.mf18 
########################################
rank.mf19 <-rankabundance(ab19)
rank.mf19 
########################################
rank.mf20 <-rankabundance(ab20)
rank.mf20 
########################################
rank.mf21 <-rankabundance(ab21)
rank.mf21 
########################################
rank.mf22 <-rankabundance(ab22)
rank.mf22 
########################################
rank.mf23 <-rankabundance(ab23)
rank.mf23 
########################################
rank.mf24 <-rankabundance(ab24)
rank.mf24 
########################################
rank.mf25 <-rankabundance(ab25)
rank.mf25 
########################################
rank.mf26 <-rankabundance(ab26)
rank.mf26 
########################################
rank.mf27 <-rankabundance(ab27)
rank.mf27 
########################################
rank.mf28 <-rankabundance(ab28)
rank.mf28 
########################################
rank.mf29 <-rankabundance(ab29)
rank.mf29 
########################################
rank.mf30 <-rankabundance(ab30)
rank.mf30 
########################################
rank.mf31 <-rankabundance(ab31)
rank.mf31 
########################################
rank.mf32 <-rankabundance(ab32)
rank.mf32 
########################################
rank.mf33 <-rankabundance(ab33)
rank.mf33 
########################################
rank.mf34 <-rankabundance(ab34)
rank.mf34 
########################################
rank.mf35 <-rankabundance(ab35)
rank.mf35 
########################################
rank.mf36 <-rankabundance(ab36)
rank.mf36 
########################################
rank.mf37 <-rankabundance(ab37)
rank.mf37 
########################################
rank.mf38 <-rankabundance(ab38)
rank.mf38 
########################################
rank.mf39 <-rankabundance(ab39)
rank.mf39 
########################################
rank.mf40 <-rankabundance(ab40)
rank.mf40 
########################################
rank.mf41 <-rankabundance(ab41)
rank.mf41 
########################################
rank.mf42 <-rankabundance(ab42)
rank.mf42 
########################################
rank.mf43 <-rankabundance(ab43)
rank.mf43 
########################################
rank.mf44 <-rankabundance(ab44)
rank.mf44 
########################################
rank.mf45 <-rankabundance(ab45)
rank.mf45 
########################################
rank.mf46 <-rankabundance(ab46)
rank.mf46 
########################################
rank.mf47 <-rankabundance(ab47)
rank.mf47 
########################################
rank.mf48 <-rankabundance(ab48)
rank.mf48 
########################################
rank.mf49 <-rankabundance(ab49)
rank.mf49 
########################################
rank.mf50 <-rankabundance(ab50)
rank.mf50 
########################################
rank.mf51 <-rankabundance(ab51)
rank.mf51 
########################################