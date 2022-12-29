#############################################################################
#############################################################################
# 0.0 IMPORTS
require(vegan)
require(ncf)
require(SoDA)
require(ggplot2)
require(cowplot)
require(ggpubr)

#############################################################################
#############################################################################
# 0.2 DATA COLLECTION 
setwd(choose.dir())
getwd()

# Correlogram spline univariate with function "spline.correlog"
ab<-read.table("abundance.txt", h=T)
xy<-read.table("xy.txt", h=T)
# data transformation by total
ab.rel<-decostand(ab, "total")

# Correlogram spline univariate with ggplot2
Fin_gg<- read.table("Fin_ggplot.txt", h=T)
Poc_gg<- read.table("Poc_ggplot.txt", h=T)
Fol_gg<- read.table("Fol_ggplot.txt", h=T)
Rre_gg<- read.table("Rre_ggplot.txt", h=T)
Prh_gg<- read.table("Prh_ggplot.txt", h=T)
Fin_ggg<- read.table("Psu_ggplot.txt", h=T)
Rmi_gg<- read.table("Rmi_ggplot.txt", h=T)
Fca_gg<- read.table("Fca_ggplot.txt", h=T)
Fst_gg <- read.table("Fst_ggplot.txt", h=T)
Rqu_gg <- read.table("Rqu_ggplot.txt", h=T)

#############################################################################
#############################################################################
# 1.0 Univariate Spline Correlograms

# Faramea involucellata (Fin)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fin")]

library(SoDA)
xy<-geoXY (x, y)

spline1 <- spline.correlog(xy[,"X"],
                            xy[,"Y"], 
                            ab.rel[,"Fin"],
                            resamp=10000, 
                            latlon=FALSE, 
                            npoints = 203)

ab_fin<-as.data.frame(spline1$real)
ac_fin() <- as.matrix(spline1$boot$boot.summary$predicted$y)

spline1[["real"]]
ab_fin <- spline1[["real"]]
ac_fin <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab_fin, "C:/Users/asus/Desktop/spl_fin1.txt")
write.table(ac_fin, "C:/Users/asus/Desktop/spl_fin2.txt")

###########################################
# Faramea oligantha (Fol)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fol")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fol"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fol<-as.data.frame(spline2$real)
ac_fol <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_fol <- spline1[["real"]]
ac_fol <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_fol1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_fol2.txt")

###########################################
# Faramea stipulacea (Fst)

x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fst")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fst"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fst<-as.data.frame(spline2$real)
ac_fst <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_fst <- spline1[["real"]]
ac_fst <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_fst1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_fst2.txt")
###########################################
# Faramea martiana (Fma)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fma")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fma"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fma<-as.data.frame(spline2$real)
ac_fma <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_fma <- spline1[["real"]]
ac_fma <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_fma1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_fma2.txt")

###########################################
# Palicourea octocuspis (Poc)

x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Poc")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Poc"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_poc<-as.data.frame(spline2$real)
ac_poc <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_poc <- spline1[["real"]]
ac_poc <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_poc1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_poc2.txt")

###########################################
# Psychotria rhytidocarpa (Prh)

x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Prh"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_prh<-as.data.frame(spline2$real)
ac_prh <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_prh <- spline1[["real"]]
ac_prh <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_prh1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_prh2.txt")

###########################################
# Psychotria subspathacea (Psu)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Psu")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Psu"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_psu<-as.data.frame(spline2$real)
ac_psu <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_psu <- spline1[["real"]]
ac_psu <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_psu1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_psu2.txt")

###########################################
# Rudgea quisquiliae (Rqu)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Rqu")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Rqu"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_rqu<-as.data.frame(spline2$real)
ac_rqu <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_rqu <- spline1[["real"]]
ac_rqu <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_rqu1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_rqu2.txt")

###########################################
# Rudgea reflexa (Rre)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Rre")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Rre"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_rre<-as.data.frame(spline2$real)
ac_rre <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_rre <- spline1[["real"]]
ac_rre <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_rre1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_rre2.txt")

###########################################
# Rudgea minutifolia (Rmi)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Rmi")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.cormilog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Rmi"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_rmi<-as.data.frame(spline2$real)
ac_rmi <- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_rmi <- spline1[["real"]]
ac_rmi <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab, "C:/Users/asus/Desktop/spl_rmi1.txt")
write.table(ac, "C:/Users/asus/Desktop/spl_rmi2.txt")


############################################################################
# Splines Correlograms with ggplot2 for the 10 most abundant species
p1<-ggplot(Fin_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 1, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b),fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p1



p1<- p1 + 
  xlab(paste("")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p1

p1 <- p1 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p1

p1<- p1+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p1

#############################################################################
p2<-ggplot(Poc_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p2


p2<- p2 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste(" ")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p2

p2 <- p2 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p2

p2<- p2+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p2

#############################################################################
p3<-ggplot(Fol_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p3



p3<- p3 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p3

p3 <- p3 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p3

p3<- p3+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p3

#############################################################################
p4<-ggplot(Rre_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p4



p4 <- p4 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste(" ")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p4

p4 <- p4 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p4

p4<- p4+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p4

#############################################################################
p5<-ggplot(Prh_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p5


p5 <- p5 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p5

p5 <- p5 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p5

p5<- p5+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p5

#############################################################################
p6<-ggplot(Fin_ggg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 1, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b),fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p6



p6<- p6 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste(" ")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p6

p6 <- p6 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p6

p6<- p6+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p6

#############################################################################
p7<-ggplot(Rmi_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p7



p7<- p7 + 
  xlab(paste("")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p7

p7 <- p7 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p7

p7<- p7+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p7


#############################################################################
p8<-ggplot(Fst_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p8



p8<- p8 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste(" ")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p8

p8 <- p8 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p8

p8<- p8+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p8


#############################################################################
p9<-ggplot(Fca_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p9



p9 <- p9 + 
  xlab(paste("Distance (m)")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p9

p9 <- p9 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p9

p9<- p9+ theme(text=element_text(size=16, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p9

#############################################################################
p10<-ggplot(Rqu_gg, aes(d,a)) +
  scale_x_continuous(limits=c(0,300)) +
  expand_limits(y=c(-1,1)) +
  scale_y_continuous(limits=c(-1,1)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 12, angle = 0),
        axis.text.y = element_text(color ="black", size = 12, angle = 0),
        plot.title = element_text(hjust = 0.5, size = 15))+ 
  geom_ribbon(aes(ymin=a, ymax=b), fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p10


p10 <- p10 + 
  xlab(paste("Distance (m)")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste(" ")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p10

p10 <- p10 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p10

p10<- p10+ theme(text=element_text(size=16, 
                                   #       family="Comic Sans MS"))
                                   #       family="CM Roman"))
                                   #       family="TT Times New Roman"))
                                   #       family="Sans"))
                                   family="TT Times New Roman"))

p10
#############################################################################
#############################################################################
#############################################################################
# Saving the 10 splines correlograms in just one figure

tiff(filename="Figure1.tiff", res=600, height=450/72*500, width=150/72*800, compression= "lzw")

# simple grid with labels and aligned plots
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
          labels = c("A", "B", "C", "D", "E",
                     "F", "G", "H", "I", "J"),
          ncol = 2)

dev.off()
