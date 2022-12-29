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

# Correlogram spline bivariate with function "spline.correlog"
ab<-read.table("abundance.txt", h=T)
xy<-read.table("xy.txt", h=T)
# data transformation by total
ab.rel<-decostand(ab, "total")

# Correlogram spline bivariate with ggplot2
Fin_Poc_gg<- read.table("1_Fin_poc_ggplot.txt", h=T)
Fin_Fol_gg<- read.table("2_Fin_Fol_ggplot.txt", h=T)
Fin_Rre_gg<- read.table("3_Fin_Rre_ggplot.txt", h=T)
Fin_Prh_gg <- read.table("4_Fin_Prh_ggplot.txt", h=T)
Fol_Poc_gg <- read.table("5_Poc_Fol_ggplot.txt", h=T)
Poc_Rre_gg <- read.table("6_Poc_Rre_ggplot.txt", h=T)
Poc_Prh_gg <- read.table("7_Poc_Prh_ggplot.txt", h=T)
Fol_Rre_gg <- read.table("8_Fol_Rre_ggplot.txt", h=T)
Fol_Prh_gg <- read.table("9_Fol_Prh_ggplot.txt", h=T)
Rre_Prh_gg <- read.table("10_Rre_Prh_ggplot.txt", h=T)

# Bubble plot
Fin_Poc <-read.table("1_Fin_Poc.txt", h=T)
Fin_Fol <-read.table("2_Fin_Fol.txt", h=T)
Fin_Rre <-read.table("3_Fin_Rre.txt", h=T)
Fin_Prh <-read.table("4_Fin_Prh.txt", h=T)
Poc_Fol <-read.table("5_Poc_Fol.txt", h=T)
Poc_Rre <-read.table("6_Poc_Rre.txt", h=T)
Poc_Prh <-read.table("7_Poc_Prh.txt", h=T)
Fol_Rre <-read.table("8_Fol_Rre.txt", h=T)
Fol_Prh <-read.table("9_Fol_Prh.txt", h=T)
Rre_Prh <-read.table("10_Rre_Prh.txt", h=T)

#############################################################################
#############################################################################
#############################################################################
#############################################################################
# 1.0 Correlogram spline bivariate with function "spline.correlog"
# Faramea involucellata (Fin) X Paclicourea octocuspis (Poc)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fin")]
w<- ab[c("Poc")]

library(SoDA)
xy<-geoXY (x, y)

spline1 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fin"],
                           ab.rel[,"Poc"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fin_poc<-as.data.frame(spline1$real)
ac_fin_poc<- as.matrix(spline1$boot$boot.summary$predicted$y)

spline1[["real"]]
ab_fin_poc <- spline1[["real"]]
ac_fin_poc <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab_fin_poc, "C:/Users/asus/Desktop/spl_fin1.txt")
write.table(ac_fin_poc, "C:/Users/asus/Desktop/spl_fin2.txt")

#############################################################################
# Faramea involucellata (Fin) X Faramea oligahtha (Fol)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fin")]
w<- ab[c("Fol")]

library(SoDA)
xy<-geoXY (x, y)

spline2 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fin"],
                           ab.rel[,"Fol"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fin_fol<-as.data.frame(spline2$real)
ac_fin_fol<- as.matrix(spline2$boot$boot.summary$predicted$y)

spline2[["real"]]
ab_fin_fol <- spline2[["real"]]
ac_fin_fol <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab_fin_fol, "C:/Users/asus/Desktop/spl_fin1.txt")
write.table(ac_fin_fol, "C:/Users/asus/Desktop/spl_fin2.txt")

#############################################################################
# Faramea involucellata (Fin) X Rudgea reflexa (Rre)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fin")]
w<- ab[c("Rre")]

library(SoDA)
xy<-geoXY (x, y)

spline3 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fin"],
                           ab.rel[,"Rre"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fin_rre<-as.data.frame(spline3$real)
ac_fin_rre<- as.matrix(spline3$boot$boot.summary$predicted$y)

spline3[["real"]]
ab_fin_rre <- spline3[["real"]]
ac_fin_rre <- as.matrix(spline3$boot$boot.summary$predicted$y)
write.table(ab_fin_rre, "C:/Users/asus/Desktop/spl_fin1.txt")
write.table(ac_fin_rre, "C:/Users/asus/Desktop/spl_fin2.txt")

#############################################################################
# Faramea involucellata (Fin) X Psychotria rhytidocarpa (Prh)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fin")]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline4 <- spline.coprhlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fin"],
                           ab.rel[,"Prh"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fin_prh<-as.data.frame(spline4$real)
ac_fin_prh<- as.matrix(spline4$boot$boot.summary$predicted$y)

spline4[["real"]]
ab_fin_prh <- spline4[["real"]]
ac_fin_prh <- as.matrix(spline4$boot$boot.summary$predicted$y)
write.table(ab_fin_prh, "C:/Users/asus/Desktop/spl_fin1.txt")
write.table(ac_fin_prh, "C:/Users/asus/Desktop/spl_fin2.txt")

#############################################################################
# Palicourea octocuspis (Poc) X Faramea oligantha (Fol)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Poc")]
w<- ab[c("Fol")]

library(SoDA)
xy<-geoXY (x, y)

spline5 <- spline.cofollog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Poc"],
                           ab.rel[,"Fol"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_poc_fol<-as.data.frame(spline5$real)
ac_poc_fol<- as.matrix(spline5$boot$boot.summary$predicted$y)

spline5[["real"]]
ab_poc_fol <- spline5[["real"]]
ac_poc_fol <- as.matrix(spline5$boot$boot.summary$predicted$y)
write.table(ab_poc_fol, "C:/Users/asus/Desktop/spl_poc_fol1.txt")
write.table(ac_poc_fol, "C:/Users/asus/Desktop/spl_poc_fol2.txt")

#############################################################################
# Palicourea octocuspis (Poc) X Rudgea reflexa (Rre)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Poc")]
w<- ab[c("Rre")]

library(SoDA)
xy<-geoXY (x, y)

spline6 <- spline.correlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Poc"],
                           ab.rel[,"Rre"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_poc_rre<-as.data.frame(spline6$real)
ac_poc_rre<- as.matrix(spline6$boot$boot.summary$predicted$y)

spline6[["real"]]
ab_poc_rre <- spline6[["real"]]
ac_poc_rre <- as.matrix(spline6$boot$boot.summary$predicted$y)
write.table(ab_poc_rre, "C:/Users/asus/Desktop/spl_poc_rre1.txt")
write.table(ac_poc_rre, "C:/Users/asus/Desktop/spl_poc_rre2.txt")

#############################################################################
# Palicourea octocuspis (Poc) X Psychotria rhytidocarpa (Prh)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Poc")]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline7 <- spline.coprhlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Poc"],
                           ab.rel[,"Prh"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_poc_prh<-as.data.frame(spline7$real)
ac_poc_prh<- as.matrix(spline7$boot$boot.summary$predicted$y)

spline7[["real"]]
ab_poc_prh <- spline7[["real"]]
ac_poc_prh <- as.matrix(spline7$boot$boot.summary$predicted$y)
write.table(ab_poc_prh, "C:/Users/asus/Desktop/spl_poc_prh1.txt")
write.table(ac_poc_prh, "C:/Users/asus/Desktop/spl_poc_prh2.txt")

#############################################################################
# Faramea oligantha (Fol) X Rudgea reflexa (Rre)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fol")]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline8 <- spline.coprhlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fol"],
                           ab.rel[,"Rre"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fol_prh<-as.data.frame(spline8$real)
ac_fol_prh<- as.matrix(spline8$boot$boot.summary$predicted$y)

spline8[["real"]]
ab_fol_prh <- spline8[["real"]]
ac_fol_prh <- as.matrix(spline8$boot$boot.summary$predicted$y)
write.table(ab_fol_prh, "C:/Users/asus/Desktop/spl_fol_rre1.txt")
write.table(ac_fol_prh, "C:/Users/asus/Desktop/spl_fol_rre2.txt")

#############################################################################
# Faramea oligantha (Fol) X Psychotria rhytidocarpa (Prh)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Fol")]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline9 <- spline.coprhlog(xy[,"X"],
                           xy[,"Y"], 
                           ab.rel[,"Fol"],
                           ab.rel[,"Prh"],
                           resamp=10000, 
                           latlon=FALSE, 
                           npoints = 203)

ab_fol_prh<-as.data.frame(spline9$real)
ac_fol_prh<- as.matrix(spline9$boot$boot.summary$predicted$y)

spline9[["real"]]
ab_fol_prh <- spline9[["real"]]
ac_fol_prh <- as.matrix(spline9$boot$boot.summary$predicted$y)
write.table(ab_fol_prh, "C:/Users/asus/Desktop/spl_fol_prh1.txt")
write.table(ac_fol_prh, "C:/Users/asus/Desktop/spl_fol_prh2.txt")

#############################################################################
# Rudgea reflexa (Rre) X Psychotria rhytidocarpa (Prh)
x<- xy[,"x"]
y<- xy[,"y"]
w<- ab[c("Rre")]
w<- ab[c("Prh")]

library(SoDA)
xy<-geoXY (x, y)

spline10 <- spline.coprhlog(xy[,"X"],
                            xy[,"Y"], 
                            ab.rel[,"Rre"],
                            ab.rel[,"Prh"],
                            resamp=10000, 
                            latlon=FALSE, 
                            npoints = 203)

ab_rre_prh<-as.data.frame(spline10$real)
ac_rre_prh<- as.matrix(spline10$boot$boot.summary$predicted$y)

spline10[["real"]]
ab_rre_prh <- spline10[["real"]]
ac_rre_prh <- as.matrix(spline10$boot$boot.summary$predicted$y)
write.table(ab_rre_prh, "C:/Users/asus/Desktop/spl_rre_prh1.txt")
write.table(ac_rre_prh, "C:/Users/asus/Desktop/spl_rre_prh2.txt")




#############################################################################
#############################################################################
#############################################################################
#############################################################################
# 2.0 Correlogram spline bivariate with ggplot2

p1<-ggplot(Fin_Poc_gg, aes(d,a)) +
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
        plot.title = element_text(hjust = 1, size = 12))+ 
  geom_ribbon(aes(ymin=a, ymax=b),fill = "#6a6a6a")+
  geom_smooth(stat="identity",fill="grey",colour="white")+
  geom_line(aes(y=y))

p1



p1<- p1 + 
  xlab(paste("Distance")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p1

p1 <- p1 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p1

p1<- p1+ theme(text=element_text(size=12, 
                                 #       family="Comic Sans MS"))
                                 #       family="CM Roman"))
                                 #       family="TT Times New Roman"))
                                 #       family="Sans"))
                                 family="TT Times New Roman"))

p1


#############################################################################
p2<-ggplot(Fin_Fol_gg, aes(d,a)) +
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
  ylab(paste("Correlation")) +
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
p3<-ggplot(Fin_Rre_gg, aes(d,a)) +
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
p4 <-ggplot(Fin_Prh_gg, aes(d,a)) +
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



p4<- p4 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
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
p5 <-ggplot(Fol_Poc_gg, aes(d,a)) +
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



p5<- p5 + 
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
p6 <-ggplot(Poc_Rre_gg, aes(d,a)) +
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

p6

p6 <- p6 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p6

p6 <- p6 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p6

p6 <- p6+ theme(text=element_text(size=16, 
                                  #       family="Comic Sans MS"))
                                  #       family="CM Roman"))
                                  #       family="TT Times New Roman"))
                                  #       family="Sans"))
                                  family="TT Times New Roman"))

p6

#############################################################################
p7 <-ggplot(Poc_Prh_gg, aes(d,a)) +
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

p7 <- p7 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p7

p7 <- p7 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p7

p7 <- p7+ theme(text=element_text(size=16, 
                                  #       family="Comic Sans MS"))
                                  #       family="CM Roman"))
                                  #       family="TT Times New Roman"))
                                  #       family="Sans"))
                                  family="TT Times New Roman"))

p7

#############################################################################

p8 <-ggplot(Fol_Rre_gg, aes(d,a)) +
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

p8 <- p8 + 
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p8

p8 <- p8 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p8

p8 <- p8+ theme(text=element_text(size=16, 
                                  #       family="Comic Sans MS"))
                                  #       family="CM Roman"))
                                  #       family="TT Times New Roman"))
                                  #       family="Sans"))
                                  family="TT Times New Roman"))

p8

#############################################################################
p9 <-ggplot(Fol_Prh_gg, aes(d,a)) +
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
  xlab(paste(" ")) +
  theme(axis.title.x=element_text(angle = 0, size = 15)) + #, face = "bold"
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p9

p9 <- p9 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p9

p9 <- p9+ theme(text=element_text(size=16, 
                                  #       family="Comic Sans MS"))
                                  #       family="CM Roman"))
                                  #       family="TT Times New Roman"))
                                  #       family="Sans"))
                                  family="TT Times New Roman"))

p9

#############################################################################
p10 <-ggplot(Rre_Prh_gg, aes(d,a)) +
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
  ylab(paste("Correlation")) +
  theme(axis.title.y=element_text(angle = 90, size = 15)) #, face = "bold"

p10

p10 <- p10 + 
  geom_hline(yintercept=0)
#labs(title = expression(paste(italic("Nome da esp?cie"))))
p10

p10 <- p10+ theme(text=element_text(size=16, 
                                    #       family="Comic Sans MS"))
                                    #       family="CM Roman"))
                                    #       family="TT Times New Roman"))
                                    #       family="Sans"))
                                    family="TT Times New Roman"))

p10


#############################################################################

tiff(filename="Bivariate.tiff", res=300, height=900/72*500, width=100/72*800, compression= "lzw")
# simple grid with labels and aligned plots
plot_grid(p1, p2, p3, p4, p5, 
          p6, p7, p8, p9, p10,
          labels = c("A", "B", "C", "D", "E",
                     "F", "G", "H", "I", "J"),
          ncol = 1)

dev.off()
#############################################################################
#############################################################################
#############################################################################
#############################################################################
# 3.0 Bubble plot for bivariate spatial correlograms of the five most abundant species

b1 <- ggplot(Fin_Poc, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  geom_count() +
  scale_size_area()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b1


angulo=-0.65
xnew_Fin_Poc <- Fin_Poc[,"x"]*cos(angulo) + Fin_Poc[,"y"]*sin(angulo)
ynew_Fin_Poc <- -Fin_Poc[,"x"]*sin(angulo) + Fin_Poc[,"y"]*cos(angulo)

b1 <- ggplot(Fin_Poc, aes(x=xnew_Fin_Poc, y=ynew_Fin_Poc, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9)+
  geom_count() +
  scale_size_area()+
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b1

b1 <- b1 + coord_fixed(ratio = -0.81)

b1


b1 <- b1 +labs(title="", x=" ", y = "Longitude")

b1

#############################################################################

b2 <- ggplot(Fin_Fol, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b2

angulo=-0.65
xnew_Fin_Fol <- Fin_Fol[,"x"]*cos(angulo) + Fin_Fol[,"y"]*sin(angulo)
ynew_Fin_Fol <- -Fin_Fol[,"x"]*sin(angulo) + Fin_Fol[,"y"]*cos(angulo)

b2 <- ggplot(Fin_Fol, aes(x=xnew_Fin_Fol, y=ynew_Fin_Fol, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b2

b2 <- b2 + coord_fixed(ratio = -0.81)

b2

b2 <- b2 +labs(title="", x=" ", y = "Longitude")

b2

b2 <- b2 + guides(color = guide_legend(order=2),
                  size = guide_legend(order=1))
#shape = guide_legend(order=3))

b2

###########################################################################

b3 <- ggplot(Fin_Rre, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b3


angulo=-0.65
xnew_Fin_Rre <- Fin_Rre[,"x"]*cos(angulo) + Fin_Rre[,"y"]*sin(angulo)
ynew_Fin_Rre <- -Fin_Rre[,"x"]*sin(angulo) + Fin_Rre[,"y"]*cos(angulo)

b3 <- ggplot(Fin_Rre, aes(x=xnew_Fin_Rre, y=ynew_Fin_Rre, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b3

b3 <- b3 + coord_fixed(ratio = -0.81)

b3

b3 <- b3 +labs(title="", x=" ", y = "Longitude")

b3

###########################################################################

b4 <- ggplot(Fin_Prh, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b4

#Alterando o ?ngulo
angulo=-0.65
xnew_Fin_Prh <- Fin_Prh[,"x"]*cos(angulo) + Fin_Prh[,"y"]*sin(angulo)
ynew_Fin_Prh <- -Fin_Prh[,"x"]*sin(angulo) + Fin_Prh[,"y"]*cos(angulo)

b4 <- ggplot(Fin_Prh, aes(x=xnew_Fin_Prh, y=ynew_Fin_Prh, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b4

b4 <- b4 + coord_fixed(ratio = -0.81)

b4

b4 <- b4 +labs(title="", x=" ", y = "Longitude")

b4
###########################################################################

b5 <- ggplot(Poc_Fol, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b5


angulo=-0.65
xnew_Poc_Fol <- Poc_Fol[,"x"]*cos(angulo) + Poc_Fol[,"y"]*sin(angulo)
ynew_Poc_Fol <- -Poc_Fol[,"x"]*sin(angulo) + Poc_Fol[,"y"]*cos(angulo)

b5 <- ggplot(Poc_Fol, aes(x=xnew_Poc_Fol, y=ynew_Poc_Fol, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b5

b5 <- b5 + coord_fixed(ratio = -0.81)

b5

b5 <- b5 +labs(title="", x=" ", y = "Longitude")

b5

###########################################################################

b6 <- ggplot(Poc_Rre, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b6

angulo=-0.65
xnew_Poc_Rre <- Poc_Rre[,"x"]*cos(angulo) + Poc_Rre[,"y"]*sin(angulo)
ynew_Poc_Rre <- -Poc_Rre[,"x"]*sin(angulo) + Poc_Rre[,"y"]*cos(angulo)

b6 <- ggplot(Poc_Rre, aes(x=xnew_Poc_Rre, y=ynew_Poc_Rre, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b6

b6 <- b6 + coord_fixed(ratio = -0.81)

b6

b6 <- b6 +labs(title="", x=" ", y = "Longitude")

b6

###########################################################################

b7 <- ggplot(Poc_Prh, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  geom_count() +
  scale_size_area()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b7


angulo=-0.65
xnew_Poc_Prh <- Poc_Prh[,"x"]*cos(angulo) + Poc_Prh[,"y"]*sin(angulo)
ynew_Poc_Prh <- -Poc_Prh[,"x"]*sin(angulo) + Poc_Prh[,"y"]*cos(angulo)

b7 <- ggplot(Poc_Prh, aes(x=xnew_Poc_Prh, y=ynew_Poc_Prh, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b7

b7 <- b7 + coord_fixed(ratio = -0.81)

b7

b7 <- b7 +labs(title="", x=" ", y = "Longitude")

b7

###########################################################################

b8 <- ggplot(Fol_Rre, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b8

angulo=-0.65
xnew_Fol_Rre <- Fol_Rre[,"x"]*cos(angulo) + Fol_Rre[,"y"]*sin(angulo)
ynew_Fol_Rre <- -Fol_Rre[,"x"]*sin(angulo) + Fol_Rre[,"y"]*cos(angulo)

b8 <- ggplot(Fol_Rre, aes(x=xnew_Fol_Rre, y=ynew_Fol_Rre, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b8

b8 <- b8 + coord_fixed(ratio = -0.81)

b8

b8 <- b8 +labs(title="", x=" ", y = "Longitude")

b8

###########################################################################


b9 <- ggplot(Fol_Prh, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b9

#Alterando o ?ngulo
angulo=-0.65
xnew_Fol_Prh <- Fol_Prh[,"x"]*cos(angulo) + Fol_Prh[,"y"]*sin(angulo)
ynew_Fol_Prh <- -Fol_Prh[,"x"]*sin(angulo) + Fol_Prh[,"y"]*cos(angulo)

b9 <- ggplot(Fol_Prh, aes(x=xnew_Fol_Prh, y=ynew_Fol_Prh, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b9

b9 <- b9 + coord_fixed(ratio = -0.81)

b9

b9 <- b9 +labs(title="", x=" ", y = "Longitude")

b9

b9 <- b9 + guides(color = guide_legend(order=2),
                  size = guide_legend(order=1))
#shape = guide_legend(order=3))

b9
###########################################################################

b10 <- ggplot(Rre_Prh, aes(x=x, y=y, size = Abundance, colour = Species)) +
  scale_fill_viridis()+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b10


angulo=-0.65
xnew_Rre_Prh <- Rre_Prh[,"x"]*cos(angulo) + Rre_Prh[,"y"]*sin(angulo)
ynew_Rre_Prh <- -Rre_Prh[,"x"]*sin(angulo) + Rre_Prh[,"y"]*cos(angulo)

b10 <- ggplot(Rre_Prh, aes(x=xnew_Rre_Prh, y=ynew_Rre_Prh, size = Abundance, colour = Species)) +
  scale_x_continuous(limits=c(-20.298, -20.292))+
  scale_y_continuous(limits=c(-41.3245, -41.3230))+
  scale_colour_manual(values = c("Black", "Grey"))+
  geom_point(alpha=0.9) +
  #ggtitle("Faramea involucellata X Faramea campanella") + 
  theme_bw() +
  theme(text=element_text(family="Times", face="plain", size=20))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

b10

b10 <- b10 + coord_fixed(ratio = -0.81)

b10

b10 <- b10 +labs(title="", x="Latitude", y = "Longitude")

b10

###########################################################################
###########################################################################

tiff(filename="Bivariado_5sp.tiff", res=600, height=1800/72*500, width=400/72*800, compression= "lzw")
# simple grid with labels and aligned plots
plot_grid(b1,b2,b3,b4,b5,
          b6,b7,b8,b9,b10,
          ncol = 1)

dev.off()