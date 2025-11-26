

######################################################
#---------------------------------------------# map

#regs1 <- read.csv("data/CoralSea/regions.csv")
#rcols1 <- regs1$col
#names(rcols1) <- regs1$name
#r.ord <- regs1$name[order(regs1$ord)]
#rownames(c1) <- c1$id

library("sf")
library("rworldmap")

sites <- unique(df1[,c("Data", "SectorRegion", "Reef", "lat", "long")])

# add colours
#rcols <- regs$col
#names(rcols) <- regs$name
# order regions
#r.ord <- regs$name[order(regs$ord)]
#sites$region <- factor(sites$region, levels=r.ord)

# world
rmap <- getMap(resolution="high") # SpatialPolygonsDataFrame
rmap.sf <- st_as_sf(rmap)

# gbr 
load("data/GBR_mapdata/reefs.RData")
summary(reefs)
gbr.sf<-st_as_sf(reefs)
head(gbr.sf)

# qld
load("data/GBR_mapdata/qld.RData")
summary(qld)
qld.sf<-st_as_sf(qld)
head(qld.sf)

length(seq(from=-85, to=85, by=10))
lat<-rep(seq(from=-85, to=85, by=10), each=18)
long<-rep(seq(from=-170, to=170, by=20), 18)
blocks<-data.frame(lat=lat, long=long, val=1)
xlim = c(142, 157)
ylim = c(-24, -10)

towns <- data.frame(name=c("Cooktown", "Cairns","Townsville", "Mackay", "Gladstone"), lat=c(-15.472405, -16.921848, -19.277811,-21.149498, -23.859962), long=c(145.250033, 145.757270, 146.769968,149.186270, 151.256776))
#towns <- st_as_sf(towns, coords = c("long", "lat")) 
#st_crs(towns) <- st_crs(qld)

map <- ggplot()+
#geom_tile(data = blocks, aes(x = long, y = lat), alpha = 0.8, fill="lightcyan2") + 
geom_sf(data=rmap.sf[!rmap.sf$NAME =="Australia",], lwd=0.01,fill="lemonchiffon", col="black") + #lemonchiffon
geom_sf(data=gbr.sf,  lwd=0.01, col="grey60", fill=NA)+ # azure
geom_sf(data=qld.sf, col="grey50", fill="lemonchiffon",lwd=0.01)+
geom_point(data=towns, aes(long, lat), shape=15, size=0.5)+
geom_text(data=towns, aes(long-0.27, lat, label=name), hjust=1, size=1.5)+
geom_point(data=sites, aes(long, lat, fill=Data), shape=21, size=1, stroke=0.1)+
coord_sf(ylim=c(ylim[1], ylim[2]), xlim=c(xlim[1], xlim[2]))+
scale_x_continuous(breaks=c(144, 148, 152, 156))+
#scale_fill_manual(values=rcols)+
scale_fill_manual(values=c("black", "red"))+
theme_minimal()+theme(legend.margin=margin(2,3,2,-1), axis.title=element_blank(), legend.title=element_blank(), 
legend.position=c(0.8, 0.92), legend.key.height=unit(1, "mm"), plot.title=element_text(size=7, hjust=0.5, face="bold"), legend.background=element_rect(colour="grey80", fill="grey97"),panel.grid.major=element_blank(), axis.text=element_text(size=6), legend.text=element_text(size=6, margin = margin(1,1,1,-3)), panel.background=element_rect(color="grey"))+ggtitle("Reef locations")
map 

######################################################
#---------------------------------------------# OVERALL PCA

# direction
#df1$PC1 <-   df1$PC1 
#vecs$PC1 <-   vecs$PC1
#df1$PC2 <-   - df1$PC2
#vecs$PC2 <-  - vecs$PC2

mag1 <- diff(range(pca1a$x[,"PC1"])) /2
mag2 <- diff(range(pca1a$x[,"PC2"])) /2

dens1 <- ggplot(df1, aes(PC1_1a))+geom_density(aes(fill=Data, col=Data), alpha=0.2, size=0.2)+scale_fill_manual(values=c("black", "red"))+scale_colour_manual(values=c("black", "red"))+guides(col="none", fill="none")

dens2 <- ggplot(df1, aes(PC2_1a))+geom_density(aes(fill=Data, col=Data), alpha=0.2, size=0.2)+scale_fill_manual(values=c("black", "red"))+scale_colour_manual(values=c("black", "red"))+guides(col="none", fill="none")

hulls <- NULL
for (d in unique(df1$Data)){
  temp <- na.omit(df1[df1$Data==d,c("PC1_1a","PC2_1a")])
  h <- temp[chull(temp$PC1_1a, temp$PC2_1a),]
  #vol <- convhulln(df,"FA")$vol / convhulln(axes,"FA")$vol *100 
  hulls <-rbind(hulls, cbind(h, Data=d))  }
head(hulls)

ppB1 <- ggplot()+
#stat_ellipse(data=df1, geom="polygon", aes(PC1_1a, PC2_1a, fill=Data), alpha=0.1)+
#geom_polygon(data=hulls, aes(PC1_1a, PC2_1a, colour=Data), fill=NA, size=0.1)+
geom_point(data=df1, aes(PC1_1a, PC2_1a, fill=Data, size=coral), shape=21, stroke=0.01)+
#scale_fill_manual(values=rcols1)+
guides(size="none")+
#stat_density_2d(data=dhw_sub, aes(PC2, PC3, col=type2, fill=type2), breaks=c(0.01), geom="polygon", alpha=0.5)+
labs(x=paste("PCA 1 (", vars1a[1], "%)", sep=""),y=paste("PCA 2 (", vars1a[2], "%)", sep=""))+
scale_radius(range=c(0.35,0.35))+
#scale_fill_viridis()+
scale_fill_manual(values=c("black", "red"))+scale_colour_manual(values=c("black", "red"))+guides(col="none")+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), legend.key.height=unit(3, "mm"), axis.text=element_text(size=7), axis.title=element_text(size=7), plot.margin=margin(0,0,5,5))
ppB1


dhw_plot <- plot_grid(
# row1
dens1+theme_void(),
ggplot()+theme_void(),
ggplot()+theme_void(),
# row2
ggplot()+theme_void(),
ggplot()+theme_void(),
ggplot()+theme_void(),
# row3
ppB1+guides(col="none", fill="none"),
ggplot()+theme_void(),
dens2+theme_void()+coord_flip(),
ncol=3, rel_heights=c(0.25,-0.16, 1), rel_widths=c(1, -0.21, 0.32), align="hv", axis="tblr")
dhw_plot


######################################################
#---------------------------------------------# PCA vectors!

vecs1a$PC2b <- vecs1a$PC2
vecs1a$PC2b[vecs1a$lab2=="complexity (1-5)"] <- vecs1a$PC2b[vecs1a$lab2=="complexity (1-5)"] -0.01
vecs1a$PC2b[vecs1a$lab2=="Herbivore %"] <- vecs1a$PC2b[vecs1a$lab2=="Herbivore %"] -0.025
vecs1a$PC2b[vecs1a$lab2=="coral %"] <- vecs1a$PC2b[vecs1a$lab2=="coral %"] +0.025
vecs1a$PC2b[vecs1a$lab2=="fish Biomass"] <- vecs1a$PC2b[vecs1a$lab2=="fish Biomass"] -0.02
vecs1a$PC2b[vecs1a$lab2=="coral Comp.2"] <- vecs1a$PC2b[vecs1a$lab2=="coral Comp.2"] -0.02

vecs1a$PC1b <- vecs1a$PC1
vecs1a$PC1b[vecs1a$lab2=="complexity (1-5)"] <- vecs1a$PC1b[vecs1a$lab2=="complexity (1-5)"] -0.25
vecs1a$PC1b[vecs1a$lab2=="fish Biomass"] <- vecs1a$PC1b[vecs1a$lab2=="fish Biomass"] -0.02
vecs1a$PC1b[vecs1a$lab2=="Morph. div (D)"] <- vecs1a$PC1b[vecs1a$lab2=="Morph. div (D)"] -0.23
vecs1a$PC1b[vecs1a$lab2=="coral Comp.2"] <- vecs1a$PC1b[vecs1a$lab2=="coral Comp.2"] -0.02

vecs1a$lab3<-vecs1a$lab2
vecs1a$lab3[vecs1a$lab3=="coral Comp.2"] <- "Comp.2"
vecs1a$lab3[vecs1a$lab3=="fish Biomass"] <- "F. Biomass"

ppA1 <- ggplot()+
geom_segment(data=vecs1a, aes(x=0, xend=PC1*mag1*0.9, y=0, yend=PC2*mag2*0.9, col=max), arrow=arrow(length=unit(1,"mm")))+
geom_text(data=vecs1a, aes(x=PC1b*mag1, y=PC2b*mag2, label=lab3), hjust=ifelse(vecs1a$PC1>0, 0, 1), size=2.3, fontface="italic")+
scale_colour_manual(values=c("slategrey","grey", "#d8b365","#998ec3", "hotpink4","lightblue3" ))+
guides(col="none")+
theme_void()+coord_cartesian(xlim=c(-3,3.5), ylim=c(-2.5, 3))
ppA1


######################################################
#---------------------------------------------# PCA loadings!


head(vecslong)
vecslong$var <- ifelse(vecslong$variable=="PC1", paste("PC1 (", vars1a[1], "%)", sep=""), 
ifelse(vecslong$variable=="PC2", paste("PC2 (", vars1a[2], "%)", sep=""), 
ifelse(vecslong$variable=="PC3", paste("PC3 (", vars1a[3], "%)", sep=""), 
ifelse(vecslong$variable=="PC4", paste("PC4 (", vars1a[4], "%)", sep=""),
ifelse(vecslong$variable=="PC5", paste("PC5 (", vars1a[5], "%)", sep=""),ifelse(vecslong$variable=="PC6", paste("PC6 (", vars1a[6], "%)", sep=""), NA))))))


labels <- data.frame(dim = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), lab=c("1. coral cover", "2. fish composition", "3. fish biomass", "4. algal composition", "5. coral composition", "6. complexity"))
vecslong$lab <- labels$lab[match(vecslong$variable, labels$dim)]


ppC1b <-ggplot()+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=value, fill=variable), stat="identity", width=0.7,  alpha=0.5)+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=selectval, fill=variable), stat="identity", width=0.7, position="stack")+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=maxval), stat="identity", fill=NA, col="black", width=0.7, size=0.2)+
geom_text(data=vecslong[vecslong$value>0,], aes(y=reorder(lab2, -ord), x=value+0.05, fill=variable, label=six), vjust=0.75)+
geom_text(data=vecslong[vecslong$value<0,], aes(y=reorder(lab2, -ord), x=value-0.05, fill=variable, label=six), vjust=0.75)+
geom_vline(xintercept=0, size=0.2)+
theme_classic()+labs(y="")+
scale_fill_manual(values=c("slategrey","grey", "#d8b365","#998ec3", "hotpink4","lightblue3" ))+
xlab("PCA axes loadings")+
#xlim(c(-0.65, 0.65))+
scale_x_continuous(breaks=c(-0.4, 0, 0.4))+
ggtitle("")+
facet_wrap(~lab+var, nrow=1)+
theme(axis.text.y=element_text(size=7), strip.text=element_text(size=7), legend.title=element_blank(), legend.position=c(1.1, 0.85), legend.text=element_text(size=7), legend.background=element_blank(), legend.key.size=unit(3, "mm"), axis.title=element_text(size=8), axis.line.y=element_blank(), strip.background=element_blank(), axis.text.x=element_text(size=6), plot.background=element_blank())
ppC1b

######################################################
#---------------------------------------------# correlation plot


colors = c("blue", "white", "red")

pcor2 <- pcor
pcor2[upper.tri(pcor2, diag=T)]<-NA

pcor3 <- melt(pcor2)
pdf <- pcor3[!is.na(pcor3$value),]
pdf$lab.x <- labs$lab2[match(pdf$Var1, labs$lab)]
pdf$lab.y <- labs$lab2[match(pdf$Var2, labs$lab)]
pdf$lab.x <- factor(pdf$lab.x, levels=c(ords$lab2))
pdf$lab.y <- factor(pdf$lab.y, levels=c(ords$lab2))
pdf$group2 <- labs$group[match(pdf$Var2, labs$lab)]
pdf$group1 <- labs$group[match(pdf$Var1, labs$lab)]
pdf$group1 <- factor(pdf$group1, levels=rev(c("complexity", "coral", "fish", "algae")))
pdf$group2 <- factor(pdf$group2, levels=c("complexity", "coral", "fish", "algae"))


tab.s<-readPNG("data/sils/tabular.png")
tab.s<-rasterGrob(tab.s, interpolate=TRUE)
fish.s<-readPNG("data/sils/fish2.png")
fish.s<-rasterGrob(fish.s, interpolate=TRUE)
comp.s<-readPNG("data/sils/comp2.png")
comp.s<-rasterGrob(comp.s, interpolate=TRUE)
alg.s<-readPNG("data/sils/algae3.png")
alg.s<-rasterGrob(alg.s, interpolate=TRUE)


pcorplot<-ggplot(pdf, aes(lab.x, lab.y))+
geom_tile(aes(fill=value), col="black", size=0.1)+
#geom_point(aes(fill=value), shape=21, col="black", size=7)+
scale_y_discrete(position = "right")+
#scale_fill_viridis()+
geom_text(aes(label=round(value, 2)), size=2, fontface="bold")+
guides(fill="none")+
#scale_fill_distiller(palette="BrBG", limits=c(-1, 1), direction=1, breaks=c(-1, 0,1))+
#scale_fill_distiller(palette="Spectral", limits=c(-1, 1), direction=1, breaks=c(-1, 0,1))+
scale_fill_gradient2(low = "red", high="blue", mid="white", midpoint=0, limit = c(-1, 1) )+
facet_grid(group1~group2, scales="free", space="free")+
theme_classic()+
#annotation_custom(tab.s, xmin=0.1, xmax=0.5, ymin=0.1, ymax=0.5)+
theme(axis.title=element_blank(), axis.line=element_blank(), 
axis.text=element_text(size=8),
legend.title=element_blank(), legend.position=c(0.1, 0.8), legend.background=element_blank(), legend.text=element_text(size=7), legend.key.width=unit(1, "mm"), legend.key.height=unit(3, "mm"), plot.background=element_blank(), panel.background=element_blank(), axis.text.y=element_text(size=7), axis.text.x=element_text(size=7, angle=90, hjust=0, vjust=0), strip.text=element_blank(), strip.background=element_blank(), panel.spacing = unit(0.5, "mm"))+coord_flip()
pcorplot




######################################################
#---------------------------------------------# fig 1


pcorplot2 <- plot_grid(NULL, plot_grid(pcorplot, NULL, nrow=1, rel_widths=c(1,-0.2), labels=c("c", ""), hjust=-5, label_size=9), NULL, ncol=1, rel_heights=c(0.1, 1, 0.25))

ppA12<- plot_grid(NULL, plot_grid(NULL, ppA1+guides(col="none"),NULL, ncol=1, rel_heights=c(0.8,1, -0.3), labels=c("", "d", ""), label_size=9, hjust=-15, vjust=7),NULL, nrow=1, rel_widths=c(-0.5, 1, -0.05))


fig1X2<- plot_grid(
plot_grid(
plot_grid(map, NULL, dhw_plot, ncol=1, labels=c("a","", "b"), label_size=9, rel_heights=c(1,0.2,1)),
plot_grid(pcorplot2, ppA12, rel_widths=c(2, 1.2)),
nrow=1, rel_widths=c(1,2)),
NULL,
ppC1b, ncol=1, rel_heights=c(2,0.05, 1), labels=c("", "", "e"), label_size=9, vjust=5, hjust=-5)+
draw_text("Metric PCA", 0.18, 0.68, size=7, fontface="bold")+
draw_text("Metric correlations", 0.62, 0.97, size=7, fontface="bold")+
draw_text("PCA dimensions", 0.75, 0.62, size=7, fontface="bold")+
draw_grob(tab.s, x=0.55, y=0.52, height=0.04, width=0.04)+
draw_grob(fish.s, x=0.7, y=0.66, height=0.04, width=0.04)+
draw_grob(comp.s, x=0.44, y=0.43, height=0.05, width=0.05)+
draw_grob(alg.s, x=0.81, y=0.765, height=0.03, width=0.05)
#draw_line(x=c(0.2, 0.31), y=c(0.11, 0.22))
#draw_line(x=c(0.2, 0.31), y=c(0.11, 0.22))
fig1X2

# ggsave("figs/Fig1.jpg", fig1X2, height=7, width=6.9)






