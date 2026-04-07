
labs <- read.csv("data/info/pcaLabs.csv")

s1coral$Reef <- s1t$Reef[match(s1coral$id, s1t$id)]
s1coral$Zone <- s1t$Zone[match(s1coral$id, s1t$id)]

siteS <- aggregate(p~size+Reef+Zone, s1coral, mean)
siteS$gp <- paste(siteS$Reef, siteS$Zone)

sizeplot <- ggplot(siteS, aes(size, p))+
geom_rect(data=NULL, aes(xmin=1, xmax=2, ymin=-Inf, ymax=Inf), fill="grey")+
geom_rect(data=NULL, aes(xmin=5.7, xmax=6.3, ymin=-Inf, ymax=Inf), fill="grey")+
geom_text(data=NULL, aes(x=1.5, y=0.9, label="small"), size=3)+
geom_text(data=NULL, aes(x=6, y=0.9, label="large"), size=3)+
geom_segment(data=NULL, aes(x=6.35, xend=6.5, y=0.9, yend=0.9), arrow=arrow(length=unit(1, "mm")))+
#scale_y_sqrt()+
geom_line(aes(group=gp), size=0.1)+ #position="stack",
theme_classic()+
labs(x="Diameter (cm)", y="Proportion of corals")+
theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8))
sizeplot 



# SIZE PCA

c1$log_fN <- log(c1$f.N)
c1$log_biomass <- log(c1$f.biomass)
c1$P.herb <- c1$N.herb / c1$f.N
c1$P.carn <- c1$N.carn / c1$f.N
c1$covsqrt <- sqrt(c1$coral)
#df1$alg_ratio <- df1$turf / (df1$cca + df1$turf)
c1$alg_ratio <- c1$cca / (c1$cca + c1$turf + c1$algae)
head(c1)

axes1 <- c("coral", "Complexity", "simpD", "log_biomass","algae", "turf", "f.richness", "P.herb", "P.carn", "alg_ratio", "FsimpD", "comp1", "comp2", "acro", "p10", "p60")
axes2 <- c("covsqrt", "Complexity", "simpD", "log_biomass", "alg_ratio", "FsimpD", "p10", "p60")
axes.list <- list(axes1, axes2)

pca <- prcomp(na.omit(c1[,axes1]), scale=T, center=T)
#biplot(pca)	
vecs <- data.frame(pca$rotation[,c("PC1", "PC2", "PC3", "PC4")], lab=rownames(pca$rotation))
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
#biplot(pca)
vecs$max <- colnames(vecs[,c("PC1", "PC2", "PC3", "PC4")])[max.col(abs(vecs[,c("PC1", "PC2", "PC3", "PC4")]), ties.method='first')]	
coord <- data.frame(pca$x)

vecs$lab2 <- labs$lab2[match(vecs$lab, labs$lab)]
vecs$size <- ifelse(vecs$lab %in% c("p10", "p60"), "size", "no")
vecs$PC2b <- vecs$PC2
vecs$PC2b[vecs$lab=="p60"] <- vecs$PC2b[vecs$lab=="p60"] +0.05
vecs$PC2b[vecs$lab=="Complexity"] <- vecs$PC2b[vecs$lab=="Complexity"] -0.025
vecs$PC2b[vecs$lab=="acro"] <- vecs$PC2b[vecs$lab=="acro"] -0.05
vecs$PC2b[vecs$lab=="coral"] <- vecs$PC2b[vecs$lab=="coral"] -0.015
vecs$PC2b[vecs$lab=="P.carn"] <- vecs$PC2b[vecs$lab=="P.carn"] -0.015


mag1 <- 6
mag2 <- 7
s.pca <- ggplot()+
geom_point(data=coord, aes(PC1, PC2),col="grey", size=0.2, shape=21, stroke=0.25)+
geom_segment(data=vecs, aes(x=0, xend=PC1*mag1, y=0, yend=PC2*mag1, col=size))+
geom_text(data=vecs, aes(PC1*mag2, PC2b*mag2, label=lab2, col=size), size=1.9, hjust=ifelse(vecs$PC1>0, 0, 1),fontface="bold", show_guide = FALSE)+ #,, direction="y", force=0.25, segment.size=0.2
scale_colour_manual(values=c("black","red"))+ ##998ec3
guides(col="none")+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), strip.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=7), legend.key.height=unit(1, "mm"), panel.margin=unit(1, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=6), strip.text=element_text(size=8))
s.pca


figS5 <- plot_grid(sizeplot, s.pca , labels=c("a","b"), label_size=9)
figS5



