

rm(list = ls())

library("ggplot2")
library("cowplot")
library("viridis")
library("reshape2")

# Calculate FD for coral sea/LTMP


######################################################
#--------------------------------------------- traits


df <- read.csv("data/traitbiogeography.csv")
colnames(df)
 df <-     df[!is.na(df$Tropical.Australia),]
df2 <- df[,c("species", "genus", "raw_growth_form", "domain", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter", "cat_skeletaldensity", "cat_colonyheight", "cat_SA_vol", "cat_spacesize", "dat_corallite", "Tropical.Australia", "reproductive_mode","raw_colonydiam","raw_growth","dat_growth","dat_colonydiameter" )]
head(df2)
nrow(df2)

cstraits <- read.csv("data/CoralSea/traits_to_CS.csv")
df2$link <- cstraits$match[match(df2$species, cstraits$species)]

#mistakes
df2$cat_growthrate[df2$link=="Acropora...Tabular"] <- 5
df2$cat_SA_vol[df2$link=="Acropora...Tabular"] <- 5
df2$cat_spacesize[df2$link=="Acropora...Tabular"] <- 5

df2$cat_SA_vol[df2$raw_growth_form=="laminar"]


# Quantify morphology
# Zawada paper: 
morph<-read.csv("data/3DLaserScannedColonies_Whole_MedtoHighQualityMeshes_Zawada_190511.csv")
head(morph)
unique(morph$GrowthForm)
morph$morph<-ifelse(morph$GrowthForm=="Corymbose","corymbose",  ifelse(morph$GrowthForm=="Massive","massive",
ifelse(morph$GrowthForm=="Tabular","tables_or_plates",
ifelse(morph$GrowthForm=="Digitate","digitate",
ifelse(morph$GrowthForm=="Submassive","massive",
ifelse(morph$GrowthForm=="Laminar","laminar",
ifelse(morph$GrowthForm=="Arborescent","branching_open",
ifelse(morph$GrowthForm=="BranchingClosed","branching_closed",
ifelse(morph$GrowthForm=="BranchingThick","columnar" ,
ifelse(morph$GrowthForm=="Encrusting","encrusting", NA))))))))))
unique(df2$raw_growth_form)
#ggplot(morph, aes(x=reorder(morph, -Vvol), y=Vvol))+geom_boxplot()+geom_point()+theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())

df2$Growth.form.typical <- df2$raw_growth_form
df2$Growth.form.typical <- ifelse(df2$Growth.form.typical %in% c( "encrusting_long_uprights"), "encrusting", ifelse(df2$Growth.form.typical=="bifacial", "laminar", ifelse(df2$Growth.form.typical=="hispidose", "branching_open", ifelse(df2$Growth.form.typical %in% c("submassive", "solitary_attached", "solitary_free"), "massive",   as.character(df2$Growth.form.typical )))))
unique(df2$Growth.form.typical)

morph$sa_vol<-  morph$ColonySA / morph$ColonyVol

s.area<-aggregate(ColonySA~morph, morph, mean)
vol<-aggregate(ColonyVol~morph, morph, mean)
cvol <- aggregate(ConvexVol~morph, morph, mean)

df2$sa<-s.area$ColonySA[match(df2$Growth.form.typical, sa_vol$morph)]
df2$vol<-vol$ColonyVol[match(df2$Growth.form.typical, sa_vol$morph)]
df2$cvol<-cvol$ConvexVol[match(df2$Growth.form.typical, sa_vol$morph)]
df2$sa_vol <- df2$sa/df2$vol
df2$spaces <- df2$cvol - df2$vol 

ggplot(df2, aes(spaces, cat_spacesize))+geom_point()+geom_smooth()
ggplot(df2, aes(sa_vol, cat_SA_vol))+geom_point()+geom_smooth()


df2$sa_vol2 <- NA
df2$sa_vol2[df2$Growth.form.typical=="branching_open"] <- 2.025
df2$sa_vol2[df2$Growth.form.typical=="tables_or_plates"] <- 1.242
df2$sa_vol2[df2$Growth.form.typical=="branching_closed"] <- 2.025
df2$sa_vol2[df2$Growth.form.typical=="digitate"] <- 2.2
df2$sa_vol2[df2$Growth.form.typical=="columnar"] <- 0.9
df2$sa_vol2[df2$Growth.form.typical=="massive"] <- 0.5
df2$sa_vol2[df2$Growth.form.typical=="corymbose"] <- 4.05
df2$sa_vol2[df2$Growth.form.typical=="encrusting"] <- 1
df2$sa_vol2[df2$Growth.form.typical=="laminar"] <- 2.176
#ggplot(df2, aes(sa_vol*10, sa_vol2))+geom_text(aes(label=Growth.form.typical))+geom_abline(slope=1)



######################################################
#--------------------------------------------- GLOBAL PCA

library("FD") 

hist(df2$dat_colonydiameter)

df2$log_growth <- log(df2$dat_growth) #log(df2$dat_growth)
df2$log_corallite <- log(df2$dat_corallite) # log(df2$dat_corallite)
df2$log_diameter <- log(df2$dat_colonydiameter) # log(df2$dat_colonydiameter)


list <- c("cat_growthrate","cat_skeletaldensity", "cat_corallitesize",
          "cat_colonydiameter","cat_colonyheight", "cat_SA_vol", "cat_spacesize")

dlist <- c("log_growth","cat_skeletaldensity", "log_corallite",
          "log_diameter","cat_colonyheight", "cat_SA_vol", "cat_spacesize")


rownames(df2) <- df2$species
pca <- prcomp(na.omit(df2[,dlist]), center=TRUE, scale.=TRUE)
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
vecs<-data.frame(varnames=rownames(pca$rotation), pca$rotation)
sum(vars[c(1:4)])

biplot(pca)

df2[,c("PC1", "PC2", "PC3", "PC4")] <- pca$x[match(rownames(df2), rownames(pca$x)),c("PC1", "PC2", "PC3", "PC4")]
head(df2)


gower<-gowdis(na.omit(df2[,c(dlist)])) # "reproductive_mode"
pco<-pcoa(gower)
df2[,c("PCo1", "PCo2", "PCo3", "PCo4")] <- pco$vectors[match(rownames(df2), rownames(pco$vectors)),c(1:4)]
head(df2)


## KERNEL DENSITY ESTIMATION ###
library("vegan")
library("ks")
library("contoureR")


df2$dim1 <- - df2$PCo1
df2$dim2 <-  df2$PCo2
df2$dim3 <- df2$PCo3
df2$dim4 <- df2$PCo4
mydata <- na.omit(df2[,c("dim1", "dim2")]) #, "dim3", "dim4"
#mydata <- df2[,list]

H <- Hpi(x=mydata)      # optimal bandwidth estimation
est<- kde(x=mydata, H=H, compute.cont=TRUE)     # kernel density estimation
cl<-contourLevels(est, prob=c(0.5, 0.05, 0.001), approx=TRUE) # set contour probabilities for drawing contour levels

coord <- data.frame(x=est$eval.points[[1]], y=est$eval.points[[2]])
est2 <- melt(est$estimate, value.name="z")
est2$x  <- coord$x[match(est2$Var1, rownames(coord))]
est2$y  <- coord$y[match(est2$Var2, rownames(coord))]
head(est2[,c("x", "y", "z")])
ggplot()+geom_tile(data=est2, aes(x, -y,fill=z))

lins <- getContourLines(est2[c(1:22000),c("x", "y", "z")], levels=cl)
propmatch <- data.frame(n=c(1:length(cl)), nom=names(cl), val=cl)
lins$nom <- propmatch$nom[match(lins$z, propmatch$val)]
head(lins)
dens <- getContourLines(est2[c(1:22000),c("x", "y", "z")], nlevel=100)
head(dens)

fit<-envfit(mydata, na.omit(df2[,dlist])) # use envfit for drawing arrows, can be also done using trait loadings
fit2 <- data.frame(fit$vectors$arrows)
fit2$lab <-  c("GR", "SD","CW", "D", "H", "SA", "IS")

fit2$dim1b <-   fit2$dim1*0.23 #3
fit2$dim2b <-  fit2$dim2*0.23 #3
fit2$dim1c <- fit2$dim1b *1.15
fit2$dim2c <- fit2$dim2b *1.15
fit2$dim2c[fit2$lab=="D"]<- fit2$dim2c[fit2$lab=="D"] +0.01
fit2$dim2c[fit2$lab=="IS"]<- fit2$dim2c[fit2$lab=="IS"] -0.01


tspace1 <- ggplot()+
geom_polygon(data=dens,aes(x, y,group=Group,fill=z), alpha=0.05)+
geom_path(data=lins,aes(x, y,group=Group, col=nom))+
#geom_point(data=df2, aes(-dim1, dim2), size=0.1, col="grey")+
geom_segment(data=fit2, aes(x=0, xend=dim1b, y=0, yend=dim2b), arrow=arrow(length=unit(2,"mm")))+ #*3
geom_segment(data=fit2, aes(x=0, xend=-dim1b, y=0, yend=-dim2b), linetype="dashed", size=0.3)+ #*3
geom_text(data=fit2, aes(dim1c, dim2c, label=lab), fontface="bold", size=3)+ # *3.5
scale_colour_manual(values=c("grey90", "grey80", "grey60"))+
scale_fill_distiller(palette="YlOrRd", direction=1)+
guides(fill="none")+labs(x="PCoA 1", y="PCoA 2")+
theme_bw()+theme(legend.title=element_blank(), legend.position=c(0.85, 0.9), legend.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(2, "mm"), legend.text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=8))
tspace1

tspace2 <- ggplot()+
geom_polygon(data=dens,aes(x, y,group=Group,fill=z), alpha=0.05)+
geom_path(data=lins,aes(x, y,group=Group, col=nom))+
#geom_point(data=df2, aes(-dim1, dim2), size=0.1, col="grey")+
geom_segment(data=fit2, aes(x=0, xend=dim1b, y=0, yend=dim2b), arrow=arrow(length=unit(2,"mm")))+ #*3
geom_segment(data=fit2, aes(x=0, xend=-dim1b, y=0, yend=-dim2b), linetype="dashed", size=0.3)+ #*3
geom_text(data=fit2, aes(dim1c, dim2c, label=lab), fontface="bold", size=2)+ # *3.5
scale_colour_manual(values=c("grey90", "grey80", "grey60"))+
scale_fill_distiller(palette="YlOrRd", direction=1)+
ggtitle("1. create trait space")+
guides(fill="none", col="none")+labs(x="PCoA 1", y="PCoA 2")+
theme_bw()+theme(legend.title=element_blank(), legend.position=c(0.85, 0.9), legend.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(2, "mm"), legend.text=element_text(size=7), axis.text=element_blank(), axis.title=element_blank(), plot.title=element_text(size=7, hjust=0.5, face="bold"))
tspace2




################################## ----- how many clusters?

#https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters

df2$dim1 <- - df2$PCo1
df2$dim2 <- df2$PCo2
df2$dim3 <- df2$PCo3
df2$dim4 <- df2$PCo4
mydata <- na.omit(df2[, c("dim1", "dim2")]) # df2[,c("dim1", "dim2") #, "dim3", "dim4" #[!is.na(df2$Tropical.Australia)]
#mydata <- df2[,list]


nrow(mydata)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
wssdat <- data.frame(n=c(1:15), wss)                                       
wssplot <- ggplot(wssdat, aes(n, wss))+
geom_rect(data=NULL, aes(xmin=5, xmax=10, ymin=-Inf, ymax=Inf), fill="grey")+
geom_line()+
geom_point()+
xlab("N clusters")+ylab("Within groups\nsum of squares")+
theme_classic() +theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8))     
wssplot  
# changes every time?    
                    
mod <- mclustBIC(mydata) # pnas method
mod

nc <- 9

library(mclust)
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
d_clust <- Mclust(as.matrix(mydata), G=1:nc) # up to 20
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
 # plot(d_clust)
m.best

mdf <- data.frame(d_clust$classification)
head(mdf)
df2$mclust <- mdf$d_clust.classification[match(rownames(df2), rownames(mdf))]

library("factoextra")
fviz_nbclust(mydata, kmeans, method = "wss") +
      geom_vline(xintercept = 3, linetype = 2)+
      labs(subtitle = "Elbow method")

#fit <- cascadeKM(scale(mydata, center = TRUE,  scale = TRUE), 1, nc, iter = 1000)
#plot(fit, sortg = TRUE, grpmts.plot = TRUE)
#calinski.best <- as.numeric(which.max(fit$results[2,]))
#cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
# 10 clusters!

#library("NbClust")
#nb <- NbClust(mydata, diss=NULL, distance = "euclidean",
#       method = "kmeans", min.nc=5, max.nc=15, 
#        index = "alllong", alphaBeale = 0.1)
#hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))
# write.csv(data.frame(nb$Best.nc[1,]), "data/optimal.all.csv")
opt <- read.csv("data/optimal.all.csv")
optplot <- ggplot(opt, aes(opt))+
geom_segment(data=NULL, aes(x=9, xend=9, y=4.5, yend=6), col="red", arrow=arrow(length=unit(1, "mm"), ends="first"))+
geom_text(data=NULL, aes(x=9, y=6.5, label="BIC"), size=3)+
geom_histogram()+
labs(x="N clusters", y="N clustering\nindices")+
scale_y_continuous(expand=c(0,0))+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8))    
optplot


k <- kmeans(mydata, centers=nc, nstart=25, iter.max=10000000)
kdf <- data.frame(k$cluster)
head(kdf)
df2$kclust <- kdf$k.cluster[match(rownames(df2), rownames(kdf))]

clust <- melt(df2[,c("dim1", "dim2", "mclust", "kclust")], id.var=c("dim1", "dim2"))
#clust$value <- as.factor(clust$value)
clust$gp <- paste(clust$variable, clust$value)
head(clust)
table(clust$gp)

library("png")
library("grid")

tab.s<-readPNG("data/sils/tab.png")
tab.s<-rasterGrob(tab.s, interpolate=TRUE)
stag.s<-readPNG("data/sils/stag.png")
stag.s<-rasterGrob(stag.s, interpolate=TRUE)
cor.s<-readPNG("data/sils/cor.png")
cor.s<-rasterGrob(cor.s, interpolate=TRUE)
mas.s<-readPNG("data/sils/pori.png")
mas.s<-rasterGrob(mas.s, interpolate=TRUE)
col.s<-readPNG("data/sils/iso.png")
col.s<-rasterGrob(col.s, interpolate=TRUE)
br.s<-readPNG("data/sils/poc.png")
br.s<-rasterGrob(br.s, interpolate=TRUE)
lam.s<-readPNG("data/sils/lam.png")
lam.s<-rasterGrob(lam.s, interpolate=TRUE)
fun.s<-readPNG("data/sils/fung.png")
fun.s<-rasterGrob(fun.s, interpolate=TRUE)
enc.s<-readPNG("data/sils/encr.png")
enc.s<-rasterGrob(enc.s, interpolate=TRUE)



cplot <- ggplot() +
  geom_point(data=df2, aes(dim1, dim2), shape=21, size=0.25, stroke=0.25, col="black")+
  #stat_ellipse(data=clust, geom="polygon", aes(dim1, dim1, group=gp, col=variable, fill=variable), alpha=0.2)+
  stat_ellipse(data=df2, aes(dim1, dim2, group=factor(kclust)), col="red", fill="red", geom="polygon", alpha=0.2)+
   stat_ellipse(data=df2, aes(dim1, dim2, group=factor(mclust)), geom="polygon", fill="cadetblue",col="cadetblue", alpha=0.2, level=0.9)+
guides(fill="none")+labs(x="PCoA 1", y="PCoA 2")+
geom_text(data=NULL, aes(x=-0.4, y=-0.25, label="BIC"), col="cadetblue", size=3, hjust=0)+
geom_text(data=NULL, aes(x=-0.4, y=-0.28, label="K-means"), col="red", size=3, hjust=0)+
theme_bw()+theme(legend.title=element_blank(), legend.position=c(0.85, 0.9), legend.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(2, "mm"), legend.text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=8))+
annotation_custom(tab.s, xmin=0.23, xmax=0.33, ymin=0.12, ymax=0.22)+
annotation_custom(cor.s, xmin=0.2, xmax=0.3, ymin=0, ymax=0.1)+
annotation_custom(stag.s, xmin=0.27, xmax=0.35, ymin=-0.1, ymax=-0.2)+
annotation_custom(br.s, xmin=0.15, xmax=0.23, ymin=-0.1, ymax=-0.2)+
annotation_custom(mas.s, xmin=-0.23, xmax=-0.35, ymin=-0.15, ymax=-0.23)+
annotation_custom(col.s, xmin=-0.1, xmax=0, ymin=-0.2, ymax=-0.28)+
annotation_custom(lam.s, xmin=-0.05, xmax=-0.15, ymin=0.2, ymax=0.3)+
annotation_custom(fun.s, xmin=-0.35, xmax=-0.44, ymin=-0.05, ymax=0.05)+
annotation_custom(enc.s, xmin=-0.25, xmax=-0.34, ymin=0.15, ymax=0.25)
cplot


cplotX <- ggplot() +
  geom_point(data=df2, aes(dim1, dim2), shape=21, size=0.25, stroke=0.25, col="black")+
  stat_ellipse(data=df2, aes(dim1, dim2, group=factor(kclust)), col="red", fill="red", geom="polygon", alpha=0.2)+
  ggtitle("2. find clusters")+
guides(fill="none")+labs(x="PCoA 1", y="PCoA 2")+
theme_bw()+theme(legend.title=element_blank(), legend.position=c(0.85, 0.9), legend.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(2, "mm"), legend.text=element_text(size=7), axis.text=element_blank(), axis.title=element_blank(), plot.title=element_text(size=7, face="bold", hjust=0.5))
cplotX


fgroups <- plot_grid(tspace1,
plot_grid(wssplot, optplot, ncol=1, align="hv", labels=c("b", "c"), label_size=10),
cplot, nrow=1, rel_widths=c(1,0.65,1), labels=c("a", "", "d"), label_size=10)
fgroups 


################################## ----- functional groups

# ADRIA uses arborescent Acropora // Tabular Acropoa // Corymbose Acropora // Corymbose other // Small massives // Large massives

 colnames(data)
#check <- cbind(df2, k=k$cluster, df2[,c("dat_growth","dat_colonydiameter","raw_colonydiam","raw_growth_form" )])
  
unique(df2$raw_growth_form)
morphcols <- c("#addd8e", "#31a354", # greens
"#a6bddb", "#1c9099", # bluegreens
"#d7b5d8","#df65b0","#dd1c77","#980043", # purples
"slategrey",
"#fdbe85","#e6550d","black",#"#a63603",    table=#fd8d3c
"#a1dab4", "#41b6c4", "#225ea8") # blues
names(morphcols) <- c("submassive", "massive", 
"solitary_attached", "solitary_free",
"encrusting", "encrusting_long_uprights", "laminar","bifacial",
"columnar", 
"digitate", "corymbose",  "tables_or_plates",
"hispidose", "branching_closed", "branching_open")
morphcols

df2$clust <- df2$kclust

avs <- aggregate(dim1 ~ clust, df2, mean)
avs$dim2 <- aggregate(dim2 ~ clust, df2, mean)$dim2

df2$raw_growth_form <- factor(df2$raw_growth_form, levels=names(morphcols))

df2$gen <- ifelse(df2$genus == "Acropora", "Acropora",ifelse( df2$genus == "Porites", "Porites", ifelse( df2$genus == "Pocillopora", "Pocillopora", "other")))
 
cplot2 <- ggplot() +
  geom_point(data=df2, aes(dim1, dim2, fill=raw_growth_form),  shape=21, stroke=0.1, size=2)+
 # geom_point(data=df2[!df2$gen=="other",], aes(dim1, dim2, col=gen),  shape=21,  size=3)+
  #geom_label(data=avs, aes(PC1, PC2, label=k), size=5)+
  stat_ellipse(data=df2, aes(dim1, dim2,,group=factor(clust)))+
  geom_label(data=avs, aes(dim1, dim2, label=clust), fontface="bold", alpha=0.6, size=3)+
  scale_fill_manual(values=morphcols)+
  labs(x="PCoA 1", y="PCoA 2")+
  theme_bw()+theme(legend.title=element_blank(), legend.key.height=unit(4, "mm"))
cplot2

df2$gen <- ifelse(df2$genus == "Acropora", "Acropora",ifelse( df2$genus == "Porites", "Porites", ifelse( df2$genus == "Pocillopora", "Pocillopora", "other")))
 unique(df2$genus)
 
freqs <- data.frame(table(df2[,c("clust", "raw_growth_form")]))
fplot <- ggplot(freqs, aes(x=as.factor(clust), y=Freq, fill=raw_growth_form))+geom_bar(stat="identity")+scale_fill_manual(values=morphcols)+theme_classic()+scale_y_continuous(expand=c(0,0))+
xlab("Cluster")+ylab("N species")
fplot


freqs2 <- data.frame(table(df2[,c("clust", "gen")]))
fplot2 <- ggplot(freqs2, aes(x=as.factor(clust), y=Freq, fill=gen))+geom_bar(stat="identity", col="grey")+theme_classic()+
scale_fill_manual(values=c("black", "white","grey", "slategrey" ))+
scale_y_continuous(expand=c(0,0))+
xlab("Cluster")+ylab("N species")+theme(legend.title=element_blank(), legend.position=c(0.8, 0.8), legend.key.size=unit(1, "mm"))
fplot2



splot <- ggplot()+
geom_rect(data=NULL, aes(xmin=-Inf, xmax=Inf, ymin=500, ymax=3000), fill="red", alpha=0.2)+
geom_rect(data=NULL, aes(xmin=-Inf, xmax=Inf, ymin=100, ymax=500), fill="orange", alpha=0.2)+
geom_rect(data=NULL, aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=100), fill="yellow", alpha=0.2)+
geom_boxplot(data=df2, aes(x=as.factor(clust), y=dat_colonydiameter), fill="grey")+
scale_y_log10()+
theme_classic()
splot



plot_grid(plot_grid(fplot+guides(fill="none"), fplot2, ncol=1), cplot2,  nrow=1,
rel_widths=c( 0.5, 1), labels=c("a", "b"))

head(df2)

######################################################
#---------------------------------------------#  Morgan



cplotMog <- ggplot() +
geom_polygon(data=dens,aes(x, y,group=Group,fill=z), alpha=0.05)+
#geom_path(data=lins,aes(x, y,group=Group, col=nom))+
 # geom_point(data=df2, aes(dim1, dim2, fill=factor(kclust)), shape=21, size=0.75, stroke=0)+
 geom_point(data=df2, aes(dim1, dim2), shape=21, size=0.25, stroke=0.25)+
  #stat_ellipse(data=clust, geom="polygon", aes(dim1, dim1, group=gp, col=variable, fill=variable), alpha=0.2)+
  stat_ellipse(data=df2, aes(dim1, dim2, group=factor(kclust)), fill=NA, col='black',geom="polygon", alpha=0.2)+
   #stat_ellipse(data=df2, aes(dim1, dim2, group=factor(mclust)), geom="polygon", fill="cadetblue",col="cadetblue", alpha=0.2, level=0.9)+
   scale_fill_distiller(palette="YlOrRd", direction=1)+
   scale_colour_manual(values=c("grey90", "grey80", "grey60"))+
labs(x="PCoA 1", y="PCoA 2")+
geom_text(data=NULL, aes(x=0.34, y=0.13, label="Tabular (Ac)"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.4, y=-0.02, label="Solitary"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.4, y=-0.24, label="Massive"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.13, y=-0.29, label="Columnar"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.15, y=0.29, label="Laminar"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.43, y=0.23, label="Encrusting"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=0.35, y=-0.15, label="Branching (Ac)"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=0.33, y=0.02, label="Corymbose"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=0.17, y=-0.18, label="Branching (other)"), col="black", size=2.5, hjust=0)+
#geom_text(data=NULL, aes(x=-0.4, y=-0.28, label="K-means"), col="red", size=3, hjust=0)+
lims(x=c(-0.55, 0.55), y=c(-0.3, 0.3))+
guides(fill="none")+
labs(x="PCoA 1", y="PCoA 2")+
theme_bw()+theme(legend.title=element_blank(), legend.position=c(0.85, 0.9), legend.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(2, "mm"), legend.text=element_text(size=7), axis.text=element_text(size=7), axis.title=element_text(size=8))+
annotation_custom(tab.s, xmin=0.35, xmax=0.45, ymin=0.12, ymax=0.22)+
annotation_custom(cor.s, xmin=0.25, xmax=0.35, ymin=0, ymax=0.1)+
annotation_custom(stag.s, xmin=0.4, xmax=0.5, ymin=0.05, ymax=-0.25)+
annotation_custom(br.s, xmin=0.15, xmax=0.23, ymin=-0.1, ymax=-0.2)+
annotation_custom(mas.s, xmin=-0.23, xmax=-0.35, ymin=-0.15, ymax=-0.23)+
annotation_custom(col.s, xmin=-0.1, xmax=0, ymin=-0.2, ymax=-0.28)+
annotation_custom(lam.s, xmin=-0.05, xmax=-0.15, ymin=0.2, ymax=0.3)+
annotation_custom(fun.s, xmin=-0.38, xmax=-0.48, ymin=-0.05, ymax=0.05)+
annotation_custom(enc.s, xmin=-0.3, xmax=-0.39, ymin=0.15, ymax=0.25)
cplotMog




######################################################
#---------------------------------------------#  Coral Sea (benthos)

# CCA.Pav/Turf.Pav a subset of pavement (Sum to Pavement)

b1og <- read.csv("data/CoralSea/CSMP-GBR_2018-2024_Coral_PIT_CLEAN.03.03.2024.csv")
b1og$Site[b1og$Site==" Mellish 6"] <- "Mellish 6"
b1og$Site <- gsub("Chilcot 2", "Chilcott 2", b1og$Site)
b1fg <- read.csv("data/CoralSea/bFG.csv")
b1og <- subset(b1og, select=-c(Pavement, X.CCA)) #Pavement, X.CCA
head(b1og) 
#unique(b1og$Site)[order(unique(b1og$Site))]

b1og$id <- paste(b1og$Site, b1og$Zone, b1og$Transect, b1og$Year)

tcols <- c("Year", "Date", "Marine.Park", "SectorRegion", "Name", "Reef", "Site","Site_notes", "Zone", "Transect", "Depth", "Complexity", "Total.Macroalgae", "Total.HARD.CORAL", "Gradient", "Total", "id")

b1t <- b1og[,tcols]
length(unique(b1t$id))
head(b1t)
nrow(b1t)

b1 <- melt(b1og, id.var=tcols, value.name="n")
b1$group <- b1fg$group[match(b1$variable, b1fg$taxon)]
head(b1)


b1$n[is.na(b1$n)] <- 0 

# total points
b1total <- aggregate(n ~ id, b1, sum)
b1t$total <- b1total$n[match(b1t$id, b1total$id)]
b1t$match <- ifelse(b1t$total==b1t$Total, "yes", "no")
table(b1t$match)
#b1t[b1t$match=="no",]
#b1[b1$id=="Boot 5 Slope 1 2023",]

# cover
b1coral <- aggregate(n ~ id, b1[b1$group=="HC",], sum)
b1t$coral <- b1coral$n[match(b1t$id, b1coral$id)]
b1t$match <- ifelse(b1t$coral==b1t$Total.HARD.CORAL, "yes", "no")
table(b1t$match)




################################## Div index
#b1$morph <- b1fg$morph7[match(b1$variable, b1fg$taxon)]
b1$morph <- b1fg$morph8[match(b1$variable, b1fg$taxon)]
b1m <- aggregate(n~morph+id, b1[!b1$morph=="none",], sum)
b1m$coral <- b1t$coral[match(b1m$id, b1t$id)]
b1m$rel_cover <-  b1m$n / b1m$coral
b1m$simpD <- b1m$rel_cover^2
b1m$shanH <- b1m$rel_cover * log(b1m$rel_cover)
coral_simp <- aggregate(simpD~id, b1m, sum)
coral_shan <- aggregate(shanH~id, b1m, sum)
b1t$simpD <- 1 - coral_simp$simpD[match(b1t$id, coral_simp$id)]
b1t$simpD[b1t$coral==0] <- min(b1t$simpD, na.rm=T) #0
b1t$shanH <- - coral_shan$shanH[match(b1t$id, coral_shan$id)]
b1t$shanH[b1t$coral==0] <- min(b1t$shanH, na.rm=T)  #0
#ggplot(b1t, aes(simpD, shanH))+geom_point()


######################################################
#---------------------------------------------#  TRAIT DIVERSITY! 

head(b1)

head(df2)

# LINK = groups in coral sea... 

coord <- aggregate(dim1~link, df2, mean)
coord$dim2 <- aggregate(dim2~link, df2, mean)$dim2
coord$dim3 <- aggregate(dim3~link, df2, mean)$dim3
coord$dim4 <- aggregate(dim4~link, df2, mean)$dim4
coord
ggplot() +
  geom_point(data=df2, aes(dim1, dim2), shape=21, size=0.5, stroke=0.25, col="black")+
geom_point(data=coord, aes(dim1, dim2), col="red")


b1$dim1 <- coord$dim1[match(b1$variable, coord$link)]
b1$dim2 <- coord$dim2[match(b1$variable, coord$link)]
b1$dim3 <- coord$dim3[match(b1$variable, coord$link)]
b1$dim4 <- coord$dim4[match(b1$variable, coord$link)]
unique(b1[b1$group=="HC",c("variable", "dim1", "dim2")])

# fdisp function
temp <- na.omit(b1[b1$group=="HC" & b1$Total.HARD.CORAL>0, c("Site", "id", "variable", "n", "dim1", "dim2", "dim3", "dim4")])
abuns <-acast(temp, id~variable, value.var="n") #n
colSums(abuns)
coord2 <- unique(temp[,c("variable", "dim1","dim2", "dim3", "dim4")])
tr.dat <- dist(data.frame(coord2[,c("dim1", "dim2", "dim3")], row.names=coord2$variable))
disp <-	data.frame(FDis=fdisp(a=abuns, tr.dat)$FDis)

b1t$FDis <- disp$FDis[match(b1t$id, rownames(disp) )]


######################################################
#---------------------------------------------#  OUTPUT

# write.csv(df2, "data/output/coraltraits.csv")

# write.csv(b1t, "data/output/coralFD.csv")


