

rm(list = ls())

library("ggplot2")
library("cowplot")
library("viridis")
library("reshape2")

# Calculate FD for coral sea/LTMP


######################################################
#--------------------------------------------- traits


df <- read.csv("data/traits/traitbiogeography.csv") # mcwilliam et al data
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


df2$Growth.form.typical <- df2$raw_growth_form
df2$Growth.form.typical <- ifelse(df2$Growth.form.typical %in% c( "encrusting_long_uprights"), "encrusting", ifelse(df2$Growth.form.typical=="bifacial", "laminar", ifelse(df2$Growth.form.typical=="hispidose", "branching_open", ifelse(df2$Growth.form.typical %in% c("submassive", "solitary_attached", "solitary_free"), "massive",   as.character(df2$Growth.form.typical )))))
unique(df2$Growth.form.typical)


######################################################
#--------------------------------------------- GLOBAL PCA

library("FD") 

hist(df2$dat_colonydiameter)

df2$log_growth <- log(df2$dat_growth) #log(df2$dat_growth)
df2$log_corallite <- log(df2$dat_corallite) # log(df2$dat_corallite)
df2$log_diameter <- log(df2$dat_colonydiameter) # log(df2$dat_colonydiameter)

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
 #      method = "kmeans", min.nc=5, max.nc=15, 
 #       index = "alllong", alphaBeale = 0.1)
# hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])))
# write.csv(data.frame(nb$Best.nc[1,]), "data/output/optimal.all.csv")

opt <- read.csv("data/output/optimal.all.csv")

optplot <- ggplot(opt, aes(nb.Best.nc.1...))+
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

tab.s<-readPNG("figs/sils/tab.png")
tab.s<-rasterGrob(tab.s, interpolate=TRUE)
stag.s<-readPNG("figs/sils/stag.png")
stag.s<-rasterGrob(stag.s, interpolate=TRUE)
cor.s<-readPNG("figs/sils/cor.png")
cor.s<-rasterGrob(cor.s, interpolate=TRUE)
mas.s<-readPNG("figs/sils/pori.png")
mas.s<-rasterGrob(mas.s, interpolate=TRUE)
col.s<-readPNG("figs/sils/iso.png")
col.s<-rasterGrob(col.s, interpolate=TRUE)
br.s<-readPNG("figs/sils/poc.png")
br.s<-rasterGrob(br.s, interpolate=TRUE)
lam.s<-readPNG("figs/sils/lam.png")
lam.s<-rasterGrob(lam.s, interpolate=TRUE)
fun.s<-readPNG("figs/sils/fung.png")
fun.s<-rasterGrob(fun.s, interpolate=TRUE)
enc.s<-readPNG("figs/sils/encr.png")
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


FIGS1 <- plot_grid(tspace1,
plot_grid(wssplot, optplot, ncol=1, align="hv", labels=c("b", "c"), label_size=10),
cplot, nrow=1, rel_widths=c(1,0.65,1), labels=c("a", "", "d"), label_size=10)
FIGS1 

# ggsave("figs/supplement/SUP1.jpg", FIGS1, height=3, width=8)


# write.csv(df2, "data/output/coraltraits.csv")




