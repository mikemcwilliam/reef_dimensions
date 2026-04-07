
rm(list = ls())

library("ggplot2")
library("cowplot")
library("viridis")
library("reshape2")
library("lubridate")
library("vegan")
library("psych")
library("png")
library("grid")

######################################################
#---------------------------------------------#  data

df1 <- read.csv("data/data.csv")
#df1 <- read.csv("data/metrics/metrics.csv")
df1$lat <- ifelse(df1$Data=="Coral Sea", -df1$lat, df1$lat)
df1 <- df1[!df1$Marine.Park=="Lord Howe Marine Park",]
head(df1)

######################################################
#---------------------------------------------#  metrics

labs <- read.csv("data/info/PCAlabs.csv")

df1$log_N <- log(df1$f.N)
df1$log_biomass <- log(df1$f.biomass)

df1$alg_ratio1 <- df1$cca / (df1$cca + df1$turf + df1$algae)
df1$alg_ratio <- sqrt(df1$alg_ratio1)

df1$P.herb1 <- df1$N.herb / df1$f.N 
df1$P.herb <- sqrt(df1$P.herb1) 

df1$P.carn1 <- df1$N.carn / df1$f.N
df1$P.carn <- sqrt(df1$P.carn1) #. df1$N.carn / df1$f.N

df1$covsqrt <- sqrt(df1$coral)
df1$richsqrt <- sqrt(df1$f.richness)
df1$acrosqrt <- sqrt(df1$acro)
df1$turfsqrt <- sqrt(df1$turf)
hist(sqrt(df1$f.richness))

axes1a <- c("covsqrt", "Complexity", "simpD", "log_biomass","algae", "turfsqrt", "richsqrt", "P.herb", "P.carn", "alg_ratio", "FsimpD", "comp1", "comp2", "acrosqrt")

######################################################
#---------------------------------------------#  full PCA

rownames(df1)<- df1$id
pca1a <- prcomp(na.omit(df1[,axes1a]), scale=T, center=T)
df1[,c("PC1_1a", "PC2_1a", "PC3_1a", "PC4_1a")] <- pca1a$x[match(rownames(df1), rownames(pca1a$x)),c("PC1", "PC2", "PC3", "PC4")]
vecs1a <- data.frame(pca1a$rotation[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")], lab=rownames(pca1a$rotation))
vecs1a$lab2 <- labs$lab2[match(vecs1a$lab, labs$lab)]
vars1a <- round((pca1a$sdev^2 / sum(pca1a$sdev^2)), 3)*100	
biplot(pca1a)
sum(vars1a[1:6])

######################################################
#---------------------------------------------#   loadings

loadings <- abs(pca1a$rotation[,1:6])
rownames(loadings)[apply(loadings,2,which.max)]

### new method
variable_scores <- apply(loadings, 1, max)
top_vars <- sort(variable_scores, decreasing=T)
top_6 <- names(top_vars)[1:6]
top_6

# select top X for PC1, remove, then select top X for PC2
n_per <- 2 # number of metrics per axis
grad <- NULL
for(i in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")){
	#i <- "PC2"
	temp <- data.frame(lab=vecs1a$lab, val = vecs1a[,i])
	temp <- temp[!temp$lab %in% c(grad$lab),]
	mets <- unique(temp$lab)
	ord <- temp[order(abs(temp$val), decreasing=T),]
	topX <- ord[c(1:n_per),]
	grad <- rbind(grad, cbind(topX, i))
	}
grad
	
# select above a certain value
cutoff <- 0.4	
grad2 <- NULL
for(i in c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")){
	#i <- "PC2"
	temp <- data.frame(lab=vecs1a$lab, val = vecs1a[,i])
	temp <- temp[!temp$lab %in% c(grad2$lab),]
	temp <- temp[abs(temp$val)>cutoff,]
	grad2 <- rbind(grad2, cbind(temp, i))
	}
grad2

vecs1a$i <- grad$i[match(vecs1a$lab, grad$lab)] # selection method

######################################################
#---------------------------------------------#   max loadings

# maximum absolute value
vecslong <- melt(vecs1a, id.var=c("i", "lab", "lab2"))
vecslong$select <- ifelse(vecslong$i==vecslong$variable, "yes", "no") 
vecslong$ord <- as.numeric(labs$ord3[match(vecslong$lab, labs$lab)])
vecslong$abs <- abs(vecslong$value)
vmax <- merge(aggregate(abs ~ lab2, max, data = vecslong), vecslong)
vecs1a$max <- vmax$variable[match(vecs1a$lab2, vmax$lab2)]
vecs1a

ggplot()+
geom_segment(data=vecs1a, aes(x=0, xend=PC1*0.9, y=0, yend=PC2*0.9, col=max), arrow=arrow(length=unit(1,"mm")))+
geom_text(data=vecs1a, aes(x=PC1, y=PC2, label=lab2), hjust=ifelse(vecs1a$PC1>0, 0, 1), size=2.3, fontface="italic")+
theme_void()


head(vecslong)
vecslong$max <- vecs1a$max[match(vecs1a$lab2, vecslong$lab2)]
vecslong$maxval <- ifelse(vecslong$max==vecslong$variable, vecslong$value, NA)
vecslong$selectval <- ifelse(vecslong$select=="yes", vecslong$value, NA)
vecslong$six <- ifelse(vecslong$lab2 %in% c("coral %", "Feeding div. (D)", "fish Biomass", "CCA ratio", "Morph. div (D)", "complexity (1-5)") & vecslong$select=="yes", "*", "")
head(vecslong)


ggplot()+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=value, fill=variable), stat="identity", width=0.7,  alpha=0.5)+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=selectval, fill=variable), stat="identity", width=0.7, position="stack")+
geom_bar(data=vecslong, aes(y=reorder(lab2, -ord), x=maxval), stat="identity", fill=NA, col="black", width=0.7, size=0.2)+
geom_text(data=vecslong[vecslong$value>0,], aes(y=reorder(lab2, -ord), x=value+0.05, fill=variable, label=six), vjust=0.75)+
geom_text(data=vecslong[vecslong$value<0,], aes(y=reorder(lab2, -ord), x=value-0.05, fill=variable, label=six), vjust=0.75)+
geom_vline(xintercept=0)+
facet_wrap(~variable, nrow=1)

######################################################
#---------------------------------------------#  correlations

library("ggcorrplot")

#library("psych")
# pairs.panels(na.omit(df1[,axes1a]), scale=T, cex.cor=3, cex=0.1)

head(labs)
ords <- data.frame(axes1a)
ords$lab2 <- labs$lab2[match(ords$axes1a, labs$lab)]
ords$ord <- labs$ord2[match(ords$axes1a, labs$lab)] # ord = pca, ord2 = fish.coral
ords <- ords[order(ords$ord),]
ords

pcor<-cor(df1[,ords$axes1a], use="pairwise.complete.obs")
#write.csv(pcor, "cor_output.csv")
ggcorrplot(pcor,type="upper")

########## new clustering method. 

pcor2 <- pcor
rownames(pcor2)<-labs$lab2[match(rownames(pcor2), labs$lab)]
dist_mat <- as.dist(1 - abs(pcor2)) # uncorrelated variables = 1
hc <- hclust(dist_mat, method="complete")
plot(hc, mean="clusters") # how do vatiables cluster based on correlation

clusters <- cutree(hc, k=6)
clust_vars <- split(names(clusters), clusters) # which variables in each cluster?
clust_vars 

cplot <- melt(clust_vars )
cplot

######################################################
#---------------------------------------------# fig1

source("figs/fig1.R")
# fig1

# ggsave("figs/Fig1.jpg", fig1X2, height=7, width=6.9)


######################################################
#---------------------------------------------#   CAN 6 METRICS PREDICT 14? 

six <- c("covsqrt", "Complexity", "simpD", "log_biomass", "alg_ratio", "FsimpD")

############# PCA RECONSTRUCTION

pc_scores <- pca1a$x[,c(1:6)]
head(pc_scores)
# Build a formula dynamically: PC1 + PC2 + ... + PCK ~ selected variables
pc_names <- paste0("PC", 1:6)
formula_str <- paste(paste(pc_names, collapse = " + "), "~", paste(six, collapse = " + "))
multi_formula <- as.formula(formula_str)
multi_formula

df1[,pc_names]<-pc_scores[match(df1$id, rownames(pc_scores)), pc_names]

# Fit multivariate linear model (each PC is a response)
mlm_model <- lm(multi_formula, data = df1)

get_r2 <- function(model) {
  sapply(pc_names, function(pc) {
    summary(lm(as.formula(paste(pc, "~", paste(six, collapse = " + "))), data = df1))$r.squared
  })
}

r2_values <- get_r2(mlm_model)
r2_values

pc_var_explained <- summary(pca1a)$importance["Proportion of Variance", names(r2_values)]
pc_var_explained

total_variance_explained <- sum(pc_var_explained * r2_values)

cat("Total variance in the original dataset explained by the selected variables:", 
    round(total_variance_explained * 100, 2), "%\n")

############# MANTEL TEST

library(vegan)  # for mantel()

# Compute distance matrices (e.g., Euclidean distance)
dist_full <- dist(scale(df1[, axes1a]))
dist_subset <- dist(scale(df1[, six]))

# Run Mantel test
mantel_result <- mantel(dist_full, dist_subset, method = "pearson", permutations = 99)

print(mantel_result)

############ COMPARE PCAS

df1x <- df1
axes2 <- axes1a
axes.list <- list(axes2, six)

coords <- NULL
vectors <- NULL
corrs <- NULL
for(i in unique(df1$Data)){
	for(j in 1:length(axes.list)){
df2 <- df1x[df1x$Data==i,]
axes.use <- axes.list[[j]]
rownames(df2) <- df2$id
pca <- prcomp(na.omit(df2[,axes.use]), scale=T, center=T)	
vecs <- data.frame(pca$rotation[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")], lab=rownames(pca$rotation))
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
vecs$max <- colnames(vecs[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")])[max.col(abs(vecs[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")]), ties.method='first')]	
simp <- ifelse(j==1, "long", "short")
pcor <- cor(df2[,axes.use], use="pairwise.complete.obs")
pcor[upper.tri(pcor, diag=T)]<-NA
corrs <- rbind(corrs, cbind(na.omit(melt(pcor)),axes=simp, data=i) )
coords <- rbind(coords, cbind(data.frame(id=rownames(pca$x), pca$x[,c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")]), axes=simp, data=i))
vectors <- rbind(vectors, cbind(data.frame(vecs), axes=simp, data=i))
}}

head(coords)
head(corrs)
head(vectors)

vectors$PC2r <- vectors$PC2
vectors$PC2r[vectors$data=="LTMP"] <-  - vectors$PC2r[vectors$data=="LTMP"]  

ggplot()+geom_segment(data=vectors, aes(x=0, xend=PC1, y=0, yend=PC2r, col=axes))+
geom_text(data=vectors, aes(x=PC1, PC2r, label=lab, col=axes))+
facet_wrap(~data)

############  ALIGN THE COORDINATES (sites) in long and short PCAs

head(coords)
clong <- melt(coords, id.var=c("axes", "data", "id"))
head(clong)
table(clong$axes)
cback <- dcast(clong, variable+id+data~axes, value.var="value")
cback$long <- as.numeric(cback$long)
cback$short <- as.numeric(cback$short)
head(cback)


# correlations between coords
lm1 <- lm(long~short, cback[cback$variable=="PC1",])
summary(lm1)

lm2 <- lm(long~short, cback[cback$variable=="PC2",])
summary(lm2)

#check by hand and reverse if needed
ggplot(data=cback, aes(short, long, col=data))+geom_point(shape=21)+facet_wrap(~variable)+geom_smooth(method="lm")


head(cback)
cback$long2 <- cback$long
#cback$long2 <- ifelse(cback$variable=="PC6" & cback$variable=="LTMP", -cback$long2, cback$long2)


# check all
preds <- NULL
for(i in unique(df1$Data)){
	for(j in unique(cback$variable)){
		temp <- cback[cback$variable==j, ]
		temp <- temp[temp$data==i, ]
		lm1 <- lm(long2~short, temp)
	slp <- coef(lm1)[2]
	low <- confint(lm1)[2,1]
	upp <- confint(lm1)[2,2]
	pval <- coef(summary(lm1))[2,"Pr(>|t|)"]
	rsq <- summary(lm1)$r.squared
	arsq <- summary(lm1)$adj.r.squared
preds <- rbind(preds, data.frame(data=i, var=j, slp, low, upp, pval, rsq, arsq))
}}
preds

preds$cod <- preds$slp * sqrt(preds$rsq) / abs(preds$slp)

plot_grid(
ggplot(preds, aes(var, abs(cod), fill=data))+geom_bar(stat="identity", position="dodge"),
ggplot(preds, aes(var, rsq, fill=data))+geom_bar(stat="identity", position="dodge"))


######################################################
#---------------------------------------------# fig2

source("figs/fig2.R")
# Fig2
# ggsave("figs/Fig2.jpg", Fig2, height=6, width=6.9)

######################################################
#---------------------------------------------# variance components

library("VCA")
library("lme4")
library("MuMIn")
library("lmerTest")

head(df1) # SectorRegion > Reef > Zone > Site
unique(df1$SectorRegion)
df1x <-  df1 #df1[!df1$SectorRegion %in% c("Elizabeth and Middleton", "Northern Great Barrier Reef", "Southern Great Barrier Reef" ,"Central Great Barrier Reef"),]
vars <- NULL
for(i in unique(df1x$Year)){
#	for(x in unique(df1x$Data)){
	for(j in  c(axes1a, "PC1_1a", "PC2_1a", "PC3_1a")){ #c("coral", "turf", "cca", "soft", "acro")
		#j <- "coral"
df1x$y <- df1x[,j]
mod1 <- lmer(y  ~ (1 | SectorRegion/Reef/Site),  data=na.omit(df1x[df1x$Zone=="Slope" & df1x$Year==i, c("y", "Site", "SectorRegion", "Reef", "Data")]))
#summary(mod1)
out <- data.frame(VarCorr(mod1))
out$sum <- sum(out$vcov)
out$p <- out$vcov / out$sum
out$Year <- i 
out$y <- j
#out$x <- x
vars <- rbind(vars, out)
}}

head(vars)
ggplot(vars[vars$y %in% c("PC1_1a"),], aes(x=grp, y=p*100,group=as.factor(Year)))+
geom_point(size=0.5)+geom_line()+facet_wrap(~y, scales="free", nrow=3)

###################################################### CORAL COMPOSITION
#---------------------------------------------# examples

maxmin <- read.csv("data/output/minmaxD.csv")
head(maxmin)
maxminT <- aggregate(n ~ id+type+t+morph8b+simpDx, maxmin[maxmin$group=="HC",], sum )
maxDs <- maxminT[maxminT$type =="high D",]
minDs <- maxminT[maxminT$type =="low D",]
maxminC <- aggregate(n ~ id+type+t+simpDx, maxminT, sum )
maxminC

ggplot(maxminT, aes(y=morph8b, x=n, fill=t))+geom_bar(stat="identity", position="dodge")+facet_grid(~type)

#---------------------------------------------# model cover vs SimpD
library("mgcv")

head(df1)

c1 <- df1[df1$Data=="Coral Sea",]
head(c1)

lm1 <- lm(simpD ~ coral, data=c1)
lm2 <- lm(simpD ~ poly(coral, 2), data=c1)
lm3 <- lm(simpD ~ poly(coral, 3), data=c1)
lm4 <- lm(simpD ~ poly(coral, 4), data=c1)
lm5 <- lm(simpD ~ poly(coral, 5), data=c1)
gam3 <- gam(simpD ~ s(coral, k=3), data=c1)
gam4 <- gam(simpD ~ s(coral, k=4), data=c1)
gam5 <- gam(simpD ~ s(coral, k=5), data=c1)
gam6 <- gam(simpD ~ s(coral, k=6), data=c1)
gam7 <- gam(simpD ~ s(coral, k=7), data=c1)
gam8 <- gam(simpD ~ s(coral, k=7), data=c1)
summary(gam8)
AIC(lm1,lm2,lm3, lm4,lm5,gam3, gam4, gam5, gam6, gam7, gam8)

fit <- data.frame(coral=seq(min(c1$coral), max(c1$coral), by=1))
fit$fit <- predict(gam7, fit)
fit$se <- predict(gam7, fit, se=T)$se.fit

ggplot()+
geom_point(data=c1, aes(coral, simpD), shape=21, size=0.5)+
geom_ribbon(data=fit, aes(x=coral, ymin=fit-se, ymax=fit+se), alpha=0.5)+
geom_line(data=fit, aes(coral, fit), col="red")

#---------------------------------------------#FDis vs Simp

c1 <- df1[df1$Data=="Coral Sea",]
head(c1)

c1$FDis[c1$coral==0]<-0

linkmod <- summary(lm(simpD~sqrt(FDis), c1))
linkmod
rsq <- round(linkmod$r.squared, 2)

ggplot(c1, aes(simpD, FDis))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
scale_y_sqrt()+
annotate(geom = 'text', x = 0.8, y = 0,    label = paste("R^2 == ", rsq), parse = TRUE, size=2.5, col="blue")


###################################################### fig3

source("figs/fig3.R")
FIG3s

# ggsave("figs/Fig3.jpg", FIG3s, height=4.8, width=6.5)

######################################################
#---------------------------------------------#  size structure

c1 <- df1[df1$Data=="Coral Sea",]
head(c1)

#c1$p10[c1$coral==0]<-0

#summary(lm(frich~coral_cover, tdf))
#summary(lm(frich~coral_cover*A_SECTOR, tdf))

bimod <- glm(p10~coral, family="quasibinomial", data=c1)
summary(bimod)
bifit <- data.frame(coral = seq(min(c1$coral, na.rm=T), max(c1$coral, na.rm=T), 1))
bifit$fit <- predict(bimod, bifit, type="response")

ggplot()+geom_line(data=bifit, aes(coral, fit))+geom_point(data=c1, aes(coral, p10), shape=21)


# lots of reefs with no small corals.. 
# reef variation trans/sites/years.. 

colz <- c("#b2182b", "#f4a582", "#92c5de", "#2166ac")
names(colz) <- c("s.neg", "ns.neg", "ns.pos", "s.pos")
head(c1)
c1$p60sqrt <- sqrt(c1$p60) #  <- ((c1$coral)/max(c1$coral)) * c1$simpD
reefs <- NULL
bifits <- NULL
Ys <- c("p10", "p60")
Xs <- c("coral")
	for(y in Ys){
		for(x in Xs){
temp <- df1[df1$Data=="Coral Sea",]
temp$x <- temp[,x]
temp$y <- temp[,y]
temp <- temp[!(is.na(temp$x) | is.na(temp$y)),]
zero.x <- aggregate(x~Reef, temp, sum) # IS THIS OK? 
temp <- temp[!temp$Reef %in% zero.x[zero.x$x==0, "Reef"],]
zero.y <- aggregate(y~Reef, temp, sum)
temp <- temp[!temp$Reef %in% zero.y[zero.y$x==0, "Reef"],]
for(i in unique(temp$Reef)){
	#i <- unique(temp$Reef)[1]
	temp2 <- temp[temp$Reef==i,]
	mod <- glm(y~x, family="quasibinomial", data=temp2)
	#summary(mod)
	slp <- coef(mod)[2]
	bifit <- data.frame(x = seq(min(temp2$x, na.rm=T), max(temp2$x, na.rm=T), 1))
	bifit$fit <- predict(mod, bifit, type="response")
	#low <- confint(mod)[2,1]
	#upp <- confint(mod)[2,2]
	pval <- coef(summary(mod))[2,"Pr(>|t|)"]
	sig <- ifelse(pval > 0.05, "ns", "s")
	dir <- ifelse(slp >0, "pos", "neg")
	sigdir <- paste(sig, dir, sep=".")
	bifits <- rbind(bifits, cbind(bifit, Reef=i, pred=x, resp=y, sig, dir, sigdir))
	reefs <- rbind(reefs, data.frame(Reef=i,  x,y, n=nrow(temp2), slp, pval, sig, dir,sigdir)) }
	}}

head(reefs)

ggplot()+
geom_line(data=bifits, aes(x, fit, group=Reef, col=sigdir))+
scale_colour_manual(values=colz)+
facet_grid(pred~resp)+
labs(x="% coral cover", y="Proportion size")


######################################################
#---------------------------------------------#  refuge vol

head(df1)

linkmod2 <- summary(lm(sqrt(shelt)~coral, df1))
linkmod2
rsq2 <- round(linkmod2$r.squared, 2)

ggplot(c1, aes(x=coral, y=shelt))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
scale_y_sqrt()+
annotate(geom = 'text', x = 75, y = 0,    label = paste("R^2 == ", rsq2), parse = TRUE, size=2.5, col="blue")

######################################################
#---------------------------------------------#  fig 4

source("figs/Fig4.R")
FIG4

# ggsave("figs/Fig4.jpg", FIG4, height=3, width=8.5)


######################################################
#---------------------------------------------#  benthos - fish

ggplot(df1, aes(coral, f.biomass, col=Data))+geom_point(size=0.5)+
#facet_wrap(~Data)+
scale_y_log10()+
geom_smooth()+
theme_classic()

reefs <- NULL
Ys <- c("f.richness", "log_biomass", "FsimpD","FsimpD2", 'P.herb')
Xs <- c("coral", "Complexity", "simpD", "alg_ratio")
#Xs <- c("soft", "turf", "simpD")
#temp <- na.omit(temp[,c("Reef","Zone","Site","id", Xs, Ys)])
	for(y in Ys){
		for(x in Xs){
		#	x <-"soft"
		#	y="turf"
temp <- df1
temp$x <- temp[,x]
temp$y <- temp[,y]
temp <- temp[!(is.na(temp$x) | is.na(temp$y)),]
for(i in unique(temp$Reef)){
	#i <- unique(temp$Reef)[1]
	temp2 <- temp[temp$Reef==i,]
	head(temp2)
	mod <- lm(y~x, temp2)
	slp <- coef(mod)[2]
	low <- confint(mod)[2,1]
	upp <- confint(mod)[2,2]
	pval <- coef(summary(mod))[2,"Pr(>|t|)"]
	rsq <- summary(mod)$r.squared
	arsq <- summary(mod)$adj.r.squared
	reefs <- rbind(reefs, data.frame(reef=i,  x,y, n=nrow(temp2),dat=unique(temp2$Data), slp, low, upp, pval, rsq, arsq)) }
	}}
	

reefs$sig <- ifelse(reefs$pval > 0.05, "ns", "s")
reefs$dir <- ifelse(reefs$slp >0, "pos", "neg")
#reefs$sec <- benth$A_SECTOR[match(reefs$reef, benth$REEF_NAME)]
reefs$sigdir2 <- paste(reefs$sig, reefs$dir, sep=".")
head(reefs)

fsimp2 <- reefs[reefs$y=="FsimpD2",]
reefs <- reefs[!reefs$y=="FsimpD2",]

colz <- c("#b2182b", "#f4a582", "#92c5de", "#2166ac")
names(colz) <- c("s.neg", "ns.neg", "ns.pos", "s.pos")

reefs$x2 <- labs$lab2[match(reefs$x, labs$lab)]
reefs$y2 <- labs$lab2[match(reefs$y, labs$lab)]

reefs$cod <- reefs$slp * sqrt(reefs$rsq) / abs(reefs$slp)

avs <- aggregate(cod~x2+y2+dat,reefs, mean)

ggplot(reefs[!reefs$slp==1,], aes(x=cod))+
geom_density(aes(fill=dat), alpha=0.5)+
geom_segment(data=avs, aes(x=cod, xend=cod, y=Inf, yend=-Inf, col=dat))+
theme_minimal()+theme(legend.title=element_blank())+
labs(x="coefficient")+
geom_vline(xintercept=0, linetype="dotted")+
facet_grid(x2~y2)

head(reefs)
slps1 <- reefs[reefs$x=="Complexity" & reefs$y=="f.richness",]
df1$sigdir1 <- slps1$sigdir[match(df1$Reef, reefs$reef)]

plot_grid( ggplot(data=data.frame(table(df1$sigdir1)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+scale_fill_manual(values=colz)+geom_bar(stat="identity"),
ggplot(df1, aes(Complexity, f.richness))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
facet_wrap(~Data)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir1), se=F, size=0.3), 
rel_widths=c(1, 2),  nrow=1)

######################################################
#---------------------------------------------#  fig 5 

source("figs/fig5.R")
FIG5
# ggsave("figs/Fig5.jpg", FIG5, height=6.5, width=6.5)


######################################################
#---------------------------------------------#  benthos - fish


slpsX <- fsimp2[fsimp2$x=="coral" & fsimp2$y=="FsimpD2",]
df1$sigdirX <- slpsX$sigdir[match(df1$Reef, fsimp2$reef)]


ggplot(df1, aes(coral, FsimpD2))+geom_point(size=0.1, col="grey")+
#geom_smooth(method="lm", aes(group=Reef), se=F, size=0.2)+
geom_smooth(method="lm", aes(group=Reef, col=sigdirX), se=F, size=0.3)+
facet_wrap(~Data)+
scale_colour_manual(values=colz)+
labs(x="% coral cover", y="Feeding diversity (3 groups)")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))



##### other
ggplot(df1, aes(coral, alg_ratio, col=turf))+geom_point(size=0.5)+
scale_colour_viridis()+
facet_wrap(~Data)+
labs(x="% coral cover", y="CCA ratio", col="% turf\ncover")+
theme_classic()+theme(strip.background=element_blank())

ggplot(df1, aes(alg_ratio, P.herb))+geom_point(size=0.5)+
scale_colour_viridis()+
#facet_wrap(~Data)+
geom_smooth()+
theme_classic()


