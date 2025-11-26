
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

c1 <- read.csv("data/output/metricsCS.csv")
c2 <- read.csv("data/output/metricsLTMP.csv")

ggplot(c2, aes(f.biomass))+geom_histogram()+scale_x_log10()

df1 <- read.csv("data/output/metrics_v1biomass.csv")
df1$lat <- ifelse(df1$Data=="Coral Sea", -df1$lat, df1$lat)

unique(df1$Site)
######################################################
#---------------------------------------------#  metrics

labs <- read.csv("data/output/PCAlabs.csv")

df1$log_N <- log(df1$f.N)
df1$log_biomass <- log(df1$f.biomass)
df1$P.herb1 <- df1$N.herb / df1$f.N 
df1$P.carn1 <- df1$N.carn / df1$f.N
df1$covsqrt <- sqrt(df1$coral)
#df1$alg_ratio <- df1$turf / (df1$cca + df1$turf)
df1$alg_ratio1 <- df1$cca / (df1$cca + df1$turf + df1$algae)

# more SQRTs? 
df1$algsqrt <- sqrt(df1$algae)
df1$richsqrt <- sqrt(df1$f.richness)
df1$log_soft <- sqrt(df1$soft)
df1$acrosqrt <- sqrt(df1$acro)
df1$turfsqrt <- sqrt(df1$turf)

df1$P.herb <- sqrt(df1$P.herb1) # df1$N.herb / df1$f.N original no SQRT
df1$P.carn <- sqrt(df1$P.carn1) #. df1$N.carn / df1$f.N
df1$alg_ratio <- sqrt(df1$alg_ratio1)
#df1$algae <- sqrt(df1$algae)

hist(sqrt(df1$f.richness))

unique(df1$Marine.Park)

df1 <- df1[!df1$Marine.Park=="Lord Howe Marine Park",]

head(df1)
axes1a <- c("covsqrt", "Complexity", "simpD", "log_biomass","algae", "turfsqrt", "richsqrt", "P.herb", "P.carn", "alg_ratio", "FsimpD", "comp1", "comp2", "acrosqrt")

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

### chat GPT method
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

########## CHAT GPT clustering method. 

dist_mat <- as.dist(1 - abs(pcor)) # uncorrelated variables = 1
hc <- hclust(dist_mat, method="complete")
plot(hc, mean="clusters") # how do vatiables cluster based on correlation

clusters <- cutree(hc, k=6)
clust_vars <- split(names(clusters), clusters) # which variables in each cluster?
clust_vars 

cplot <- melt(clust_vars )
cplot

representatives <- sapply(clust_vars, function(group_vars){
	if(length(group_vars)==1) return(group_vars)
	sub_cor <- abs(pcor[group_vars, group_vars])
	avg_cor <- rowMeans(sub_cor) # select var with highest average corellation
	group_vars[which.max(avg_cor)]
	})
representatives

######################################################
#---------------------------------------------# fig1

source("figs/fig1.R")
# fig1X2

# ggsave("figs/Fig1.jpg", fig1X2, height=7, width=6.9)



######################################################
#---------------------------------------------#   CAN 6 METRICS PREDICT 14? 
# And one extra CS with colony size? 

six <- c("covsqrt", "Complexity", "simpD", "log_biomass", "alg_ratio", "FsimpD")

# CHAT GPT recommend

# PCA reconstruction. 

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

head(df1)


# MULTIVARIATE REGRESSION

# # Assume df is your original data frame with all 14 numeric variables
# and selected_vars is a character vector of your 6 selected variable names

# Define predictor and response variables
predictors <- df1[, six]
responses <- df1[, setdiff(axes1a, six)]  # The remaining 8 variables

# Fit multivariate linear regression
mlm_model2 <- lm(as.matrix(responses) ~ as.matrix(predictors))

# Summary per response variable
summary_list <- summary(mlm_model2)

# Extract R-squared values for each predicted variable
r_squared <- summary_list$r.squared
r_squared <- sapply(1:ncol(responses), function(i) {
  summary(lm(responses[, i] ~ as.matrix(predictors)))$r.squared
})
names(r_squared) <- colnames(responses)

# Print R² for each predicted variable
print(r_squared)


r2_df <- data.frame(Variable = names(r_squared), R2 = r_squared)

ggplot(r2_df, aes(x = reorder(Variable, R2), y = R2)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "R² of predicting variables from 6 selected predictors",
       x = "Predicted Variable",
       y = "R²") +
  theme_minimal()


## MANTEL TEST

library(vegan)  # for mantel()

# Compute distance matrices (e.g., Euclidean distance)
dist_full <- dist(scale(df1[, axes1a]))
dist_subset <- dist(scale(df1[, six]))

# Run Mantel test
mantel_result <- mantel(dist_full, dist_subset, method = "pearson", permutations = 99)

print(mantel_result)


# suggested plots 

df_pca <- data.frame(
  PC = factor(paste0("PC", 1:6)),
  VarianceExplained = pca1a$sdev[1:6]^2 / sum(pca1a$sdev^2),  # proportion of variance per PC
  R2 = r2_values  # your PCA reconstruction R2 vector
)

ggplot(df_pca, aes(x = PC)) +
  geom_col(aes(y = VarianceExplained), fill = "lightblue") +
  geom_line(aes(y = R2), group = 1, color = "red", size = 1) +
  geom_point(aes(y = R2), color = "red", size = 2) +
  scale_y_continuous(
    name = "Variance Explained (Proportion)",
    sec.axis = sec_axis(~., name = "R² of Reconstruction")
  ) +
  labs(title = "Variance Explained and Reconstruction R² per Principal Component") +
  theme_minimal()

varplot <- ggplot(df_pca, aes(x = PC)) +
  geom_col(aes(y = R2), fill = "lightblue", col="black", size=0.1, width=0.7) +
  #geom_line(aes(y = R2), group = 1, color = "red", size = 1) +
  #geom_point(aes(y = R2), color = "red", size = 2) +
  scale_y_continuous( name = "R² of Reconstruction", expand=c(0,0))  +
  labs(title = "Original PCA ~ 6 metrics", x="Original PCA axis") +
  theme_classic()+theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.text.y=element_text(size=8),axis.text.x=element_text(size=8, angle=45, hjust=1),axis.title=element_text(size=8), axis.line=element_line(size=0.1))
varplot

dist_full_vec <- as.vector(dist_full)
dist_subset_vec <- as.vector(dist_subset)

df_dist <- data.frame(
  Dist_Full = dist_full_vec,
  Dist_Subset = dist_subset_vec
)
nrow(df_dist)


my_breaks <- c(2, 4, 6)
my_labels <- sapply(my_breaks,function(i)as.expression(bquote(10^ .(i))))

distplot <- ggplot(df_dist, aes(x =  Dist_Subset, y =Dist_Full)) +
  geom_hex(data=df_dist,aes(y = Dist_Full, x = Dist_Subset)) +
  geom_smooth(method = "lm", color = "red", size=0.25) +
 geom_abline()+
 scale_fill_viridis( trans="log",breaks = c(1, 10^my_breaks), labels = c("1", my_labels))+ 
  labs( y = "Distance - 14 metrics",
       x = "Distance - 6 metrics") +
       ggtitle("Pairwise distances")+
  theme_classic()+theme(legend.key.width=unit(0.75, "mm"), legend.key.height=unit(2.5, "mm"), legend.position=c(0.92, 0.23), legend.title=element_text(size=8), legend.text=element_text(size=7), legend.background=element_blank(), axis.text=element_text(size=8),axis.title=element_text(size=8), axis.line=element_line(size=0.1),plot.title=element_text(size=8, hjust=0.5, face="bold"))
distplot

distGrob <- ggplotGrob(distplot)



 axes2 <- c("coral", "Complexity", "simpD", "log_biomass","algae", "turf", "f.richness", "P.herb", "P.carn", "alg_ratio", "FsimpD", "comp1", "comp2", "acro") # ORIGINAL REPORT
#axes1a
axes2 <- axes1a


axes.list <- list(axes2, six)

#df1$comp2 <- ifelse(df1$Data=="LTMP", -df1$comp2, df1$comp2)

df1x <- df1#[!(df1$Data=="Coral Sea" & df1$Marine.Park=="GBR"),]

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

vectors$lab2 <- labs$lab2[match(vectors$lab, labs$lab)]
vectors$data2 <- ifelse(vectors$data=="LTMP","LTMP", "CSMP")

mag1 <- diff(range(pca$x[,"PC1"])) /2
mag2 <- diff(range(pca$x[,"PC2"])) /2


unique(vectors$lab)
vectors$PC2b <- vectors$PC2r
vectors$PC2b[vectors$data=="LTMP" &  vectors$axes=="long" & vectors$lab=="P.herb"] <- vectors$PC2b[vectors$data=="LTMP"  &  vectors$axes=="long" & vectors$lab=="P.herb"]  +0.05
vectors$PC2b[vectors$data=="LTMP" &  vectors$axes=="long" & vectors$lab=="alg_ratio"] <- vectors$PC2b[vectors$data=="LTMP"  &  vectors$axes=="long" & vectors$lab=="alg_ratio"]  +0.05
vectors$PC2b[vectors$data=="LTMP" &  vectors$axes=="long" & vectors$lab=="covsqrt"] <- vectors$PC2b[vectors$data=="LTMP"  &  vectors$axes=="long" & vectors$lab=="covsqrt"]  -0.025
vectors$PC2b[vectors$data=="LTMP" &  vectors$axes=="long" & vectors$lab=="log_biomass"] <- vectors$PC2b[vectors$data=="LTMP"  &  vectors$axes=="long" & vectors$lab=="log_biomass"]  -0.05
vectors$PC2b[vectors$data2=="CSMP" &  vectors$axes=="long" & vectors$lab=="log_biomass"] <- vectors$PC2b[vectors$data2=="CSMP"  &  vectors$axes=="long" & vectors$lab=="log_biomass"]  -0.05
vectors$PC2b[vectors$data2=="CSMP" &  vectors$axes=="long" & vectors$lab=="Complexity"] <- vectors$PC2b[vectors$data2=="CSMP"  &  vectors$axes=="long" & vectors$lab=="Complexity"]  +0.05
vectors$PC2b[vectors$data2=="CSMP" &  vectors$axes=="long" & vectors$lab=="acrosqrt"] <- vectors$PC2b[vectors$data2=="CSMP"  &  vectors$axes=="long" & vectors$lab=="acrosqrt"]  +0.05
vectors$PC2b[vectors$data2=="CSMP" &  vectors$axes=="long" & vectors$lab=="simpD"] <- vectors$PC2b[vectors$data2=="CSMP"  &  vectors$axes=="long" & vectors$lab=="simpD"]  -0.05
vectors$PC2b[vectors$data2=="CSMP" &  vectors$axes=="long" & vectors$lab=="P.carn"] <- vectors$PC2b[vectors$data2=="CSMP"  &  vectors$axes=="long" & vectors$lab=="P.carn"]  -0.05

twovecs <- ggplot()+geom_segment(data=vectors, aes(x=0, xend=PC1*mag1*0.9, y=0, yend=PC2r*mag2*0.9, col=axes), arrow=arrow(length=unit(1,"mm")))+
geom_text(data=vectors, aes(x=PC1*mag1, PC2b*mag2, label=lab2, col=axes), hjust=ifelse(vectors$PC1>0, 0, 1), size=2.3, fontface="italic")+
facet_wrap(~data2)+
scale_colour_manual(values=c("grey", "black"))+
guides(col="none")+
theme_void()+
coord_cartesian(xlim=c(-3.5,4), ylim=c(-3, 3.5))
twovecs


vecslong2 <- melt(subset(vectors, select=-c(PC2r, PC2b)), id.var=c("lab", "max", "axes", "data", "lab2", "data2"))
head(vecslong2)
unique(vecslong2$variable)
vecslong2$abs <- abs(vecslong2$value)
vmax2 <- merge(aggregate(abs ~ lab+axes+data, max, data = vecslong2), vecslong2)
vecslong2$max2 <- vmax2$variable[match(vecslong2$lab, vmax2$lab)]
vecslong2$maxval <- ifelse(vecslong2$max2==vecslong2$variable, vecslong2$value, 0)
head(vecslong2)

ggplot()+
geom_bar(data=vecslong2[vecslong2$axes=="long",], aes(y=lab, x=value, fill=data), stat="identity", position="dodge", alpha=0.5)+
#geom_bar(data=vecslong2[vecslong2$axes=="long",], aes(y=lab, x=maxval, group=data), stat="identity", position="dodge", fill=NA, col="black")+
facet_grid(data~variable, scales="free_y")


# ALIGN THE COORDINATES (sites) in long and short PCAs
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

cback$data2 <- ifelse(cback$data=="LTMP", "LTMP", "CSMP")

predplot <- ggplot(data=cback[cback$variable %in% c("PC1", "PC2"),], aes(short, long, col=data2))+
geom_point(shape=21, size=1, stroke=0.25)+
#geom_point(shape=21)+
facet_wrap(~variable, scales="free")+
geom_smooth(method="lm", size=0.8)+
labs(x="PC coordinates - 6 metrics", y="PC coordinates - 14 metrics")+
scale_colour_manual(values=c("black", "red"))+
lims(x=c(-6, 6), y=c(-7, 7))+
theme_classic()+theme(axis.text=element_text(size=8),axis.title=element_text(size=8), axis.line=element_line(size=0.1), legend.title=element_blank(), legend.text=element_text(size=8),legend.position=c(0.1, 0.85), legend.key.size=unit(1, "mm"),strip.background=element_blank(), legend.background=element_blank())
predplot

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


predplot2 <- ggplot(preds, aes(var, abs(slp), fill=data, ymin=abs(low), ymax=abs(upp)))+
geom_bar(stat="identity", position=position_dodge())+
geom_errorbar( position=position_dodge(width=0.9), width=0.1)+
#coord_flip()+
scale_fill_manual(values=c("black", "red"))+
theme_classic()
predplot2

predplot3 <- ggplot(preds, aes(var, abs(slp), fill=data, ymin=abs(low), ymax=abs(upp)))+
geom_errorbar( width=0.1, stroke=0.5)+
geom_point(stat="identity", shape=21)+
coord_flip()+
lims(y=c(0,1.1))+
guides(fill="none")+
labs(y="Slope of PCA axes (6 vs 14 metrics)",x="")+
scale_fill_manual(values=c("black", "red"))+
#theme_minimal()+
theme_classic()+
theme(axis.text=element_text(size=8),axis.title=element_text(size=8), axis.line.x=element_line(size=0.1), axis.line.y=element_blank(), legend.title=element_blank(), legend.text=element_text(size=8),legend.position=c(0.1, 0.85), legend.key.size=unit(1, "mm"),strip.background=element_blank(), legend.background=element_blank())
predplot3


plot_grid(
plot_grid( distGrob, varplot, ncol=1),
plot_grid(twovecs, plot_grid(predplot+guides(col='none'), predplot2+guides(fill='none'), rel_widths=c(2, 1)), ncol=1),
rel_widths=c(0.5, 1))

plot_grid(
plot_grid(
plot_grid( distGrob, varplot, nrow=1, rel_widths=c(1, 0.7), labels=c("a", "b"), label_size=9),
NULL,
plot_grid(predplot, predplot3+guides(fill='none')+coord_flip(), rel_heights=c(2, 1), ncol=1, labels=c("d", ""), label_size=9),
ncol=1, rel_heights=c(1, 0.1, 1.4)),
plot_grid(NULL, twovecs+facet_wrap(~data2, ncol=1), ncol=1, rel_heights=c(0.1, 1),labels=c("c", ""), label_size=9, hjust=-10), rel_widths=c(1.5, 1))+
draw_text("PCA vectors\n(6 vs 14 metrics)", 0.8, 0.97, size=8, fontface="bold")+
draw_text("PCA coordinates (6 vs 14 metrics)", 0.33, 0.57, size=8, fontface="bold")

######################################################
#---------------------------------------------# variance components

# https://stats.stackexchange.com/questions/460247/nested-anova-and-assigning-proportion-of-variation


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
vars$lab <- labs$lab2[match(vars$y, labs$lab)]
vars$group <- labs$group[match(vars$y, labs$lab)]

#ggplot(vars, aes(grp, p, fill=as.factor(Year)))+geom_bar(stat="identity", position="dodge")
vars$group <- factor(vars$group, levels=c("coral", "fish", "algae","complexity"))
vars$scale <- ifelse(vars$grp=="SectorRegion", "Region", ifelse(vars$grp=="Reef:SectorRegion", "Reef", ifelse(vars$grp=="Site:(Reef:SectorRegion)", "Site", ifelse(vars$grp=="Residual", "Transect",NA))))

vars$scale  <- factor(vars$scale , levels=c("Region", "Reef", "Site","Transect"))
vars$line <- paste( vars$Year) #paste(vars$x, vars$Year)
ggplot(vars[!vars$y %in% c("PC1_1a", "PC2_1a", "PC3_1a"),], aes(x=scale, y=p*100,col=group, group=as.factor(Year)))+
geom_point(size=0.5)+
geom_line(aes(group=line))+
theme_classic()+guides(col="none")+
scale_colour_manual(values=c("black", "slategrey", "darkgreen", "darkred"))+
labs(y="% variance explained")+
ylim(c(0, 70))+
theme(axis.text.x=element_text(angle=45, hjust=1), strip.background=element_blank(), axis.line=element_line(size=0.1), axis.title.x=element_blank(), legend.title=element_blank())+
facet_wrap(group~lab, scales="free", nrow=3)

vars$y2 <- gsub("_1a","",vars$y)
ggplot(vars[vars$y %in% c("PC1_1a", "PC2_1a", "PC3_1a"),], aes(x=scale, y=p*100,col=group, group=as.factor(Year)))+
geom_point(size=0.5)+
geom_line(aes(group=line))+
theme_classic()+guides(col="none")+
scale_colour_manual(values=c("black", "slategrey", "darkgreen", "darkred"))+
labs(y="% variance explained")+
facet_wrap(~y2, scales="free", nrow=1)+
ylim(c(0, 70))+
theme(axis.text.x=element_text(angle=45, hjust=1), strip.background=element_blank(), axis.line=element_line(size=0.1), axis.title.x=element_blank(), legend.title=element_blank())


# other methods

#df1$SectorRegion <- as.factor(df1$SectorRegion)
#df1$Reef <- as.factor(df1$Reef)
#df1$Site <- as.factor(df1$Site)

#tab1 <- fitVCA(coral  ~ SectorRegion/Reef/Site, na.omit(df1[df1$Zone=="Slope" & df1$Year==2019, c("coral", "Site", "SectorRegion", "Reef")]))
#tab1
# varPlot(form=coral~(SectorRegion + Reef)/Zone/Site, Data=df1)

#aov1 <- aov(coral~Site+Zone+Reef+SectorRegion, data=df1)
#summary(aov1)



######################################################
#---------------------------------------------#  DEMO FIG FOR SIMPS-D
library("mgcv")


b1og <- read.csv("data/CoralSea/CSMP-GBR_2018-2024_Coral_PIT_CLEAN.03.03.2024.csv")
b1fg <- read.csv("data/CoralSea/bFG.csv")
b1og$id <- paste(b1og$Site, b1og$Zone, b1og$Transect, b1og$Year)
tcols <- c("Year", "Date", "Marine.Park", "SectorRegion", "Name", "Reef", "Site","Site_notes", "Zone", "Transect", "Depth", "Complexity", "Total.Macroalgae", "Total.HARD.CORAL", "Gradient", "Total", "id")
b1t <- b1og[,tcols]
b1 <- melt(b1og, id.var=tcols, value.name="n")
head(b1fg)
b1$group <- b1fg$group[match(b1$variable, b1fg$taxon)]
head(b1)

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

#morphs
b1$morph <- b1fg$morph8[match(b1$variable, b1fg$taxon)]
b1m <- aggregate(n~morph+id, b1[!b1$morph=="none",], sum)
head(b1m)


b1m[b1m$id %in% c1[c1$simpD == max(c1$simpD),"id"],]

maxD <- b1m[b1m$id %in% c1[c1$simpD == max(c1$simpD[c1$coral>0]),"id"][1],]
maxDlow <- b1m[b1m$id %in% c1[c1$simpD == max(c1$simpD[c1$coral <20]),"id"][1],]
maxDhigh <- b1m[b1m$id %in% c1[c1$simpD == max(c1$simpD[c1$coral > 60]),"id"][1],]
minD1 <- b1m[b1m$id %in% c1[c1$simpD == min(c1$simpD[c1$coral>10]),"id"][1],]
minD2 <- b1m[b1m$id %in% c1[c1$simpD == min(c1$simpD[c1$coral>60]),"id"][1],]

maxDs <- rbind(cbind(maxD,t="ii"), cbind(maxDlow,t="i"), cbind(maxDhigh,t="iii"))
minDs <- rbind(cbind(minD1,t="iv"), cbind(minD2,t="v"))


##############
# model cover vs SimpD

lm1 <- lm(simpD ~ coral, data=c1)
summary(lm1)
lm2 <- lm(simpD ~ poly(coral, 2), data=c1)
summary(lm2)
lm3 <- lm(simpD ~ poly(coral, 3), data=c1)
summary(lm3)
lm4 <- lm(simpD ~ poly(coral, 4), data=c1)
summary(lm4)
lm5 <- lm(simpD ~ poly(coral, 5), data=c1)
summary(lm5)
gam3 <- gam(simpD ~ s(coral, k=3), data=c1)
summary(gam3)
gam4 <- gam(simpD ~ s(coral, k=4), data=c1)
summary(gam4)
gam5 <- gam(simpD ~ s(coral, k=5), data=c1)
summary(gam5)
gam6 <- gam(simpD ~ s(coral, k=6), data=c1)
summary(gam6)
gam7 <- gam(simpD ~ s(coral, k=7), data=c1)
summary(gam7)
gam8 <- gam(simpD ~ s(coral, k=7), data=c1)
summary(gam8)
AIC(lm1,lm2,lm3, lm4,lm5,gam3, gam4, gam5, gam6, gam7, gam8)

fit <- data.frame(coral=seq(min(c1$coral), max(c1$coral), by=1))
fit$fit <- predict(gam7, fit)
fit$se <- predict(gam7, fit, se=T)$se.fit

covD <- ggplot()+
geom_point(data=c1, aes(coral, simpD), shape=21, size=0.5)+
geom_ribbon(data=fit, aes(x=coral, ymin=fit-se, ymax=fit+se), alpha=0.5)+
geom_line(data=fit, aes(coral, fit), col="red")+
geom_point(data=c1[c1$id %in% unique(maxD$id),],aes(coral, simpD), fill="#91bfdb", shape=21, size=2.5)+
geom_point(data=c1[c1$id %in% unique(maxDlow$id),],aes(coral, simpD), fill="#e0f3f8", shape=21, size=2.5)+
geom_point(data=c1[c1$id %in% unique(maxDhigh$id),],aes(coral, simpD), fill="#4575b4",, shape=21, size=2.5)+
geom_point(data=c1[c1$id %in% unique(minD1$id),],aes(coral, simpD), fill="#fee0b6", shape=21, size=2.5)+
geom_point(data=c1[c1$id %in% unique(minD2$id),],aes(coral, simpD), fill="#b35806", shape=21, size=2.5)+
lims(y=c(-0.15,1.05))+
#geom_smooth(data=c1, aes(coral, simpD))+
labs(x="% coral cover", y="Morphological\ndiversity (D)")+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8)) 
covD


head(maxDs)
head(minDs)
maxmin <- rbind(cbind(maxDs, type="high D"), cbind(minDs, type="low D"))

highlow <- ggplot(maxmin, aes(x=morph, y=n, fill=t))+
geom_bar(stat="identity", position="dodge", col="black", size=0.05)+
facet_wrap(~type, ncol=1, scales="free_y")+
scale_fill_manual(values=c("#e0f3f8", "#91bfdb", "#4575b4","#fee0b6", "#b35806"))+
scale_y_continuous(expand=c(0,0))+
ylab("% cover")+theme_classic()+theme(axis.line=element_line(size=0.1),axis.title.x=element_blank(), axis.title.y=element_text(size=8), axis.text.y=element_text(size=8), axis.text.x=element_text(size=7, angle=45, hjust=1), plot.title=element_text(hjust=0.5, size=8), strip.background=element_blank())
highlow

highD <- ggplot(maxDs[!maxDs$morph %in% c("Acr. Bushy"),] , aes(y=morph, x=n, fill=t))+
geom_bar(stat="identity", position="dodge", col="black", size=0.05)+
ggtitle("high D")+
scale_x_continuous(expand=c(0,0), breaks=c(0, 5, 10, 15), limits=c(0, 20))+
facet_wrap(~t)+
#scale_fill_manual(values=c("#d8daeb", "#998ec3", "#542788"))+
scale_fill_manual(values=c("#e0f3f8", "#91bfdb", "#4575b4"))+
xlab("% cover")+theme_classic()+theme(axis.line=element_line(size=0.1),axis.title.y=element_blank(), axis.title.x=element_text(size=8), axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), 
plot.title=element_text(hjust=0.5, size=8), strip.text=element_blank(), strip.background=element_blank())
highD


lowD <- ggplot(minDs[!minDs$morph %in% c("Acr. Bushy"),] , aes(y=morph, x=n, fill=t))+
geom_bar(stat="identity", position="dodge", col="black", size=0.05)+
ggtitle("low D")+
facet_wrap(~t)+
scale_x_continuous(expand=c(0,0))+
scale_fill_manual(values=c("#fee0b6", "#b35806"))+
xlab("% cover")+theme_classic()+theme(axis.line=element_line(size=0.1),axis.title.y=element_blank(),axis.title.x=element_text(size=8), axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), 
plot.title=element_text(hjust=0.5, size=8), strip.text=element_blank(), strip.background=element_blank())
lowD 


demoplot2 <- plot_grid(covD, highlow+guides(fill="none"), rel_widths=c(1,0.5))
demoplot2

demoplot <- plot_grid(
highD+guides(fill="none"),
covD, 
plot_grid(lowD+guides(fill="none"), NULL, rel_widths=c(1, 0.25)), 
ncol=1, rel_heights=c(1, 2, 1), labels=c("a", "b", "c"), label_size=10)
demoplot


######################################################
#---------------------------------------------#  hard coral FD

df2 <- read.csv("data/output/coraltraits.csv")
b1t <- read.csv("data/output/coralFD.csv")

unique(b1t$Marine.Park)


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

maxFD <- b1[b1$id %in% b1t[b1t$FDis %in% max(b1t$FDis, na.rm=T),"id"][1],]
head(maxFD)

siteFD <- aggregate(FDis~Site, b1t, mean)
siteFD <- siteFD[order(siteFD$FDis),]
head(siteFD)
tail(siteFD)
max(siteFD$FDis, na.rm=T)

maxFDs <-  b1[b1$Site %in% siteFD[nrow(siteFD)-1,"Site"],] #temp[temp$Site %in% siteFD[siteFD$FDis %in% max(siteFD$FDis, na.rm=T),"Site"][1],]
maxplot <- aggregate(n~Site+dim1+dim2+variable, maxFDs, sum )
minFDs <- b1[b1$Site %in% siteFD[1,"Site"],]  #temp[temp$Site %in% siteFD[siteFD$FDis %in% min(siteFD$FDis, na.rm=T),"Site"][1],]
minplot <- aggregate(n~Site+dim1+dim2+variable, minFDs, sum )
head(minplot)


FDdemo <- ggplot() +
  geom_point(data=df2, aes(dim1, dim2), size=0.15,  col="grey")+
   stat_ellipse(data=df2, aes(dim1, dim2, group=factor(kclust)), col="grey", fill="grey", geom="polygon", alpha=0.2, size=0.2)+
  geom_segment(data=maxplot[maxplot$n > 0,], aes(x=0, xend=dim1, y=0, yend=dim2), col="black", size=0.3)+ #"#5ab4ac"
geom_point(data=maxplot[maxplot$n > 0,], aes(dim1, dim2, size=n), fill="black", shape=21, stroke=0.5, col="white")+
geom_segment(data=minplot[minplot$n > 0,], aes(x=0, xend=dim1, y=0, yend=dim2), col="red", linetype="dotted")+
  geom_point(data=minplot[minplot$n > 0,], aes(dim1, dim2, size=n), col="red", shape=21)+
scale_radius(range=c(1,6))+guides(size="none")+
geom_text(data=NULL, aes(x=-0.4, y=-0.24, label="High FDis"), col="black", size=2.5, hjust=0)+
geom_text(data=NULL, aes(x=-0.4, y=-0.28, label="Low FDis"), col="red", size=2.5, hjust=0)+
labs(x="PCoA 1", y="PCoA 2")+
theme_bw()+theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=6), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
FDdemo



linkmod <- summary(lm(simpD~sqrt(FDis), b1t))
linkmod
rsq <- round(linkmod$r.squared, 2)



FDtoD <- ggplot(b1t, aes(simpD, FDis))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
labs(x="Morph. diversity (D)", y="Functional dispersion")+
annotate(geom = 'text', x = 0.8, y = 0,    label = paste("R^2 == ", rsq), parse = TRUE, size=2.5, col="blue")+
scale_y_sqrt()+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8)) 


bquote(italic(R)^2==.(rsq))

FDdemoplot <- plot_grid(
plot_grid(FDdemo, NULL, ncol=1, rel_heights=c(1, 0.7)),
FDtoD, rel_widths=c(1.3, 2), labels=c("", ""), label_size=10)
FDdemoplot


tspace <- readPNG("groups.png")
tspace<-rasterGrob(tspace, interpolate=TRUE)

FDdemoplot2 <- plot_grid(
plot_grid(NULL, FDdemo, NULL, nrow=1, rel_widths=c(0.5, 1,0.5)),
FDtoD, ncol=1, rel_heights=c(1, 2), labels=c("d", "e"), label_size=10)
FDdemoplot2

demo <- plot_grid(demoplot2, FDdemoplot, ncol=1,labels=c("a", "b"), label_size=10)+
draw_line(x=c(0.25, 0.25), y=c(0.19, 0.22))+
draw_line(x=c(0.25, 0.32), y=c(0.19, 0.19), arrow=arrow(length=unit(1, "mm")))+
draw_line(x=c(0.65, 0.65), y=c(0.69, 0.65))+
draw_line(x=c(0.7, 0.65), y=c(0.67, 0.67))+
draw_line(x=c(0.6, 0.6), y=c(0.86, 0.90))+
draw_line(x=c(0.6, 0.65), y=c(0.88, 0.88))
demo

# other script


demo2 <- plot_grid(demoplot, NULL, FDdemoplot2, nrow=1, rel_widths=c(1, 0.1, 1))+
draw_line(x=c(0.2, 0.17), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.27, 0.255), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.37, 0.365), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.13, 0.17), y=c(0.36, 0.2), size=0.1)+
draw_line(x=c(0.4, 0.33), y=c(0.36, 0.2), size=0.1)+
draw_line(x=c(0.65, 0.58), y=c(0.8, 0.8), size=0.1)+
draw_line(x=c(0.58, 0.58), y=c(0.8, 0.7), size=0.1, arrow=arrow(length=unit(1, "mm")))
demo2


######################################################
#---------------------------------------------#  size structure
colz <- c("#b2182b", "#f4a582", "#92c5de", "#2166ac")
names(colz) <- c("s.neg", "ns.neg", "ns.pos", "s.pos")
s1$Date[s1$id=="Saumarez 5 Crest 3 2020"] <- "18/02/2020" # check
s1$Date[s1$id=="Kenn 3 Slope 1 2018"] <- "10/12/2018" # check

s1$Reef[s1$Reef=="Sweetlips"] <- "Sweetlip"
s1$Site <- gsub("Sweetlips", "Sweetlip", s1$Site)
s1$Reef[s1$Reef=="Whitetip"] <- "White tip"
s1$Site <- gsub("Whitetip", "White tip", s1$Site)
s1$Site <- gsub("5A", "5a", s1$Site)
s1$Site <- gsub("3A", "3a", s1$Site)
s1$Site <- gsub("Herlad 6", "Herald 6", s1$Site)
s1$Site <- gsub("Flinders  6", "Flinders 6", s1$Site)
s1$Site <- gsub("Chilcot 2", "Chilcott 2", s1$Site)
s1$Zone <- gsub("crest", "Crest", s1$Zone)
s1$Zone <- gsub("slope", "Slope", s1$Zone)
unique(s1$Site)[order(unique(s1$Site))]

s1 <- read.csv("data/CoralSea/CSMP-GBR_2018-2024_Coral_Sea_Health_and_Recruitment_Surveys.csv")
s1$id <- paste(s1$Site, s1$Zone,s1$Transect, s1$Year)
length(unique(s1$id))
scols <- c("Year", "Date",  "SectorRegion", "Reef", "Site", "Zone", "Transect", "id")
s1t <- unique(s1[,scols])
s1FG <- read.csv("data/CoralSea/sFG.csv")
s1$pit <- s1FG$pit_name[match(s1$Genus, s1FG$taxon)]
s1$group <- s1FG$group[match(s1$Genus, s1FG$taxon)]
head(s1)

s1s <- melt(s1[,c("id","Reef","Site", "Zone","pit", "group", "Total.recruits", "X6.10cm", "X11.20cm","X21.40cm", "X41.60cm", "X.60cm")], id=c("id","Reef","Site","Zone", "pit", "group"))
head(s1)
s1s$size <- ifelse(s1s$variable=="Total.recruits", "00 - 05", ifelse(s1s$variable=="X6.10cm", "06 - 10", ifelse(s1s$variable=="X11.20cm", "11 - 20", ifelse(s1s$variable=="X21.40cm", "21 - 40", ifelse(s1s$variable=="X41.60cm", "41 - 60", ifelse(s1s$variable=="X.60cm", "60+", NA))))))
head(s1s)
s1coral <- s1s[s1s$group=="HC",]
#s1coral <- aggregate(value~., s1coral, sum)
s1coral <- aggregate(value~., subset(s1coral, select=-c(pit)), sum) # totals..
head(s1coral)

cover <- aggregate(value~id, s1coral, sum)
s1coral$total <- cover$value[match(s1coral$id, cover$id)]
s1coral$p <- s1coral$value/s1coral$total

siteS <- aggregate(p~size+Reef+Zone, s1coral, mean)
siteS$gp <- paste(siteS$Reef, siteS$Zone)

sizeplot <- ggplot(siteS, aes(size, p))+
geom_rect(data=NULL, aes(xmin=1, xmax=2, ymin=-Inf, ymax=Inf), fill="grey")+
geom_rect(data=NULL, aes(xmin=5.7, xmax=6.3, ymin=-Inf, ymax=Inf), fill="grey")+
geom_text(data=NULL, aes(x=1.5, y=0.9, label="small"), size=3)+
geom_text(data=NULL, aes(x=6, y=0.9, label="large"), size=3)+
geom_segment(data=NULL, aes(x=6.35, xend=6.5, y=0.9, yend=0.9), arrow=arrow(length=unit(1, "mm")))+
geom_line(aes(group=gp), size=0.1)+ #position="stack",
theme_classic()+
labs(x="Diameter (cm)", y="Proportion of corals")+
theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8))
sizeplot 

p10x <- s1coral[s1coral$size %in% c("00 - 05", "06 -10"),]
head(p10x)
p60x <- s1coral[s1coral$size %in% c("60+"),]
head(p10x)

#---------------------------------------------# reef-level relationships

#summary(lm(frich~coral_cover, tdf))
#summary(lm(frich~coral_cover*A_SECTOR, tdf))


bimod <- glm(p10~coral, family="quasibinomial", data=c1)
summary(bimod)
bifit <- data.frame(coral = seq(min(c1$coral, na.rm=T), max(c1$coral, na.rm=T), 1))
bifit$fit <- predict(bimod, bifit, type="response")

ggplot()+geom_line(data=bifit, aes(coral, fit))+geom_point(data=c1, aes(coral, p10), shape=21)


# lots of reefs with no small corals.. 

# reef variation trans/sites/years.. 

hist(((c1$coral)/max(c1$coral)))
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
temp <- c1
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

ggplot(reefs, aes(x=slp))+
geom_density(alpha=0.5, fill="grey")+
#geom_segment(data=avs, aes(x=slp, xend=cod, y=Inf, yend=-Inf))+
theme_minimal()+
lims(x=c(-0.5, 0.7))+
labs(x="coefficient")+
geom_vline(xintercept=0, linetype="dotted")+facet_grid(x~y)


p10dat <- c1[,c("p10", "coral", "Reef", "Zone")]
slpdat <- reefs[reefs$y=="p10",]
p10dat[,c("slp", "dir", "sig", "sigdir")] <- slpdat[match(p10dat$Reef, slpdat$Reef), c("slp", "dir", "sig", "sigdir")]
p10fit <- bifits[bifits$resp=="p10",]
p10plot <- ggplot(p10dat[p10dat$coral > 0 & p10dat$Zone %in% c("Crest", "Slope"),], aes( coral, p10))+
geom_point(shape=21, aes(col=sigdir), size=0.1)+
#scale_y_sqrt()+scale_x_sqrt()+
#geom_smooth(se=F, method="lm", aes(group=Reef, col=sigdir), size=0.5)+
geom_line(data=p10fit, aes(x, fit, group=Reef, col=sigdir))+
#geom_smooth(col="black")+
#facet_wrap(~Zone)+
labs(x="% coral cover", y="Proportion corals < 10cm")+
scale_colour_manual(values=colz)+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.2), axis.text=element_text(size=8), axis.title=element_text(size=8))
p10plot

p60dat <- c1[,c("p60", "coral", "Reef", "Zone")]
slpdat <- reefs[reefs$y=="p60",]
p60dat[,c("slp", "dir", "sig", "sigdir")] <- slpdat[match(p60dat$Reef, slpdat$Reef), c("slp", "dir", "sig", "sigdir")]
p60fit <- bifits[bifits$resp=="p60",]
p60plot <- ggplot(p60dat[p60dat$coral > 0 & p10dat$Zone %in% c("Crest", "Slope"),], aes(coral, p60))+
geom_point(shape=21, aes(col=sigdir), size=0.1)+
#scale_y_sqrt()+scale_x_sqrt()+
#geom_smooth(se=F, method="lm", aes(group=Reef, col=sigdir), size=0.5)+
#geom_smooth(col="black")+
geom_line(data=p60fit, aes(x, fit, group=Reef, col=sigdir))+
#facet_wrap(~Reef)+
labs(x="% coral cover", y="Proportion corals > 60cm")+
scale_colour_manual(values=colz)+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.2), axis.text=element_text(size=8), axis.title=element_text(size=8))
p60plot 


sizecov <- plot_grid(NULL,p10plot+guides(col="none"), p60plot+guides(col="none"),  nrow=1, rel_widths=c(0.4, 1,1), labels=c("", "a", 'b'), label_size=10)+
draw_label("Negative (p < 0.05)", x=0.01, y=0.8, colour=colz[1], size=7, hjust=0)+
draw_label("Negative (p > 0.05)", x=0.01, y=0.75, colour=colz[2], size=7, hjust=0)+
draw_label("Positive (p > 0.05)", x=0.01, y=0.7, colour=colz[3], size=7, hjust=0)+
draw_label("Positive (p < 0.05)", x=0.01, y=0.65, colour=colz[4], size=7, hjust=0)
sizecov


sizecov2 <- plot_grid(p10plot+guides(col="none"), p60plot+guides(col="none"),ncol=1,  labels=c("c", 'd'), label_size=10)+
draw_label("Negative (p < 0.05)", x=0.66, y=0.95, colour=colz[1], size=7, hjust=0)+
draw_label("Negative (p > 0.05)", x=0.66, y=0.92, colour=colz[2], size=7, hjust=0)+
draw_label("Positive (p > 0.05)", x=0.66, y=0.89, colour=colz[3], size=7, hjust=0)+
draw_label("Positive (p < 0.05)", x=0.66, y=0.86, colour=colz[4], size=7, hjust=0)
sizecov2


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

plot_grid(plot_grid(sizeplot, s.pca, NULL, labels=c("a", "b", ""), label_size=10,rel_widths=c(1,0.9,0.2), nrow=1),
sizecov, ncol=1, labels=c("", "c"), label_size=10)



plot_grid(s.pca, sizecov, labels=c("a", "", ""), label_size=10, rel_widths=c(1,2))

plot_grid(demo, sizecov2, rel_widths=c(1.5, 1))


######################################################
#---------------------------------------------#  refuge vol

head(c1)

linkmod2 <- summary(lm(sqrt(shelt)~coral, c1))
linkmod2
rsq2 <- round(linkmod2$r.squared, 2)

FDtoD <- ggplot(b1t, aes(simpD, FDis))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
labs(x="Morph. diversity (D)", y="Functional dispersion")+
annotate(geom = 'text', x = 0.8, y = 0,    label = paste("R^2 == ", rsq), parse = TRUE, size=2.5, col="blue")+
scale_y_sqrt()+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8)) 

refplot <- ggplot(c1, aes(x=coral, y=shelt))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
scale_y_sqrt()+
annotate(geom = 'text', x = 75, y = 0,    label = paste("R^2 == ", rsq2), parse = TRUE, size=2.5, col="blue")+
labs(x="% coral cover", y=expression(Shelter~Volume~(dm^3)))+
scale_colour_manual(values=colz)+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.2), axis.text=element_text(size=8), axis.title=element_text(size=8))
refplot


plot_grid(NULL,
plot_grid(demo, NULL, plot_grid(sizecov2, refplot, ncol=1, rel_heights=c(1.8, 1), labels=c("", "e"), label_size=10), rel_widths=c(1.8,0.15, 1), nrow=1),
rel_heights=c(0.05,1), ncol=1)+
draw_label("Morphological and functional diversity", x=0.21, y=0.97,  size=9, hjust=0, fontface="bold")+
draw_label("Size structure", x=0.8, y=0.97,  size=9, hjust=0, fontface="bold")

plot_grid(sizecov, refplot, rel_widths=c(2.23, 1), labels=c("", "c"), label_size=10)




######################################################
#---------------------------------------------#  benthos - fish


ggplot(df1, aes(coral, FsimpD2))+geom_point()+
geom_smooth(method="lm", aes(group=Reef), se=F, size=0.2)+
facet_wrap(~Data)+
labs(x="% coral cover", y="Feeding diversity (3 groups)")+
theme_classic()

######################################################
#---------------------------------------------#  benthos - fish


df2$coral1 <- df2$coral + 1
mod16 <- nls(f.richness ~ b*log(coral1) + c, data=na.omit(df2), start=list(b=48.6, c=-21.6))
summary(mod16)



new.data <- data.frame(coral1=seq(1, 100, 1))
new.data$fit <- predict(mod16, new.data)
head(new.data)

fit16 <- cbind(new.data, data.frame(predFit(mod16, newdata = new.data, interval = "confidence", level= 0.95)))


ggplot(df2, aes(coral, f.richness))+geom_point(size=0.5)+
#facet_wrap(~Data)+
geom_line(data=new.data, aes(coral1, fit), col="red")+
geom_smooth(method="lm", formula=y~poly(x, 3), col="green")+
geom_smooth()+
theme_classic()



head(df1)
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



ggplot(df1, aes(coral, f.biomass, col=Data))+geom_point(size=0.5)+
#facet_wrap(~Data)+
scale_y_log10()+
geom_smooth()+
theme_classic()

head(df1)

hist(df1$f.richness)

exp(14)/1000
exp(11.8)/1000

reefs <- NULL
Ys <- c("f.richness", "log_biomass", "FsimpD", 'P.herb')
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



colz <- c("#b2182b", "#f4a582", "#92c5de", "#2166ac")
names(colz) <- c("s.neg", "ns.neg", "ns.pos", "s.pos")

head(reefs)
reefs$x2 <- labs$lab2[match(reefs$x, labs$lab)]
reefs$y2 <- labs$lab2[match(reefs$y, labs$lab)]

reefs$cod <- reefs$slp * sqrt(reefs$rsq) / abs(reefs$slp)

hist(reefs$cod)

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

bars1 <- ggplot(data=data.frame(table(df1$sigdir1)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass1 <- ggplot(df1, aes(Complexity, f.richness))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
guides(col="none")+
facet_wrap(~Data)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir1), se=F, size=0.3)+
#annotation_custom( ggplotGrob(bars1),   xmin = 5.1, xmax = 6, ymin = 0, ymax =20)+
labs(x="Complexity (1-5)", y="Fish richness")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))

plot_grid(ass1, bars1, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1)

slps2 <- reefs[reefs$x=="Complexity" & reefs$y=="log_biomass",]
df1$sigdir2 <- slps2$sigdir[match(df1$Reef, reefs$reef)]
bars2 <- ggplot(data=data.frame(table(df1$sigdir2)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass2 <- ggplot(df1, aes(Complexity, log(f.biomass)))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
#scale_y_log10()+
guides(col="none")+
facet_wrap(~Data)+
#annotation_custom( ggplotGrob(bars2),   xmin = 4, xmax = 5, ymin = 3, ymax = 5)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir2), se=F, size=0.3)+
labs(x="Complexity (1-5)", y="Fish biomass")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))
ass2

plot_grid(ass2, bars2, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1)


slps3 <- reefs[reefs$x=="simpD" & reefs$y=="f.richness",]
df1$sigdir3 <- slps3$sigdir[match(df1$Reef, reefs$reef)]
bars3 <- ggplot(data=data.frame(table(df1$sigdir3)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass3 <- ggplot(df1, aes(simpD, f.richness))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
#scale_y_log10()+
facet_wrap(~Data)+
guides(col="none")+
#annotation_custom( ggplotGrob(bars3),   xmin = 0.75, xmax = 1, ymin = 0, ymax = 20)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir3), se=F, size=0.3)+
labs(x="Morphological diversity (D)", y="Fish richness")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))
ass3

plot_grid(ass3, bars3, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1)

slps4 <- reefs[reefs$x=="coral" & reefs$y=="log_biomass",]
df1$sigdir4 <- slps4$sigdir[match(df1$Reef, reefs$reef)]
bars4 <- ggplot(data=data.frame(table(df1$sigdir4)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass4 <- ggplot(df1, aes(coral, log(f.biomass)))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
#scale_y_log10()+
facet_wrap(~Data)+
guides(col="none")+
#annotation_custom( ggplotGrob(bars3),   xmin = 0.75, xmax = 1, ymin = 0, ymax = 20)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir4), se=F, size=0.3)+
labs(x="% coral cover", y="Fish biomass")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))
ass4


slps5 <- reefs[reefs$x=="coral" & reefs$y=="f.richness",]
df1$sigdir5 <- slps5$sigdir[match(df1$Reef, reefs$reef)]
bars5 <- ggplot(data=data.frame(table(df1$sigdir5)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass5 <- ggplot(df1, aes(coral, f.richness))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
#scale_y_log10()+
facet_wrap(~Data)+
guides(col="none")+
#annotation_custom( ggplotGrob(bars3),   xmin = 0.75, xmax = 1, ymin = 0, ymax = 20)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir5), se=F, size=0.3)+
labs(x="% coral cover", y="Fish richness")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))
ass5

slps6 <- reefs[reefs$x=="coral" & reefs$y=="FsimpD",]
df1$sigdir6 <- slps6$sigdir[match(df1$Reef, reefs$reef)]
bars6 <- ggplot(data=data.frame(table(df1$sigdir6)), aes(reorder(Var1, -Freq), Freq, fill=Var1))+
scale_fill_manual(values=colz)+geom_bar(stat="identity")+theme_void()+guides(fill="none")
ass6 <- ggplot(df1, aes(coral, FsimpD))+
geom_point(size=0.1, col="grey")+
scale_colour_manual(values=colz)+
#scale_y_log10()+
facet_wrap(~Data)+
guides(col="none")+
#annotation_custom( ggplotGrob(bars3),   xmin = 0.75, xmax = 1, ymin = 0, ymax = 20)+
geom_smooth(method="lm", aes(group=Reef, col=sigdir6), se=F, size=0.3)+
labs(x="% coral cover", y="Feeding diversity (D)")+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.1), axis.text=element_text(size=7), axis.title=element_text(size=8))
ass6


plot_grid(
plot_grid(ass1, bars1, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass2, bars2, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass3, bars3, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass4, bars4, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass5, bars5, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass6, bars6, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
labels=c("a", "b", "c", "d", "e", "f"), label_size=9, nrow=3)



	#slp.lin <- fixef(mod.lin)["scalelogtval"]
	#conf.lin <- confint(mod.lin)["scalelogtval",]
	#AICs <- AICc(mod.sq, mod.lin)[,2]
	#aic.sq <- AICs[1]
	#aic.lin <- AICs[2]
	best <- ifelse(aic.sq > aic.lin, "lin", "quad")
	p.lin <- coef(summary(mod.lin))["scalelogtval","Pr(>|t|)"]
	p.sq <- coef(summary(mod.sq))["scalelogtval","Pr(>|t|)"]
	p.sq.quad <- coef(summary(mod.sq))["scalelogtval_sq","Pr(>|t|)"]
	best <- ifelse(aic.sq > aic.lin, "lin", "quad")
	# warn <- paste( warnings(), collapse = '') 
lmer_mods <- rbind(lmer_mods, data.frame(modID=lmer_ids[i], slp.sq, slp.upp.sq = slp.conf.sq[1], slp.low.sq=slp.conf.sq[2], quad.sq,  quad.upp.sq = quad.conf.sq[1], quad.low.sq = quad.conf.sq[2], slp.lin, slp.upp.lin = conf.lin[1], slp.low.lin=conf.lin[2], aic.sq, aic.lin, best, p.lin, p.sq, p.sq.quad))












