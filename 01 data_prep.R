
rm(list = ls())

library("ggplot2")
library("cowplot")
library("viridis")
library("reshape2")
library("lubridate")
library("vegan")
library("psych")

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
nrow(b1t)

b1 <- melt(b1og, id.var=tcols, value.name="n")
b1$group <- b1fg$group[match(b1$variable, b1fg$taxon)]
head(b1)
b1fg
unique(b1fg$taxon)
unique(b1fg$morph7)

length(unique(b1$id))

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
###
b1alg <- aggregate(n ~ id, b1[b1$group=="MA",], sum)
b1t$algae <- b1alg$n[match(b1t$id, b1alg$id)]
b1t$match <- ifelse(b1t$algae==b1t$Total.Macroalgae, "yes", "no")
table(b1t$match)
###
b1soc <- aggregate(n ~ id, b1[b1$group=="SoC",], sum)
b1t$soft <- b1soc$n[match(b1t$id, b1soc$id)]
###
b1turf <- aggregate(n ~ id, b1[b1$group=="Turf",], sum)
b1t$turf <- b1turf$n[match(b1t$id, b1turf$id)]
###
b1cca <- aggregate(n ~ id, b1[b1$group=="CCA",], sum)
b1t$cca <- b1cca$n[match(b1t$id, b1cca$id)]
###
b1sand <- aggregate(n ~ id, b1[b1$group=="Sand",], sum)
b1t$sand <- b1sand$n[match(b1t$id, b1sand$id)]
###
b1tother <- aggregate(n ~ id, b1[b1$group=="Other",], sum)
b1t$other <- b1tother$n[match(b1t$id, b1tother$id)]
head(b1t)

# divide by total to get cov
b1comp <- melt(b1t[,c("coral", "algae", "turf", "cca", "other", "sand","soft", "id", "total")], id.var=c("id","total"))
b1comp$cov <- (b1comp$value/b1comp$total)*100
#ggplot(b1comp, aes(id, cov, fill=variable))+geom_bar(stat="identity")+theme(axis.text.x=element_blank(), axis.ticks=element_blank())

# Div index
b1$morph <- b1fg$morph7[match(b1$variable, b1fg$taxon)]
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

# Acro 
b1$acro <- b1fg$acro[match(b1$variable, b1fg$taxon)]
b1ac <- aggregate(n ~ id, b1[b1$acro=="yes",], sum)
b1t$acro <- b1ac$n[match(b1t$id, b1ac$id)] / b1t$coral


######################################################
#---------------------------------------------#  MULTIDIMENSIONAL?

## multidimensional composition...
mdsdat0 <- acast(b1m[b1m$coral>0,], id~morph, value.var="rel_cover") #n
pca <- prcomp(sqrt(mdsdat0), scale=T, center=T)
vecs <- data.frame(pca$rotation[,c("PC1", "PC2", "PC3", "PC4")], lab=rownames(pca$rotation))
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100	
#biplot(pca)
b1t[,c("comp1", "comp2")] <- pca$x[match(b1t$id, rownames(pca$x)), c("PC1", "PC2")]

ggplot()+geom_point(data=b1t, aes(comp1, comp2, col=simpD))+
scale_colour_viridis()+
labs(x=paste("PCA 1 (", vars[1], "%)", sep=""),y=paste("PCA 2 (", vars[2], "%)", sep=""))+
geom_segment(data=vecs, aes(x=0, xend=PC1*3, y=0, yend=PC2*3))+
geom_text(data=vecs, aes(x=PC1*4, y=PC2*4, label=lab))

ggplot(b1t, aes(comp1, simpD))+geom_point()


######################################################
#---------------------------------------------#  Coral Sea (size)

s1 <- read.csv("data/CoralSea/CSMP-GBR_2018-2024_Coral_Sea_Health_and_Recruitment_Surveys.csv")
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

s1$id <- paste(s1$Site, s1$Zone,s1$Transect, s1$Year)
length(unique(s1$id))

unique(s1$Reef)

s1FG <- read.csv("data/CoralSea/sFG.csv")
s1$pit <- s1FG$pit_name[match(s1$Genus, s1FG$taxon)]
s1$group <- s1FG$group[match(s1$Genus, s1FG$taxon)]
head(s1)



unique(s1$pit)

scols <- c("Year", "Date",  "SectorRegion", "Reef", "Site", "Zone", "Transect", "id")
s1t <- unique(s1[,scols])
length(unique(s1t$id))
nrow(s1t)
#enns <- data.frame(table(s1t$id))

head(s1)
s1s <- melt(s1[,c("id", "pit", "group", "Total.recruits", "X6.10cm", "X11.20cm","X21.40cm", "X41.60cm", "X.60cm")], id=c("id", "pit", "group"))
s1s$size <- ifelse(s1s$variable=="Total.recruits", "00 - 05", ifelse(s1s$variable=="X6.10cm", "06 - 10", ifelse(s1s$variable=="X11.20cm", "11 - 20", ifelse(s1s$variable=="X21.40cm", "21 - 40", ifelse(s1s$variable=="X41.60cm", "41 - 60", ifelse(s1s$variable=="X.60cm", "60+", NA))))))
head(s1s)


s1coral <- s1s[s1s$group=="HC",]
s1coral <- aggregate(value~., s1coral, sum)
s1coral <- aggregate(value~., subset(s1coral, select=-c(pit)), sum)

ggplot(s1coral[s1coral$id %in% unique(s1coral$id)[100:110],], aes(size, value))+geom_point()+geom_line(aes(group=id), size=0.1)

# proportions
stot <- aggregate(value~., subset(s1coral, select=-c(size, variable)), sum)
s1coral$tot <- stot$value[match(s1coral$id, stot$id)]
s1coral$p <- s1coral$value / s1coral$tot
head(s1coral)

ggplot(s1coral[s1coral$id %in% unique(s1coral$id)[100:110],], aes(size, p))+geom_point()+geom_line(aes(group=id), size=0.1)

# proportion of corals under 10
p10 <- s1coral[s1coral$size %in% c("00 - 05", "06 - 10"),]
p10 <- aggregate(p ~ id, p10, sum)
nrow(p10)
s1t$p10 <- p10$p[match(s1t$id, p10$id)]

# proportion of corals over 20 (top 3 groups)
p20 <- s1coral[s1coral$size %in% c("21 - 40", "41 - 40", "60+"),]
p20 <- aggregate(p ~ id, p20, sum)
nrow(p20)
s1t$p20 <- p20$p[match(s1t$id, p20$id)]

# proportion of corals over 60
p60 <- s1coral[s1coral$size %in% c("60+"),]
nrow(p60)
s1t$p60 <- p60$p[match(s1t$id, p60$id)]

ggplot(s1t, aes(Year, p60))+geom_point()


######################################################
#---------------------------------------------#  shelter vol

head(s1s)

# recruits as zero?
s1s$max <- ifelse(s1s$variable=="Total.recruits",0, ifelse(s1s$variable=="X6.10cm", 10, ifelse(s1s$variable=="X11.20cm", 20, ifelse(s1s$variable=="X21.40cm", 40, ifelse(s1s$variable=="X41.60cm", 60, ifelse(s1s$variable=="X.60cm", 100, NA))))))

s1s$min <- ifelse(s1s$variable=="Total.recruits",0, ifelse(s1s$variable=="X6.10cm", 6, ifelse(s1s$variable=="X11.20cm", 11, ifelse(s1s$variable=="X21.40cm", 21, ifelse(s1s$variable=="X41.60cm", 41, ifelse(s1s$variable=="X.60cm", 60, NA))))))


s1s$morph <- b1fg$morph8[match(s1s$pit, b1fg$taxon)]
unique(s1s$morph)
s1s$morph[s1s$morph=="Foliose"] <- "Acr. tabular" # assume similar
s1s$morph[s1s$morph=="Branching (other)"] <- "Branching" 
s1s$morph[s1s$morph=="Acr. branching"] <- "Branching" 
s1s$morph[s1s$morph=="Acr. bushy"] <- "Acr. tabular" # assume similar
s1s$morph[s1s$morph=="Acr. Bushy"] <- "Acr. tabular" # assume similar
s1s$morph[s1s$morph=="Bushy (other)"] <- "Columnar" # assume similar
unique(s1s[,c("pit", "morph")])


shelt <- na.omit(s1s[!s1s$morph %in% c("none", "Encrusting"),]) # assume zero
unique(shelt$morph)
head(shelt)

# log(S) = a + log(D)b # ()a = intercept, b = slope_
shelt$a <- ifelse(shelt$morph=="Massive", -10.2, ifelse(shelt$morph=="Columnar", -8.5, 
ifelse(shelt$morph=="Acr. tabular", -8.66, ifelse(shelt$morph=="Branching", -9.41, NA))))
shelt$b <- ifelse(shelt$morph=="Massive", 2.91, ifelse(shelt$morph=="Columnar", 2.74, 
ifelse(shelt$morph=="Acr. tabular", 2.83, ifelse(shelt$morph=="Branching", 3, NA))))

shelt$shelt <- exp(shelt$a + (shelt$b * log(shelt$max))) # min / max
shelt$shelt2 <- shelt$shelt * shelt$value
head(shelt)

shelt$planar <- pi*(shelt$max/2)^2 # min / max
shelt$planar2 <- shelt$planar * shelt$value

ggplot(shelt, aes(x=morph, y=shelt))+geom_boxplot()+scale_y_sqrt()

svol <- aggregate(shelt2~id, shelt, sum)
svol$planar2 <- aggregate(planar2~id, shelt, sum)$planar2


ggplot(svol, aes(shelt2))+geom_histogram()+scale_x_log10()

ggplot(svol, aes(((planar2*0.0001)/10)*100, shelt2))+geom_point()



######################################################
#---------------------------------------------#  Coral Sea (fish)

# each transect is dplit into 3 sizes (50x5, 50x4, 50*2) and depth can vary within a transect
# sometimes big/small have different days!

f1 <- read.csv("data/CoralSea/CSMP-GBR_2018-2024_Fish_Surveys_combined.csv")
f1$Species[f1$Species=="Stegostoma fasciatum "] <- "Stegostoma fasciatum"
#f1$Date <- parse_date_time(f1$Date, orders = c('dmy', 'ymd'))
f1$Reef[f1$Reef=="Diane Banks"] <- "Diane"
f1$Site <- gsub("Diane Banks", "Diane", f1$Site)
unique(f1$Site)

head(f1) # already long

f1$id <- paste(f1$Site, f1$Zone, f1$Transect, f1$Year)
length(unique(f1$id))
length(unique(f1$Species))
#f1[f1$id=="Milln 1 Crest 2 2022",]

f1[is.na(f1$N),] # two NAs. Remove?
f1 <- f1[!is.na(f1$N),]

fcols <- c("Year", "Marine.Park", "SectorRegion", "Reef", "Site", "Zone", "Transect", "id") 
# removing transect size and depth? sites notes? 

unique(f1$Site_notes)

f1t <- unique(f1[,fcols])
length(unique(f1t$id))
nrow(f1t)
#enns <- data.frame(table(f1t$id))

f1fg <- read.csv("data/CoralSea/fFGs.csv")
f1fg$Species[f1fg$Species=="Chromis retrofasciata "] <- "Chromis retrofasciata"
head(f1fg)
nrow(f1fg)

unique(f1fg$Functional.group)

# biomass = a * L^b
# https://www.fishbase.se/manual/english/FishBaseThe_LENGTH_WEIGHT_Table.htm
f1$FG <- f1fg$Functional.group[match(f1$Species, f1fg$Species)]
f1$a <- f1fg$a[match(f1$Species, f1fg$Species)]
f1$b <- f1fg$b[match(f1$Species, f1fg$Species)]
f1$biomass <- f1$a  * (f1$Size)^f1$b
#f1[f1$biomass==max(f1$biomass),]
head(f1)

# total number/mass
f1N <- aggregate(N ~ id, f1, sum)
f1t$N <- f1N$N[match(f1t$id, f1N$id)]

f1$biomassT <- f1$biomass * f1$N
f1B <- aggregate(biomass ~ id, f1, sum)
f1BT <- aggregate(biomassT ~ id, f1, sum)
f1t$biomass <- f1B$biomass[match(f1t$id, f1B$id)]
f1t$biomassT <- f1BT$biomassT[match(f1t$id, f1BT$id)]
f1$pres <- 1
f1R <- aggregate(pres~id, unique(f1[,c("id", "Species", "pres")]), sum) 
f1t$richness <- f1R$pres[match(f1t$id, f1R$id)]
head(f1t)

# max fish biomass? 
fmax <- f1t[which.max(f1t$biomass),]
fmax
f1[f1$id %in% fmax$id,c("id", "Species","Size", "N", "biomass", "biomassT")]

ggplot(f1t, aes(N, biomassT))+geom_point()+scale_y_log10()+scale_x_log10()

# hist(f1t$richness)

# herbivores?
unique(f1$FG)
f1$FG2 <- ifelse(f1$FG %in% c("Carnivore", "Piscivore", "Benthic invertivore", "Corallivore"), "Carnivore", ifelse(f1$FG %in% c("Grazer", "Scraper", "Excavator", "Browser"), "Herbivore", "Other")) # excavators = parrotfish
f1herb <- aggregate(N ~ id, f1[f1$FG2=="Herbivore",], sum)
f1t$Nherb <- f1herb$N[match(f1t$id, f1herb$id)]
f1t$Nherb[is.na(f1t$Nherb)] <- 0 
f1carn <- aggregate(N ~ id, f1[f1$FG2=="Carnivore",], sum)
f1t$Ncarn <- f1carn$N[match(f1t$id, f1carn$id)]
f1t$Ncarn[is.na(f1t$Ncarn)] <- 0 

# D for feeding groups?
head(f1)
# Div index
f1tot <- aggregate(N~id, f1, sum)
f1f <- aggregate(N~FG+id, f1, sum)
f1f$tot <- f1tot$N[match(f1f$id, f1tot$id)]
f1f$rel <-  f1f$N / f1f$tot
f1f$simpD <- f1f$rel^2
f1f$shanH <- f1f$rel * log(f1f$rel)
fish_simp <- aggregate(simpD~id, f1f, sum)
fish_shan <- aggregate(shanH~id, f1f, sum)
f1t$FsimpD <- 1 - fish_simp$simpD[match(f1t$id, fish_simp$id)]
#f1t$simpD[f1t$coral==0] <- min(b1t$simpD, na.rm=T) #0
f1t$FshanH <- - fish_shan$shanH[match(f1t$id, fish_shan$id)]
#b1t$shanH[b1t$coral==0] <- min(b1t$shanH, na.rm=T)  #0
#ggplot(f1t, aes(FsimpD, FshanH))+geom_point()

#unique(f1[f1$FG=="Browser","Species"])

# simple D for feeding groups? collapse herbivores... 

head(f1)
unique(f1$FG)
f1$FG3 <- ifelse(f1$FG %in% c("Grazer", "Browser", "Scraper", "Excavator", "Farmer"), "Herbivore", ifelse(f1$FG %in% c("Benthic invertivore", "Piscivore", "Carnivore", "Corallivore"), "Carnivore", f1$FG))
unique(f1$FG3)

f1f2 <- aggregate(N~FG3+id, f1, sum)
f1f2$tot <- f1tot$N[match(f1f2$id, f1tot$id)]
f1f2$rel <-  f1f2$N / f1f2$tot
f1f2$simpD <- f1f2$rel^2
f1f2$shanH <- f1f2$rel * log(f1f2$rel)
fish_simp2 <- aggregate(simpD~id, f1f2, sum)
fish_shan2 <- aggregate(shanH~id, f1f2, sum)
f1t$FsimpD2 <- 1 - fish_simp2$simpD[match(f1t$id, fish_simp2$id)]
#f1t$simpD[f1t$coral==0] <- min(b1t$simpD, na.rm=T) #0
f1t$FshanH2 <- - fish_shan2$shanH[match(f1t$id, fish_shan2$id)]

ggplot(f1t, aes(FsimpD, FsimpD2))+geom_point(shape=21)


######################################################
#---------------------------------------------#  MULTIDIMENSIONAL FISH!?

head(f1f)

## multidimensional composition...
fmdsdat0 <- acast(f1f, id~FG, value.var="N") #n
head(fmdsdat0)
fmdsdat0[is.na(fmdsdat0)] <- 0
pca <- prcomp(sqrt(fmdsdat0), scale=T, center=T)
vecs <- data.frame(pca$rotation[,c("PC1", "PC2", "PC3", "PC4")], lab=rownames(pca$rotation))
vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100	
#biplot(pca)
f1t[,c("fcomp1", "fcomp2")] <- pca$x[match(f1t$id, rownames(pca$x)), c("PC1", "PC2")]


ggplot()+geom_point(data=f1t, aes(fcomp1, fcomp2, col=FsimpD))+
scale_colour_viridis()+
labs(x=paste("PCA 1 (", vars[1], "%)", sep=""),y=paste("PCA 2 (", vars[2], "%)", sep=""))+
geom_segment(data=vecs, aes(x=0, xend=PC1*3, y=0, yend=PC2*3))+
geom_text(data=vecs, aes(x=PC1*4, y=PC2*4, label=lab))


ggplot(f1t, aes(fcomp1, FsimpD))+geom_point()


######################################################
#---------------------------------------------#  Coral Sea (commbine)

unique(s1t$id)
unique(b1t$id)

nrow(f1t)
nrow(b1t)
nrow(s1t)
check <- data.frame(id = unique(c(f1t$id, b1t$id, s1t$id)))
nrow(check)
check$pit <- b1t$id[match(check$id, b1t$id)]
check$fish <- f1t$id[match(check$id, f1t$id)]
check$size <- s1t$id[match(check$id, s1t$id)]
nrow(na.omit(check))

check2 <- data.frame(Site = unique(c(f1t$Site, b1t$Site, s1t$Site)))
nrow(check2)
check2$pit <- b1t$Site[match(check2$Site, b1t$Site)]
check2$fish <- f1t$Site[match(check2$Site, f1t$Site)]
check2$size <- s1t$Site[match(check2$Site, s1t$Site)]
check2
nrow(na.omit(check2))

check[is.na(check$fish),] # no Fish
lengthcheck[is.na(check$fish),]
# check[is.na(check$pit),] # all PIT present in fish 

c1 <- b1t
c1[,c("f.biomass","f.biomassT", "f.N", "f.richness", "N.herb", "N.carn", "FshanH", "FsimpD","FsimpD2", "fcomp1", "fcomp2")] <- f1t[match(c1$id, f1t$id), c("biomass","biomassT", "N", "richness", "Nherb", "Ncarn", "FshanH", "FsimpD","FsimpD2", "fcomp1", "fcomp2")]
head(c1)

c1[,c("p10", "p20", "p60")] <- s1t[match(c1$id, s1t$id), c("p10", "p20", "p60")]
head(c1)

c1$shelt <- svol$shelt2[match(c1$id, svol$id)]
c1$planar <- svol$planar2[match(c1$id, svol$id)]

unique(c1$SectorRegion)
unique(c1$Reef)
unique(c1$Year)
length(unique(c1$id))
nrow(unique(na.omit(c1[,c("id","coral", "f.N", "p10", "Complexity")])))

head(c1)

ggplot(c1, aes(coral, shelt))+geom_point()+geom_smooth()+scale_y_sqrt()
summary(lm(sqrt(shelt)~coral, c1))

#ggplot(c1, aes(algae, N.herb/f.N))+geom_point()+geom_smooth(se=F, method="lm")+scale_x_sqrt()+scale_y_sqrt()
#ggplot(c1, aes(p60, coral))+geom_point()+geom_smooth(se=F, method="lm")+scale_x_sqrt()+scale_y_sqrt()

ggplot(c1, aes(f.N, f.biomass))+geom_point()+geom_smooth(se=F, method="lm")+scale_x_log10()+scale_y_log10()

######################################################
#---------------------------------------------#  Coral Sea LAT/LONGs - GET FROM MORGAN

ll1 <- read.csv("data/CoralSea/latlongs.csv")
c1$lat <- ll1$GPS.Smm[match(c1$Reef, ll1$ReefMM)]
c1$long <- ll1$GPS.Emm[match(c1$Reef, ll1$ReefMM)]
head(c1)

######################################################
#---------------------------------------------#  LTMP benthos

load(file = "data/LTMP/all_benthic_data.RData")
ls()
b2 <- data.frame(all_benthic_data)
head(b2)
unique(b2$cREPORT_YEAR)
nrow(b2)
unique(b2$COMP_2021)

# create site IDs
b2$id <- paste(b2$A_SECTOR, b2$SHELF, b2$REEF_NAME, b2$SITE_NO,  b2$TRANSECT_NO, b2$REPORT_YEAR,  sep="_")
head(b2)

# add region names
regs2 <- read.csv("data/LTMP/regions.csv")
b2$region <- regs2$name[match(b2$A_SECTOR, regs2$abb)]
head(b2)

# add lat longs
#ll <- read.csv("data/LTMP/latlongs.csv")
#lluse <- c("REEF_LAT", "REEF_LONG")
#b2[,lluse] <- ll[match(paste(b2$REEF_NAME, b2$SITE_NO), paste(ll$REEF_NAME, ll$SITE_NO)), lluse]

# missing lat/long
#b2$REEF_LAT[b2$REEF_NAME=="HELIX REEF"] <- -18.623812
#b2$REEF_LONG[b2$REEF_NAME=="HELIX REEF"] <- 147.295934

tax <- read.csv("data/LTMP/benthicTaxa.csv")
head(tax)

# taxon ID (lowest resolution)
unique(b2$COMP_2021) #68
tax$benthID <- paste(tax$GROUP_CODE, tax$FAMILY_2021, tax$COMP_2021) # 128
b2$benthID <- paste(b2$GROUP_CODE, b2$FAMILY_2021, b2$COMP_2021) # 128

# broad taxon groups
b2$morph7 <- tax$morph7[match(b2$benthID, tax$benthID)]
b2$group <- tax$group[match(b2$benthID, tax$benthID)]

unique(b2[b2$cREPORT_YEAR %in% c(2006, 2023), c("morph7", "cREPORT_YEAR")])
unique(b2[b2$cREPORT_YEAR %in% c(2006, 2023), c("morph7", "cREPORT_YEAR")])

# SUBSET TO ALIGN WITH CORAL SEA TIMINGS!?

b2 <- b2[b2$REPORT_YEAR %in% c("2018", "2019","2020", "2021", "2022", "2023", "2024"),]
head(b2)

######################################################
#---------------------------------------------#  LTMP complexity

# started in 2014

# complexity data
load(file = "data/LTMP/complexity.rdata") 
ls()
head(complexity)
head(b2)
unique(complexity$REPORT_YEAR)

# create site IDs
complexity$id <- paste(complexity$A_SECTOR,  complexity$SHELF, complexity$REEF_NAME, complexity$SITE_NO, complexity$TRANSECT_NO, complexity$REPORT_YEAR,sep="_")

allc <- data.frame(id=unique(c(b2$id[b2$REPORT_YEAR %in% c(2022, 2023, 2024)], complexity$id)))
allc$b2 <- unique(b2$id)[match(allc$id, unique(b2$id))]
allc$comp <- unique(complexity$id)[match(allc$id, unique(complexity$id))]
head(allc)
nrow(allc)

nomatch <- allc[is.na(allc$comp) | is.na(allc$b2),]  
nomatch  #just 1 missing

b2$Complexity <- complexity$HABITAT_COMPLEXITY[match(b2$id, complexity$id)]
unique(b2$Complexity)
head(b2)

head(b2[b2$REPORT_YEAR==2022,])

######################################################
#---------------------------------------------#  LTMP fish

# benthos/fish sent 6th Cept 2024.. 
# update on fish sent 26th Sept 2024?

#load(file = "data/LTMP/latest.fish.RData") # incorrect (shorter).. use the one sent later. 
#ls()
#f2x <- data.frame(latest.fish)
#head(f2x)
#hist(f2x$LENGTH)

load(file = "data/LTMP/260924_data_for_mike.RData")
ls()
f2 <- data.frame(data_for_Mike)
unique(f2$REPORT_YEAR)
nrow(f2)
head(f2)

f2 <- f2[f2$ABUNDANCE > 0, ]
nrow(f2)

f2$id <- paste(f2$A_SECTOR, f2$SHELF, f2$REEF_NAME, f2$SITE_NO, f2$TRANSECT_NO,f2$REPORT_YEAR, sep="_")
head(f2)

# species spelling (just for aligning FGs)
f2$spp <- paste(f2$GENUS, f2$SPECIES)
f2f1t <- read.csv("data/LTMP/LTMP_to_CS.csv")
f2$spp2 <- f2f1t$cs[match(f2$spp, f2f1t$ltmp)] # if spp = ltmp, replace spp with cs
f2$spp2 <- ifelse(is.na(f2$spp2), f2$spp, f2$spp2)
head(f2)

head(f1fg)
sp <- data.frame(all = unique(c(f2$spp2,f1fg$Species)))
sp$f2 <- unique(f2$spp2)[match(sp$all, unique(f2$spp2))]
sp$cs <- f1fg$Species[match(sp$all, f1fg$Species)]
head(sp)
sp[order(sp$all),]
#sp[is.na(sp$cs) | is.na(sp$f2),]
nrow(sp)
length(unique(f2$spp2))

f2fg <- read.csv("data/LTMP/ltmpFG.csv") # replace eventually
nrow(na.omit(f2fg))
nrow(f2fg)

# fill gaps
nodat <- f2fg$spp2[is.na(f2fg$FG)]
gaps <- f2[f2$spp2 %in% nodat,]
gapN <- aggregate(ABUNDANCE~spp2, gaps, sum)
gapN
sum(gapN$ABUNDANCE)/sum(f2$ABUNDANCE)*100 # 0.1-0.3% lacking
head(f2)

# add
f2$a <- f2fg$a[match(f2$spp2, f2fg$spp2)]
f2$b <- f2fg$b[match(f2$spp2, f2fg$spp2)]
f2$FG <- f2fg$FG[match(f2$spp2, f2fg$spp2)]


# SUBSET TO ALIGN WITH CORAL SEA TIMINGS!?

f2 <- f2[f2$REPORT_YEAR %in% c("2018", "2019","2020", "2021", "2022", "2023", "2024"),]
head(f2)



######################################################
#---------------------------------------------#  combine transects

all <- data.frame(id=unique(c(b2$id, f2$id)))
all$b2 <- unique(b2$id)[match(all$id, unique(b2$id))]
all$f2 <- unique(f2$id)[match(all$id, unique(f2$id))]
head(all)
nrow(all)

nomatch <- all[is.na(all$f2) | is.na(all$b2),]  
nomatch  # no name reef missing for benthos - only wreck island 2021 missing benthos

head(b2)
lcols <- c("A_SECTOR", "SHELF", "AIMS_REEF_NAME", "REEF_NAME", "FULLREEF_ID", "REPORT_YEAR", "SITE_NO" ,"TRANSECT_NO", "region", "id", "Complexity")
c2 <- unique(b2[,lcols])
head(c2)
length(unique(b2$id))
length(unique(b2$region))
length(unique(b2$REEF_NAME))



######################################################
#---------------------------------------------#  add metrics

# cover
unique(b2$group)
b2coral <- aggregate(cover~id, b2[b2$GROUP_CODE=="HC", ], sum)
c2$coral <- b2coral$cover[match(c2$id, b2coral$id)]
b2alg <- aggregate(cover ~ id, b2[b2$group=="Macroalgae",], sum)
c2$algae <- b2alg$cover[match(c2$id, b2alg$id)]
b2turf <- aggregate(cover ~ id, b2[b2$group=="Turf Algae",], sum)
c2$turf <- b2turf$cover[match(c2$id, b2turf$id)]
b2cca <- aggregate(cover ~ id, b2[b2$group=="CCA",], sum)
c2$cca <- b2cca$cover[match(c2$id, b2cca$id)]
b2soft <- aggregate(cover ~ id, b2[b2$group=="Soft Coral",], sum)
c2$soft <- b2soft$cover[match(c2$id, b2soft$id)]
# sand?
# other
head(c2)


# Acro 
b2$acro <- tax$acro[match(b2$COMP_2021, tax$COMP_2021)]
b2ac <- aggregate(cover ~ id, b2[b2$acro=="yes",], sum)
c2$acro <- b2ac$cover[match(c2$id, b2ac$id)] / c2$coral

covplot <- ggplot(c2, aes(x=coral))+geom_histogram(aes(fill=region))+
geom_vline(xintercept=mean(c2$coral ), linetype="dotted")+
scale_y_continuous(expand=c(0,0))+theme_minimal()
covplot

# fish 
f2$pres <- ifelse(f2$ABUNDANCE==0, 0, 1)
f2rich <- aggregate(pres~id, unique(f2[,c("spp2", "id", "pres")]), sum)  #needed due to multiple sizes
c2$f.richness <- f2rich$pres[match(c2$id, f2rich$id)]
f2abun <- aggregate(ABUNDANCE~id, f2, sum) 
c2$f.N <- f2abun$ABUNDANCE[match(c2$id, f2abun$id)]
f2$f.biomass <- f2$a  * (f2$LENGTH)^f2$b # are there sufficient length values???

f2$f.biomassT <- f2$f.biomass * f2$ABUNDANCE

head(f2) 
f2B <- aggregate(f.biomass ~ id, f2, sum)
f2BT <- aggregate(f.biomassT ~ id, f2, sum)
c2$f.biomass <- f2B$f.biomass[match(c2$id, f2B$id)]
c2$f.biomassT <- f2BT$f.biomassT[match(c2$id, f2BT$id)]


head(c2)

# max fish biomass? 
fmax <- c2[which.max(c2$f.biomass),]
fmax
f2[f2$id %in% fmax$id,c("id", "spp","LENGTH", "ABUNDANCE", "f.biomass","f.biomassT")]

ggplot(c2, aes(f.N, f.biomass))+geom_point()+scale_y_log10()+scale_x_log10()


# coral DIV indices
head(b2)
b2m <- aggregate(cover~morph7+id, b2[b2$group=="Hard Coral",], sum) # check it's getting all groups in each year
b2m$coral <- c2$coral[match(b2m$id, c2$id)]
b2m$rel_cover <-  b2m$cover / b2m$coral
hist(b2m$rel_cover)
b2m$simpD <- b2m$rel_cover^2
b2m$shanH <- b2m$rel_cover * log(b2m$rel_cover)
coral_simp2 <- aggregate(simpD~id, b2m, sum)
coral_shan2 <- aggregate(shanH~id, b2m, sum)
c2$simpD <- 1 - coral_simp2$simpD[match(c2$id, coral_simp2$id)]
c2$shanH <- - coral_shan2$shanH[match(c2$id, coral_shan2$id)]

c2$simpD[c2$coral==0] <- min(c2$simpD, na.rm=T) #0 ????
c2$shanH[c2$coral==0] <- min(c2$shanH, na.rm=T)  #0 
ggplot(c2, aes(simpD, shanH))+geom_point()

# fish comp
unique(f2$FG)
f2$FG2 <- ifelse(f2$FG %in% c("Carnivore", "Piscivore", "Benthic invertivore", "Corallivore"), "Carnivore", ifelse(f2$FG %in% c("Grazer", "Scraper", "Excavator", "Browser"), "Herbivore", "Other")) # excavators = parrotfish
f2herb <- aggregate(ABUNDANCE ~ id, f2[f2$FG2=="Herbivore",], sum)
c2$N.herb <- f2herb$ABUNDANCE[match(c2$id, f2herb$id)]
c2$N.herb[is.na(c2$Nherb)] <- 0 
f2carn <- aggregate(ABUNDANCE ~ id, f2[f2$FG2=="Carnivore",], sum)
c2$N.carn <- f2carn$ABUNDANCE[match(c2$id, f2carn$id)]
c2$N.carn[is.na(c2$Ncarn)] <- 0

f2$Species <- f2$spp
FGcheck <- unique(rbind(f2[,c("Species", "FG", "FG2")], f1[,c("Species", "FG", "FG2")]))
nrow(FGcheck[is.na(FGcheck$FG),])
FGcheck$FG[is.na(FGcheck$FG)]<-""
FGcheck$FG2[is.na(FGcheck$FG)]<-""
nrow(FGcheck)
head(FGcheck)
unique(FGcheck$FG)

#write.csv(FGcheck, "fish_groups.csv")

# fish DIV indices
head(f2)
# Div index
f2f <- aggregate(ABUNDANCE~FG+id, f2, sum)
f2f$tot <- c2$f.N[match(f2f$id, c2$id)]
f2f$rel <-  f2f$ABUNDANCE / f2f$tot
#hist(f2f$rel)
f2f$simpD <- f2f$rel^2
f2f$shanH <- f2f$rel * log(f2f$rel)
fish_simp2 <- aggregate(simpD~id, f2f, sum)
fish_shan2 <- aggregate(shanH~id, f2f, sum)
c2$FsimpD <- 1 - fish_simp2$simpD[match(c2$id, fish_simp2$id)]
c2$FshanH <- - fish_shan2$shanH[match(c2$id, fish_shan2$id)]
#b1t$shanH[b1t$coral==0] <- min(b1t$shanH, na.rm=T)  #0
ggplot(c2, aes(FsimpD, FshanH))+geom_point()


# simple D for feeding groups?
head(f2)
unique(f2$FG)
f2$FG3 <- ifelse(f2$FG %in% c("Carnivore", "Piscivore", "Benthic invertivore", "Corallivore"), "Carnivore", ifelse(f2$FG %in% c("Grazer", "Scraper", "Excavator", "Browser", "Farmer"), "Herbivore", f2$FG)) # excavators = parrotfish
unique(f2$FG3)

f2f2 <- aggregate(ABUNDANCE~FG3+id, f2, sum)
f2f2$tot <- c2$f.N[match(f2f2$id, c2$id)]
f2f2$rel <-  f2f2$ABUNDANCE / f2f2$tot
#hist(f2f$rel)
f2f2$simpD <- f2f2$rel^2
f2f2$shanH <- f2f2$rel * log(f2f2$rel)
fish_simp2b <- aggregate(simpD~id, f2f2, sum)
fish_shan2b <- aggregate(shanH~id, f2f2, sum)
c2$FsimpD2 <- 1 - fish_simp2b$simpD[match(c2$id, fish_simp2b$id)]
c2$FshanH2 <- - fish_shan2b$shanH[match(c2$id, fish_shan2b$id)]
#b1t$shanH[b1t$coral==0] <- min(b1t$shanH, na.rm=T)  #0

ggplot(c2, aes(FsimpD, FsimpD2))+geom_point()



ggplot(c2, aes(f.biomass, f.richness))+geom_point()+scale_x_log10()+scale_y_log10()


nrow(unique(na.omit(c2[,c("id","coral", "f.N", "Complexity")])))


######################################################
#---------------------------------------------#  MULTIDIMENSIONAL?

## multidimensional composition...

head(b2m)

mdsdat2 <- acast(b2m[b2m$coral>0,], id~morph7, value.var="rel_cover") #n
pca2 <- prcomp(sqrt(mdsdat2), scale=T, center=T)
vecs2 <- data.frame(pca2$rotation[,c("PC1", "PC2", "PC3", "PC4")], lab=rownames(pca2$rotation))
vars2 <- round((pca2$sdev^2 / sum(pca2$sdev^2)), 3)*100	
#biplot(pca)
c2[,c("comp1", "comp2")] <- pca2$x[match(c2$id, rownames(pca2$x)), c("PC1", "PC2")]

ggplot()+geom_point(data=c2, aes(comp1, comp2, col=simpD))+
scale_colour_viridis()+
labs(x=paste("PCA 1 (", vars[1], "%)", sep=""),y=paste("PCA 2 (", vars[2], "%)", sep=""))+
geom_segment(data=vecs2, aes(x=0, xend=PC1*3, y=0, yend=PC2*3))+
geom_text(data=vecs2, aes(x=PC1*4, y=PC2*4, label=lab))


######################################################
#---------------------------------------------#  MULTIDIMENSIONAL fish?

## multidimensional composition...

head(f2f)

fmdsdat2 <- acast(f2f, id~FG, value.var="ABUNDANCE") #n
fmdsdat2[is.na(fmdsdat2)]<-0
pca2 <- prcomp(sqrt(fmdsdat2), scale=T, center=T)
vecs2 <- data.frame(pca2$rotation[,c("PC1", "PC2", "PC3", "PC4")], lab=rownames(pca2$rotation))
vars2 <- round((pca2$sdev^2 / sum(pca2$sdev^2)), 3)*100	
#biplot(pca2)
c2[,c("fcomp1", "fcomp2")] <- pca2$x[match(c2$id, rownames(pca2$x)), c("PC1", "PC2")]

fcomp1 <- ggplot()+geom_point(data=c2, aes(fcomp1, fcomp2, col=FsimpD))+
scale_colour_viridis()+
labs(x=paste("PCA 1 (", vars[1], "%)", sep=""),y=paste("PCA 2 (", vars[2], "%)", sep=""))+
geom_segment(data=vecs2, aes(x=0, xend=PC1*3, y=0, yend=PC2*3))+
geom_text(data=vecs2, aes(x=PC1*4, y=PC2*4, label=lab))

fcomp2 <- ggplot()+geom_point(data=c1, aes(fcomp1, fcomp2, col=FsimpD))+
scale_colour_viridis()+
labs(x=paste("PCA 1 (", vars[1], "%)", sep=""),y=paste("PCA 2 (", vars[2], "%)", sep=""))+
geom_segment(data=vecs2, aes(x=0, xend=PC1*3, y=0, yend=PC2*3))+
geom_text(data=vecs2, aes(x=PC1*4, y=PC2*4, label=lab))

plot_grid(fcomp1, fcomp2)


######################################################
#---------------------------------------------# latlong - need more

head(c2)
ll2 <- read.csv("data/LTMP/latlongs.csv")
head(ll2)
c2$lat <- ll2$REEF_LAT[match(c2$REEF_NAME, ll2$REEF_NAME)]
c2$long <- ll2$REEF_LONG[match(c2$REEF_NAME, ll2$REEF_NAME)]

######################################################
#---------------------------------------------#  combine/save

metrics <- c("id", "coral", "algae", "turf", "soft", "cca", "simpD", "shanH", "acro", "comp1", "comp2", "f.biomass","f.biomassT", "f.N", "f.richness", "N.herb", "N.carn", "FshanH", "FsimpD", "FsimpD2","fcomp1", "fcomp2","Complexity", "lat", "long")

head(c1)
head(c2)

comb2 <- cbind(c2[,metrics,], data.frame(Data="LTMP", Marine.Park="GBR", SectorRegion=c2$region, Reef=c2$REEF_NAME, Site=paste(c2$REEF_NAME, c2$SITE_NO), Transect=c2$TRANSECT_NO, Zone="Slope", Year=c2$REPORT_YEAR))
head(comb2)

comb1 <- cbind(c1[,metrics], Data="Coral Sea", c1[,c("Marine.Park","SectorRegion", "Reef", "Site", "Transect","Zone", "Year")])

colnames(comb2)==colnames(comb1)

comb <- rbind(comb2, comb1)

write.csv(c2, "data/output/metricsLTMP.csv")
write.csv(c1, "data/output/metricsCS.csv")
write.csv(comb, "data/output/metrics.csv")




