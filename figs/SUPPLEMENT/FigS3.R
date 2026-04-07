

# Supplemental  figures

head(vars)
vars$lab <- labs$lab2[match(vars$y, labs$lab)]
vars$group <- labs$group[match(vars$y, labs$lab)]

vars$group <- factor(vars$group, levels=c("coral", "fish", "algae","complexity"))
vars$scale <- ifelse(vars$grp=="SectorRegion", "Region", ifelse(vars$grp=="Reef:SectorRegion", "Reef", ifelse(vars$grp=="Site:(Reef:SectorRegion)", "Site", ifelse(vars$grp=="Residual", "Transect",NA))))

vars$scale  <- factor(vars$scale , levels=c("Region", "Reef", "Site","Transect"))
vars$line <- paste( vars$Year) #paste(vars$x, vars$Year)

SUP3 <- ggplot(vars[!vars$y %in% c("PC1_1a", "PC2_1a", "PC3_1a"),], aes(x=scale, y=p*100,col=group, group=as.factor(Year)))+
geom_point(size=0.5)+
geom_line(aes(group=line))+
theme_classic()+guides(col="none")+
scale_colour_manual(values=c("black", "slategrey", "darkgreen", "darkred"))+
labs(y="% variance explained")+
ylim(c(0, 70))+
theme(axis.text.x=element_text(angle=45, hjust=1), strip.background=element_blank(), axis.line=element_line(size=0.1), axis.title.x=element_blank(), legend.title=element_blank())+
facet_wrap(group~lab, scales="free", nrow=3)
SUP3


# ggsave("figs/SUPPLEMENT/SUP3.jpg", SUP3, height=7, width=7)
