


p10dat <- df1[,c("p10", "coral", "Reef", "Zone")]
slpdat <- reefs[reefs$y=="p10",]
p10dat[,c("slp", "dir", "sig", "sigdir")] <- slpdat[match(p10dat$Reef, slpdat$Reef), c("slp", "dir", "sig", "sigdir")]
p10fit <- bifits[bifits$resp=="p10",]

p10plot <- ggplot(p10dat[p10dat$coral > 0 & p10dat$Zone %in% c("Crest", "Slope"),], aes( coral, p10))+
geom_point(shape=21, aes(col=sigdir), size=0.1)+
geom_line(data=p10fit, aes(x, fit, group=Reef, col=sigdir))+
labs(x="% coral cover", y="Proportion corals < 10cm")+
scale_colour_manual(values=colz)+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.2), axis.text=element_text(size=8), axis.title=element_text(size=8))
p10plot

p60dat <- df1[,c("p60", "coral", "Reef", "Zone")]
slpdat <- reefs[reefs$y=="p60",]
p60dat[,c("slp", "dir", "sig", "sigdir")] <- slpdat[match(p60dat$Reef, slpdat$Reef), c("slp", "dir", "sig", "sigdir")]
p60fit <- bifits[bifits$resp=="p60",]
p60plot <- ggplot(p60dat[p60dat$coral > 0 & p10dat$Zone %in% c("Crest", "Slope"),], aes(coral, p60))+
geom_point(shape=21, aes(col=sigdir), size=0.1)+
geom_line(data=p60fit, aes(x, fit, group=Reef, col=sigdir))+
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


refplot <- ggplot(c1, aes(x=coral, y=shelt))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
scale_y_sqrt()+
annotate(geom = 'text', x = 75, y = 0,    label = paste("R^2 == ", rsq2), parse = TRUE, size=2.5, col="blue")+
labs(x="% coral cover", y=expression(Shelter~Volume~(dm^3)))+
scale_colour_manual(values=colz)+
theme_classic()+theme(strip.background=element_blank(), axis.line=element_line(size=0.2), axis.text=element_text(size=8), axis.title=element_text(size=8))
refplot


FIG4 <- plot_grid(sizecov, refplot, rel_widths=c(2.23, 1), labels=c("", "c"), label_size=10)
FIG4
