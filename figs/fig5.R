



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


FIG5 <- plot_grid(
plot_grid(ass1, bars1, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass2, bars2, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass3, bars3, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass4, bars4, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass5, bars5, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
plot_grid(ass6, bars6, NULL, rel_widths=c(1,0.1, 0.05), align="hv", axis="tb", nrow=1),
labels=c("a", "b", "c", "d", "e", "f"), label_size=9, nrow=3)
FIG5

