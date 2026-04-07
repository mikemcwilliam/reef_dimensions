

covD <- ggplot()+
geom_point(data=c1, aes(coral, simpD), shape=21, size=0.5)+
geom_ribbon(data=fit, aes(x=coral, ymin=fit-se, ymax=fit+se), alpha=0.5)+
geom_line(data=fit, aes(coral, fit), col="red")+
geom_point(data=maxminC,aes(n, simpDx, fill=id), shape=21, size=2.5)+
scale_fill_manual(values=c("#e0f3f8","#4575b4","#91bfdb", "#fee0b6", "#b35806"))+
lims(y=c(-0.15,1.05))+
guides(fill="none")+
#geom_smooth(data=c1, aes(coral, simpD))+
labs(x="% coral cover", y="Morphological\ndiversity (D)")+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8)) 
covD



######################################################
#---------------------------------------------#  DEMO FIG FOR SIMPS-D

highD <- ggplot(maxDs[!maxDs$morph8 %in% c("Acr. Bushy"),] , aes(y=morph8, x=n, fill=t))+
geom_bar(stat="identity", position="dodge", col="black", size=0.05)+
ggtitle("high D")+
scale_x_continuous(expand=c(0,0), breaks=c(0, 5, 10, 15), limits=c(0, 20))+
facet_wrap(~t)+
#scale_fill_manual(values=c("#d8daeb", "#998ec3", "#542788"))+
scale_fill_manual(values=c("#e0f3f8", "#91bfdb", "#4575b4"))+
xlab("% cover")+theme_classic()+theme(axis.line=element_line(size=0.1),axis.title.y=element_blank(), axis.title.x=element_text(size=8), axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), 
plot.title=element_text(hjust=0.5, size=8), strip.text=element_blank(), strip.background=element_blank())
highD


lowD <- ggplot(minDs[!minDs$morph8 %in% c("Acr. Bushy"),] , aes(y=morph8, x=n, fill=t))+
geom_bar(stat="identity", position="dodge", col="black", size=0.05)+
ggtitle("low D")+
facet_wrap(~t)+
scale_x_continuous(expand=c(0,0))+
scale_fill_manual(values=c("#fee0b6", "#b35806"))+
xlab("% cover")+theme_classic()+theme(axis.line=element_line(size=0.1),axis.title.y=element_blank(),axis.title.x=element_text(size=8), axis.text.x=element_text(size=8), axis.text.y=element_text(size=7), 
plot.title=element_text(hjust=0.5, size=8), strip.text=element_blank(), strip.background=element_blank())
lowD 


demoplot <- plot_grid(
highD+guides(fill="none"),
covD, 
plot_grid(lowD+guides(fill="none"), NULL, rel_widths=c(1, 0.25)), 
ncol=1, rel_heights=c(1, 2, 1), labels=c("a", "b", "c"), label_size=10)
demoplot


######################################################
#---------------------------------------------#  hard coral FD

df2 <- read.csv("data/output/coraltraits.csv")
exampleFD <- read.csv("data/output/exampleFD.csv")
minplot <- exampleFD[exampleFD$type=="LowFD",]
maxplot <- exampleFD[exampleFD$type=="HighFD",]


coord <- aggregate(dim1~link, df2, mean)
coord$dim2 <- aggregate(dim2~link, df2, mean)$dim2
coord$dim3 <- aggregate(dim3~link, df2, mean)$dim3
coord$dim4 <- aggregate(dim4~link, df2, mean)$dim4
coord

ggplot() +
  geom_point(data=df2, aes(dim1, dim2), shape=21, size=0.5, stroke=0.25, col="black")+
geom_point(data=coord, aes(dim1, dim2), col="red")

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


FDtoD <- ggplot(df1, aes(simpD, FDis))+
geom_point(shape=21, size=0.5)+
geom_smooth(method="lm", col="blue", size=0.3)+
labs(x="Morph. diversity (D)", y="Functional dispersion")+
annotate(geom = 'text', x = 0.8, y = 0,    label = paste("R^2 == ", rsq), parse = TRUE, size=2.5, col="blue")+
scale_y_sqrt()+
theme_classic()+theme(axis.line=element_line(size=0.1), axis.text=element_text(size=8), axis.title=element_text(size=8)) 
FDtoD 

FDdemoplot2 <- plot_grid(
plot_grid(NULL, FDdemo, NULL, nrow=1, rel_widths=c(0.5, 1,0.5)),
FDtoD, ncol=1, rel_heights=c(1, 2), labels=c("d", "e"), label_size=10)
FDdemoplot2

FIG3 <- plot_grid(demoplot, NULL, FDdemoplot2, nrow=1, rel_widths=c(1, 0.1, 1))+
draw_line(x=c(0.2, 0.17), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.27, 0.255), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.37, 0.365), y=c(0.78, 0.68), size=0.1)+
draw_line(x=c(0.13, 0.17), y=c(0.36, 0.2), size=0.1)+
draw_line(x=c(0.4, 0.33), y=c(0.36, 0.2), size=0.1)+
draw_line(x=c(0.65, 0.58), y=c(0.8, 0.8), size=0.1)+
draw_line(x=c(0.58, 0.58), y=c(0.8, 0.7), size=0.1, arrow=arrow(length=unit(1, "mm")))
FIG3


