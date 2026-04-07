

################################################

#  MANTEL PLOT 


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

################################################

# PCA reconstruction

df_pca <- data.frame(
  PC = factor(paste0("PC", 1:6)),
  VarianceExplained = pca1a$sdev[1:6]^2 / sum(pca1a$sdev^2),  # proportion of variance per PC
  R2 = r2_values  # your PCA reconstruction R2 vector
)


varplot <- ggplot(df_pca, aes(x = PC)) +
  geom_col(aes(y = R2), fill = "lightblue", col="black", size=0.1, width=0.7) +
  #geom_line(aes(y = R2), group = 1, color = "red", size = 1) +
  #geom_point(aes(y = R2), color = "red", size = 2) +
  scale_y_continuous( name = "R² of Reconstruction", expand=c(0,0))  +
  labs(title = "Original PCA ~ 6 metrics", x="Original PCA axis") +
  theme_classic()+theme(plot.title=element_text(size=8, hjust=0.5, face="bold"), axis.text.y=element_text(size=8),axis.text.x=element_text(size=8, angle=45, hjust=1),axis.title=element_text(size=8), axis.line=element_line(size=0.1))
varplot


################################################ COMPARE PCAS

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

################################################ COMPARE PCAS2

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





################################################

# FIGURE 2

Fig2 <- plot_grid(
plot_grid(
plot_grid( distGrob, varplot, nrow=1, rel_widths=c(1, 0.7), labels=c("a", "b"), label_size=9),
NULL,
plot_grid(predplot, predplot3+guides(fill='none')+coord_flip(), rel_heights=c(2, 1), ncol=1, labels=c("d", ""), label_size=9),
ncol=1, rel_heights=c(1, 0.1, 1.4)),
plot_grid(NULL, twovecs+facet_wrap(~data2, ncol=1), ncol=1, rel_heights=c(0.1, 1),labels=c("c", ""), label_size=9, hjust=-10), rel_widths=c(1.5, 1))+
draw_text("PCA vectors\n(6 vs 14 metrics)", 0.8, 0.97, size=8, fontface="bold")+
draw_text("PCA coordinates (6 vs 14 metrics)", 0.33, 0.57, size=8, fontface="bold")
Fig2
