

# 

vectors$group <- labs$group[match(vectors$lab, labs$lab)]
vectors$lab2 <- labs$lab2[match(vectors$lab, labs$lab)]

coords$PC2 <- ifelse(coords$data=="LTMP", -coords$PC2, coords$PC2)
vectors$PC2 <- ifelse(vectors$data=="LTMP", -vectors$PC2, vectors$PC2)

coords$axes2 <- ifelse(coords$axes=="short", "6 metrics", "14 metrics")
vectors$axes2 <- ifelse(vectors$axes=="short", "6 metrics", "14 metrics")

head(vectors)
vectors$PC1b <- vectors$PC1
vectors$PC2b <- vectors$PC2
vectors$PC2b[vectors$lab2=="CCA ratio" & vectors$data=="LTMP" ] <- vectors$PC2b[vectors$lab2=="CCA ratio" & vectors$data=="LTMP" ] +0.08
vectors$PC2b[vectors$lab2=="acropora ratio" & vectors$axes=="long" ] <- vectors$PC2b[vectors$lab2=="acropora ratio" & vectors$axes=="long" ] -0.08
vectors$PC2b[vectors$lab2=="complexity (1-5)" & vectors$data=="LTMP" & vectors$axes=="long"] <- vectors$PC2b[vectors$lab2=="complexity (1-5)" & vectors$data=="LTMP" & vectors$axes=="long"] +0.08
vectors$PC2b[vectors$lab2=="complexity (1-5)" & vectors$data=="Coral Sea" & vectors$axes=="long"] <- vectors$PC2b[vectors$lab2=="complexity (1-5)" & vectors$data=="Coral Sea" & vectors$axes=="long"] -0.05
vectors$PC2b[vectors$lab2=="Morph. div (D)" & vectors$data=="Coral Sea" & vectors$axes=="long"] <- vectors$PC2b[vectors$lab2=="Morph. div (D)" & vectors$data=="Coral Sea" & vectors$axes=="long"] +0.05
vectors$PC2b[vectors$lab2=="Morph. div (D)" & vectors$data=="LTMP" & vectors$axes=="long"] <- vectors$PC2b[vectors$lab2=="Morph. div (D)" & vectors$data=="LTMP" & vectors$axes=="long"] -0.015
vectors$PC1b[vectors$lab2=="complexity (1-5)"  & vectors$axes=="short"] <- vectors$PC1b[vectors$lab2=="complexity (1-5)" & vectors$axes=="short"] -0.1

mag1<- 5
mag2<- 6
pcas <- ggplot()+
geom_point(data=coords, aes(PC1, PC2),col="grey", size=0.1, shape=21, stroke=0.25)+
geom_segment(data=vectors, aes(x=0, xend=PC1*mag1, y=0, yend=PC2*mag1, col=group))+
geom_text(data=vectors, aes(PC1b*mag2, PC2b*mag2, label=lab2, col=group), size=1.9, hjust=ifelse(vectors$PC1b>0, 0, 1),fontface="bold", show_guide = FALSE)+ #,, direction="y", force=0.25, segment.size=0.2
scale_colour_manual(values=c("darkgreen", "darkorchid", "darkred","black"))+ ##998ec3
facet_grid(axes2~data)+
theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), strip.background=element_blank(), legend.title=element_blank(), legend.text=element_text(size=7), legend.key.height=unit(1, "mm"), panel.margin=unit(1, "mm"), axis.text=element_text(size=6), axis.title=element_text(size=6), strip.text=element_text(size=8))
pcas
