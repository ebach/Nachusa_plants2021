#Elizabeth Bach
#Nachusa Grassland long-term vegetation monitoring
#Multivariate analysis
#January 2021

rm(list=ls())
library(reshape2)
library(plyr)
library(ggplot2)
library(vegan)
library(gridExtra)
library(indicspecies)
library(RColorBrewer)

#read in data, "Nachusa_veg_transect_data1994_2020vFinal.csv"
Nachusa.data<-read.csv(file.choose(),header=T, na.strings="NA")
head(Nachusa.data[,1:11])

#cast wide
Nachusa.wide<-dcast(Nachusa.data, Transect + Sample_Year + Quadrat ~ newcol, sum, value.var = "cover_value")
#484 rows, 429 cols
head(Nachusa.wide)

#reduce data to transects >3 sampling points, P24 now included (data error corrected)
Nachusa.wide2<-droplevels(subset(Nachusa.wide, Nachusa.wide$Transect=="P08"|Nachusa.wide$Transect=="P13"|Nachusa.wide$Transect=="P14"|Nachusa.wide$Transect=="P15"|Nachusa.wide$Transect=="P17"|Nachusa.wide$Transect=="P19"|Nachusa.wide$Transect=="P23"|Nachusa.wide$Transect=="P24"|
                                   Nachusa.wide$Transect=="P26"|Nachusa.wide$Transect=="P31"|Nachusa.wide$Transect=="P31"|Nachusa.wide$Transect=="P32"|Nachusa.wide$Transect=="P39"|Nachusa.wide$Transect=="P41"))
#585 observations

#add habitat variable to wide data
habitat.class<-"initialize"
for(i in 1:nrow(Nachusa.wide2)){
  if(Nachusa.wide2[i,1]=="P08"){
    habitat.class[i]<-paste0("savanna")
  }else if(Nachusa.wide2[i,1]=="P13"){
    habitat.class[i]<-paste0("wetland_restored")
  }else if(Nachusa.wide2[i,1]=="P14"){
    habitat.class[i]<-paste0("native_prairie")
  }else if(Nachusa.wide2[i,1]=="P15"){
    habitat.class[i]<-paste0("native_prairie")
  }else if(Nachusa.wide2[i,1]=="P17"){
    habitat.class[i]<-paste0("native_prairie")
  }else if(Nachusa.wide2[i,1]=="P19"){
    habitat.class[i]<-paste0("native_prairie")
  }else if(Nachusa.wide2[i,1]=="P23"){
    habitat.class[i]<-paste0("prairie_restored")
  }else if(Nachusa.wide2[i,1]=="P26"){
    habitat.class[i]<-paste0("prairie_restored")
  }else if(Nachusa.wide2[i,1]=="P31"){
    habitat.class[i]<-paste0("prairie_restored")
  }else if(Nachusa.wide2[i,1]=="P32"){
    habitat.class[i]<-paste0("prairie_restored")
  }else if(Nachusa.wide2[i,1]=="P39"){
    habitat.class[i]<-paste0("savanna")
  }else if(Nachusa.wide2[i,1]=="P41"){
    habitat.class[i]<-paste0("savanna")
  }else if(Nachusa.wide2[i,1]=="P24"){
    habitat.class[i]<-paste0("savanna")
  }
}
Nachusa.wide2$habitat<-as.factor(habitat.class)

#add site names for Nachusa graphs
site.names<-"initialize"
for(i in 1:nrow(Nachusa.wide2)){
  if(Nachusa.wide2[i,1]=="P08"){
    site.names[i]<-paste0("W.Heinkel")
  }else if(Nachusa.wide2[i,1]=="P13"){
    site.names[i]<-paste0("Prairie Potholes")
  }else if(Nachusa.wide2[i,1]=="P14"){
    site.names[i]<-paste0("Schaffer Knob")
  }else if(Nachusa.wide2[i,1]=="P15"){
    site.names[i]<-paste0("Main Unit Knobs")
  }else if(Nachusa.wide2[i,1]=="P17"){
    site.names[i]<-paste0("E.Heinkel Clearing")
  }else if(Nachusa.wide2[i,1]=="P19"){
    site.names[i]<-paste0("Doug's Knob")
  }else if(Nachusa.wide2[i,1]=="P23"){
    site.names[i]<-paste0("1991 Planting")
  }else if(Nachusa.wide2[i,1]=="P26"){
    site.names[i]<-paste0("Doug/Dot Planting")
  }else if(Nachusa.wide2[i,1]=="P31"){
    site.names[i]<-paste0("1986 Oak Island")
  }else if(Nachusa.wide2[i,1]=="P32"){
    site.names[i]<-paste0("1987 Oak Island")
  }else if(Nachusa.wide2[i,1]=="P39"){
    site.names[i]<-paste0("W.Heinkel-Jay")
  }else if(Nachusa.wide2[i,1]=="P41"){
    site.names[i]<-paste0("E.Heinkel")
  }else if(Nachusa.wide2[i,1]=="P24"){
    site.names[i]<-paste0("Shabbona")
  }
}
Nachusa.wide2$site<-as.factor(site.names)

#convert to presence/absence
data.pa<-cbind(Nachusa.wide2[,c(1:3,430:431)],decostand(Nachusa.wide2[,-c(1:3,430:431)],"pa"))
data.pa$Sample_Year<-as.factor(data.pa$Sample_Year)

#drop wetland site, doesn't make sense to include here, even for anaylsis within the transect
data.pa2<-droplevels(subset(data.pa, data.pa$Transect!="P13"))

#Use the adonis() function to perform PERMANOVA (non-parametric, multivariate analysis of variance) on the communities of each sample
#Again, exclude any columns of metadata (here, columns 1-8)
#First run the full model, look for any interactions with factors
adonis(data.pa2[,-c(1:5)]~data.pa2$habitat*data.pa2$Transect*data.pa2$Sample_Year, permutations=999)
#interactions are strongly significant, we want to look at changes across time in individual transcets
#look at each habitat individually

#modify NMDS.ggplot code to group by transect rather than year
ggplot.NMDS<-function(XX,ZZ,COLORS){
  library(ggplot2)
  MDS1<-data.frame(scores(XX))$NMDS1
  MDS2<-data.frame(scores(XX))$NMDS2
  Site<-ZZ
  
  NMDS<-data.frame(MDS1,MDS2,Site)
  
  NMDS.mean=aggregate(NMDS[,1:2],list(group=Site),mean)
  
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  for(g in levels(NMDS$Site)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Site==g,],
                                                     veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                  ,group=g))
  }
  
  X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Site),size=5) +
    theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))+theme(legend.title=element_text(size=11),legend.text=element_text(size=11), axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X1    
}


#native prairie (I would expect these to differ the most from each other)
native.prairie<-droplevels(subset(data.pa2, data.pa2$habitat=="native_prairie"))
native.prairie$Sample_Year<-as.factor(native.prairie$Sample_Year)
native.prairie$TransYear<-as.factor(paste(native.prairie$Transect, native.prairie$Sample_Year, sep=":"))
native.prairie$SiteYear<-as.factor(paste(native.prairie$site, native.prairie$Sample_Year, sep=":"))
#All native prairie transects have 15 quadrats, so keep as is

adonis(native.prairie[,-c(1:5,432:433)]~native.prairie$Transect*native.prairie$Sample_Year, permutations=999)
#strong interaction between transect and year
#TukeyHSD to untangle
mds.dist<-metaMDSdist(decostand(native.prairie[,-c(1:5,432:433)], "pa"),k=2, index="jaccard", autotransform=FALSE,na.rm=TRUE)
np.stat<-betadisper(mds.dist, native.prairie$TransYear, type="median")
np.stat

#Now use TukeysHSD to contrast the median dispersion of the groups and find differences
TukeyHSD(np.stat)

mds.pa2<-metaMDS(decostand(native.prairie[,-c(1:5,432:433)],"pa" ),distance="jaccard", k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#stress=0.185
mds.pa3<-metaMDS(decostand(native.prairie[,-c(1:5,432:433)],"pa" ),distance="jaccard", k=3,autotransform=FALSE, na.rm=TRUE)
mds.pa3
#stress=0.134

#will use 2D, not that much advantage of 3D

#graph

#colors to add in sampling year dimension
colors.14<-c(rev(brewer.pal(n=3, name = "YlOrRd")),rev(brewer.pal(n=3, name="Greens")),rev(brewer.pal(n=5, name="PuRd")),rev(brewer.pal(n=3,name="Blues")))
 
#interaction graph
np.graph2<-ggplot.NMDS(mds.pa2, native.prairie$TransYear, colors.14)+annotate("text", x=1.1, y=-1.0, label="stress=0.185", size=6)+labs(title = "Native Prairie")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())
np.graph2

#Nachusa labels
#rearrange colors to reflect factors levels (which are re-ordered alphabetically)
colors.14N<-c(rev(brewer.pal(n=3, name = "YlOrRd")),rev(brewer.pal(n=5, name="Greens")),rev(brewer.pal(n=3, name="PuRd")),rev(brewer.pal(n=3,name="Blues")))

np.graph2N<-ggplot.NMDS(mds.pa2, native.prairie$SiteYear, colors.14N)+annotate("text", x=1, y=-1.0, label="stress=0.185", size=6)+labs(title = "Native Prairie")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())
np.graph2N

#Look statistically at year changes within each transect, as vector analysis
#This tests for differences in centroids between transects
envectors<-envfit(mds.pa2 ~ native.prairie$Transect, na.rm=TRUE)
head(envectors)
centroids.trans<-data.frame(envectors$factors[1:4])
centroids.trans
#highly significant, as expected from the ADONIS

#Now layer this to run Sample_Yea vectors from these centroids, sample year needs to be a continuous variable, not a factor
native.prairie$Sample_Year<-as.numeric(paste(native.prairie$Sample_Year))

envectors2<-envfit(mds.pa2 ~ native.prairie$Transect + native.prairie$Sample_Year, na.rm=TRUE)
envectors2
#overall Sample_Year vector is not significant, transects likely shifting in diff directions

#Run the adonis within each transect, but use scores from full ordination
#subset native prairie dataset and run adonis to see if Sample_Year significant
P14<-droplevels(subset(native.prairie, native.prairie$Transect=="P14"))
adonis(P14[,-c(1:5,432:433)]~P14$Sample_Year, permutations=999)
#sample year strongly significant

P15<-droplevels(subset(native.prairie, native.prairie$Transect=="P15"))
adonis(P15[,-c(1:5,432:433)]~P15$Sample_Year, permutations=999)
#Sample year strongly significant

P17<-droplevels(subset(native.prairie, native.prairie$Transect=="P17"))
adonis(P17[,-c(1:5,432:433)]~P17$Sample_Year, permutations=999)
#Sample year strongly significant

P19<-droplevels(subset(native.prairie, native.prairie$Transect=="P19"))
adonis(P19[,-c(1:5,432:433)]~P19$Sample_Year, permutations=999)
#Sample year strongly significant

#run envfit regressions, keep original scores from full NMDS (e.g. not a new dispersal for each transect)
MDS1<-data.frame(scores(mds.pa2))$NMDS1
MDS2<-data.frame(scores(mds.pa2))$NMDS2
Sample_Year<-as.numeric(paste(native.prairie$Sample_Year))
Transect<-native.prairie$Transect
np.scores<-data.frame(MDS1,MDS2,Sample_Year, Transect)

P14.scores<-droplevels(subset(np.scores, np.scores$Transect=="P14"))

envectors14<-envfit(P14.scores[,1:2] ~ P14.scores$Sample_Year, na.rm=TRUE)
head(envectors14)
#Sample_Year P=0.231

P15.scores<-droplevels(subset(np.scores, np.scores$Transect=="P15"))
envectors15<-envfit(P15.scores[,1:2] ~ P15.scores$Sample_Year, na.rm=TRUE)
head(envectors15)
#Sample_Year P=0.001

P17.scores<-droplevels(subset(np.scores, np.scores$Transect=="P17"))
envectors17<-envfit(P17.scores[,1:2] ~ P17.scores$Sample_Year, na.rm=TRUE)
head(envectors17)
#Sample_Year P=0.001

P19.scores<-droplevels(subset(np.scores, np.scores$Transect=="P19"))
envectors19<-envfit(P19.scores[,1:2] ~ P19.scores$Sample_Year, na.rm=TRUE)
head(envectors19)
#Sample_Year, P=0.006


#extract vectors, center them at origin, or centroid for each transect?
vectors.15<-data.frame(envectors15$vectors[1:4])
vectors.17<-data.frame(envectors17$vectors[1:4])
vectors.19<-data.frame(envectors19$vectors[1:4])
vectors.np<-rbind(vectors.15, vectors.17, vectors.19)
vectors.np

#add vector to graph
np.graph3<-ggplot.NMDS(mds.pa2, native.prairie$TransYear, colors.14)+annotate("text", x=1.1, y=-1.0, label="stress=0.185", size=6)+labs(title = "Native Prairie")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())+
  geom_segment(data=vectors.np, aes(x=0,xend=arrows.MDS1, y=0,yend=arrows.MDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)
np.graph3

#Nachusa labels
np.graph5<-ggplot.NMDS(mds.pa2, native.prairie$SiteYear, colors.14N)+annotate("text", x=1.1, y=-1.0, label="stress=0.185", size=6)+labs(title = "Native Prairie")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())+
  geom_segment(data=vectors.np, aes(x=0,xend=arrows.MDS1, y=0,yend=arrows.MDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)
np.graph5

#indicator species analysis
#summarize by species, use wide data
transect.np<-multipatt(native.prairie[,-c(1:5,432:433)], cluster=native.prairie$Transect, func = "IndVal", duleg=TRUE)
summary(transect.np)
#Be clear in methods this is using frequency, not cover!

year.np<-multipatt(native.prairie[,-c(1:5,432:433)], cluster=native.prairie$Sample_Year, func = "IndVal", duleg=TRUE)
summary(year.np)

transyear.np<-multipatt(native.prairie[,-c(1:5,432:433)], cluster=native.prairie$TransYear, func = "IndVal", duleg=TRUE)
summary(transyear.np)


#planted prairie
planted.prairie<-droplevels(subset(data.pa2, data.pa2$habitat=="prairie_restored" & data.pa2$Quadrat!="Q11"&data.pa2$Quadrat!="Q12"&data.pa2$Quadrat!="Q13"&data.pa2$Quadrat!="Q14"&data.pa2$Quadrat!="Q15"))
planted.prairie$Sample_Year<-as.factor(planted.prairie$Sample_Year)
planted.prairie$TransYear<-as.factor(paste(planted.prairie$Transect, planted.prairie$Sample_Year, sep=":"))
planted.prairie$SiteYear<-as.factor(paste(planted.prairie$site, planted.prairie$Sample_Year, sep=":"))

mds.pa2<-metaMDS(decostand(planted.prairie[,-c(1:5,432:433)],"pa" ),distance="jaccard", k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#stress for 2D=0.246
mds.pa3<-metaMDS(decostand(planted.prairie[,-c(1:5,432:433)],"pa" ),distance="jaccard", k=3,autotransform=FALSE, na.rm=TRUE)
mds.pa3
#stress for 3D=0.165, this is substantially lower than 2D, definitely use 3D
#scores into dataframe
MDS1<-data.frame(scores(mds.pa3))$NMDS1
MDS2<-data.frame(scores(mds.pa3))$NMDS2
MDS3<-data.frame(scores(mds.pa3))$NMDS3
Transect<-planted.prairie$Transect
Sample_Year<-planted.prairie$Sample_Year
TransYear<-planted.prairie$TransYear
SiteYear<-planted.prairie$SiteYear
rest.prairie.NMDS<-data.frame(MDS1, MDS2, MDS3, Transect, Sample_Year,TransYear, SiteYear)

#ggplot2 doesn't do 3d, use plotly
library(plotly)

#graph with the interaction
rest.prairie3D<-plot_ly(x=rest.prairie.NMDS$MDS1, y=rest.prairie.NMDS$MDS2, z=rest.prairie.NMDS$MDS3)
add_markers(rest.prairie3D, color= ~TransYear, colors=c(rev(brewer.pal(n=4, name = "YlOrRd")),rev(brewer.pal(n=3, name="Greens")),rev(brewer.pal(n=5, name="PuRd")),rev(brewer.pal(n=4,name="Blues"))))

#Nachusa labels
rest.prairie3DN<-plot_ly(x=rest.prairie.NMDS$MDS1, y=rest.prairie.NMDS$MDS2, z=rest.prairie.NMDS$MDS3)
add_markers(rest.prairie3D, color= ~SiteYear, colors=c(rev(brewer.pal(n=5, name = "YlOrRd")),rev(brewer.pal(n=4, name="Greens")),rev(brewer.pal(n=4, name="PuRd")),rev(brewer.pal(n=3,name="Blues"))))

#stats
adonis(planted.prairie[,-c(1:5,432:433)]~planted.prairie$Transect*planted.prairie$Sample_Year, permutations=999)
#strong interaction between transect and year

#Tukey's to pull apart differences
mds.dist3<-metaMDSdist(decostand(planted.prairie[,-c(1:5,432:433)], "pa"),k=3, index="jaccard", autotransform=FALSE,na.rm=TRUE)
rp.stat<-betadisper(mds.dist3, planted.prairie$TransYear, type="median")
rp.stat

#Now use TukeysHSD to contrast the median dispersion of the groups and find differences
TukeyHSD(rp.stat)
#some of these are statistically different, but the categorical comparisons are really meaningful, take the envfit approach

#Now layer this to run Sample_Year vectors for each transect, pull scores from full ordination by transect
P23<-droplevels(subset(planted.prairie, planted.prairie$Transect=="P23"))
adonis(P23[,-c(1:5,432:433)]~P23$Sample_Year, permutations=999)
#sample year strongly significant

#run envfit regressions with original NMDS scores, use rest.prairie.NMDS data frame
P23.scores<-droplevels(subset(rest.prairie.NMDS, rest.prairie.NMDS$Transect=="P23"))
P23.scores$Sample_Year<-as.numeric(paste(P23.scores$Sample_Year))
envectors23<-envfit(P23.scores[,1:3] ~ P23.scores$Sample_Year, choices=c(1:3), na.rm=TRUE)
head(envectors23)
#Sample_Year P=0.001

P26<-droplevels(subset(planted.prairie, planted.prairie$Transect=="P26"))
adonis(P26[,-c(1:5,432:433)]~P26$Sample_Year, permutations=999)
#Sample year highly significant

#run envfit regressions with original NMDS scores, use rest.prairie.NMDS data frame
P26.scores<-droplevels(subset(rest.prairie.NMDS, rest.prairie.NMDS$Transect=="P26"))
P26.scores$Sample_Year<-as.numeric(paste(P26.scores$Sample_Year))
envectors26<-envfit(P26.scores[,1:3] ~ P26.scores$Sample_Year, choices=c(1:3), na.rm=TRUE)
head(envectors26)
#Sample_Year P=0.001

P31<-droplevels(subset(planted.prairie, planted.prairie$Transect=="P31"))
adonis(P31[,-c(1:5,432:433)]~P31$Sample_Year, permutations=999)
#Sample_Year highly significant

P31.scores<-droplevels(subset(rest.prairie.NMDS, rest.prairie.NMDS$Transect=="P31"))
P31.scores$Sample_Year<-as.numeric(paste(P31.scores$Sample_Year))
envectors31<-envfit(P31.scores[,1:3] ~ P31.scores$Sample_Year, choices=c(1:3), na.rm=TRUE)
head(envectors31)
#Sample year P=0.001

P32<-droplevels(subset(planted.prairie, planted.prairie$Transect=="P32"))
adonis(P32[,-c(1:5,432:433)]~P32$Sample_Year, permutations=999)
#Sample_year highly significant

P32.scores<-droplevels(subset(rest.prairie.NMDS, rest.prairie.NMDS$Transect=="P32"))
P32.scores$Sample_Year<-as.numeric(paste(P32.scores$Sample_Year))
envectors32<-envfit(P32.scores[,1:3] ~ P32.scores$Sample_Year, choices=c(1:3), na.rm=TRUE)
head(envectors31)
#Sample_Year, P=0.001

#Sample year is a significant correlation for each transect, but likely going different directions
#extract vectors, center them at origin, or centroid for each transect?
origin<-c(0,0,0,0,0,0)
vectors.23<-data.frame(envectors23$vectors[1:4])
vectors.23c<-data.frame(rbind(vectors.23, origin))
vectors.26<-data.frame(envectors26$vectors[1:4])
vectors.26c<-data.frame(rbind(vectors.26, origin))
vectors.31<-data.frame(envectors31$vectors[1:4])
vectors.31c<-data.frame(rbind(vectors.31, origin))
vectors.32<-data.frame(envectors32$vectors[1:4])
vectors.32c<-data.frame(rbind(vectors.32, origin))

plot_ly() %>%
  add_markers(data=rest.prairie.NMDS, x=rest.prairie.NMDS$MDS1, y=rest.prairie.NMDS$MDS2, z=rest.prairie.NMDS$MDS3, color= ~TransYear, colors=c(rev(brewer.pal(n=4, name = "YlOrRd")),rev(brewer.pal(n=3, name="Greens")),rev(brewer.pal(n=5, name="PuRd")),rev(brewer.pal(n=4,name="Blues"))), mode="markers")%>%
  add_trace(data=vectors.23c, x=vectors.23c$arrows.MDS1, y=vectors.23c$arrows.MDS2, z=vectors.23c$arrows.MDS3,mode="line",type="scatter3d", color=I("red"), name="P23")%>%
  add_trace(data=vectors.26c, x=vectors.26c$arrows.MDS1, y=vectors.26c$arrows.MDS2, z=vectors.26c$arrows.MDS3,mode="line",type="scatter3d", color=I("green"), name="P26")%>%
  add_trace(data=vectors.31c, x=vectors.31c$arrows.MDS1, y=vectors.31c$arrows.MDS2, z=vectors.31c$arrows.MDS3,mode="line",type="scatter3d", color=I("purple"),name="P31")%>%
  add_trace(data=vectors.32c, x=vectors.32c$arrows.MDS1, y=vectors.32c$arrows.MDS2, z=vectors.32c$arrows.MDS3,mode="line",type="scatter3d", color=I("blue"), name="P32")

#indicator species analysis

#summarize by species, use wide data
transect.sp<-multipatt(planted.prairie[,-c(1:5,432:433)], cluster=planted.prairie$Transect, func = "IndVal", duleg=TRUE)
summary(transect.sp)

year.sp<-multipatt(planted.prairie[,-c(1:5,432:433)], cluster=planted.prairie$Sample_Year, func = "IndVal", duleg=TRUE)
summary(year.sp)

transyear.sp<-multipatt(planted.prairie[,-c(1:5,432:433)], cluster=planted.prairie$TransYear, func = "IndVal", duleg=TRUE)
summary(transyear.sp)

#savannas
savanna<-droplevels(subset(data.pa2, data.pa2$habitat=="savanna" & data.pa2$Quadrat!="Q11"& data.pa2$Quadrat!="Q12"& data.pa2$Quadrat!="Q13"& data.pa2$Quadrat!="Q14"& data.pa2$Quadrat!="Q15"& data.pa2$Quadrat!="Q16"))
savanna$Sample_Year<-as.factor(paste(savanna$Sample_Year))
savanna$TransYear<-as.factor(paste(savanna$Transect, savanna$Sample_Year, sep=":"))
savanna$SiteYear<-as.factor(paste(savanna$site, savanna$Sample_Year, sep=":"))

adonis(savanna[,-c(1:5,432:433)]~savanna$Transect*savanna$Sample_Year, permutations=999)
#strong interaction between transect and year
adonis(savanna[,-c(1:5,432:433)]~savanna$Transect, permutations=999)
#See if P41 groups out really weird
mds.dist<-metaMDSdist(decostand(savanna[,-c(1:5,432:433)], "pa"),k=2, index="jaccard", autotransform=FALSE,na.rm=TRUE)
savanna.stat<-betadisper(mds.dist, savanna$Transect, type="median")
savanna.stat

#Now use TukeysHSD to contrast the median dispersion of the groups and find differences
TukeyHSD(savanna.stat)
#P39-P08 are different
#P41 different from P08 and P24, but not P39
#keep it in the mix

mds.pa2<-metaMDS(decostand(savanna[,-c(1:5,432:433)],"pa" ),distance="jaccard", k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#stress=0.133

#graph with the interaction
colors.13<-c(rev(brewer.pal(n=3, name = "YlOrRd")),rev(brewer.pal(n=3, name="Greens")),rev(brewer.pal(n=4, name="PuRd")),rev(brewer.pal(n=3,name="Blues")))

savanna.graph<-ggplot.NMDS(mds.pa2, savanna$TransYear, colors.13)+annotate("text", x=-0.80, y=1.35, label="stress=0.133", size=5)+labs(title = "Savanna")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid = element_blank())
savanna.graph
#roughly x-axis is transect difference, y-axis is the sampling year difference

#Nachusa labels
savanna.graphN<-ggplot.NMDS(mds.pa2, savanna$SiteYear, colors.13)+annotate("text", x=-0.80, y=1.35, label="stress=0.133", size=5)+labs(title = "Savanna")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid = element_blank())
savanna.graphN

#stats to understand interaction
#TukeyHSD to untangle
mds.dist<-metaMDSdist(decostand(savanna[,-c(1:5,432:433)], "pa"),k=2, index="jaccard", autotransform=FALSE,na.rm=TRUE)
savanna.stat<-betadisper(mds.dist, savanna$TransYear, type="median")
savanna.stat

#Now use TukeysHSD to contrast the median dispersion of the groups and find differences
TukeyHSD(savanna.stat)
#some of these contrasts are significant, but not sure they mean anything, also, correcting for multiple comparisons likely elimnates most of them
#take the envfit approach

savanna$Sample_Year<-as.numeric(paste(savanna$Sample_Year))
envectors2<-envfit(mds.dist ~ savanna$Transect + savanna$Sample_Year,na.rm=TRUE)
envectors2
#Sample Year vector P=0.21! interaction is real, so break this down
#Transect centroids are strongly significant, P=0.001

#Pull out transect level data
P08<-droplevels(subset(savanna, savanna$Transect=="P08"))
adonis(P08[,-c(1:5,432:433)]~P08$Sample_Year, permutations=999)
#Sample_Year highly significant

P24<-droplevels(subset(savanna, savanna$Transect=="P24"))
adonis(P24[,-c(1:5,432:433)]~P24$Sample_Year, permutations=999)
#Sample_Year highly significant

P39<-droplevels(subset(savanna, savanna$Transect=="P39"))
adonis(P39[,-c(1:5,432:433)]~P39$Sample_Year, permutations=999)
#Sample_Year highly significant

P41<-droplevels(subset(savanna, savanna$Transect=="P41"))
adonis(P41[,-c(1:5,432:433)]~P41$Sample_Year, permutations=999)
#Sample_Year highly significant

#pull out NMDS scores and run regressions with envfit
MDS1<-data.frame(scores(mds.pa2))$NMDS1
MDS2<-data.frame(scores(mds.pa2))$NMDS2
Transect<-savanna$Transect
Sample_Year<-as.numeric(paste(savanna$Sample_Year))
TransYear<-savanna$TransYear
SiteYear<-savanna$SiteYear
savanna.NMDS<-data.frame(MDS1, MDS2, Transect, Sample_Year,TransYear,SiteYear)

P08.scores<-droplevels(subset(savanna.NMDS, savanna.NMDS$Transect=="P08"))
envectors08<-envfit(P08.scores[,1:2] ~ P08.scores$Sample_Year, na.rm=TRUE)
head(envectors08)
#Sample_Year P=0.001

P24.scores<-droplevels(subset(savanna.NMDS, savanna.NMDS$Transect=="P24"))
envectors24<-envfit(P24.scores[,1:2] ~ P24.scores$Sample_Year, na.rm=TRUE)
head(envectors24)
#Sample_Year P=0.001

P39.scores<-droplevels(subset(savanna.NMDS, savanna.NMDS$Transect=="P39"))
envectors39<-envfit(P39.scores[,1:2] ~ P39.scores$Sample_Year, na.rm=TRUE)
head(envectors39)
#Sample_Year P=0.001

P41.scores<-droplevels(subset(savanna.NMDS, savanna.NMDS$Transect=="P41"))
envectors41<-envfit(P41.scores[,1:2] ~ P41.scores$Sample_Year, na.rm=TRUE)
head(envectors41)
#Sample_Year P=0.016

#interaction may be because P41 is not significant
#others may also be going different directions, so follow-up on that

#extract vectors, center them at origin, or centroid for each transect?
vectors.08<-data.frame(envectors08$vectors[1:4])
vectors.24<-data.frame(envectors24$vectors[1:4])
vectors.39<-data.frame(envectors39$vectors[1:4])
vectors.41<-data.frame(envectors41$vectors[1:4])
vectors.sav<-rbind(vectors.08, vectors.24, vectors.39, vectors.41)
vectors.sav

#add vector to graph
savanna.graph3<-ggplot.NMDS(mds.pa2, savanna$TransYear, colors.13)+annotate("text", x=1.1, y=-1.0, label="stress=0.133", size=6)+labs(title = "Savanna")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())+
  geom_segment(data=vectors.sav, aes(x=0,xend=arrows.MDS1, y=0,yend=arrows.MDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)
savanna.graph3

#P08 and P39 are moving in a similar direction, and P24 & P41 are different

#Nachusa labels
savanna.graph3N<-ggplot.NMDS(mds.pa2, savanna$SiteYear, colors.13)+annotate("text", x=1.1, y=-1.0, label="stress=0.133", size=6)+labs(title = "Savanna")+theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), panel.grid= element_blank())+
  geom_segment(data=vectors.sav, aes(x=0,xend=arrows.MDS1, y=0,yend=arrows.MDS2),arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",size=1,inherit_aes=FALSE)
savanna.graph3N

#indicator species analysis
#summarize by species, use wide data
transect.sav<-multipatt(savanna[,-c(1:5,432:433)], cluster=savanna$Transect, func = "IndVal", duleg=TRUE)
summary(transect.sav)

year.sav<-multipatt(savanna[,-c(1:5,432:433)], cluster=savanna$Sample_Year, func = "IndVal", duleg=TRUE)
summary(year.sav)

transyear.sav<-multipatt(savanna[,-c(1:5,432:433)], cluster=savanna$TransYear, func = "IndVal", duleg=TRUE)
summary(transyear.sav)
