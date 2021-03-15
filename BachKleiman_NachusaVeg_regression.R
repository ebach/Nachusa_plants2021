#Elizabeth Bach
#Nachusa Grassland long-term vegetation monitoring
#Regression analysis
#March 2021

rm(list=ls())
library(reshape2)
library(plyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(data.table)

#read in data, "Nachusa_veg_transect_data1994_2020_widev2.csv"
Nachusa.data<-read.csv(file.choose(),header=T, na.strings="NA")

#diversity metrics, FQI, 
#calculating richness, %native, mean C
#diversity metrics, per transect

Nachusa.data2<-unite(Nachusa.data, Trans_Year, c(6,4))
Nachusa.data2$Trans_Year<-as.factor(Nachusa.data2$Trans_Year)

richness.trans<-as.data.frame(specnumber(Nachusa.data2[,-c(1:9)], groups=Nachusa.data2$Trans_Year))
range(richness.trans)
max(richness.trans)
#max is 83 species in P86 in 2009, this was the first year of growth for the Holland 2008 planting, original data shows 79 for total species
#most of these are spot on, a few are off by 1-3 species, nomenclature thing?

rownames<-rownames(richness.trans)
transect.richness2<-cbind(rownames, data.frame(richness.trans, row.names=NULL))
colnames(transect.richness2)[2]<-"species.num"
quadrat.size<-unique(Nachusa.data2[,c(4,8)])
transect.richness2.5<-merge(transect.richness2, quadrat.size, by.x="rownames", by.y = "Trans_Year", all.y=FALSE)
transect.richness3<-separate(transect.richness2.5, 1, sep="_", into=c("transect","year"), extra="merge")
transect.richness3$species.m<-transect.richness3$species.num/transect.richness3$quadrat_size_m
transect.richness3$year<-as.factor(transect.richness3$year)
transect.richness3$transect<-as.factor(transect.richness3$transect)

#translate data into pa only to double-check
Nachusa.data.pa<-Nachusa.data2[,-c(1:9)]
Nachusa.data.pa[Nachusa.data.pa>0]<-1
Nachusa.data.pa2<-cbind(Nachusa.data2[,2:8], Nachusa.data.pa)

richness.trans.pa<-as.data.frame(specnumber(Nachusa.data.pa2[,-c(1:7)], groups=Nachusa.data.pa2$Trans_Year))
range(richness.trans.pa)
max(richness.trans.pa)
rownames<-rownames(richness.trans.pa)
transect.richness.pa2<-cbind(rownames, data.frame(richness.trans.pa, row.names=NULL))
colnames(transect.richness.pa2)[2]<-"species.num"
richness.trans.pa3<-separate(transect.richness.pa2, 1, sep="_", into=c("year","transect"), extra="merge")
#give exactly the same numbers as full dataset, so not a coding or data type issue

#filter out to look at transect with 3 or more observations
sample.n<-ddply(transect.richness3, .(transect), summarize, N=length(year))
sample.replicated<-droplevels(subset(sample.n, sample.n$N>=3))
sample.replicated2<-droplevels(subset(sample.replicated, sample.replicated$transect!="Didier"))

richness.trans4<-merge(sample.replicated2, transect.richness3, by="transect", all.x=TRUE, all.y = FALSE)
richness.trans4$transect<-as.factor(richness.trans4$transect)

#graph
#use species.num, NOT species.m, which extrapolates 0.5m plots to 1m plots
total.richness2<-ggplot(richness.trans4, aes(year,species.num))+geom_point()+facet_wrap(~transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
total.richness2

#regression analysis
#not using mixed models, as there isn't really a random effect within the transects
richness.dt<-data.table(richness.trans4)
#year can't be a factor here
richness.dt$year<-as.numeric(paste(richness.dt$year))
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}
regression.results = richness.dt[,linreg(year ~ species.num),by=transect]
regression.results

#only P19 is significant, P=0.036
#P41 is marginal, P=0.057
#see publication graph at end of code

#%native cover, mean C
#use long form data structure
#"Nachusa_vegdata_long_final.csv"
Nachusa.data.long<-read.csv(file.choose(),header=T, na.strings="NA")
Nachusa.data.long[Nachusa.data.long=="NA"]<-NA

#native
#for the full dataset
unique(Nachusa.data.long$species)
native<-ddply(Nachusa.data.long, .(Transect, Sample_Year,Quadrat, Native.), summarize, n=sum(value))

#skim out 0's, this works by skimming out all the unobserved values. Pay attention with downstream analysis, this is ok for these calcs, but might be an issue elsewhere
Nachusa.data.long.pa3<-droplevels(subset(Nachusa.data.long, Nachusa.data.long$value>0))
native4<-ddply(Nachusa.data.long.pa3, .(Transect, Sample_Year, Native.), summarize, n=length(unique(species)))
p.native<-dcast(native4, Transect + Sample_Year ~ Native., fill = 0)
p.native$p.native<-p.native$native/(p.native$native+p.native$`non-native`)

#graph p.native
p.natives<-ggplot(p.native, aes(Sample_Year,p.native))+geom_point()+facet_wrap(~Transect)
p.natives
#publication graph at end of code

#regression analysis
#proportion native
native.dt<-data.table(p.native)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = native.dt[,linreg(Sample_Year~p.native),by=Transect]
regression.results
#P31 is only statistically significant, P=0.002

#mean C
C.calc<-ddply(Nachusa.data.long.pa3, .(Transect, Sample_Year, Quadrat), summarize, n=length(unique(species)), m.C=sum(C/n))
mean.C<-ddply(C.calc,.(Transect, Sample_Year), summarize, mean.C=mean(m.C), SE=sd(m.C)/sqrt(length(n)), max=mean.C+SE, min=mean.C-SE)

#graph mean C
mean.C.graph<-ggplot(mean.C, aes(Sample_Year,mean.C))+geom_pointrange(aes(x=Sample_Year, ymin=min, ymax=max), stat="identity")+facet_wrap(~Transect)
mean.C.graph
#publication graph at end of code

#regression analysis
meanC.dt<-data.table(mean.C)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = meanC.dt[,linreg(Sample_Year~mean.C),by=Transect]
regression.results
#P17, P=0.035
#P31, P=0.004
#P32, P=0.033

#Graphs
#graph titles for transects
library(stringr)
plot.titles<-c("Savanna","Wetland restored 1995","Native Prairie","Native Prairie","Native Prairie","Native Prairie","Prairie planted 1991", "Prairie planted 1994","Prairie planted 1986","Prairie planted 1987", "Savanna","Savanna")
plot.titles2<-str_wrap(plot.titles, width=15)
names(plot.titles2)<-c("P08","P13","P14","P15","P17","P19","P23", "P26","P31","P32","P39","P41")

#for publication

#total richness
richness.dt$transect<-factor(richness.dt$transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P39","P41","P13"))
total.rich.graph3<-ggplot(richness.dt, aes(year,species.num))+geom_point(size=3)+facet_wrap(~transect, labeller = labeller(transect=plot.titles2))+ylab("total plant species")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(richness.dt, richness.dt$transect=="P19"), aes(year,species.num),method = "lm", se=FALSE, color="darkgreen", size=2)
total.rich.graph3

#prop native
p.native$Transect<-factor(p.native$Transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P39","P41","P13"))
p.native2<-as.data.frame(p.native)
NativeP.graph4<-ggplot(p.native2, aes(Sample_Year,p.native))+geom_point(size=3)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("proportion native plant species")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(p.native2, p.native2$Transect=="P31"), aes(Sample_Year, p.native),method = "lm", se=FALSE, color="darkgreen", size=2)
NativeP.graph4

#mean C
mean.C$Transect<-factor(mean.C$Transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P39","P41","P13"))
mean.C2<-as.data.frame(mean.C)
meanC.graph3<-ggplot(mean.C2, aes(Sample_Year,mean.C))+geom_point(size=3)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("mean coefficient of conservatism")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P17"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P31"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P32"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)
meanC.graph3


