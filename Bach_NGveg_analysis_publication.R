#Elizabeth Bach
#Nachusa Grasslands long-term veg data
#Analysis for ESE Decade on Restoration manuscript
#Dec 2020
#revisions July 2021
#double-checked/cleaned September 2021

rm(list=ls())
library(plyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(reshape2)

#for total species richness, read in wide-form data, "Nachusa_veg_transect_data1994_2016.csv"
Nachusa.data<-read.csv(file.choose(),header=T, na.strings="NA")

native<-ddply(Nachusa.data, .(Transect, Sample_Year,Quadrat, native2016), summarize, n=sum(cover_value))

#skim out 0's, this works by skimming out all the unobserved values. Pay attention with downstream analysis, this is ok for these calcs, but might be an issue elsewhere
Nachusa.data.long2<-droplevels(subset(Nachusa.data, Nachusa.data$cover_value>0))
native4<-ddply(Nachusa.data.long2, .(Transect, Sample_Year, native2016), summarize, n=length(unique(Scientific_Name)))
p.native<-dcast(native4, Transect + Sample_Year ~ native2016, fill = 0, value.var = "n")
p.native$p.native<-p.native$native/(p.native$native+p.native$`non-native`)

#graph changes in total native species (number)
total.natives<-ggplot(p.native, aes(Sample_Year,native))+geom_point()+facet_wrap(~Transect)
total.natives

#filter out to look at transect with 3 or more observations
#this has already been done, but for some reason, merging this with p.native categorizes the columns as recursive vectors, which can be graphed
#also adds "N" value to data frame
sample.n<-ddply(p.native, .(Transect), summarize, N=length(native))
sample.replicated<-droplevels(subset(sample.n, sample.n$N>=3))
#filter out wetland, P13
sample.replicated.pub<-droplevels(subset(sample.replicated, sample.replicated$Transect!="P13"))


#graph titles for transects
library(stringr)
plot.titles<-c("Savanna","Native Prairie","Native Prairie","Native Prairie","Native Prairie","Prairie planted 1991", "Savanna","Prairie planted 1994","Prairie planted 1986","Prairie planted 1987", "Savanna","Savanna")
plot.titles2<-str_wrap(plot.titles, width=15)
names(plot.titles2)<-c("P08","P14","P15","P17","P19","P23","P24", "P26","P31","P32","P39","P41")

p.native2<-merge(sample.replicated.pub, p.native, by="Transect", all.x=TRUE, all.y = FALSE)
p.native2$Transect<-factor(p.native2$Transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P24","P39","P41"))

natives.graph<-ggplot(p.native2, aes(Sample_Year,native))+geom_point()+facet_wrap(~Transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
natives.graph

#regression analysis
library(data.table)
#total native
native.dt<-data.table(p.native2)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], intercept=coefficients(m)[1], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = native.dt[,linreg(Sample_Year~native),by=Transect]
regression.results

#only those with Pvals >0.05
regression.results2<-droplevels(subset(regression.results, regression.results$p.val<=0.06))
regression.results2

#proportion native in community
#proportion native
regression.results.p = native.dt[,linreg(Sample_Year~p.native),by=Transect]
regression.results.p

#add regression lines to graph
#for poster
NativeP.graph3<-ggplot(p.native2, aes(Sample_Year,p.native))+geom_point(size=32)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("proportion native plant species")+
  theme(axis.line=element_line(colour="black", size=24), axis.ticks=element_line(colour="black", size=15), aspect.ratio=1, axis.text.y=element_text(colour="black", size=200), axis.text.x=element_text(colour="black", size=200, angle=90, hjust=1), axis.ticks.length=unit(3,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=10), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=250, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=150, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(p.native2, p.native2$Transect=="P31"), aes(Sample_Year, p.native),method = "lm", se=FALSE, color="darkgreen", size=20)
NativeP.graph3

#for publication
NativeP.graph4<-ggplot(p.native2, aes(Sample_Year,p.native))+geom_point(size=3)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("proportion native plant species")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(p.native2, p.native2$Transect=="P31"), aes(Sample_Year, p.native),method = "lm", se=FALSE, color="darkgreen", size=2)
NativeP.graph4


#C-values
#filter out unidentified species, species with missing Cvalues
Nachusa.data.long.C<-droplevels(subset(Nachusa.data, !is.na(Nachusa.data$C2016)))
#filter out 0 values, so only present species counted
Nachusa.data.long.C2<-droplevels(subset(Nachusa.data.long.C, Nachusa.data.long.C$cover_value>0))

C.calc<-data.frame(ddply(Nachusa.data.long.C2, .(Transect, Sample_Year, Quadrat), summarize, n=length(unique(Scientific_Name)), m.C=sum(C2016/n)))
mean.C<-data.frame(ddply(C.calc,.(Transect, Sample_Year), summarize, mean.C=mean(m.C), SE=sd(m.C)/sqrt(length(n)), max=mean.C+SE, min=mean.C-SE))


#filter out to look at transect with 3 or more observations
#this has already been done, but for some reason, merging this with p.native categorizes the columns as recursive vectors, which can be graphed
#also adds "N" value to data frame
C.n<-ddply(mean.C, .(Transect), summarize, N=length(mean.C))
C.replicated<-droplevels(subset(C.n, sample.n$N>=3))
#filter out wetland, P13
C.replicated.pub<-droplevels(subset(C.replicated, C.replicated$Transect!="P13"))

mean.C2<-merge(C.replicated.pub, mean.C, by="Transect", all.x=TRUE, all.y = FALSE)
mean.C2$Transect<-factor(mean.C2$Transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P24","P39","P41"))

#regression analysis
meanC.dt<-data.table(mean.C2)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results.C = meanC.dt[,linreg(Sample_Year~mean.C),by=Transect]
regression.results.C

#graph for poster
meanC.graph2<-ggplot(mean.C2, aes(Sample_Year,mean.C))+geom_point(size=32)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("mean coefficient of conservatism")+
  theme(axis.line=element_line(colour="black", size=24), axis.ticks=element_line(colour="black", size=15), aspect.ratio=1, axis.text.y=element_text(colour="black", size=200), axis.text.x=element_text(colour="black", size=200, angle=90, hjust=1), axis.ticks.length=unit(3,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=10), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=250, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=150, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P17"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=20)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P31"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=20)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P32"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=20)
meanC.graph2

#graphfor publication
meanC.graph3<-ggplot(mean.C2, aes(Sample_Year,mean.C))+geom_point(size=3)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("mean coefficient of conservatism")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P17"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P31"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)+geom_smooth(data=subset(mean.C2, mean.C2$Transect=="P32"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)
meanC.graph3


#total richness from long form data
richness<-ddply(Nachusa.data.long2, .(Transect, Sample_Year), summarize, rich=length(unique(Scientific_Name)))
richness$rich<-as.numeric(richness$rich)

#skim down to >3 observations
sample.rich<-ddply(richness, .(Transect), summarize, N=length(rich))
sample.rich.rep<-droplevels(subset(sample.rich, sample.rich$N>=3))
#filter out wetland, P13
sample.rich.rep.pub<-droplevels(subset(sample.rich.rep, sample.rich.rep$Transect!="P13"))

richness2<-merge(sample.rich.rep.pub, richness, by="Transect", all.x=TRUE, all.y = FALSE)
richness2$Transect<-factor(richness2$Transect, levels = c("P14","P15","P17","P19","P31","P32","P23","P26","P08","P24","P39","P41"))

#total plant species
richness.dt<-data.table(richness2)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], intercept=coefficients(m)[1], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1], p.val=summary(m)$coefficients[2,4])
}

regression.results.rich = richness.dt[,linreg(Sample_Year~rich),by=Transect]
regression.results.rich

#subset for significant lines
rich.sig<-droplevels(subset(richness2, richness2$Transect=="P19"))


#graph for poster
total.rich.graph<-ggplot(richness2, aes(Sample_Year,rich))+geom_point(size=32)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("total number of plant species")+
  theme(axis.line=element_line(colour="black", size=24), axis.ticks=element_line(colour="black", size=15), aspect.ratio=1, axis.text.y=element_text(colour="black", size=200), axis.text.x=element_text(colour="black", size=200, angle=90, hjust=1), axis.ticks.length=unit(3,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=10), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=250, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=150, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(richness2, richness2$Transect=="P19"), aes(Sample_Year, rich),method = "lm", se=FALSE, color="darkgreen", size=20)+geom_smooth(data=subset(richness2, richness2$Transect=="P41"), aes(Sample_Year, rich),method = "lm", se=FALSE, color="darkgreen", size=20)
total.rich.graph

#graph for publication
total.rich.graph3<-ggplot(richness2, aes(Sample_Year,rich))+geom_point(size=3)+facet_wrap(~Transect, labeller = labeller(Transect=plot.titles2))+ylab("total plant species")+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(richness2, richness2$Transect=="P19"), aes(Sample_Year, rich),method = "lm", se=FALSE, color="darkgreen", size=2)
total.rich.graph3

#Look more closely at Doug's Knob, is there really this big increase in species?
P.19<-droplevels(subset(Nachusa.data, Nachusa.data$Transect=="P19"))
P.19.observed<-droplevels(subset(P.19, P.19$value>0))
#summarize by year
P.19.sum<-ddply(P.19.observed, .(Sample_Year), summarize, rich=length(unique(species)))
#1994 is 17, 2013 is 43fq
#This matches the raw data, so accept it

#Revise to create figure grouped by habitat type, not response variable
#bind response variables into single data frame per habitat
np.richness<-droplevels(subset(richness2, richness2$Transect=="P14"|richness2$Transect=="P15"|richness2$Transect=="P17"|richness2$Transect=="P19"))
np.propN<-droplevels(subset(p.native2, p.native2$Transect=="P14"|p.native2$Transect=="P15"|p.native2$Transect=="P17"|p.native2$Transect=="P19"))
np.meanC<-droplevels(subset(mean.C2, mean.C2$Transect=="P14"|mean.C2$Transect=="P15"|mean.C2$Transect=="P17"|mean.C2$Transect=="P19"))

pp.richness<-droplevels(subset(richness2, richness2$Transect=="P23"|richness2$Transect=="P26"|richness2$Transect=="P31"|richness2$Transect=="P32"))
pp.propN<-droplevels(subset(p.native2, p.native2$Transect=="P23"|p.native2$Transect=="P26"|p.native2$Transect=="P31"|p.native2$Transect=="P32"))
pp.meanC<-droplevels(subset(mean.C2, mean.C2$Transect=="P23"|mean.C2$Transect=="P26"|mean.C2$Transect=="P31"|mean.C2$Transect=="P32"))

sav.richness<-droplevels(subset(richness2, richness2$Transect=="P08"|richness2$Transect=="P24"|richness2$Transect=="P39"|richness2$Transect=="P41"))
sav.propN<-droplevels(subset(p.native2, p.native2$Transect=="P08"|p.native2$Transect=="P24"|p.native2$Transect=="P39"|p.native2$Transect=="P41"))
sav.meanC<-droplevels(subset(mean.C2, mean.C2$Transect=="P08"|mean.C2$Transect=="P24"|mean.C2$Transect=="P39"|mean.C2$Transect=="P41"))


#Native prairie
prairie.rich<-ggplot(np.richness, aes(Sample_Year,rich))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("total plant species")+labs(title="Native Prairie")+
  ylim(0,50)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(np.richness, np.richness$Transect=="P19"), aes(Sample_Year, rich),method = "lm", se=FALSE, color="darkgreen", size=2)
prairie.rich

prairie.propN<-ggplot(np.propN, aes(Sample_Year,p.native))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("p(native species)")+
  ylim(0,1)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
prairie.propN

prairie.meanC<-ggplot(np.meanC, aes(Sample_Year,mean.C))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("mean C")+
  ylim(0,8)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(np.meanC, np.meanC$Transect=="P17"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)
prairie.meanC

library(gridExtra)
gRich<-ggplotGrob(prairie.rich)
gProp<-ggplotGrob(prairie.propN)
gMeanC<-ggplotGrob(prairie.meanC)
grid::grid.newpage()
grid::grid.draw(rbind(gRich, gProp, gMeanC))

#planted prairie
restprairie.rich<-ggplot(pp.richness, aes(Sample_Year,rich))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("total plant species")+labs(title="Planted Prairie")+
  ylim(0,50)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
restprairie.rich

restprairie.propN<-ggplot(pp.propN, aes(Sample_Year,p.native))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("p(native species)")+
  ylim(0,1)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(pp.propN, pp.propN$Transect=="P31"), aes(Sample_Year, p.native),method = "lm", se=FALSE, color="darkgreen", size=2)
restprairie.propN

restprairie.meanC<-ggplot(pp.meanC, aes(Sample_Year,mean.C))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("mean C")+
  ylim(0,8)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))+
  geom_smooth(data=subset(pp.meanC, pp.meanC$Transect=="P31"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)+geom_smooth(data=subset(pp.meanC, pp.meanC$Transect=="P32"), aes(Sample_Year, mean.C),method = "lm", se=FALSE, color="darkgreen", size=2)
restprairie.meanC

gRichP<-ggplotGrob(restprairie.rich)
gPropP<-ggplotGrob(restprairie.propN)
gMeanCP<-ggplotGrob(restprairie.meanC)
grid::grid.newpage()
grid::grid.draw(rbind(gRichP, gPropP, gMeanCP))


#Savanna
sav.rich<-ggplot(sav.richness, aes(Sample_Year,rich))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("total plant species")+labs(title="Savanna")+
  ylim(0,50)+xlim(1994, 2015)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
sav.rich

sav.propN2<-ggplot(sav.propN, aes(Sample_Year,p.native))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("p(native species)")+
  ylim(0,1)+xlim(1994, 2015)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_blank(), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
sav.propN2

sav.meanC2<-ggplot(sav.meanC, aes(Sample_Year,mean.C))+geom_point(size=3)+facet_wrap(~Transect, nrow=1)+ylab("mean C")+
  ylim(0,8)+xlim(1994, 2015)+theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=18, face="bold"), axis.title.x=element_blank(), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
sav.meanC2

gRichS<-ggplotGrob(sav.rich)
gPropS<-ggplotGrob(sav.propN2)
gMeanCS<-ggplotGrob(sav.meanC2)
grid::grid.newpage()
grid::grid.draw(rbind(gRichS, gPropS, gMeanCS))
