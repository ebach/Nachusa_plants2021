#Elizabeth Bach
#Nachusa Grassland long-term vegetation monitoring
#Revising figures for publication
#28 June 2021

rm(list=ls())
library(reshape2)
library(plyr)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(data.table)

#read in data, "Nachusa_veg_transect_data1994_2020_wide.csv"
Nachusa.data<-read.csv(file.choose(),header=T, na.strings="NA")

#diversity metrics, FQI, 
#calculating richness, shannons, and evenness

#diversity metrics, per quadrat
richness<-rowSums(Nachusa.data[,-c(1:10)]>0)
head(richness)
hist(richness)



#do not have abundance data, cannot do Shannon's, Simpson's should work
simpson<-diversity(Nachusa.data[,-c(1:10)], "invsimpson")
hist(simpson)
#using reciprical Simpsons, which is 1/D, this produces essentially a look at total species, max being the total community
#These values range from 1-25
#evenness is not useful without abundance data

#diversity metrics, per m2
#Nachusa.data$richness.m<-(richness/Nachusa.data$quadrat_size_m)
#Nachusa.data$simpson.m<-simpson/Nachusa.data$quadrat_size_m

#diversity metrics, per transect
##Pick up here!!##
Nachusa.data2<-unite(Nachusa.data, Trans_Year, c(4,6,7))
Nachusa.data2$Trans_Year<-as.factor(Nachusa.data2$Trans_Year)

richness.trans<-as.data.frame(specnumber(Nachusa.data2[,-c(1:8)], groups=Nachusa.data2$Trans_Year))
range(richness.trans)
max(richness.trans)
#max is 83 species in P86 in 2009, this was the first year of growth for the Holland 2008 planting, original data shows 79 for total species
#most of these are spot on, a few are off by 1-3 species, nomenclature thing?

rownames<-rownames(richness.trans)
transect.richness2<-cbind(rownames, data.frame(richness.trans, row.names=NULL))
colnames(transect.richness2)[2]<-"species.num"
quadrat.size<-unique(Nachusa.data2[,c(4,7)])
transect.richness2.5<-merge(transect.richness2, quadrat.size, by.x="rownames", by.y = "Trans_Year", all.y=FALSE)
transect.richness3<-separate(transect.richness2.5, 1, sep="_", into=c("year","transect"), extra="merge")
transect.richness3$species.m<-transect.richness3$species.num/transect.richness3$quadrat_size_m
transect.richness3$year<-as.factor(transect.richness3$year)
transect.richness3$transect<-as.factor(transect.richness3$transect)

#translate data into pa only
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

#graph changes in richness
total.richness<-ggplot(transect.richness3, aes(year,species.num))+geom_point()+facet_wrap(~transect)
total.richness

#filter out to look at transect with 3 or more observations
sample.n<-ddply(transect.richness3, .(transect), summarize, N=length(year))
sample.replicated<-droplevels(subset(sample.n, sample.n$N>=3))

richness.trans4<-merge(sample.replicated, transect.richness3, by="transect", all.x=TRUE, all.y = FALSE)
richness.trans4$transect<-as.factor(richness.trans4$transect)

#graph
total.richness2<-ggplot(richness.trans4, aes(year,species.num))+geom_point()+facet_wrap(~transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
total.richness2

#regression analysis
#not using mixed models, as there isn't really a random effect within the transects, may revise for final poster
richness.dt<-data.table(richness.trans4)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = richness.dt[,linreg(year ~ species.num),by=transect]
regression.results


#?Indicator species?

#%native cover, FQI, mean C
#use long form data structure
#"Nachusa_veg_transect_data1994-2020.csv"
Nachusa.data.long<-read.csv(file.choose(),header=T, na.strings="NA")
Nachusa.data.long[Nachusa.data.long=="NA"]<-NA
#transform so all is presence/absence
Nachusa.data.long.pa<-Nachusa.data.long
Nachusa.data.long.pa$cover_value[Nachusa.data.long.pa$cover_value>0]<-1
range(Nachusa.data.long.pa$cover_value)

#figure out where the NA values for native are coming from
#levels(Nachusa.data.long.pa$native.2016)
#Didier.23<-droplevels(subset(Nachusa.data.long.pa, Nachusa.data.long.pa$Transect=="Didier_23"))

#may need to go in a handle some subspecies, trying to pull out those not hitting 2016 database from plant list
#na.native<-droplevels(subset(Nachusa.data.long.pa, is.na(Nachusa.data.long.pa$native.2016)))
#na.native<-droplevels(subset(data.full7, is.na(data.full7$native.2016)))
#species<-levels(na.native$newcol)
#species
#the only species with NAs for both 1994 an 2016 are truly NAs. I think I've fixed all the others


#physiognomy, native, and duration can use the 1994 data without revision, check if any NAs in those coulums
#is.na(na.native2$Native.1994)
#hmm, subset these and see what's going on, may be a mis-spelling
#NAs.1994<-droplevels(subset(na.native2, is.na(na.native2$Native.1994)))

#correct misspellings, do this to correct Nachusa.data.long in wrangling code
#'Anenome virginiana','Anemone virginiana'
#'Artemisia campestris','Artemisia campestris subsp. caudata'
#Plant list flips Botrypus virginianus back to Botrychium, which is not in 2016 file, so keep that
#see also Baptisia lactea, 
#Plant list reverting seveal Dichantheliums back to panicum, I guess fix this 
#not sure what's up with that one observation of Gaura biennis?? Granger_08? Maybe there's a space


#Unknown species will be NAs: Aster NA, Carex NA, Bromus NA, 

#filter out true NAs, as they can't count in these metrics
#pick up here, I've cleaned up the underlying data, need to dig in and think through what's happening here
Nachusa.data.long.pa2<-droplevels(subset(Nachusa.data.long.pa, !is.na(Nachusa.data.long.pa$native.2016)|!is.na(Nachusa.data.long.pa$C.1994)))
range(Nachusa.data.long.pa2$cover_value)

#subset out P14 - Schafer to see how/if numbers match
Nachusa.P14<-droplevels(subset(Nachusa.data.long.pa2, Nachusa.data.long.pa2$Transect=="P14"))
range(Nachusa.P14$cover_value)
#native<-ddply(Nachusa.P14, .(Transect, Sample_Year,Quadrat,newcol, native.2016), summarize, n=sum(cover_value))
native2<-ddply(Nachusa.P14, .(Transect, Sample_Year,native.2016), summarize, n=length(unique(newcol)))
#this seems to be working, had to prevent double counting between the quadrats
#but not quite, some species missing, at least from 1994
P14.2011<-droplevels(subset(Nachusa.data.long.pa2, Nachusa.data.long.pa2$Sample_Year==2011&Nachusa.data.long.pa2$Transect=="P14"))
levels(Nachusa.P14$newcol)

D.lind<-droplevels(subset(Nachusa.data.long, Nachusa.data.long$Taxon=="Dichanthelium lindheimeri"))

#see if there are others Plant List is not updating to current taxonomy for 2016 values
NA.2016<-droplevels(subset(Nachusa.data.long.pa,is.na(Nachusa.data.long.pa$native.2016)))
levels(NA.2016$newcol)


#OK, back to the %native
p.native<-dcast(native2, Transect + Sample_Year ~ native.2016, fill = 0)
p.native$p.native<-p.native$native/(p.native$native+p.native$`non-native`)

#native
#for the full dataset
unique(Nachusa.data.long.pa2$newcol)
native<-ddply(Nachusa.data.long.pa2, .(Transect, Sample_Year,Quadrat, native.2016), summarize, n=sum(cover_value))
native3<-ddply(Nachusa.data.long.pa2, .(Transect, Sample_Year,Quadrat, native.2016), summarize, n=length(unique(newcol)))

#skim out 0's, this works by skimming out all the unobserved values. Pay attention with downstream analysis, this is ok for these calcs, but might be an issue elsewhere
Nachusa.data.long.pa3<-droplevels(subset(Nachusa.data.long.pa2, Nachusa.data.long.pa2$cover_value>0))
native4<-ddply(Nachusa.data.long.pa3, .(Transect, Sample_Year, native.2016), summarize, n=length(unique(newcol)))
p.native<-dcast(native4, Transect + Sample_Year ~ native.2016, fill = 0)
p.native$p.native<-p.native$native/(p.native$native+p.native$`non-native`)

#graph changes in total native species (number)
total.natives<-ggplot(p.native, aes(Sample_Year,native))+geom_point()+facet_wrap(~Transect)
total.natives

#filter out to look at transect with 3 or more observations
sample.n<-ddply(p.native, .(Transect), summarize, N=length(native))
sample.replicated<-droplevels(subset(sample.n, sample.n$N>=3))

p.native2<-merge(sample.replicated, p.native, by="Transect", all.x=TRUE, all.y = FALSE)

natives.graph<-ggplot(p.native2, aes(Sample_Year,native))+geom_point()+facet_wrap(~Transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
natives.graph

#graph p.native
p.natives<-ggplot(p.native, aes(Sample_Year,p.native))+geom_point()+facet_wrap(~Transect)
p.natives

#multiple observations
p.native.graph<-ggplot(p.native2, aes(Sample_Year,p.native))+geom_point()+facet_wrap(~Transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
p.native.graph

#regression analysis
#total native
native.dt<-data.table(p.native2)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = native.dt[,linreg(Sample_Year~native),by=Transect]
regression.results

#proportion native
regression.results = native.dt[,linreg(Sample_Year~p.native),by=Transect]
regression.results

#mean C
C.calc<-ddply(Nachusa.data.long.pa3, .(Transect, Sample_Year, Quadrat), summarize, n=length(unique(newcol)), m.C=sum(C.2016/n))
mean.C<-ddply(C.calc,.(Transect, Sample_Year), summarize, mean.C=mean(m.C), SE=sd(m.C)/sqrt(length(n)), max=mean.C+SE, min=mean.C-SE)

#graph mean C
mean.C.graph<-ggplot(mean.C, aes(Sample_Year,mean.C))+geom_pointrange(aes(x=Sample_Year, ymin=min, ymax=max), stat="identity")+facet_wrap(~Transect)
mean.C.graph

#filter out to look at transect with 3 or more observations
sample.n<-ddply(mean.C, .(Transect), summarize, N=length(mean.C))
sample.replicated<-droplevels(subset(sample.n, sample.n$N>=3))

mean.C2<-merge(sample.replicated, mean.C, by="Transect", all.x=TRUE, all.y = FALSE)

meanC.graph<-ggplot(mean.C2, aes(Sample_Year,mean.C))+geom_pointrange(aes(x=Sample_Year, ymin=min, ymax=max), stat="identity")+facet_wrap(~Transect)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, axis.text.y=element_text(colour="black", size=12), axis.text.x=element_text(colour="black", size=10, angle=90, hjust=1), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=0.5), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title.y=element_text(colour="black", size=12), axis.title.x=element_blank(), plot.title = element_text(colour="black", size=12, hjust=0.5))
meanC.graph

#regression analysis
meanC.dt<-data.table(mean.C2)
linreg = function (formula) {
  m=lm(formula)
  list(slope=coefficients(m)[2], adj.r2=summary(m)$adj.r.squared, f=summary(m)$fstatistic[1],p.val=summary(m)$coefficients[2,4])
}

regression.results = meanC.dt[,linreg(Sample_Year~mean.C),by=Transect]
regression.results

#FQI





####
#something off about 2013
data.2013<-droplevels(subset(Nachusa.data.long.pa2, Nachusa.data.long.pa2$Sample_Year==2013))
levels(data.2013$newcol)
#235 species
range(data.2013$cover_value)
data.2013.P08<-droplevels(subset(data.2013, data.2013$Transect=="P08"))
levels(data.2013.P08$newcol)
range(data.2013.P08$cover_value)
P08.native<-ddply(data.2013.P08, .(Transect, Sample_Year,Quadrat, native.2016), summarize, n=sum(cover_value))
P08.native2<-ddply(data.2013.P08, .(Transect, Sample_Year,Quadrat, native.2016), summarize, n=length(unique(levels(newcol))))

#0-30, so on the right track, but not 
n.native<-ddply(native, .(Transect, Sample_Year,native.2016), summarize, n=length(n))
range(n.native$n)
sum.total<-ddply(n.native, .(Transect, Sample_Year), summarize, total=sum(n))
sum.native<-droplevels(subset(n.native, n.native$native.2016=="native"))
p.native<-merge(sum.total, sum.native, by=c("Transect","Sample_Year"), all=TRUE)
p.native$p_native<-(p.native$total.y/p.native$total.x)

range(p.native$total.x)
range(p.native$total.y)

#graph number of native species


total.natives<-ggplot(p.native, aes(Sample_Year,total.x))+geom_point()+facet_wrap(~Transect)
total.natives

#graph proportion of native species



#FQI