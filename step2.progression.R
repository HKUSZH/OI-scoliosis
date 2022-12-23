library(tidyverse)
library(scales)
library(ggplot2)
library(gplots)


load("SCOLIOSIS.RData")

SCOLIOSISdate<-SCOLIOSIS[, grep("^Date", colnames(SCOLIOSIS))]


longiCobb<-sapply(LCobbAll, function(x){	apply(x[-1,], 1, max, na.rm=T)})

longiDate<-sapply(1:290, function(i){ as.Date(as.matrix(SCOLIOSISdate)[i,])})
longiAge<-sapply(1:290, function(i){
		DATEi<-as.Date(as.matrix(SCOLIOSISdate)[i,])
		as.numeric(difftime(DATEi, as.Date(SCOLIOSIS$dob[i]), units="day"))/365
	})
longiSurgAge<-sapply(1:290, function(i){
		DATEi<-as.Date(as.matrix(SCOLIOSISdate)[i,])
		as.numeric(difftime(DATEi, SCOLIOSIS$surgeryDate[i], units="day"))/365
	})
longiCobb[longiCobb==-Inf]<-NA

table(colSums(!is.na(longiDate)),  colSums(!is.na(longiAge )))
table(colSums(!is.na(longiDate)),  colSums(!is.na(longiCobb)))

CobbDiff<-apply(longiCobb, 2, function(x)diff(c(0,x)))
AgeDiff<-apply(longiAge, 2, function(x)diff(c(0,x)))

table(colSums(!is.na(CobbDiff)),  colSums(!is.na(AgeDiff)))

RateProg<-CobbDiff/AgeDiff

table(colSums(!is.na(CobbDiff)),  colSums(!is.na(RateProg)))

############################
indSurg<-which(!is.na(SCOLIOSIS$surgeryDate))
cobbsurg<-CobbDiff[, !is.na(SCOLIOSIS$surgeryDate)]
cobbnonsurg<-CobbDiff[, is.na(SCOLIOSIS$surgeryDate)]

inddrop<-which(apply(CobbDiff, 2, min, na.rm=T) < -10 & SCOLIOSIS$surgery=="noOp")
as.matrix(SCOLIOSIS[inddrop, 1:4])
CobbDiff[, inddrop]




matNA<-matrix(seq(prod(dim(longiCobb))), nrow=10)
matNA[is.na(longiCobb)]<-NA
matNA[, indSurg]
table(is.na(matNA))
LindnonNA<-sapply(apply(cobbsurg, 2, which.min),seq)

LindNA<-sapply(apply(cobbsurg, 2, which.min), function(i) seq(i, 10))

sapply(seq(length(indSurg)), function(i){
	matNA[LindNA[[i]], indSurg[i]]<<-NA
})
table(is.na(matNA))
matNA[, indSurg]
indSetNA<-which(is.na(matNA))
indGetNonNA<-which(!is.na(matNA))
############################
apply(longiCobb[, !is.na(SCOLIOSIS$surgeryDate)], 2, max, na.rm=T)
summary(apply(cobbsurg, 2, min, na.rm=T))
############################
indRisser<-grep("Risser", colnames(SCOLIOSIS))
longiRisser<-t(SCOLIOSIS[, indRisser])

table(colSums(!is.na(longiRisser)),  colSums(!is.na(RateProg)))

vecMajor<-as.vector(longiCobb)
vecAge<-as.vector(longiAge)
vecRisser<-as.vector(longiRisser)
vecRateProg<-as.vector(RateProg)

matInd<-rep(seq(nrow(SCOLIOSIS)), rep(10, nrow(SCOLIOSIS)))

matrix(vecLLD, nrow=10)[,1:5]


vecGeno<-rep(SCOLIOSIS$geno, rep(10, nrow(SCOLIOSIS)))
vecGeno2<-rep(SCOLIOSIS$geno2, rep(10, nrow(SCOLIOSIS)))
vecGeno3<-rep(SCOLIOSIS$geno3, rep(10, nrow(SCOLIOSIS)))
vecZscores<-rep(SCOLIOSIS$Zscores, rep(10, nrow(SCOLIOSIS)))
vecLLD<-rep(SCOLIOSIS$LLD, rep(10, nrow(SCOLIOSIS)))
vecQuantQual<-rep(SCOLIOSIS$quantQual, rep(10, nrow(SCOLIOSIS)))
vecGender<-rep(SCOLIOSIS$gender, rep(10, nrow(SCOLIOSIS)))
vecDrugs2<-rep(SCOLIOSIS$drugs2, rep(10, nrow(SCOLIOSIS)))
vecName<-rep(SCOLIOSIS$nameChi, rep(10, nrow(SCOLIOSIS)))
vecGenotype<-rep(SCOLIOSIS$Genotype, rep(10, nrow(SCOLIOSIS)))
vecSurgery<-rep(SCOLIOSIS$surgery, rep(10, nrow(SCOLIOSIS)))
vecMaxCobb2<-rep(maxCobb2, rep(10, nrow(SCOLIOSIS)))

############################
ProgAgeGroup<-ceiling(vecAge/5)
ProgAgeGroup[ProgAgeGroup>7]<-8

indRates<-intersect(which(vecRateProg>=0 & vecRateProg<20), indGetNonNA)

boxplot(vecRateProg[indRates]~ProgAgeGroup[indRates], outline=F, xlab="age (yrs)", ylab="degree/year")

vecGeno4<-vecGeno3
vecGeno4[vecGenotype=="COL1A1"]<-"aaCOL1A1"
vecGeno4[vecGenotype=="COL1A2"]<-"COL1A2"

COBB<-data.frame(vecName, vecRateProg, vecMajor, vecAge, vecGender, vecSurgery, vecMaxCobb2, vecRisser,
	ProgAgeGroup, vecGenotype, vecGeno, vecGeno2, vecGeno3, vecGeno4,
	vecZscores, vecLLD, vecQuantQual, vecDrugs2)

COBB1<-droplevels(COBB[indGetNonNA, ])
COBB2<-droplevels(COBB[indRates, ])
COBB3<-droplevels(COBB[vecGeno2=="WNT1", ])
COBB4<-droplevels(COBB2[COBB2$vecMaxCobb2>=10, ])
COBB5<-droplevels(COBB[intersect(indGetNonNA, which(vecGenotype%in%c("COL1A1", "COL1A2"))), ])

boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB1)
boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB2)
boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB2[COBB2$vecGeno4 =="aaCOL1A1", ], ylim=c(0, 12))
boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB2[COBB2$vecGeno4 =="COL1A2", ], ylim=c(0, 12))

boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB3)

boxplot(vecRateProg~ProgAgeGroup, outline=F, data=COBB2[COBB2$vecGeno=="COL", ])

ggplot(COBB, 
			aes(x=vecAge,y=vecRateProg,group=vecName,label=vecName)) +
		geom_line(aes(color=vecGeno), size=1, alpha=0.5) + 
		geom_point(aes(color=vecGeno, shape=vecGender),size=3)  + ylim(c(0, 10))
pdf("Cobb.Curve.pdf")
	ggplot(COBB5[COBB5$vecGenotype=="COL1A1", ], 
			aes(x=vecAge,y=vecMajor)) +
		geom_point(size=3) +
		geom_line(aes(group=vecName), size=1, alpha=0.5) + 
		geom_smooth() + xlim(c(0, 25)) +  ylim(c(0, 125)) + ggtitle("COL1A1")

	ggplot(COBB5[COBB5$vecGenotype=="COL1A2", ], 
			aes(x=vecAge,y=vecMajor)) +
		geom_point(size=3) + 
		#geom_text(aes(label=vecName)) +
		geom_line(aes(group=vecName), size=1, alpha=0.5) + 
		geom_smooth() + xlim(c(0, 25)) +  ylim(c(0, 125)) + ggtitle("COL1A2")
dev.off()


data=COBB5[COBB5$vecGenotype=="COL1A1", ], aes(x=vecAge,y=vecMajor))

###########################
source("helpers.R")

fit_1 <- lm(vecRateProg ~ ProgAgeGroup + vecGender + vecGeno2+
				vecDrugs2 + vecZscores + vecLLD, data=COBB2)
summary(fit_1)

fit_2 <- lm(vecRateProg ~ vecGender , data=COBB2)
summary(fit_2)

fit_2 <- lm(vecRateProg ~ as.character(ProgAgeGroup)  , data=COBB2)
summary(fit_2)

fit_2 <- lm(vecRateProg ~ as.character(vecRisser)  , data=COBB2)
summary(fit_2)

fit_2 <- lm(vecRateProg ~ as.character(ProgAgeGroup) + as.character(vecRisser) , data=COBB2)
summary(fit_2)


fit_2 <- lm(vecRateProg ~  as.character(ProgAgeGroup) + vecGeno4, data=COBB2)
summary(fit_2)
confint(fit_2)
MAKEREPORTSanova(fit_2, PVAL="Pr(>|t|)")

fit_2 <- lm(vecRateProg ~ vecGeno4, data=COBB2[COBB2$vecGeno4 %in% c("aaCOL1A1", "COL1A2"), ])
summary(fit_2)
confint(fit_2)+2.91975


fit_2 <- lm(vecRateProg ~ factor(ProgAgeGroup), data=COBB2)
summary(fit_2)
confint(fit_2)

fit_2 <- lm(vecRateProg ~ vecDrugs2, data=COBB2)
summary(fit_2)
confint(fit_2)

