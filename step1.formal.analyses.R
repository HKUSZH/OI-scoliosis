library(tidyverse)
library(scales)
library(ggplot2)
library(gplots)


genetics<-readxl::read_excel("OI-seq-summary-Nov-29-2022.xlsx", sheet = "SZH-OI")


patientInfo<-readxl::read_excel("patientInfo.xlsx", sheet = "basicInfo")
patientInfo2<-patientInfo[match(SCOLIOSIS$hospital_id, as.matrix(patientInfo)[,20]), ]
patientInfo2[["DOB"]]<-substr(as.matrix(patientInfo2)[,6], 1, 10)

SCOLIOSIS0<-readxl::read_excel("Dec19-OIwithscoliosisandgenetype_aug30V2.0.LLD.xlsx", sheet = "pkchen")
SCOLIOSIS<-SCOLIOSIS0[-which(SCOLIOSIS0$Genotype%in%c("PLOD2", "P4HB")), ]

tabGeno<-table(SCOLIOSIS$Genotype, SCOLIOSIS$Daniel, useNA="always")[,-5]
tabGeno<-tabGeno[order(rowSums(tabGeno), decreasing=T), ]
tabGeno1<-tabGeno[setdiff(rownames(tabGeno), c(NA, "No mutation")), ]
tabGeno2<-tabGeno[c(6, 3), ]
tabGeno<-rbind(tabGeno1, tabGeno2)

 percent(rowSums(tabGeno)/290, 0.1)
 percent(as.numeric(table(SCOLIOSIS$gender))/290, 0.1)
 percent(as.numeric(table(SCOLIOSIS$Daniel))/290, 0.1)

tabSurg<-table(SCOLIOSIS$surgery, SCOLIOSIS$Daniel)

tabGender<-table(SCOLIOSIS$gender, SCOLIOSIS$Daniel)

tabDrugs<-table(SCOLIOSIS$drugs, SCOLIOSIS$Daniel, useNA="always")[,-5]


SCOLIOSIS[which(grepl("COL", SCOLIOSIS$Genotype) & !is.na(SCOLIOSIS$quantQual)), ]


apply(table(SCOLIOSIS$quantQual[SCOLIOSIS$Genotype=="COL1A1"], SCOLIOSIS$Daniel[SCOLIOSIS$Genotype=="COL1A1"]), 2, function(x)paste0(x, collapse=" v. "))
apply(table(SCOLIOSIS$quantQual[SCOLIOSIS$Genotype=="COL1A2"], SCOLIOSIS$Daniel[SCOLIOSIS$Genotype=="COL1A2"]), 2, function(x)paste0(x, collapse=" v. "))

################################################
SCOLIOSIS.col1<-SCOLIOSIS[grep("COL1", SCOLIOSIS$Genotype), ]

matAge<-round(sapply(split(SCOLIOSIS$currentAge, SCOLIOSIS$Daniel), summary), 1)
strAge<-paste0(round(matAge["Mean",], 1), " (", matAge["1st Qu.",], " to ",  matAge["3rd Qu.",], ")")

 as.matrix(percent(rowSums(tabGeno)/221,0.1))

 as.matrix(percent(rowSums(tabGeno)/290,0.1))
 as.matrix(percent(rowSums(tabDrugs[order(rowSums(tabDrugs), decreasing=T), ])/290,0.1))
 as.matrix(percent(rowSums(tabSurg[order(rowSums(tabSurg), decreasing=T), ])/290,0.1))


MatFig1<-rbind( 
	tabGender, strAge, 
	tabGeno,
	tabSurg[order(rowSums(tabSurg), decreasing=T), ], 
	tabDrugs[order(rowSums(tabDrugs), decreasing=T), ])


xlsx::write.xlsx(MatFig1, "data-out/Table-1.xlsx", 
	sheetName = "Table1", 
	col.names = TRUE, row.names = TRUE, append = F)
###############################################

SCOLIOSISdate<-SCOLIOSIS[, grep("^Date", colnames(SCOLIOSIS))]

lastVisit<-as.Date(apply(SCOLIOSISdate, 1 , function(x)rev(x[!is.na(x)])[1]))
AgelastVisit<-as.numeric(difftime(lastVisit, as.Date(SCOLIOSIS$dob), units ="days"))/365

###############################################

matAge2<-round(sapply(split(AgelastVisit, SCOLIOSIS$Daniel), summary), 1)
strAge2<-paste0(round(matAge2["Mean",], 1), " (", matAge2["1st Qu.",], " to ",  matAge2["3rd Qu.",], ")")


ageDiff<-as.numeric(difftime(as.Date(as.matrix(patientInfo2)[,6]), as.Date(SCOLIOSIS$dob), units ="days"))

SCOLIOSIS[["newDOB"]]<-patientInfo2[,6]
SCOLIOSIS[which(ageDiff!=0), c("nameChi", "dob", "newDOB")]


summary(lm(SCOLIOSIS$currentAge ~ as.character(SCOLIOSIS$Daniel)))
summary(lm(AgelastVisit ~ as.character(SCOLIOSIS$Daniel)))


chisq.test(rbind(c(29, 117), c(16, 28)))
###############################################
###############################################
allVisit<-apply(SCOLIOSISdate, 1 , function(x)x[!is.na(x)])
tabXray<-table(sapply(allVisit, length), SCOLIOSIS$Daniel)

indDates<-grep("Date", colnames(SCOLIOSIS))
LCobbAll<-lapply(seq(nrow(SCOLIOSIS)), function(i){
	x<-as.matrix(SCOLIOSIS)[i, ]
	Tho<-x[indDates+1]
	TL<-x[indDates+2]
	Lumb<-x[indDates+3]
	tmp<-cbind(as.numeric(gsub("\\(.*|supine", "", Tho)),
		as.numeric(gsub("\\(.*|supine", "", TL)),
		as.numeric(gsub("\\(.*|supine", "", Lumb)))
	return(tmp)
})


maxCobb<-t(sapply(LCobbAll, function(mat){
	apply(mat, 1, max, na.rm=T)
}))
maxCobb[maxCobb==-Inf]<- NA
maxCobb2<-apply(maxCobb[,-1], 1, max, na.rm=T)

MATmaxCobb<-maxCobb[,-1]

#####################################
CobbSites<-sapply(seq(nrow(SCOLIOSIS)), function(i){
	cobbi<-MATmaxCobb[i,]
	ind.max<-which.max(cobbi)

	x<-as.matrix(SCOLIOSIS)[i, ]
	Tho<-x[indDates[ind.max]+1]
	TL<-x[indDates[ind.max]+2]
	Lumb<-x[indDates[ind.max]+3]

	cobbs<-cbind(as.numeric(gsub("\\(.*|supine", "", Tho)),
		as.numeric(gsub("\\(.*|supine", "", TL)),
		as.numeric(gsub("\\(.*|supine", "", Lumb)))

	sitei<-cbind(gsub("^.*\\(|supine|\\)", "", Tho),
		gsub("^.*\\(|supine|\\)", "", TL),
		gsub("^.*\\(|supine|\\)", "", Lumb))
	return(sitei[which.max(cobbs)])
})
EndSite<-sapply(strsplit(CobbSites, "-"), function(x)x[1])
StartSite<-sapply(strsplit(CobbSites, "-"), function(x)x[2])

IVDs<-c(paste0("T",seq(12)), paste0("L",seq(5)))
IVDs<-factor(IVDs, levels=IVDs)

numStart<-match(StartSite, IVDs)
numEnd<-match(EndSite, IVDs)
table(is.na(numStart), is.na(numStart))

apex<-IVDs[(numStart+numEnd)/2]

tabApex<-table(apex[SCOLIOSIS$severity!="non-scoliotic"], SCOLIOSIS$severity[SCOLIOSIS$severity!="non-scoliotic"])
matApex<-t(sapply(seq(8), function(i){
	colSums(tabApex[seq(2*i-1, 2*i), ])
}))
#####################################

DateMaxCobb<-as.vector(sapply(seq(nrow(MATmaxCobb)), function(i){
	indi<-which.max(MATmaxCobb[i,])
	as.matrix(SCOLIOSISdate)[i,indi]
}))
AgeMaxCobb<-as.numeric(difftime(DateMaxCobb, as.Date(SCOLIOSIS$dob), units ="days"))/365

sum(!is.na(maxCobb[,-1]))

SCOLIOSIS[maxCobb2==-Inf, ]

table(rowSums(!is.na(maxCobb[,-1])), sapply(allVisit, length))

rbind(table(maxCobb2<10, SCOLIOSIS$Daniel)[2, ],
	table(maxCobb2>=10 & maxCobb2 <25, SCOLIOSIS$Daniel)[2, ],
	table(maxCobb2>=25 & maxCobb2 <50, SCOLIOSIS$Daniel)[2, ],
	table(maxCobb2>50, SCOLIOSIS$Daniel)[2, ])

#####################################################
#####################################################

SCOLIOSIS[["severity"]]<-"--"
SCOLIOSIS[["severity"]][maxCobb2<10]<-"non-scoliotic"
SCOLIOSIS[["severity"]][maxCobb2>=10 & maxCobb2 <25]<-"mild"
SCOLIOSIS[["severity"]][maxCobb2>=25 & maxCobb2 <50]<-"moderate"
SCOLIOSIS[["severity"]][maxCobb2>50]<-"severe"

#####################################################

tableGenoSeve<-table(SCOLIOSIS$Genotype, SCOLIOSIS[["severity"]], useNA="always")
tableGenoSeve1<-tableGenoSeve[setdiff(rownames(tableGenoSeve), c(NA, "No mutation")), ]
tableGenoSeve1<-tableGenoSeve1[order(rowSums(tableGenoSeve1), decreasing=T), ]

tableGenoSeve2<-tableGenoSeve[c(10, 16), ]
tableGenoSeve<-rbind(tableGenoSeve1, tableGenoSeve2)

chisq.test(tableGenoSeve[c("WNT1", "IFITM5"), 1:4])

chisq.test(rbind(tableGenoSeve[c("WNT1", "IFITM5"), ],
	colSums(tableGenoSeve[c("SERPINF1", "FKBP10", "P3H1", "BMP1", "SERPINH1"), ]))[, 1:4])

chisq.test(rbind(tableGenoSeve[c("IFITM5"), ],
	colSums(tableGenoSeve[c("SERPINF1", "FKBP10", "P3H1", "BMP1", "SERPINH1"), ]))[, 1:4])

chisq.test(rbind(colSums(tableGenoSeve[1:15, ]), tableGenoSeve[16, ])[,1:4])





tabQualquntCol1<-table(SCOLIOSIS$quantQual[SCOLIOSIS$Genotype=="COL1A1"], SCOLIOSIS$severity[SCOLIOSIS$Genotype=="COL1A1"], useNA="always")
tabQualquntCol2<-table(SCOLIOSIS$quantQual[SCOLIOSIS$Genotype=="COL1A2"], SCOLIOSIS$severity[SCOLIOSIS$Genotype=="COL1A2"], useNA="always")

chisq.test(tabQualquntCol1[1:2, c(3,1,2,4)])
chisq.test(tabQualquntCol2[1:2, c(3,1,2,4)])

matMaxcobb<-sapply(split(maxCobb2, SCOLIOSIS$severity)[c(3,1,2,4)], summary)
strMaxcobb<-paste0(round(matMaxcobb["Mean",], 1), " (", matMaxcobb["1st Qu.",], "~",  matMaxcobb["3rd Qu.",], ")")

#####################################################
triradiate<-gsub("-.*", "", SCOLIOSIS$triradiate_cartilage_closure)
triradiateDate<-gsub("[a-z]-*", "", SCOLIOSIS$triradiate_cartilage_closure)
triradiateAge<-difftime(as.Date(triradiateDate), as.Date(SCOLIOSIS$dob))/365

table(triradiate, is.na(triradiateAge))


tabTriSev<-table(round(triradiateAge/2), SCOLIOSIS$severity)[, c(3,1,2,4)]
tabTriSev2<-table(triradiate, SCOLIOSIS$severity)[, c(3,1,2,4)]

tabTriSev2<-table(triradiate, SCOLIOSIS$severity)[, c(3,1,2,4)]

tabTriLLD<-table(SCOLIOSIS$LLD, SCOLIOSIS$severity)
percent(as.numeric(tabTriLLD[2,]/table(SCOLIOSIS$severity)))[ c(3,1,2,4)]

#####################################################
tabDrugsSev<-table(SCOLIOSIS$drugs, SCOLIOSIS$severity)[,c(3,1,2,4)]
	tabDrugsSev[order(rowSums(tabDrugsSev), decreasing=T), ]

tabGenderSev<-table(SCOLIOSIS$gender, SCOLIOSIS$severity)[,c(3,1,2,4)]

tabSurgSev<-table(SCOLIOSIS$surgery, SCOLIOSIS$severity)[c(2,3,1),c(3,1,2,4)]
tabDrugsSev<-table(SCOLIOSIS$drugs, SCOLIOSIS$severity, useNA="always")[,c(3,1,2,4)]
tabDrugsSev[order(rowSums(tabDrugsSev), decreasing=T), ]

matAge3<-round(sapply(split(SCOLIOSIS$firstBPAge[SCOLIOSIS$drugs!="never"], 
	SCOLIOSIS$severity[SCOLIOSIS$drugs!="never"]), summary), 1)[, c(3,1,2,4)]
strAge3<-paste0(round(matAge3["Mean",], 1), " (", matAge3["1st Qu.",], "~",  matAge3["3rd Qu.",], ")")

summary(SCOLIOSIS$firstBPAge[SCOLIOSIS$drugs!="never"])

summary(lm(SCOLIOSIS$firstBPAge [SCOLIOSIS$firstBPAge!=-1] ~ SCOLIOSIS$severity[SCOLIOSIS$firstBPAge!=-1]))
#####################################################
matAge4<-round(sapply(split(AgeMaxCobb, 
	SCOLIOSIS$severity), summary), 1)[, c(3,1,2,4)]
strAge4<-paste0(round(matAge4["Mean",], 1), " (", matAge4["1st Qu.",], "~",  matAge4["3rd Qu.",], ")")

summary(AgeMaxCobb)
summary(lm(AgeMaxCobb~SCOLIOSIS$severity))

#####################################################
SCOLIOSIS[["Zscores"]]<-as.numeric(SCOLIOSIS$latestBMDZscores)
matBMD<-sapply(split(SCOLIOSIS$Zscores, SCOLIOSIS$severity), function(x)summary(x[!is.na(x)]))[,c(3,1,2,4)]
strBMDz<-paste0(round(matBMD["Mean",], 1), " (", matBMD["1st Qu.",], "~",  matBMD["3rd Qu.",], ")")

summary(SCOLIOSIS[["Zscores"]], na.rm=T)

summary(lm(SCOLIOSIS$Zscores ~ SCOLIOSIS$severity))
TukeyHSD(aov(SCOLIOSIS$Zscores ~ SCOLIOSIS$severity))


#####################################################
indDates<-grep("^Date", colnames(SCOLIOSIS))
indRisser<-grep("Risser", colnames(SCOLIOSIS))
maxRisser<-apply(as.matrix(SCOLIOSIS[,indRisser]), 1, function(x){
	max(x[which(!is.na(x))])
})

chisq.test(table(maxRisser, triradiate))

sum(!is.na(SCOLIOSIS[,indDates]) & is.na(SCOLIOSIS[,indRisser]))

indinconst1<-which(rowSums(!is.na(SCOLIOSIS[,indDates]) & is.na(SCOLIOSIS[,indRisser]))>0)
indinconst2<-which(rowSums(is.na(SCOLIOSIS[,indDates]) & !is.na(SCOLIOSIS[,indRisser]))>0)
t(SCOLIOSIS[c(indinconst1, indinconst2), ])

 t(sapply(split(AgelastVisit[SCOLIOSIS$gender=="Male"],  maxRisser[SCOLIOSIS$gender=="Male"]), range))
 t(sapply(split(AgelastVisit[SCOLIOSIS$gender=="Female"],  maxRisser[SCOLIOSIS$gender=="Female"]), range))


tabRisSev<-table(maxRisser, SCOLIOSIS$severity)[,c(3,1,2,4)]
rbind(tabRisSev[1,],
	rowSums(tabRisSev[2:3,]),
	rowSums(tabRisSev[3:4,]),
	tabRisSev[6,])

summary(lm(maxCobb2 ~ AgelastVisit + SCOLIOSIS$gender))
summary(lm(maxCobb2 ~ AgelastVisit + as.factor(SCOLIOSIS$Daniel)+factor(SCOLIOSIS$LLD)))

#####################################################
chisq.test(table(SCOLIOSIS$LLD, SCOLIOSIS$severity))
#####################################################
firstBPAge<-SCOLIOSIS$firstBPAge
firstBPAge[firstBPAge==-1]<-NA

gender<-factor(SCOLIOSIS$gender, levels=c("Male", "Female"))

firstBPAgeGroup<-rep("", length(firstBPAge))
firstBPAgeGroup[is.na(firstBPAge)]<-"anever"
firstBPAgeGroup[!is.na(firstBPAge)]<-paste0("Grp_",formatC(ceiling(firstBPAge[!is.na(firstBPAge)]/6), flag=0, width=2))
firstBPAgeGroup[firstBPAge>6]<-"Grp_02"

geno<-SCOLIOSIS$Genotype
geno[grepl("COL", SCOLIOSIS$Genotype)]<- "COL"
geno[!grepl("COL", SCOLIOSIS$Genotype)]<- "nonCOL"
geno[is.na(SCOLIOSIS$Genotype)]<- "notTested"
geno[SCOLIOSIS$Genotype == "No mutation"]<- "noMut"

geno2<-geno
#geno2[grepl("COL1A1", SCOLIOSIS$Genotype)]<- "aaCOL1A1"
#geno2[grepl("COL1A2", SCOLIOSIS$Genotype)]<- "COL1A2"
geno2[geno=="nonCOL" & grepl("IFIT", SCOLIOSIS$Genotype)]<- "IFITM5"
geno2[geno=="nonCOL" & grepl("WNT", SCOLIOSIS$Genotype)]<- "WNT1"


geno3<-geno2
geno3[grepl("FKBP10", SCOLIOSIS$Genotype)]<- "FKBP10"
geno3[grepl("SERPINF1", SCOLIOSIS$Genotype)]<- "SERPINF1"
geno3[grepl(",", SCOLIOSIS$Genotype)]<- "cmp"


inheritance<-rep("AR", nrow(SCOLIOSIS))
inheritance[grepl("IFIT", SCOLIOSIS$Genotype)]<- "AD"

drugs <-SCOLIOSIS$drugs 
drugs[! drugs %in%c("Pamidronate", "Zoledronate")]<- NA 
drugs[firstBPAgeGroup=="never"]<-"never"

summary(lm(maxCobb2 ~ SCOLIOSIS$Sillence))

gMod1<-glm(maxCobb2 ~ AgeMaxCobb + firstBPAge + as.factor(SCOLIOSIS$LLD) + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser, family="quasipoisson")
gMod2<-glm(maxCobb2 ~ AgeMaxCobb + firstBPAge  + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser, family="quasipoisson")

Mod1<-lm(log(maxCobb2+1) ~ AgeMaxCobb + firstBPAgeGroup + as.factor(SCOLIOSIS$LLD) + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser)
Mod2<-lm(log(maxCobb2+1) ~ AgeMaxCobb + firstBPAgeGroup  + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser)
anova(Mod1, Mod2)

Mod3<-lm(log(maxCobb2+1) ~ AgeMaxCobb + firstBPAgeGroup  + SCOLIOSIS$gender + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser)
Mod4<-lm((maxCobb2+1) ~ AgeMaxCobb + SCOLIOSIS$Genotype + drugs +SCOLIOSIS$Zscores + maxRisser)
Mod5<-lm((maxCobb2+1) ~ AgeMaxCobb + SCOLIOSIS$Genotype + firstBPAgeGroup  +SCOLIOSIS$Zscores + maxRisser)

Mod5<-lm(maxCobb2 ~ AgeMaxCobb + SCOLIOSIS$gender)
Mod5<-lm(maxCobb2 ~ AgeMaxCobb + SCOLIOSIS$gender + drugs)
Mod5<-lm(maxCobb2 ~ AgeMaxCobb  + drugs + )
Mod5<-lm(maxCobb2 ~ AgeMaxCobb  + drugs + SCOLIOSIS$Zscores)
Mod5<-lm(maxCobb2 ~ AgeMaxCobb  + drugs + SCOLIOSIS$Zscores + maxRisser + inheritance)
Mod5<-lm(maxCobb2 ~ AgeMaxCobb  + drugs + SCOLIOSIS$Zscores + maxRisser + Genotype)

anova(Mod3, Mod4)


summary(lm(maxCobb2 ~ AgeMaxCobb + firstBPAge + as.factor(SCOLIOSIS$LLD) + SCOLIOSIS$Genotype + drugs + SCOLIOSIS$Zscores + maxRisser))

###########################
hasScoliosis<- SCOLIOSIS$severity != "non-scoliotic"
isSevere<- SCOLIOSIS$severity == "severe"
isModerSevere<- SCOLIOSIS$severity %in% c("moderate", "severe")


MAKEREPORTS<-function(fit1){
	s1<-summary(fit1)
	mat1<-cbind(confint(fit1), s1$coefficients)
	#print(mat1)
	CI<-apply(round(mat1[,1:2], 2), 1, paste0, collapse= " ~ ")
	str1<-paste0(round(mat1[, "Estimate"], 2), " (", CI, ")", sep="")

	mat2<-round(exp(cbind(confint(fit1), OR = coef(fit1))), 2)
	#print(mat2)
	CI2<-apply(round(mat2[,1:2], 2), 1, paste0, collapse= " ~ ")
	str2<-paste0(round(mat2[, "OR"], 2), " (", CI2, ")", sep="")

	pval<-round(s1$coefficients[, "Pr(>|z|)"], 3)
	cbind("log(OR) (CI)"=str1, 
		"OR (CI)"=str2,
		pval)
}

logistic_1 <- glm(hasScoliosis ~ AgeMaxCobb  + gender + geno2 +
				drugs +  SCOLIOSIS$Zscores + factor(SCOLIOSIS$LLD), 
                      family = "binomial")
rep1<-MAKEREPORTS(logistic_1)


logistic_2 <- glm(isModerSevere~ AgeMaxCobb  + gender +geno2 +
				drugs +  SCOLIOSIS$Zscores + factor(SCOLIOSIS$LLD), 
                      family = "binomial")
rep2<-MAKEREPORTS(logistic_2)


logistic_3 <- glm(isSevere~ AgeMaxCobb  + gender + geno2 +
				drugs +  SCOLIOSIS$Zscores + factor(SCOLIOSIS$LLD), 
                      family = "binomial")
rep3<-MAKEREPORTS(logistic_3)

xlsx::write.xlsx(cbind(rep1, rep2, rep3), "data-out/Table-5.logistic.xlsx", 
	sheetName = "Table5", 
	col.names = TRUE, row.names = TRUE, append = F)





