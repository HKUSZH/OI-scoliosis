library(tidyverse)
library(scales)
library(ggplot2)
library(gplots)
library(patchwork)

load("COBB.RData")

indRates1<-intersect(which(COBB$vecRateProg> -5 & COBB$vecRateProg<30), indGetNonNA)
indRates2<-intersect(which(COBB$vecRateProg>0 & COBB$vecRateProg<30), indGetNonNA)

COBB666<-COBB[!is.na(COBB$vecMajor), ]
COBB606<-COBB[intersect(which(!is.na(COBB$vecMajor)), indGetNonNA), ]

COBB553<-COBB[indRates1, ]
COBB361<-COBB[indRates2, ]

getNewProgRate<-function(){
pred2$fit
}

calDev<-function(AGESEG, forecastk, COBBfuture){
	deviateCobb<-rep(NA, nrow(COBBfuture))
	for(i in 1:nrow(COBBfuture)){
		agei<-COBBfuture$vecAge[i]
		indFrom<-max(which(AGESEG<agei))
		indTo<-min(which(AGESEG>agei))
		deltaAge<-AGESEG[indTo]-AGESEG[indFrom]
		deltaCobb<-forecastk[indTo]-forecastk[indFrom]
		SLOPE<-(deltaCobb/deltaAge)

		estCobbi<-forecastk[indFrom]+(agei-AGESEG[indFrom])* SLOPE
		deviateCobb[i]<- estCobbi-COBBfuture$vecMajor[i]
	}
	names(deviateCobb)<-COBBfuture$vecAge
	return(deviateCobb)
}

Laccuracy<-list()
for(m in 1:11){
Laccuracy[[m]]<-list()
INDEX<-0.9+m/10

for(n in 1:10){
ratio<-n/10

pdf("Kalman-like.leave.one.out.pdf")
	LDEV<-list()
	for(i in c(10, 14, 21, 151, 155)){
		iPat<-unique(COBB553$vecName)[seq(290)[-i]]
		indTrain<-which(COBB553$vecName %in% iPat)
		indTest<-which(!COBB553$vecName %in% iPat)
		COBB.Train<-COBB553[indTrain, ]
		COBB.Test<-COBB553[indTest, ]

		ggplot(COBB.Train, aes(vecAgeMidPoint, vecRateProg)) +
			geom_point() + 
			geom_smooth()


		fit1<-loess(vecRateProg~vecAgeMidPoint, data=COBB.Train)
		pred1<-predict(fit1, newdata=data.frame(vecAgeMidPoint=seq(1,25)), se = TRUE)

		DRAW<-function(PRED, PLOT=1, MAIN=""){
			Cobb.fit<-cumsum(PRED$fit)
			se.cum<-cumsum(PRED$se.fit^2)
			if(PLOT==1)
				plot(Cobb.fit, type='n', pch=16, ylim=c(0, 150), main=MAIN, xlab="Age (years)")
			if(PLOT!=1)
				lines(seq_along(Cobb.fit),Cobb.fit, type='b')
			lines(seq_along(Cobb.fit),Cobb.fit-se.cum*1.96, lty=2, col=2)
			lines(seq_along(Cobb.fit),Cobb.fit+se.cum*1.96, lty=2, col=2)
		}

		lenu<-length(unique(COBB.Test$vecName))
		for(j in 1:lenu){
			COBB.j<-COBB.Test[COBB.Test$vecName %in% unique(COBB.Test$vecName)[j], ]
			vecAgeMPj<-COBB.j$vecAgeMidPoint
			vecAgeMPj<-c(vecAgeMPj, max(vecAgeMPj)+5)
			pred2<-predict(fit1, newdata=data.frame(vecAgeMidPoint=vecAgeMPj), se = TRUE)
			progj<-pred1$fit
			progsej<-pred1$se.fit
			MAIN1<-COBB.j$vecName[1]

			DRAW(pred1, MAIN=i)
			#DRAW(pred2, 2)
			LDEV[[i]]<-list()
			for(k in 1:nrow(COBB.j)){		
				if(k==1){
					ratek<-pred2$fit
				}else{
					ratek<-pred2$fit #getNewProgRate(forecastk, seforecastk)
				}
				OFFSETk<-(COBB.j$vecRateProg[k]-ratek[k])	
				OFFSETk<- sign(OFFSETk)*abs(OFFSETk)^(1.2)
				cat("OFFSETk", OFFSETk, "\n")	
				#progj<-progj+OFFSETk*0.1
				points(COBB.j$vecAge[1:k], COBB.j$vecMajor[1:k], col=k+1, pch=16)
				points(COBB.j$vecAge[-seq(k)], COBB.j$vecMajor[-seq(k)], col=k+1, pch=16,cex=2)

				lines(c(0, COBB.j$vecAge[1:k]), c(0, COBB.j$vecMajor[1:k]), col="red")

				indK<-seq(floor(COBB.j$vecAge[k]))
				forecastk<- cumsum(c(0,progj[-indK]))+ COBB.j$vecMajor[k]
				#forecastk<-forecastk[-2]

				seforecastk<-c(0, cumsum(progsej^2)[-indK])
				#seforecastk<-seforecastk[-2]
				ageOffset<-1#abs(COBB.j$vecAge[k]-round(COBB.j$vecAge[k]))

				agek<-c(COBB.j$vecAge[k], ageOffset+as.integer(names(forecastk)[-1]))
				#points(agek, forecastk, lty=2, col=k+1)
				lines(agek, forecastk, lty=2, lwd=2, col=k+1)
				polygon(c(agek, rev(agek)), 
					y = c((forecastk-seforecastk*1.96), rev(forecastk+seforecastk*1.96)), 
					border = NA, col="#80808040")

				#lines(agek, (forecastk-seforecastk*1.96), lty=2, col=k+1)
				#lines(agek, (forecastk+seforecastk*1.96), lty=2, col=k+1)
	############################
				COBB.jfuture<-COBB.j[-seq(k),]

				LDEV[[i]][[k]]<-calDev( agek, forecastk, COBB.jfuture)
			}
		}
	}

dev.off()

deviate1<-abs(unlist(sapply(LDEV, function(x)x[[1]])))

m1<-mean(deviate1, na.rm=T)
sd1<-sd(deviate1, na.rm=T)

deviate2<-abs(sapply(LDEV, function(x){y<-x[[1]]; rev(y)[1]}))
summary(deviate2)
m2<-mean(deviate2, na.rm=T)
sd2<-sd(deviate2, na.rm=T)

deviate3<-abs(unlist(sapply(LDEV, function(x){x[length(x)-1][1]})))
summary(deviate3)
m3<-mean(deviate3, na.rm=T)
sd3<-sd(deviate3, na.rm=T)

accu<-round(cbind(c(m1, sd1),
	c(m2, sd2),
	c(m3, sd3)), 2)

cat("m=", m, "\tn=", n, "\tINDEX=", INDEX, "\tratio=",ratio, "\t",
	paste(apply(accu, 2, paste0, collapse="+-"), collapse="\t"),
	"\t", date(), "\n")


Laccuracy[[m]][[n]]<-accu


}}




