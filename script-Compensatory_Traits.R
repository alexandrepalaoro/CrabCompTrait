#################################################################################
## These are the analyses used to test whether aeglids have compensatory       ##
## traits, or if the water acts a compensatory trait. First, we test           ##
## if male and female crabs have different walking performances under two      ##
## conditions: with water flow, or no water flow. Afterwards, we test if       ##
## heavier males walk faster or slower than lighter males. For dimorphism,     ##
## first we test if males have longer pereopods than females. Second, we test  ##
## if males with heavier claws have proportionally larger pereopods.           ##
## And that's about what you will see here :D                                  ##
#################################################################################

rm(list=ls()) # Cleaning up the workspace

## Packages needed
library(lme4)
library(phia)
library(scales)
library(lmerTest)

comptrait<-read.csv("compensatory_data.csv",h=T)
head(comptrait)

comptrait2 <- read.csv("Razao_Dados.csv", header = T)

comptrait2<- comptrait2[which(comptrait2$side=="left"),]

###
# First test: do males walker faster or slower than females?
##

plot(vm~scale(cc),data=comptrait2,las=1,bty='l',cex=1.3,pch=21,
     bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     xlab="Cephalothorax length (scaled and centered)",
     ylab="Maximum speed (cm/s)")

# There is one female that is walking waaay faster than everyone.
# Let's check the analysis. 

m1<-lmer(vm~scale(cc)*sex+fluxo + (1 | id),data=comptrait2)
plot(m1)

summary(m1)

testInteractions(m1,pairwise="sex",slope="scale(cc)")
testInteractions(m1,pairwise="sex")

# Residuals seem fine and females walk, on average, faster than females
# Now, let's plot the regression lines calculated by the lmer
(coef.m1<-m1@beta)

plot(vm~scale(cc),data=comptrait2,las=1,bty='l',cex=1.3,pch=21,
     bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     xlab="Cephalothorax length (scaled and centered)",
     ylab="Maximum speed (cm/s)")

#Female's regression line without flow AND with flow
curve((coef.m1[1]+(coef.m1[2]*x)),add=T,lwd=2,col='grey',
      from=-1.5,to=0.5)
curve((coef.m1[1]+coef.m1[4])+(coef.m1[2]*x),add=T,lwd=2,col='grey',
      lty=2,from=-1.5,to=0.5) #flow

#Male's regression line
curve(((coef.m1[1]+coef.m1[3])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=2,col='black',from=-0.5,to=1.6)
curve(((coef.m1[1]+coef.m1[3]+coef.m1[4])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=2,col='black',from=-0.5,to=1.6,lty=2)

## Now, we will perform what you could call a sensibility analysis.
## We will perform the exact same test above, but this time without 
## the fastest female and te slowest male. 

m1.s<-lmer(vm~scale(cc)*sex+fluxo + (1 | id),
         data=comptrait2[comptrait2$vm!=max(comptrait2$vm)&comptrait2$vm!=min(comptrait2$vm),])
plot(m1.s)
summary(m1.s)

testInteractions(m1.s,pairwise='sex',slope='scale(cc)')
testInteractions(m1.s,pairwise='sex')

## No difference in the results. So, we can keep the first analysis.

tiff(file="Figure2.tiff",units="mm",width=180,height=150,res=600,
     compression="lzw")

par(mar=c(5,5,2,2)+0.1)

plot(vm~scale(cc),data=comptrait2,las=1,bty='l',cex=1.3,pch=21,
     bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     xlab="Cephalothorax length (scaled and centered)",
     ylab="Maximum speed (cm/s)",ann=F)
mtext(side=2,text="Maximum speed (cm/s)",line=4)
mtext(side=1,text="Cephalothorax length (scaled and centered)",line=3)
#Female's regression line
curve((coef.m1[1]+(coef.m1[2]*x)),add=T,lwd=2,col='grey',
      from=-1.5,to=0.5) #no-flow
curve((coef.m1[1]+coef.m1[4])+(coef.m1[2]*x),add=T,lwd=2,col='grey',
      lty=2,from=-1.5,to=0.5) #flow

#Male's regression line
curve(((coef.m1[1]+coef.m1[3])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=2,col='black',from=-0.5,to=1.6) #no-flow
curve(((coef.m1[1]+coef.m1[3]+coef.m1[4])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=2,col='black',from=-0.5,to=1.6,lty=2) #flow
legend("topright",c("Female","Male"),cex=1.3,bty='n',
       pch=21,pt.bg= c(alpha("gray",0.75), alpha("black",0.75)))

dev.off()

####===== Finding a good proxy for claw weight ====####

## We would have to test whether individuals with heavier claws had
## longer pereiopods. But we cannot compare weight, which has a cubic
## scaling pattern, to length, which scales linearly - we could find
## a relation simply because the scaling patterns differ. 
## To avoid that, we need to find a good proxy for weight that also 
## scales linearly. We have two candidates: claw length and weight.
## To select between them, we will use the R-square of a linear 
## regression: the highest R-square will be used as the proxy in the
## next analysis.

proxy.test<-read.csv("dados_pere.csv", header = T) # Carregar os dados 

## Lefties only
proxy.test<-proxy.test[which(proxy.test$flank=="canhoto" & 
                               proxy.test$side == "left"),] 

## Now, to the tests

m1<-lm(log(peso_propodo)~scale(log(cp)), data = proxy.test)
plot(m1,which=1)

m2<-lm(log(peso_propodo)~scale(log(ap)),data=proxy.test)
plot(m2,which=1)

summary(m1) #Adjusted R-squared: 0.9861
summary(m2) #Adjusted R-squared: 0.9935

## Claw height wins! Now, to the real test

data1<-read.csv("Esquerdo_Dados.csv", header = T)

## Let's look at how the data scales before jumping on the tests

par(mar=c(5,4,4,2))
plot(p~scale(ap), data = data1, pch = 21,ann=T,  
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)], 
     xlab = "Claw height (scaled and centered)", ylab = "Pereopod length (cm)",
     cex = 1.3, bty = "l",las =1)

## Seems fine, although a little bit hump-shaped. Let's perform
## the LMMs and check the residuals

m3<-lmer(p~scale(ap)*sex + (1 | id), data = data1)
plot(m3,col=data1$sex)

## Males are a bit more dispersed than females... Let's check the
## log version

plot(log(p)~scale(log(ap)), data = data1, pch = 21,ann=T,  
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)], 
     xlab = "Claw height (scaled and centered)", ylab = "Pereopod length (cm)",
     cex = 1.3, bty = "l",las =1)
## More linear than the previous analysis... Now, the residuals

m3.log<-lmer(log(p)~scale(log(ap))*sex+(1|id),data=data1)
plot(m3.log,col=data1$sex)

## Much better! Let's check the results!

summary(m3)
summary(m3.log)

## Pretty much identical...

testInteractions(m3,pairwise = 'sex',slope='scale(ap)')
testInteractions(m3,pairwise = 'sex')

testInteractions(m3.log,pairwise = 'sex',slope='scale(log(ap))')
testInteractions(m3.log,pairwise = 'sex')

## Identical, it really does not matter. Thus, we will keep the log-version because 
## it seems better fitted. 

(coef.m3<-m3.log@beta)

plot(log(p)~scale(log(ap)), data = data1, pch = 21,ann=T,  
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)], 
     xlab = "Claw height (scaled and centered)", ylab = "Pereopod length (cm)",
     cex = 1.3, bty = "l",las =1)

#Female
curve(coef.m3[1]+(coef.m3[2]*x),add=T,col='grey',lwd=2,from=-2,to=0.2)

#Male
curve((coef.m3[1]+coef.m3[3])+((coef.m3[2]+coef.m3[4])*x),add=T,col='black',
      lwd=2,from=-1,to=2.1)

## That may be why females walk faster. Females have longer walking legs on average,
## and females also have proportionally longer legs. Cool.

## Now, for our last test: do females have more muscles on their legs than males?

plot(musc1~scale(peso_propodo), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (mg)",
     cex = 1.2, bty = "l",las =1)

## Uff, highly exponential. That is a bit tricky to handle without using complex
## models. Let's try logging it.

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (log)",
     cex = 1.2, bty = "l",las =1)

## Not *THAT* linear, but we could use a Gamma distribution... But first,
## let's check how a linear regression handles it.

m4<-lmer(log(musc1+1)~scale(log(peso_propodo))*sex + (1 | id), 
         data = data1)
plot(m4,col=data1$sex)

## Okay, residuals increase with the mean value. Time to bring out the big
## guns: the Gamma distribution.

m4.gam<-glmer(log(musc1+1)~scale(log(peso_propodo))*sex + (1|id),data=data1,
          family="Gamma"(link="log"))
plot(m4.gam,col=data1$sex)

## Worked as a charm :D Now, let's check the results

summary(m4.gam)

testInteractions(m4.gam,pairwise="sex",slope="scale(log(peso_propodo))")
testInteractions(m4.gam,pairwise="sex")

## Not difference, but the mean effect is borderline significant (p=0.06):
## Females have more muscle on average. Let's check the plot...

(coef.m4<-m4.gam@beta)

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (log)",
     cex = 1.2, bty = "l",las =1)

#Females
curve(exp(coef.m4[1]+(coef.m4[2]*x)),add=T,col='grey',lwd=2,
      from=-2,to=0.2)
#Males
curve(exp((coef.m4[1]+coef.m4[3])+((coef.m4[2]+coef.m4[4])*x)),add=T,col='black',
      lwd=2,from=-1,to=2.1)

## Now, I will plot the figure andn that should be it for the day


tiff(file="Figure3-PANEL.tiff",units="mm",width=220,height=120,res=600,
     compression="lzw")
par(mfrow=c(1,2),mar=c(4,4,2,3))
plot(log(p)~scale(log(ap)), data = data1, pch = 21,ann=T,  
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)], 
     xlab = "Claw height (scaled and centered)", ylab = "Pereopod length (cm)",
     cex = 1.3, bty = "l",las =1)
mtext("(a)",4,las=1,padj=-12,adj=-1)
legend("topleft",c("Female","Male"),cex=1.3,bty='n',
       pch=21,pt.bg= c(alpha("gray",0.75), alpha("black",0.75)))

#Female
curve(coef.m3[1]+(coef.m3[2]*x),add=T,col='grey',lwd=2,from=-2,to=0.2)

#Male
curve((coef.m3[1]+coef.m3[3])+((coef.m3[2]+coef.m3[4])*x),add=T,col='black',
      lwd=2,from=-1,to=2.1)

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (log)",
     cex = 1.2, bty = "l",las =1)
mtext("(b)",4,las=1,padj=-12,adj=-1)

#Females
curve(exp(coef.m4[1]+(coef.m4[2]*x)),add=T,col='grey',lwd=2,
      from=-2,to=0.2)
#Males
curve(exp((coef.m4[1]+coef.m4[3])+((coef.m4[2]+coef.m4[4])*x)),add=T,col='black',
      lwd=2,from=-1,to=2.1)
dev.off()
