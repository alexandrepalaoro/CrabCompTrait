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

png(file="Figure2.png", units='mm',width=200,height=150,res=600)
tiff(file="Figure2.tiff",units="mm",width=200,height=150,res=600,
     compression="lzw")

par(mar=c(5,5,2,2)+0.1)

plot(vm~scale(cc),data=comptrait2,las=1,bty='l',cex=1.5,pch=21,
     bg=c(alpha("grey",0.75),alpha("black",0.75))[as.numeric(sex)],
     xlab="Cephalothorax length (scaled and centered)",
     ylab="Maximum speed (cm/s)",ann=F)
mtext(side=2,text="Maximum speed (cm/s)",line=4)
mtext(side=1,text="Cephalothorax length (scaled and centered)",line=3)
#Female's regression line
curve((coef.m1[1]+(coef.m1[2]*x)),add=T,lwd=3,col='grey',
      from=-1.5,to=0.5) #no-flow
curve((coef.m1[1]+coef.m1[4])+(coef.m1[2]*x),add=T,lwd=3,col='grey',
      lty=2,from=-1.5,to=0.5) #flow

#Male's regression line
curve(((coef.m1[1]+coef.m1[3])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=3,col='black',from=-0.5,to=1.6) #no-flow
curve(((coef.m1[1]+coef.m1[3]+coef.m1[4])+((coef.m1[2]+coef.m1[5])*x)),add=T,
      lwd=3,col='black',from=-0.5,to=1.6,lty=2) #flow
legend("topright",c("Female","Male"),cex=1.3,bty='n',
       pch=21,pt.bg= c(alpha("gray",0.75), alpha("black",0.75)))

dev.off()

####===== Finding a good proxy for claw weight ====####

## We would have to test whether individuals with heavier claws had
## longer pereiopods. But we cannot compare weight, which has a cubic
## scaling pattern, to length, which scales linearly - we could find
## a relation simply because the scaling patterns differ. 
## To avoid that, we need to find a good proxy for weight that also 
## scales linearly. We have two candidates: claw length and height.
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

## Claw height wins! Now, to the real tests

data1<-read.csv("Esquerdo_Dados.csv", header = T)

## First, we have to test if males and females have asymmetric legs. However,
## since legs are correlated to body size, we regress leg length against body size
## and then use the residuals to test.

resid.leg<-lmer(p~sex*scale(cc)+(1|id),data=data1)
plot(resid.pereop)
summary(resid.pereop)

## Seems fine, so let's add to the main data.frame

data1$resid<-residuals(resid.pereop)

## Now, we divide in two data.frames. One for the pereopods on the left side of the body, 
## and another for the pereopods on the right side

left<-data1[data1$side=="left",]
right<-data1[data1$side=="right",]

## Plotting...

plot(left$resid~right$resid,bty='l',las=1,cex=1.5,
     pch=21,bg=c(alpha("lightgrey",0.75),alpha("grey30",0.75))[as.numeric(left$sex)])

## There are two possible female outliers. So we have to perform a sensibility test 
## once again. But first let's test without messing around.

asym<-lm(left$resid~right$resid*left$sex)
plot(asym,which=1,col=left$sex)
summary(asym)
summary.aov(asym)

## No difference for the sex factor; males are not more asymmetrical (or symmetrical)
## than females. The next step is checking the impact of those outliers in our analysis
## So, we will remove the outliers and rerun the analysis.

left1<-left[left$resid!=max(left$resid)&left$resid!=min(left$resid),]
right1<-right[right$resid!=max(right$resid)&right$resid!=min(right$resid),]

asym1<-lm(left1$resid~right1$resid*left1$sex)
plot(asym1,which=1,col=left1$sex)
summary(asym1)

# Results did not change much... Let's see the plot

plot(left1$resid~right1$resid,bty='l',las=1,cex=1.5,
     pch=21,bg=c(alpha("lightgrey",0.75),alpha("grey30",0.75))[as.numeric(left1$sex)])
asym1.p<-lm(left1$resid~right1$resid*left1$sex-1)

curve(coef(asym1.p)[2]+(coef(asym1.p)[1]*x),add=T,lwd=2,lty=2,col=alpha("lightgrey",0.75))
curve(coef(asym1.p)[3]+((coef(asym1.p)[1]+coef(asym1.p)[4])*x),add=T,lwd=2,lty=2,col=alpha("grey30",0.75))

## I will leave this here to make plotting the curves easier later on
asym.p<-lm(left$resid~right$resid*left$sex-1)
coef(asym.p)[1]+coef(asym.p)[4]
confint(asym.p)[1,1]+confint(asym.p)[4,1]
confint(asym.p)[1,2]+confint(asym.p)[4,2]

## Now, moving to the other tests

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
     cex = 1.5, bty = "l",las =1)

#Female
curve(coef.m3[1]+(coef.m3[2]*x),add=T,col='grey',lwd=2,from=-2,to=0.2)

#Male
curve((coef.m3[1]+coef.m3[3])+((coef.m3[2]+coef.m3[4])*x),add=T,col='black',
      lwd=2,from=-1,to=2.1)

## That may be why females walk faster. Females have longer walking legs on average,
## and females also have proportionally longer legs. Cool.

## Okay, but now we need to know whether bigger males have relatively larger legs than
## smaller males. To do this, we ran a similar model as ran above. The key differences are:
# 1. We used only males;
# 2. Both variables were scaled and centered.
## We did this to compare the male's slope with an isometric slope (beta = 1). If the 
## confidence interval of the males' slope is higher than 1, than bigger males have 
## relatively larger legs than small males. 

## Making a data frame with males only
males<-data1[data1$sex=="male",]

## Running the same model while scaling everything and excluding the sex factor
m3.males<-lmer(scale(log(p))~scale(log(ap))+(1|id),data=males)
summary(m3.males)

## Looking at the estimates and their confidence intervals
round(m3.males@beta,digits=5)
confint(m3.males)

## The estimated slope is around 1. Thus, bigger males do not have relatively
## longer legs than smaller males.

## Checking the plot..

plot(scale(log(p))~scale(log(ap)),data=males,
     las=1,bty='l',cex=1.5,pch=21,bg=alpha('black',0.75),
     xlab="Claw height (scaled and centered)",
     ylab="Pereopod length (scaled and centered)")

## This is the isometric line
abline(0,1,lwd=5,lty=2,col=alpha('grey',0.75))
## This is the estimated line
curve(m3.males@beta[1]+(m3.males@beta[2]*x),add=T,
      col=alpha('indianred1',0.75),lwd=5)


## Now, do females have more muscles on their legs than males?

plot(musc1~scale(peso_propodo), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (mg)",
     cex = 1.5, bty = "l",las =1)

## Uff, highly exponential. That is a bit tricky to handle without using complex
## models. Let's try logging it.

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Pereopod weight (scaled and centered)", 
     ylab = "Pereopod muscle's weight (log)",
     cex = 1.5, bty = "l",las =1)

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

## No difference, but the mean effect is borderline significant (p=0.06):
## Females have more muscle on average. Let's check the plot...

(coef.m4<-m4.gam@beta)

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Claw weight (scaled and centered)", 
     ylab = "Pereopod muscle weight (log)",
     cex = 1.5, bty = "l",las =1)

#Females
curve(exp(coef.m4[1]+(coef.m4[2]*x)),add=T,col='grey',lwd=3,
      from=-2,to=0.2)
#Males
curve(exp((coef.m4[1]+coef.m4[3])+((coef.m4[2]+coef.m4[4])*x)),add=T,col='black',
      lwd=3,from=-1,to=2.1)

## Our last test is the test between males. Do bigger males have relatively more muscles than
## smaller males?
## This analysis is a bit trickier. Since we used a Gamma distribution above
## we wanted to keep the same type of analysis. However, in this analysis
## we need to standardize our variables to compare them to an isometric line.
## When we standardize variables (scale() in R), we automatically generate 
## negative values. And Gamma distribution does not accept negative values
## (in the dependent variable). We thus needed to work around that.
## First, we were using GLMMs is to be able to use the pereopods of both sides
## of the individual. If it weren't for that, we could use "simpler" analysis. 
## Second, we were using a Gamma distribution because the variance increase with 
## the mean, which hindered using distributions that rely on "homokedasticity".
## The first step was to use only one side of the individual, that way we would not
## have random effects. The next step would be to use an analysis that did not rely
## on homokedasticity of our errors. That is why we opted for a GLS. In GLS we 
## can choose a weight to modify the variance-covariance matrix, thus relaxing 
## the need for homokedasticity in the errors. And that is what we did. We performed
## a GLS using only the right side of the individual (we performed we both sides
## results were identical).
## So, here goes the analysis.

males1<-males[males$side=="right",]
library(nlme)

peso.scale<-scale(log(males1$peso_propodo))
m4.males1<-gls(scale(log(musc1+1))~peso.scale,data=males1,
               weights = varExp(form=~as.vector(peso.scale)))
plot(m4.males1)
summary(m4.males1)
confint(m4.males1)

## The plot:
plot(scale(log(musc1+1))~scale(log(peso_propodo)),data=males1,
     las=1,bty='l',cex=1.5,pch=21,bg=alpha('black',0.75),
     xlab="Claw weight (scaled and centered)",
     ylab="Pereopod muscle weight (scaled and centered)")
abline(0,1,lwd=5,lty=2,col=alpha('grey',0.75))
curve(coef(m4.males1)[1]+(coef(m4.males1)[2]*x),add=T,
      col=alpha('indianred1',0.75),lwd=5)

## Now, I will plot both pannel figures and that should be it for the day
# Sexual dimorphism first
png(file="Figure3-PANEL.png", units='mm',width=300,height=90,res=600)

tiff(file="Figure3-PANEL.tiff",units="mm",width=300,height=90,res=600,
     compression="lzw")

par(mfrow=c(1,3),mar=c(4,4,2,3))

plot(left$resid~right$resid,bty='l',las=1,cex=1.5,
     pch=21,bg=c(alpha("gray",0.75),alpha("black",0.75))[as.numeric(left$sex)],
     ylab="Residuals of the left pereopod length",
     xlab="Residuals of the left pereopod length")
curve(coef(asym.p)[2]+(coef(asym.p)[1]*x),add=T,lwd=3,col=alpha("gray",0.75))
curve(coef(asym.p)[3]+((coef(asym.p)[1]+coef(asym.p)[4])*x),add=T,lwd=3,
      col=alpha("black",0.75),from=-0.7,to=0.7)
mtext("(a)",4,las=1,padj=-12,adj=-1)

plot(log(p)~scale(log(ap)), data = data1, pch = 21,ann=T,  
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)], 
     xlab = "Claw height (scaled and centered)", ylab = "Pereopod length (cm)",
     cex = 1.5, bty = "l",las =1)
mtext("(b)",4,las=1,padj=-12,adj=-1)
legend("topleft",c("Female","Male"),cex=1.3,bty='n',
       pch=21,pt.bg= c(alpha("gray",0.75), alpha("black",0.75)))

#Female
curve(coef.m3[1]+(coef.m3[2]*x),add=T,col='grey',lwd=3,from=-2,to=0.2)

#Male
curve((coef.m3[1]+coef.m3[3])+((coef.m3[2]+coef.m3[4])*x),add=T,col='black',
      lwd=3,from=-1,to=2.1)

plot(log(musc1+1)~scale(log(peso_propodo)), data = data1, pch = 21, 
     bg = c(alpha("gray",0.75), alpha("black",0.75))[as.numeric(sex)],
     xlab = "Claw weight (scaled and centered)", 
     ylab = "Pereopod muscle weight (log)",
     cex = 1.5, bty = "l",las =1)
mtext("(c)",4,las=1,padj=-12,adj=-1)

#Females
curve(exp(coef.m4[1]+(coef.m4[2]*x)),add=T,col='grey',lwd=3,
      from=-2,to=0.2)
#Males
curve(exp((coef.m4[1]+coef.m4[3])+((coef.m4[2]+coef.m4[4])*x)),add=T,col='black',
      lwd=3,from=-1,to=2.1)
dev.off()

# Analysis within males

png(file="Figure4-PANEL.png", units='mm',width=250,height=120,res=600)

tiff(file="Figure4-PANEL.tiff",units="mm",width=250,height=120,res=600,
     compression="lzw")

par(mfrow=c(1,2))

plot(scale(log(p))~scale(log(ap)),data=males,
     las=1,bty='l',cex=1.5,pch=21,bg=alpha('black',0.75),
     xlab="Claw height (scaled and centered)",
     ylab="Pereopod length (scaled and centered)")
abline(0,1,lwd=5,lty=2,col=alpha('grey',0.75))
curve(m3.males@beta[1]+(m3.males@beta[2]*x),add=T,
      col=alpha('indianred1',0.75),lwd=5)
mtext("(a)",4,las=1,padj=-12,adj=-1)

plot(scale(log(musc1+1))~scale(log(peso_propodo)),data=males1,
     las=1,bty='l',cex=1.5,pch=21,bg=alpha('black',0.75),
     xlab="Claw weight (scaled and centered)",
     ylab="Pereopod muscle weight (scaled and centered)")
abline(0,1,lwd=5,lty=2,col=alpha('grey',0.75))
curve(coef(m4.males1)[1]+(coef(m4.males1)[2]*x),add=T,
      col=alpha('indianred1',0.75),lwd=5)
mtext("(b)",4,las=1,padj=-12,adj=-1)
dev.off()
