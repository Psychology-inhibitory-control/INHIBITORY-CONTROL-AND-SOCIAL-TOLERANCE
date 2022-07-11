setwd("C:/Users/Marine/Desktop/ARTSPECIESCOMPA")


#######LOAD PACKAGES
library(languageR)
library(nlme)
library(car)
library(reshape)
library(ggplot2)
library(lme4)
library (pastecs)
library(nlme)
library(reshape)
library(tidyverse)
library(dplyr)
library(lmerTest)
library(stats)
library(multcomp)
library(tidyr)
library(stats)
library(forcats)
library(data.table)
library(dataPreparation)
library(rstatix)
library(RVAideMemoir)
#############
########REVERSAL LEARNING

DATA_TRIAL_RL<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/NB_TRIAL.csv", header=TRUE)
summary(DATA_TRIAL_RL)
DATA_TRIAL_RL$rule<-as.factor(DATA_TRIAL_RL$rule)

Modelnbtrials<-lmer(perf~ rule+ age +sex+ tolerance +(1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
summary(Modelnbtrials)

Modelnbtrials0<-lmer(perf~ 1 +(1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
Modelnbtrialsrule<-lmer(perf~ rule +(1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
anova(Modelnbtrials0, Modelnbtrialsrule)

##############
DATA_TRIAL_RL$rule<-as.factor(DATA_TRIAL_RL$rule)

DATA_TRIAL_RL%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
ggplot( aes(x=age, y=perf, group=tolerance, color = tolerance)) + 
  
  geom_point(aes(shape=tolerance, color=tolerance)) + 
  geom_smooth(aes(fill=tolerance),method=lm,  fullrange=TRUE) +
  theme_classic(base_size = 22)+
  labs( x = "Age", y = "Number of trials to learn the rules")


###GRAPHS

p<-DATA_TRIAL_RL%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot( aes(x=rule, y=perf, fill=tolerance)) +
  scale_fill_manual(
    values = c(
      "low" = "#722090",
      "medium" = "#339999", "high"="#FFE933"), labels = c("Low", "Medium", "High"))+

  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black", size=0.6, alpha=0.3) +

  labs(x = "\n Type of rule", y = "Number of trials to learn the rules \n") + labs(fill="Degree of Tolerance") + 

  scale_x_discrete(labels=c("Acquisition rule", "Reversed rule")) +
  
  theme_classic(base_size = 22)+

  theme(legend.position =  c(0.8, 0.8))+
  theme(
    axis.title.x = element_text(vjust=4.5),
    axis.title.y = element_text(vjust=0.01))

p
#########
########
DATA_TRIAL_RL_1<-subset(DATA_TRIAL_RL, rule==1)
psych::describeBy(DATA_TRIAL_RL_1$perf, group = DATA_TRIAL_RL_1$tolerance)


Modelnbtrials<-lmer(perf~ rule +tolerance +  (1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
summary(Modelnbtrials)
##RULE1
leveneTest(DATA_TRIAL_RL_1$perf, DATA_TRIAL_RL_1$tolerance ,center = median)
DiffModel<-kruskal.test(perf~tolerance, data=DATA_TRIAL_RL_1)
DiffModel
pairwise.wilcox.test(DATA_TRIAL_RL_1$perf, DATA_TRIAL_RL_1$tolerance,
                     p.adjust.method = "BH")

psych::describeBy(DATA_TRIAL_RL$perf, group = DATA_TRIAL_RL$tolerance)

###############

###RULE2
DATA_TRIAL_RL_2<-subset(DATA_TRIAL_RL, rule==2)
leveneTest(DATA_TRIAL_RL_2$perf, DATA_TRIAL_RL_2$tolerance ,center = median)

DiffModel<-aov(perf~tolerance, data=DATA_TRIAL_RL_2)
DiffModel0<-aov(perf~1, data=DATA_TRIAL_RL_2)
summary(DiffModel)

anova(DiffModel,DiffModel0)

summary.lm(DiffModel)

##POST HOC TUKEY
postHocs<-glht(DiffModel, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

psych::describeBy(DATA_TRIAL_RL_2$perf, group = DATA_TRIAL_RL_2$tolerance)

########################################################################################
####PER SPECIES

#####RHESUS
DATA_TRIAL_RL_rhesus<-subset(DATA_TRIAL_RL, tolerance=="low")
Modelnbtrials_rhesus<-lmer(perf~ location + age +sex+ rule +(1|subj_id),data=DATA_TRIAL_RL_rhesus, REML=FALSE)

Modelnbtrials_rhesus_loca<-lmer(perf~ location  *age +(1|subj_id),data=DATA_TRIAL_RL_rhesus, REML=FALSE)
summary(Modelnbtrials_rhesus_loca)

Modelnbtrials_rhesus<-lmer(perf~ 1 +(1|subj_id),data=DATA_TRIAL_RL_rhesus, REML=FALSE)
anova(Modelnbtrials_rhesus_loca, Modelnbtrials_rhesus)


Modelruler<-lmer(perf~location + (1|subj_id),data=DATA_TRIAL_RL_rhesus, REML=FALSE)
Modelrule0r<-lmer(perf~  1+ (1|subj_id),data=DATA_TRIAL_RL_rhesus, REML=FALSE)


anova(Modelrule0r, Modelruler)

#


#########
########RANK EFFECT
DATA_TRIAL_RL_rhesusF<-subset(DATA_TRIAL_RL_rhesus, sex=="         F")
Modelnbtrials_rhesusF<-lmer(perf~ location + age + RANK+ (1|subj_id),data=DATA_TRIAL_RL_rhesusF, REML=FALSE)
summary(Modelnbtrials_rhesusF)
Modelnbtrials_rhesusF0<-lmer(perf~ location + age +(1|subj_id),data=DATA_TRIAL_RL_rhesusF, REML=FALSE)
anova(Modelnbtrials_rhesusF,Modelnbtrials_rhesusF0)
t<-coef(summary(Modelnbtrials_rhesus))

write.table(t, file = "tablenbtrial.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
########

 #####FASCICULARIS
DATA_TRIAL_RL_fasci<-subset(DATA_TRIAL_RL, tolerance=="medium")
Modelnbtrials_fasci<-lmer(perf~  age +sex+rule + (1|subj_id),data=DATA_TRIAL_RL_fasci, REML=FALSE)
Modelrulef<-lmer(perf~ rule + (1|subj_id),data=DATA_TRIAL_RL_fasci, REML=FALSE)
Modelrule0f<-lmer(perf~  1 + (1|subj_id),data=DATA_TRIAL_RL_fasci, REML=FALSE)


anova(Modelrule0f, Modelrulef)
summary(Modelnbtrials_fasci)

t<-coef(summary(Modelnbtrials_fasci))

write.table(t, file = "tablenbtrial.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
########

#########
########RANK EFFECT
DATA_TRIAL_RL_fasciF<-subset(DATA_TRIAL_RL_fasci, sex=="         F")
Modelnbtrials_fasciF<-lmer(perf~ sex+ age + RANK+ (1|subj_id),data=DATA_TRIAL_RL_fasciF, REML=FALSE)
summary(Modelnbtrials_fasciF)
Modelnbtrials_fasciF0<-lmer(perf~ age +(1|subj_id),data=DATA_TRIAL_RL_fasciF, REML=FALSE)
anova(Modelnbtrials_fasciF,Modelnbtrials_fasciF0)
#
#####TONKEANS
DATA_TRIAL_RL_tonk<-subset(DATA_TRIAL_RL, tolerance=="high")
Modelnbtrials_tonk<-lmer(perf~  age +sex+rule + (1|subj_id),data=DATA_TRIAL_RL_tonk, REML=FALSE)
summary(Modelnbtrials_tonk)


Modelrulet<-lmer(perf~ rule + (1|subj_id),data=DATA_TRIAL_RL_tonk, REML=FALSE)
Modelrule0t<-lmer(perf~  1 + (1|subj_id),data=DATA_TRIAL_RL_tonk, REML=FALSE)


anova(Modelrule0t, Modelrulet)
summary(Modelnbtrials_fasci)

t<-coef(summary(Modelnbtrials_tonk))

write.table(t, file = "tablenbtrial.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


####
########RANK EFFECT
DATA_TRIAL_RL_tonkM<-subset(DATA_TRIAL_RL_tonk, sex=="         M")
Modelnbtrials_tonkM<-lmer(perf~  age + RANK+ (1|subj_id),data=DATA_TRIAL_RL_tonkM, REML=FALSE)
summary(Modelnbtrials_tonkM)
Modelnbtrials_tonkM0<-lmer(perf~ age +(1|subj_id),data=DATA_TRIAL_RL_tonkM, REML=FALSE)
anova(Modelnbtrials_tonkM,Modelnbtrials_tonkM0)

#

psych::describeBy(DATA_TRIAL_RL$perf, group = DATA_TRIAL_RL$tolerance)
#

postHocs<-glht(Modelnbtrials, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

##
##BACKWARD SELECTION
Modelsex<-lmer(perf~ tolerance + age +sex+rule + (1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
summary(Modelsex)
#
Modelage<-lmer(perf~ tolerance + age +rule + (1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
summary(Modelage)
#
Modeloptimum<-lmer(perf~ tolerance +rule + (1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
summary(Modeloptimum)
#
Modeltolerance<-lmer(perf~ rule + (1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)
MOdelrule<-lmer(perf~ tolerance +(1|subj_id),data=DATA_TRIAL_RL, REML=FALSE)

anova(Modeloptimum,Modeltolerance)

##
t<-coef(summary(Modelnbtrials))

write.table(t, file = "tablenbtrial.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
###GRAPHS


##

p<-DATA_TRIAL_RL%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  ggplot(aes(x=rule, y=perf, fill=tolerance))+
  
 # stat_boxplot(geom = "errorbar",
#               width = 0.15) +
  geom_boxplot() +

 # geom_jitter(color="black", size=0.6, alpha=0.3) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(x = "\nRègle", y = "Nombre d'essais pour apprendre\n") +
  scale_x_discrete(labels=c("Acquisition", "Inversée")) +
  theme_classic()
p


p <- DATA_ALL_STROOP_S1 %>%
  mutate(type_picture = fct_relevel(type_picture, 
                                    "None", "Object", "Neutral", "Threat"
  )) %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  ggplot(aes(x=type_picture, y=perf)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() +
  geom_jitter(color="black", size=0.6, alpha=0.3) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(x = "\nType of picture", y = "Distraction control score \n") +
  scale_x_discrete(labels=c("Control", "Object", "Neutral", "Threat")) +
  theme_classic(base_size = 20)+
  facet_grid(. ~ tolerance)


















##
DATA_TRIAL_RL_1<-subset(DATA_TRIAL_RL,rule==1)
DATA_TRIAL_RL_2<-subset(DATA_TRIAL_RL,rule==2)
psych::describeBy(DATA_TRIAL_RL$perf, group =DATA_TRIAL_RL$tolerance)

psych::describeBy(DATA_TRIAL_RL_2$perf, group =DATA_TRIAL_RL_2$tolerance)


##########################################################################################################
##################DIFFERENCE TRIALS
###########################################################################################################

DATA_DIFF_TRIAL<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/DIFF_TRIAL.csv", header=TRUE)
summary(DATA_DIFF_TRIAL)
###
DATA_DIFF_TRIAL_rhesus<-subset(DATA_DIFF_TRIAL, tolerance=="low" )

psych::describeBy(DATA_DIFF_TRIAL$perf, group =DATA_DIFF_TRIAL$tolerance)
##

p<-DATA_DIFF_TRIAL%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  ggplot( aes(x=tolerance, y=perf, fill=tolerance)) +
  scale_fill_viridis(discrete = T, name = "Level of tolerance", labels = c("Less", "Intermediate", "Highly"))+
  geom_boxplot() +
  labs(x = "\nDegree of tolerance", y = " Number of trials: Acquisition-Reversed rule\n") +
  scale_x_discrete(labels=c("Low", "Medium", "High")) +
  theme_classic(base_size = 22)

p 


leveneTest(DATA_DIFF_TRIAL$perf, DATA_DIFF_TRIAL$tolerance ,center = median)
leveneTest(DATA_DIFF_TRIAL_rhesus$perf, DATA_DIFF_TRIAL_rhesus$location ,center = median)

#Levene's Test for Homogeneity of Variance (center = median)
#Df F value Pr(>F)
#group  2  0.5158 0.6001

DiffModel<-aov(perf~tolerance, data=DATA_DIFF_TRIAL)
DiffModel<-aov(perf~age, data=DATA_DIFF_TRIAL_rhesus)
DiffModel<-aov(perf~sex, data=DATA_DIFF_TRIAL_rhesus)
DiffModel<-aov(perf~location, data=DATA_DIFF_TRIAL_rhesus)
summary(DiffModel)

summary.lm(DiffModel)

##POST HOC TUKEY
postHocs<-glht(DiffModel, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

#######################################################################################################
#########ACCURACY ON A TRIAL


####
DATA_ACCU_RL<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/ACCU_TRIAL.csv", header=TRUE)
DATA_ACCU_RL_rhesus<-subset(DATA_ACCU_RL, tolerance=="low")
DATA_ACCU_RL_fasci<-subset(DATA_ACCU_RL, tolerance=="medium")
DATA_ACCU_RL_tonk<-subset(DATA_ACCU_RL, tolerance=="high")


wilcox.test
DATA_ACCU_RL<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/ACCU_TRIAL.csv", header=TRUE)
DATATASK4<-DATA_ACCU_RL

DATATASK4_Wilcox<-subset(DATATASK4,tolerance=='high' &  session==1 & rule ==1)

DATATASK4_Wilcox= aggregate((DATATASK4_Wilcox$success), by=list(DATATASK4_Wilcox$subj_id), mean, na.rm=TRUE)
summary (DATATASK4_Wilcox)
colnames(DATATASK4_Wilcox) = c('subj_id','perf' )
wilcox.test(DATATASK4_Wilcox$perf, mu = 0.5)
#######################################################






summary(DATA_ACCU_RL)


psych::describeBy(DATA_ACCU_RL$success, group =DATA_ACCU_RL$tolerance)

##FASCi
Modelaccufasci<- glmer(success~ rule+ session +sex+age + (1|subj_id), data = DATA_ACCU_RL_fasci,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelaccufasci)

##TONK
Modelaccutonk<- glmer(success~ rule+ sex+ age+session+ (1|subj_id), data = DATA_ACCU_RL_tonk,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelaccutonk)

###RHESUS
Modelaccurhesus<- glmer(success~ rule+ session +sex+age +location+ (1|subj_id), data = DATA_ACCU_RL_rhesus,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelaccurhesus0<- glmer(success~ rule+ session +sex+age + (1|subj_id), data = DATA_ACCU_RL_rhesus,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Modelaccurhesus,Modelaccurhesus0)

summary(Modelaccurhesus)
##
t<-coef(summary(Modelaccurhesus))

write.table(t, file = "accur.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
###



DATA_ACCU_RL$tolerance <- factor(DATA_ACCU_RL$tolerance, levels = c("low", "medium","high"))

Modelaccu<- glmer(success~ rule+ session +sex+age + tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelaccu)
#

DATA_ACCU_RLna <- DATA_ACCU_RL %>%
  drop_na(success)
#
postHocs<-glht(Modelaccu, linfct=mcp(tolerance="Tukey"))
summary(postHocs)
##
t<-coef(summary(Modelaccu))

write.table(t, file = "accur.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#
Modelsex<- glmer(success~ rule+ session +sex+age + tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelaccu)
#
Modelopti<- glmer(success~ rule+ session+age + tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelopti2<- glmer(success~ rule+ session+age + tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova (Modelopti,Modelsex)
nobs(Modelopti)
nobs(Modeltol)
nobs(Modelrule)
nobs(Modelsession)
nobs(Modelage)
DATA_ACCU_RL %>% drop_na(success)

row.has.na <- apply(DATA_ACCU_RL, 1, function(rule){any(is.na(rule))})
sum(row.has.na)
#
Modelage<- glmer(success~ rule +  session+tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Modelopti,Modelage)
#

Modelrule<- glmer(success~ age +  session+tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Modelopti,Modelrule)
Modelopti<- glmer(success~ rule+ session+age + tolerance +(1|subj_id), data = data2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#
Modelsession<- glmer(success~ age + rule+tolerance +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Modelopti,Modelsession)
#
Modeltol<- glmer(success~ age + rule+session +(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)) )
anova (Modelopti,Modeltol)

############


###################################################################
no_missing_2 <-DATA_ACCU_RL %>%
  drop_na(rule, session, tolerance)

R1<-subset(DATA_ACCU_RL, rule==1) 
R2<-subset(DATA_ACCU_RL, rule==2)
 
##############RULE 1

ModelR1<- glmer(success~  session + tolerance+(1|subj_id), data = R1,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(ModelR1)
#
ModelR2<- glmer(success~  session +  tolerance+(1|subj_id), data = R2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(ModelR2)

postHocs<-glht(ModelR1, linfct=mcp(tolerance="Tukey"))
summary(postHocs)
postHocs<-glht(ModelR2, linfct=mcp(tolerance="Tukey"))
summary(postHocs)
######################
Modelsex<- glmer(success~ rule+ session +sex+age + tolerance+(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelsex)
#
Modelage<- glmer(success~ rule+ session +age + tolerance+(1|subj_id), data = DATA_ACCU_RL,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelage)

anova(Modelsex,Modelage)
anova(Modelopti,Modelage)

###################################
Modelopti<- glmer(success~tolerance+ session+(1|subj_id), data = no_missing_R2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelopti)
Modelrule<- glmer(success~ session +tolerance+(1|subj_id), data = no_missing_2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelsession<- glmer(success~ rule+tolerance+(1|subj_id), data = no_missing_2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

Modeltol<- glmer(success~ session+ (1|subj_id), data =no_missing_R2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modeltol)
Model0<- glmer(success~ 1+(1|subj_id), data =no_missing_2,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))



anova(Modelopti,Modelrule)
anova(Modelopti,Modelsession)
anova(Modelopti,Modeltol)
anova(Model0,Modeltol)




######################################################################################################
###NUMBER OF TAPS
DATA_ACCU_RL<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/ACCU_TRIAL.csv", header=TRUE)

DATA_ACCU_RL_wrong<-subset(DATA_ACCU_RL, success==0)
##

#DATA_ACCU_RL_wrong$session<-as.factor(DATA_ACCU_RL_wrong$session)
####

DATA_ACCU_RL_wrong$tolerance <- factor(DATA_ACCU_RL_wrong$tolerance, levels = c("low", "medium","high"))
############

psych::describeBy(DATA_ACCU_RL_wrong$nb_taps, group =DATA_ACCU_RL_wrong$tolerance)

#look at distribution of residals
hist(log(DATA_ACCU_RL_wrong$nb_taps))
par(mfrow=c(1,2))
qqnorm(residuals(Modeltaps))
hist(residuals(Modeltaps))
################
psych::describeBy(DATA_ACCU_RL_wrong$nb_taps, group =DATA_ACCU_RL_wrong$rule)
###################################

Modeltaps<-lmer(log(nb_taps)~ age + sex + tolerance+ rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong, REML=FALSE)
summary(Modeltaps)
postHocs<-glht(Modeltaps, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

DATA_ACCU_RL_wrong$tolerance <- factor(DATA_ACCU_RL_wrong$tolerance, levels = c("low", "medium","high"))
t<-coef(summary(Modeltaps))

write.table(t, file = "taps.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

######Sans les 6 rhesus

DATA_ACCU_RL_2<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/react2.csv", header=TRUE)
DATA_ACCU_RL_wrong_witout<-subset(DATA_ACCU_RL_2, success==0)

Modeltaps<-lmer(log(nb_taps)~ age + sex + tolerance+ rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_witout, REML=FALSE)
summary(Modeltaps)
postHocs<-glht(Modeltaps, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

###REPEATABILITY

rpt(log(nb_taps)~ (1 | subj_id), grname = "subj_id", data = DATA_ACCU_RL_wrong_tonk, datatype = "Gaussian", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(log(nb_taps)~ session +rule + sex+(1 | subj_id), grname = "subj_id", data = DATA_ACCU_RL_wrong_tonk, datatype = "Gaussian", nboot =1000, npermut =  1000)


###########Backward

Modeltapsopti<-lmer(log(nb_taps)~  tolerance+ rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong, REML=FALSE)
summary(Modeltapsopti)
Modeltapstol<-lmer(log(nb_taps)~   rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong, REML=FALSE)
Modeltapsrule<-lmer(log(nb_taps)~  tolerance+  session +(1|subj_id),data=DATA_ACCU_RL_wrong, REML=FALSE)
Modeltapssession<-lmer(log(nb_taps)~  tolerance+  rule +(1|subj_id),data=DATA_ACCU_RL_wrong, REML=FALSE)

anova(Modeltapsopti, Modeltapsrule)
anova(Modeltapsopti, Modeltapssession)
anova(Modeltapsopti, Modeltapstol)
###
##

##############RHESUS
DATA_ACCU_RL_wrong_rhe<-subset(DATA_ACCU_RL_wrong, tolerance=="low")


Modeltapsrhe<-lmer(log(nb_taps)~ age + sex + location+ rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)
summary(Modeltapsrhe)

########RANK EFFECT

DATA_ACCU_RL_wrong_rhe_F<-subset(DATA_ACCU_RL_wrong_rhe, sex=="F" )
ModeltapsrheF<-lmer(log(nb_taps)~ age + sex +rule+ session + location+ RANK+(1|subj_id),data=DATA_ACCU_RL_wrong_rhe_F, REML=FALSE)
summary(ModeltapsrheF)
ModeltapsrheF0<-lmer(log(nb_taps)~ age +location+ rule+ session + (1|subj_id),data=DATA_ACCU_RL_wrong_rhe_F, REML=FALSE)

anova(ModeltapsrheF, ModeltapsrheF0)

########

Modeltapsrheopti<-lmer(log(nb_taps)~   location +rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)

Modeltapsrheloca<-lmer(log(nb_taps)~   rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)
Modeltapsrherule<-lmer(log(nb_taps)~  location+session +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)
Modeltapsrhesession<-lmer(log(nb_taps)~  location+ rule+(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)

anova(Modeltapsrheopti,Modeltapsrheloca)
anova(Modeltapsrheopti,Modeltapsrherule)
anova(Modeltapsrheopti,Modeltapsrhesession)

Modeltapsrhe0<-lmer(log(nb_taps)~  1 +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)

Modeltapsrhetest<-lmer(log(nb_taps)~  location +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)
anova(Modeltapsrhe0,Modeltapsrhetest )
##############
Modeltapsrhe0<-lmer(log(nb_taps)~  sex +  session +(1|subj_id),data=DATA_ACCU_RL_wrong_rhe, REML=FALSE)

anova(Modeltapsrhe0, Modeltapsrhe)

t<-coef(summary(Modeltapsrhe))

write.table(t, file = "tapsrhe.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

psych::describeBy(DATA_ACCU_RL_wrong_rhe$nb_taps, group =DATA_ACCU_RL_wrong_rhe$location)





##############FASCI
DATA_ACCU_RL_wrong_fasci<-subset(DATA_ACCU_RL_wrong, tolerance=="medium")

Modeltapsfasci<-lmer(log(nb_taps)~ age + sex + rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_fasci, REML=FALSE)
summary(Modeltapsfasci)

Modeltapsfascim<-lmer(log(nb_taps)~  rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_fasci, REML=FALSE)
Modeltapsfasci0<-lmer(log(nb_taps)~ session +(1|subj_id),data=DATA_ACCU_RL_wrong_fasci, REML=FALSE)
Modeltapsfascirule<-lmer(log(nb_taps)~ rule +(1|subj_id),data=DATA_ACCU_RL_wrong_fasci, REML=FALSE)
anova(Modeltapsfascim,Modeltapsfasci0)
anova(Modeltapsfascim,Modeltapsfascirule)

t<-coef(summary(Modeltapsfasci))

write.table(t, file = "tapsfasci.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

##
########RANK EFFECT

DATA_ACCU_RL_wrong_fasci_F<-subset(DATA_ACCU_RL_wrong_fasci, sex=="F")
ModeltapsfasciF<-lmer(log(nb_taps)~ age +rule+ session + RANK+(1|subj_id),data=DATA_ACCU_RL_wrong_fasci_F, REML=FALSE)
summary(ModeltapsfasciF)
ModeltapsfasciF0<-lmer(log(nb_taps)~ age + rule+ session + (1|subj_id),data=DATA_ACCU_RL_wrong_fasci_F, REML=FALSE)

anova(ModeltapsfasciF, ModeltapsfasciF0)

######TONK
DATA_ACCU_RL_wrong_tonk<-subset(DATA_ACCU_RL_wrong, tolerance=="high")

Modeltapstonk<-lmer(log(nb_taps)~ age + sex + rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_tonk, REML=FALSE)
summary(Modeltapstonk)


Modeltapstonk<-lmer(log(nb_taps)~  sex + rule+ session +(1|subj_id),data=DATA_ACCU_RL_wrong_tonk, REML=FALSE)

Modeltapstonk0<-lmer(log(nb_taps)~  session+  sex+(1|subj_id),data=DATA_ACCU_RL_wrong_tonk, REML=FALSE)

anova(Modeltapstonk0, Modeltapstonk)

t<-coef(summary(Modeltapstonk))

write.table(t, file = "tapstnk.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

psych::describeBy(DATA_ACCU_RL_wrong_tonk$nb_taps, group =DATA_ACCU_RL_wrong_tonk$sex)
######
########RANK EFFECT

DATA_ACCU_RL_wrong_tonk_M<-subset(DATA_ACCU_RL_wrong_tonk, sex=="M")
ModeltapstonkM<-lmer(log(nb_taps)~ age +rule+ session + RANK+(1|subj_id),data=DATA_ACCU_RL_wrong_tonk_M, REML=FALSE)
summary(ModeltapstonkM)
ModeltapstonkM0<-lmer(log(nb_taps)~ age + rule+ session + (1|subj_id),data=DATA_ACCU_RL_wrong_tonk_M, REML=FALSE)

anova(ModeltapstonkM, ModeltapstonkM0)



##GRAPHS


#######################################
DATA_ACCU_RL_wrong$rule<-as.factor(DATA_ACCU_RL_wrong$rule)
p<- DATA_ACCU_RL_wrong %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot(aes(x=tolerance, y=nb_taps)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(x = "\nDegree of Tolerance", y = "Number of taps on the wrong stimulus \n") +
  
  theme_classic(base_size = 20)

p

p<- DATA_ACCU_RL_wrong %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot(aes(x=session, y=nb_taps)) +
  
  geom_jitter(aes(colour = rule))+
  geom_boxplot(scale = "area")+
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red") +
  labs(x = "\nSession number", y = "Number of perseverative errors\n") +  
  
  
  
  theme_classic(base_size = 24)+
  facet_grid(tolerance~ rule,
  labeller = labeller(
    rule = c(`1` = "acquisition rule", `2` = "reversed rule"),
    tolerance = c(`low` = "low tolerance", `medium` = "medium tolerance", `high` = "high tolerance")))
  
  

p 

#####
 
###################

##########


psych::describeBy(DATA_ACCU_RL_rhesus$good_resp, group =DATA_ACCU_RL_rhesus$location)




bar <- DATA_ACCU_RL_wrong%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot() + 
  stat_summary(aes(tolerance, good_resp), 
               fun.y = mean, geom = "bar", position="dodge",width = 0.6) + 
  stat_summary(aes(tolerance,good_resp), 
               fun.data = mean_cl_normal, geom = "errorbar", 
               position=position_dodge(width=0.9), width = 0.2, size=1.1) + # position dodge separates the bars
  labs(x = "\nTolerance Degree", y = "Proportion of success \n") + 
  #scale_x_discrete(labels=c("Go", "No-go")) +
  theme(axis.text.x=element_text(size=16))+
  
  
  theme(axis.text=element_text(size=14), axis.title = element_text(size=28,family="serif", face="bold"))+
  
  theme(text=element_text(family="serif",  size=28))+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_grid(. ~ session)

bar
