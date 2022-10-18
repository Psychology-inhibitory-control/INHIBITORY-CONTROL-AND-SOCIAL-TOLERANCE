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


##################################GONOGO
#############################################################################################################

DATA_ALL_GONO<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/DATAALLGONOGO.csv", header=TRUE)
DATA_ALL_GONO = droplevels(DATA_ALL_GONO)

DATA_ALL_GONO$success<-as.factor(DATA_ALL_GONO$success)
summary(DATA_ALL_GONO)

DATA_ALL_GONO<-subset(DATA_ALL_GONO,type=="NoGo")


DATA_ALL_GONO<-subset(DATA_ALL_GONO,att==1)
DATA_ALL_GONO = droplevels(DATA_ALL_GONO)
summary(DATA_ALL_GONO)
##
##RFEPEATABILITY

rpt(success~ (1 | subj_id), grname = "subj_id", data = DATA_ALL_GONO, datatype = "Binary", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(success~  session +tolerance+  (1 | subj_id), grname = "subj_id", data = DATA_ALL_GONO, datatype = "Binary", nboot =1000, npermut =  1000)
#########
####TEST CHANCE LEVEL FOR SESSION 5

DATAGONO5<-subset(DATA_ALL_GONO,session==5)
psych::describeBy(DATAGONO5$success, group = DATAGONO5$subj_id)
DATA1= aggregate(DATAGONO5$success, by=list(DATAGONO5$subj_id,DATAGONO5$tolerance), mean, na.rm=TRUE)
colnames(DATA1) = c('subj_id','tolerance'  ,'perf2' )
DATA1_low<-subset(DATA1,tolerance=='low')
DATA1_medium<-subset(DATA1,tolerance=='medium')
DATA1_high<-subset(DATA1,tolerance=='high')


wilcox.test(DATA1_low$perf2, mu=0.5)
wilcox.test(DATA1_medium$perf2, mu=0.5)
wilcox.test(DATA1_high$perf2, mu=0.5)


summary(DATA1_low)
summary(DATA1_medium)
summary(DATA1_high)

sd(DATA1_low$perf2)
mean(DATA1_low$perf2)

sd(DATA1_medium$perf2)
mean(DATA1_medium$perf2)

sd(DATA1_high$perf2)
mean(DATA1_high$perf2)



###GRAPH

###

#DATA_ALL_GONO$session<-as.factor(DATA_ALL_GONO$session)
#ici
bar <- DATA_ALL_GONO%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot() + 
  scale_fill_manual(
    values = c(
      "low" = "#fcffa4",
      "medium" = "#f98e09", "high"="#bc3754"), labels = c("Low", "Medium", "High"))+


  
stat_summary(aes(session, success, fill=tolerance), 
               fun.y = mean, geom = "bar", position="dodge",width = 0.9) + 
  stat_summary(aes(session, success, fill=tolerance), 
               fun.data = mean_cl_normal, geom = "errorbar", 
               position=position_dodge(width=0.9), width = 0.2, size=0.7) + # position dodge separates the bars
  labs(x = "\nSession number", y = "Proportion of success \n")+labs(fill= "Degree of Tolerance") + 

  

  theme_classic(base_size = 22)+
  theme(legend.position =  c(0.3, 0.9),  legend.title = element_text(size=12),   legend.text = element_text(size=11))

bar+
  theme(
    axis.title.x = element_text(vjust=4.5),
    axis.title.y = element_text(vjust=0.01))

help(theme)

###NEW GAPH
DATA_g=aggregate(DATA_ALL_GONO$success, by=list(DATA_ALL_GONO$subj_id,DATA_ALL_GONO$session, DATA_ALL_GONO$tol), mean, na.rm=TRUE)

colnames(DATA_g) = c('subj_id','session','tolerance', 'perf' )


bar <- DATA_g%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot() +
  
  geom_point(aes(session, perf, 
                 group = tolerance, color =tolerance, shape = tolerance),position = position_jitterdodge(jitter.width = .4, 
                                                                                                         dodge.width = .7)) +
  
  scale_color_manual(values =  c(
    "low" = "#f5ea1b",
    "medium" = "#f98e09", "high"="#bc3754"), labels = c("Low", "Medium", "High"))+
  
  stat_summary(aes(session, perf, 
                   group = tolerance, color =tolerance),fun.y = mean, na.rm = TRUE, 
               geom = "point", 
               size = 4, color = "red", 
               position = position_dodge(width = .6)) +
  stat_summary(aes(session, perf, 
                   group = tolerance, color =tolerance),fun.data = mean_cl_normal, na.rm = TRUE, 
               geom = "errorbar", width = .2, color = "black",
               position = position_dodge(width = .6)) +
  
  theme_classic(base_size = 22)+
  theme(legend.position =  c(0.3, 0.9),  legend.title = element_text(size=12),   legend.text = element_text(size=11))+
  
  labs(x = "\n Session", y = "Proportion of success \n") + labs(fill="Degree of Tolerance")  +
  theme(legend.position = "none")+
  theme(
    axis.title.x = element_text(vjust=6.5),
    axis.title.y = element_text(vjust=0.01))


bar

###################################################################################################
###GLMER SUCCESS GONOGO
DATA_ALL_GONO$tolerance <- factor(DATA_ALL_GONO$tolerance, levels = c("low", "medium","high"))
##
Modelgono<- glmer(success~ sex + age + tolerance+ trial +session+ (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelgono)
###
t<-coef(summary(Modelgono))

write.table(t, file = "gono.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

###NORMAL DISTRIBUTION RESIDUALS

library(DHARMa)
res = simulateResiduals(Modelgono)
plot(res)

###


Modelgono0<- glmer(success~ 1+ (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Modelgono0,Modelgono)
summary(Modelgono)
##

postHocs<-glht(Modelgono, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

######
DATA_ALL_GONO_S5<-subset(DATA_ALL_GONO, session==5)
ModelgonoS5<- glmer(success~ sex + age + tolerance+ trial +(1|subj_id), data = DATA_ALL_GONO_S5,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(ModelgonoS5)

postHocs<-glht(ModelgonoS5, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

t<-coef(summary(ModelgonoS5))

write.table(t, file = "gono.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


##############PER SPECIES DIFFEENCES
#RFEPEATABILITY
DATA_ALL_GONO_rhesus$success<-as.factor(DATA_ALL_GONO_rhesus$success)
DATA_ALL_GONO_fasci$success<-as.factor(DATA_ALL_GONO_fasci$success)
DATA_ALL_GONO_tonk$success<-as.factor(DATA_ALL_GONO_tonk$success)


rpt(success~ (1 | subj_id), grname = "subj_id", data = DATA_ALL_GONO_tonk, datatype = "Binary", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(success~  session + (1 | subj_id), grname = "subj_id", data = DATA_ALL_GONO_fasci, datatype = "Binary", nboot =1000, npermut =  1000)
#
###TONKEAN#######################################################################################################

DATA_ALL_GONO_tonk<-subset(DATA_ALL_GONO,tolerance=="high")

Modelgonotonk<- glmer(success~  age + sex +  session+ trial+ (1|subj_id), data = DATA_ALL_GONO_tonk,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelgonotonk)

##############

Baseline <-glmer(success~ 1 + (1|subj_id), data = DATA_ALL_GONO_tonk,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelgonotonksession<- glmer(success~    session+   (1|subj_id), data = DATA_ALL_GONO_tonk,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Baseline, Modelgonotonksession)

t<-coef(summary(Modelgonotonk))

write.table(t, file = "gono.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

####RANK EFFECT
DATA_ALL_GONO_tonk_M<-subset(DATA_ALL_GONO_tonk,sex=="M")
Modelgonotonk_M<- glmer(success~  age + RANK+  session+ trial+ (1|subj_id), data = DATA_ALL_GONO_tonk_M,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelgonotonk_M0<- glmer(success~  age + session+ trial+ (1|subj_id), data = DATA_ALL_GONO_tonk_M,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(Modelgonotonk_M,Modelgonotonk_M0)
summary(Modelgonotonk_M)

###########################################################################################################
#FASCICULARIS
########

DATA_ALL_GONO_fasci<-subset(DATA_ALL_GONO,tolerance=="medium")
Modelgonofasci<- glmer(success~  age + sex + session+ trial+ (1|subj_id), data = DATA_ALL_GONO_fasci,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelgonofasci)


Baseline <-glmer(success~ 1+ (1|subj_id), data = DATA_ALL_GONO_fasci,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelgonofascisession<- glmer(success~  session+ (1|subj_id), data = DATA_ALL_GONO_fasci,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Baseline, Modelgonofascisession)

t<-coef(summary(Modelgonofasci))

write.table(t, file = "gono.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

####RANK EFFECT
DATA_ALL_GONO_fasci_F<-subset(DATA_ALL_GONO_fasci,sex=="F")
Modelgonofasci_F<- glmer(success~  age + RANK+  session+ trial+ (1|subj_id), data = DATA_ALL_GONO_fasci_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelgonofasci_F0<- glmer(success~  age + session+ trial+ (1|subj_id), data = DATA_ALL_GONO_fasci_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(Modelgonofasci_F, Modelgonofasci_F0)
summary(Modelgonofasci_F)


################################################################################
#######RHESUS

DATA_ALL_GONO_rhesus<-subset(DATA_ALL_GONO,tolerance=="low")
Modelgonorhesus<- glmer(success~  age + sex + trial+ session +location+ (1|subj_id), data = DATA_ALL_GONO_rhesus,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelgonorhesus)


Baseline <-glmer(success~ 1+ (1|subj_id), data = DATA_ALL_GONO_rhesus,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelgonorhesussession<- glmer(success~  session+ (1|subj_id), data = DATA_ALL_GONO_rhesus,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Baseline, Modelgonorhesussession)
##
t<-coef(summary(Modelgonorhesus))
write.table(t, file = "gono.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

####RANK EFFECT
DATA_ALL_GONO_rhesus_F<-subset(DATA_ALL_GONO_rhesus,sex=="F")
Modelgonorhesus_F<- glmer(success~  age + RANK+  session+ trial+ (1|subj_id), data = DATA_ALL_GONO_rhesus_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

Modelgonorhesus_F0<- glmer(success~  age +   session+ trial+ (1|subj_id), data = DATA_ALL_GONO_rhesus_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Modelgonorhesus_F,Modelgonorhesus_F0)

##
Modelgonorhesus_F<- glmer(success~    session+  (1|subj_id), data = DATA_ALL_GONO_rhesus_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

Modelgonorhesus_F0<- glmer(success~  1+ (1|subj_id), data = DATA_ALL_GONO_rhesus_F,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Modelgonorhesus_F,Modelgonorhesus_F0)



summary(Modelgonorhesus_F)

###########################
############


##
###BACKWARD SELECTION
Modelsex<- glmer(success~ sex+ age + tolerance+ session+ trial+ (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelsex)
##
Modelage<- glmer(success~  age + tolerance+ session+ trial+ (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelage)
#
Modeltrial<- glmer(success~  tolerance+ session+ trial+ (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modeltrial)
#
Modelopti<- glmer(success~  tolerance+ session+  (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modeltolerance)
#
Modelsession<- glmer(success~ tolerance+  (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

Modeltolerance<- glmer(success~ session+  (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova(Modelsession, Modelopti)

anova(Modeltolerance, Modelopti)

Baseline<- glmer(success~ 1+  (1|subj_id), data = DATA_ALL_GONO,family=binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Baseline, Modelsession)


###

##
psych::describeBy(DATA_ALL_GONO_rhesus$success, group =DATA_ALL_GONO_rhesus$location)
psych::describeBy(DATA_ALL_GONO$success, group =DATA_ALL_GONO$tolerance)

##



