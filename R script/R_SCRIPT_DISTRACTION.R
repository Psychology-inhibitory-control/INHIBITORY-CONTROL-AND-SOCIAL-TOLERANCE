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
library(rptR)
        
#########################################################
#########################################################
####################DISTRACTION TASK#####################

#########################################################
##########LOAD DATA

DATA_ALL_STROOP<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/STROOPALL.csv", header=TRUE)
DATA_ALL_STROOP<-subset(DATA_ALL_STROOP,react_time>200 & react_time<35000)
DATA_ALL_STROOP = droplevels(DATA_ALL_STROOP)
summary(DATA_ALL_STROOP)

#######################################################
##TANSFORM DATA to HAVE NORMALLY DISTIBUTED RESIDUALS
max<-max(DATA_ALL_STROOP$perf, na.rm = TRUE)
DATA_ALL_STROOP$perf2<-sqrt((max+1)-DATA_ALL_STROOP$perf)
##

####################
###########MEAN RESPONSE LATENCY FOR CONTROL TRIALS
DATA_ALL_STROOP_NONE<-subset(DATA_ALL_STROOP, stim=="None")
psych::describeBy(DATA_ALL_STROOP_NONE$react_time, group = DATA_ALL_STROOP_NONE$tolerance)
psych::describeBy(DATA_ALL_STROOP_NONE$react_time)

Modelreact<-lmer(react_time~ age +sex+RANK +tolerance + trial+session+ (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modelreact)
##

########NORMAL DISTRIBUTION RESIDUALS
par(mfrow=c(1,2))
qqnorm(residuals(Modelreact))
hist(residuals(Modelreact))

#PRINT TABLE
##############
t<-coef(summary(Modelreact))

write.table(t, file = "REACT.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
#BACKWARD SELECTION
##NO EFFECT 
Modelrank<-lmer(react_time~ age +sex+RANK +tolerance + trial+session+ (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modelrank)
Modeltrial<-lmer(react_time~ age +sex +tolerance + trial+session+ (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modeltrial)
Modelsex<-lmer(react_time~ age +sex +tolerance +session+  (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modelsex)
Modeltolerance<-lmer(react_time~ age +tolerance +session+  (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modeltolerance)

Modeloptimum<-lmer(react_time~ age +session+ (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)
summary(Modeloptimum)
#######
Modelage<-lmer(react_time~  session+ (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)

Modelsession<-lmer(react_time~ age + (1|subj_id),data=DATA_ALL_STROOP_NONE, REML=FALSE)

anova(Modeloptimum, Modelage)

anova(Modeloptimum, Modelsession) 

###############################################################################################################
##############DISTACTION CONTROL SCORE
DATA_ALL_STROOP<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/STROOPALL.csv", header=TRUE)
DATA_ALL_STROOP<-subset(DATA_ALL_STROOP,react_time>200 & react_time<35000)
DATA_ALL_STROOP = droplevels(DATA_ALL_STROOP)
summary(DATA_ALL_STROOP)

#######################################################
##TANSFORM DATA to HAVE NORMALLY DISTIBUTED RESIDUALS
max<-max(DATA_ALL_STROOP$perf, na.rm = TRUE)
DATA_ALL_STROOP$perf2<-sqrt((max+1)-DATA_ALL_STROOP$perf)
###

######SET CONTRASTS

.none_vs_pict<-c(1,-1,0,0)
.none_vs_object<-c(0, -1, 1, 0)
.none_vs_threat<-c(0,-1,0,1)
m<-cbind(.none_vs_neutral,.none_vs_object,.none_vs_threat)
contrasts(DATA_ALL_STROOP$type_picture)<-m
##
.none_vs_pic<-c(1, -3, 1, 1)
.object_vs_face<-c(1, 0, -2, 1)
.threat_vs_neutral<-c(1, 0, 0, -1)
m<-cbind(.none_vs_pic,.object_vs_face,.threat_vs_neutral)

contrasts(DATA_ALL_STROOP$type_picture)<-m
#################################################################################

################################################################################
##ALL SESSIONS
################################################################################
psych::describeBy(DATA_ALL_STROOP$perf, group = DATA_ALL_STROOP$tolerance)    

DATA_ALL_STROOP$tolerance <- factor(DATA_ALL_STROOP$tolerance, levels = c("low", "medium","high"))

Modelstroop<-lmer(perf2~ age +sex+tolerance + type_picture+ trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
summary(Modelstroop)
#######
Modelstroop0<-lmer(perf2~ 1+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)

anova(Modelstroop,Modelstroop0)
##
summary(DATA_ALL_STROOP_face)
sd(DATA_ALL_STROOP_face$perf)
psych::describeBy(DATA_ALL_STROOP$perf2, group = DATA_ALL_STROOP$tolerance) 

DATA_ALL_STROOP_face<-subset(DATA_ALL_STROOP, type_picture=="Threat" | type_picture=="Neutral"| type_picture=="Object")
psych::describeBy(DATA_ALL_STROOP_face$perf, group = DATA_ALL_STROOP_face$type_picture) 
##CHECKING ASSUMPTIONS NORMLITY RESIDUALS

par(mfrow=c(1,2))
qqnorm(residuals(Modelstroop))
hist(residuals(Modelstroop))
#####################

####TUKEY TASK

postHocs<-glht(Modelstroop, linfct=mcp(tolerance="Tukey"))
summary(postHocs)
postHocs<-glht(Modelstroop, linfct=mcp(type_picture="Tukey"))
summary(postHocs)

t<-coef(summary(postHocs))

write.table(t, file = "STROOPALL.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
#######################
#PRINT TABLE
##############
t<-coef(summary(Modelstroop))

write.table(t, file = "STROOPALL.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#BACKWARD SELECTION
##NO EFFECT

Modelstroopage<-lmer(perf2~ age +sex+tolerance + type_picture+ trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
summary(Modelstroopage)
Modelstroopsex<-lmer(perf2~sex+tolerance + type_picture+ trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
summary(Modelstroopsex)
##
Modelstroopopti<-lmer(perf2~ tolerance + type_picture+ trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
summary(Modelstroopopti)
#
Modelstroopotoler<-lmer(perf2~ type_picture+ trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
Modelstroopotpict<-lmer(perf2~ tolerance +trial+session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
Modelstrooptrial<-lmer(perf2~ type_picture+ tolerance +session+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)
Modelstroopsession<-lmer(perf2~ type_picture+ tolerance +trial+ (1|subj_id),data=DATA_ALL_STROOP, REML=FALSE)

anova(Modelstroopopti,Modelstroopotoler)
anova(Modelstroopopti,Modelstroopotpict)
      anova(Modelstroopopti,Modelstroopsession)
            anova(Modelstroopopti,Modelstrooptrial)
            
#####################################################
####SESSION 1
            
            
                       
DATA_ALL_STROOP_S1<-subset(DATA_ALL_STROOP, session==1 )
            DATA_ALL_STROOP_S1_pic<-subset(DATA_ALL_STROOP, stim=="pic" )
       summary(DATA_ALL_STROOP_S1_pic)     
       sd(DATA_ALL_STROOP_S1_pic$perf)         
            
psych::describeBy(DATA_ALL_STROOP_S1$perf, group = DATA_ALL_STROOP_S1$type_picture)                        

DATA_ALL_STROOP_face_S1<-subset(DATA_ALL_STROOP_S1, type_picture=="Threat" | type_picture=="Neutral"| type_picture=="Object")
psych::describeBy(DATA_ALL_STROOP_face_S1$perf, group = DATA_ALL_STROOP_face_S1$type_picture) 
summary(DATA_ALL_STROOP_face_S1)
sd(DATA_ALL_STROOP_face_S1$perf)
#############################SESSION 1

#########################################

Modelstroop1<-lmer(perf2~ age +sex +tolerance + type_picture+ trial+(1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
summary(Modelstroop1)
#######
Modelstroop0<-lmer(perf2~ 1+ (1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)

anova(Modelstroop1,Modelstroop0)

##CHECKING ASSUMPTIONS NORMLITY RESIDUALS

par(mfrow=c(1,2))
qqnorm(residuals(Modelstroop1))
hist(residuals(Modelstroop1))
#####################
psych::describeBy(DATA_ALL_STROOP_S1$perf2, group = DATA_ALL_STROOP_S1$tolerance)
####TUKEY TASK

postHocs<-glht(Modelstroop1, linfct=mcp(tolerance="Tukey"))
summary(postHocs)
postHocs<-glht(Modelstroop1, linfct=mcp(type_picture="Tukey"))
summary(postHocs)
#######################
#PRINT TABLE
##############
t<-coef(summary(Modelstroop1))

write.table(t, file = "STROOPALL.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#####BACKWARD SELECTION STROOP 1


Modelstroop1age<-lmer(perf2~ age +sex+tolerance + type_picture+ trial+(1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
summary(Modelstroop1age)
#
Modelstroop1sex<-lmer(perf2~ sex +tolerance + type_picture+ trial+(1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
summary(Modelstroop1sex)
#
Modelstroop1trial<-lmer(perf2~ tolerance + type_picture+ trial+(1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
summary(Modelstroop1trial)
#

#
#
Modelstroop1opti<-lmer(perf2~ tolerance + type_picture+ (1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
summary(Modelstroop1opti)
##
Modelstroop1tolerance<-lmer(perf2~  type_picture+ (1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)
Modelstroop1typepict<-lmer(perf2~  tolerance+ (1|subj_id),data=DATA_ALL_STROOP_S1, REML=FALSE)

anova(Modelstroop1opti,Modelstroop1tolerance)
anova(Modelstroop1opti,Modelstroop1typepict)

###################################################################################################################
###EFFECT SPECIES

#####RHESUS
################
DATA_ALL_STROOPL<-subset(DATA_ALL_STROOP,tolerance=="low" & session==1)
##
Modelstroop_rhesus<-lmer(perf2~ age + sex +type_picture+ location +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
summary(Modelstroop_rhesus)

postHocsR<-glht(Modelstroop_rhesus, linfct=mcp(type_picture="Tukey"))
summary(postHocsR)
##
psych::describeBy(DATA_ALL_STROOPL$perf, group = DATA_ALL_STROOPL$sex)    
psych::describeBy(DATA_ALL_STROOPL$perf, group = DATA_ALL_STROOPL$type_picture)  
psych::describeBy(DATA_ALL_STROOPL$perf, group = DATA_ALL_STROOPL$stim)  
###
t<-coef(summary(Modelstroop_rhesus))

write.table(t, file = "STROOPALL.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
###########################################
#BACKWARD SELECTION
Modelstroop_rhesus_trial<-lmer(perf2~ age + sex +type_picture+ location +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
summary(Modelstroop_rhesus_trial)
Modelstroop_rhesus_location<-lmer(perf2~ age + sex +type_picture+ location +  (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
summary(Modelstroop_rhesus_location)
Modelstroop_rhesus_age<-lmer(perf2~ age + sex +type_picture+   (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
summary(Modelstroop_rhesus_age)
Modelstroop_rhesus_opti<-lmer(perf2~ sex +type_picture+   (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
summary(Modelstroop_rhesus_opti)
##
anova(Modelstroop_rhesus_location,Modelstroop_rhesus_age)

Modelstroop_rhesus_sex<-lmer(perf2~ type_picture+   (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)

Modelstroop_rhesus_pict<-lmer(perf2~ sex +  (1|subj_id/type_picture),data=DATA_ALL_STROOPL, REML=FALSE)
###
anova(Modelstroop_rhesus_opti, Modelstroop_rhesus_sex)
anova(Modelstroop_rhesus_opti, Modelstroop_rhesus_pict)
##########################
#EFFECT OF THE RANK
DATA_ALL_STROOPL_rank<-subset(DATA_ALL_STROOP,tolerance=="low"& sex=="F")

Modelstroop_rhesus_rank<-lmer(perf2~ age +  RANK +type_picture+ location +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPL_rank, REML=FALSE)

summary(Modelstroop_rhesus_rank)
Modelstroop_rhesus_opti<-lmer(perf2~ age  +type_picture+ location +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPL_rank, REML=FALSE)
anova(Modelstroop_rhesus_rank,Modelstroop_rhesus_opti)

t<-coef(summary(Modelstroop_rhesus_rank))

write.table(t, file = "rank.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


#####FASCICULARIS
################
##
################
DATA_ALL_STROOPM<-subset(DATA_ALL_STROOP,tolerance=="medium" & session==1)
##
Modelstroop_fasci<-lmer(perf2~ age + sex + type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)
summary(Modelstroop_fasci)

postHocsR<-glht(Modelstroop_fasci, linfct=mcp(type_picture="Tukey"))
summary(postHocsR)
##
psych::describeBy(DATA_ALL_STROOPM$perf, group = DATA_ALL_STROOPM$sex)    
psych::describeBy(DATA_ALL_STROOPM$perf, group = DATA_ALL_STROOPM$type_picture) 
###
t<-coef(summary(Modelstroop_fasci))

write.table(t, file = "STROOPALLM.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


DATA_ALL_STROOPM_face<-subset (DATA_ALL_STROOPM, type_picture=="Threat"| type_picture=="Neutral")
summary(DATA_ALL_STROOPM_face)
sd(DATA_ALL_STROOPM_face$perf)

#####
###REPEATABILITY
rpt(perf2~ (1 | subj_id), grname = "subj_id", data = DATA_ALL_STROOPL, datatype = "Gaussian", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(perf2~ sex+ type_picture+  (1 | subj_id), grname = "subj_id", data = DATA_ALL_STROOPL, datatype = "Gaussian", nboot =1000, npermut =  1000)



##########################


#BACKWARD SELECTION
Modelstroop_fasci_age<-lmer(perf2~ age + sex +type_picture+trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)
summary(Modelstroop_fasci_age)
Modelstroop_fasci_trial<-lmer(perf2~  sex +type_picture+ trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)
summary(Modelstroop_fasci_trial)
Modelstroop_fasci_opti<-lmer(perf2~ sex +type_picture+   (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)
summary(Modelstroop_fasci_opti)
##
Modelstroop_fasci_sex<-lmer(perf2~ type_picture+   (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)

Modelstroop_fasci_pict<-lmer(perf2~ sex +  (1|subj_id/type_picture),data=DATA_ALL_STROOPM, REML=FALSE)
###
anova(Modelstroop_fasci_opti, Modelstroop_fasci_sex)
anova(Modelstroop_fasci_opti, Modelstroop_fasci_pict)




###########################
#EFFECT OF THE RANK
DATA_ALL_STROOPM_rank<-subset(DATA_ALL_STROOPH, sex=="M")

Modelstroop_fasci_rank<-lmer(perf2~ age +  RANK +type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPM_rank, REML=FALSE)
summary(Modelstroop_fasci_rank)
Modelstroop_fasci_opti<-lmer(perf2~ age  +type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPM_rank, REML=FALSE)
anova(Modelstroop_fasci_rank,Modelstroop_fasci_opti)

summary(Modelstroop_fasci_rank)
t<-coef(summary(Modelstroop_fasci_rank))

write.table(t, file = "rank.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

postHocsR<-glht(Modelstroop_fasci_rank, linfct=mcp(type_picture="Tukey"))
summary(postHocsR)
psych::describeBy(DATA_ALL_STROOPM_rank$perf, group = DATA_ALL_STROOPM_rank$type_picture)
#####TONKEAN
################
################
DATA_ALL_STROOPH<-subset(DATA_ALL_STROOP,tolerance=="high" & session==1)
##
Modelstroop_tonk<-lmer(perf2~ experience +age + sex + type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPH, REML=FALSE)
summary(Modelstroop_tonk)

Modelstroop0tonk<-lmer(perf2~age + sex + type_picture +trial+   (1|subj_id/type_picture),data=DATA_ALL_STROOPH, REML=FALSE)
Modelstroop_tonk_exp<-lmer(perf2~ experience +age + sex + type_picture +trial+    (1|subj_id/type_picture),data=DATA_ALL_STROOPH, REML=FALSE)

Modelstroop0tonk_pict<-lmer(perf2~experience +age + sex + type_picture +trial+    (1|subj_id/type_picture),data=DATA_ALL_STROOPH, REML=FALSE)
Modelstroop_tonk_pict<-lmer(perf2~ experience +age + sex +trial+    (1|subj_id/type_picture),data=DATA_ALL_STROOPH, REML=FALSE)
anova(Modelstroop0tonk_pict,Modelstroop_tonk_pict)

anova(Modelstroop_tonk_exp,Modelstroop0tonk)

postHocsR<-glht(Modelstroop_tonk, linfct=mcp(type_picture="Tukey"))
summary(postHocsR)
##
psych::describeBy(DATA_ALL_STROOPH$perf, group = DATA_ALL_STROOPH$sex)    

###
t<-coef(summary(Modelstroop_tonk))

write.table(t, file = "STROOPALLH.txt", sep = "\t",
            row.names = TRUE, col.names = NA)
##########################
#EFFECT OF THE RANK
DATA_ALL_STROOPH_rank<-subset(DATA_ALL_STROOP,tolerance=="high"& sex=="M")
Modelstroop_tonk_opti<-lmer(perf2~ age  +type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPH_rank, REML=FALSE)
anova(Modelstroop_tonk_rank,Modelstroop_tonk_opti)

Modelstroop_tonk_rank<-lmer(perf2~ age +  RANK +type_picture +trial+  (1|subj_id/type_picture),data=DATA_ALL_STROOPH_rank, REML=FALSE)
summary(Modelstroop_tonk_rank)
t<-coef(summary(Modelstroop_tonk_rank))

write.table(t, file = "rank.txt", sep = "\t",
            row.names = TRUE, col.names = NA)




###
########################################################################################################
##GRAPHS ici

####
b<-boxplot(DATA_ALL_STROOP_S1_T$perf ~ DATA_ALL_STROOP_S1_T$tolerance , col=terrain.colors(4) )

# Add DATA_ALL_STROOP points
ggplot(DATA_ALL_STROOP_S1, aes(x= tolerance,y = perf, fill=tolerance)) + 
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot()+
geom_jitter(color="black", size=0.6, alpha=0.3) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  

##Effect tolerance
   

p <-DATA_ALL_STROOP_S1%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  ggplot(aes(x=tolerance, y=perf, fill=tolerance)) +
  
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot(outlier.shape = NA)  +
  scale_fill_manual(
    values = c(
      "low" = "#722090",
      "medium" = "#339999", "high"="#FFE933"))+

  geom_jitter(color="black", size=0.6, alpha=0.3) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=6, color="red", fill="red") +
  labs(x = "\nDegree of Tolerance", y = "Distraction control score \n") +
 scale_x_discrete(labels=c("Low", "Medium", "High")) +
  theme_classic(base_size = 22)+
  theme(legend.position = "none")+
  theme(
    axis.title.x = element_text(vjust=4.5),
    axis.title.y = element_text(vjust=0.01))


p 

install.packages("viridis")
library(viridis)
########
##Effect type pict

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

p

###
##Effect type picture tolerance
p<- DATA_ALL_STROOP %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  ggplot( aes(x=sex, y=perf2, fill=type_picture)) +
  geom_boxplot() +
  theme_classic(base_size = 20)+
  facet_grid(. ~ tolerance)
p 
##
###Effect sex
p<- DATA_ALL_STROOP_S1 %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  ggplot(aes(x=sex, y=perf, fill=sex)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) + 
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.3) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=8, color="red", fill="red") +
  labs(x = "\nSex of the subjects", y = "Distraction control score \n") +
  scale_x_discrete(labels=c("Female", "Male"))+
  
  theme_classic(base_size = 20)+
  
  facet_grid(. ~ tolerance)

p 



###########################################


#######ABANDON
library(pscl)
library(rcompanion)
library(multcompView)
library(emmeans)


DATA_ALL_STROOP_abdn<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/abd.csv", header=TRUE)


DATA_ALL_STROOP_abdn = droplevels(DATA_ALL_STROOP_abdn)
summary(DATA_ALL_STROOP_abdn)

.none_vs_neutral<-c(1,-1,0,0)
.none_vs_object<-c(0, -1, 1, 0)
.none_vs_threat<-c(0,-1,0,1)
m<-cbind(.none_vs_neutral,.none_vs_object,.none_vs_threat)
contrasts(DATA_ALL_STROOP_abdn$type_picture)<-m
###########
.none_vs_pic<-c(1, -3, 1, 1)
.object_vs_face<-c(1, 0, -2, 1)
.threat_vs_neutral<-c(1, 0, 0, -1)
m<-cbind(.none_vs_pic,.object_vs_face,.threat_vs_neutral)

contrasts(DATA_ALL_STROOP_abdn$type_picture)<-m
################################################

DATA_ALL_STROOP_abdn$tolerance <- factor(DATA_ALL_STROOP_abdn$tolerance, levels = c("low", "medium","high"))

Modelpoiss<- glmer(abdn~ age + sex + tolerance+ type_picture + (1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
##
#
psych::describeBy(DATA_ALL_STROOP_abdn$abdn, group = DATA_ALL_STROOP_abdn$type_picture) 

##RFEPEATABILITY
rpt(abdn~ (1 | subj_id), grname = "subj_id", data = DATA_ALL_STROOP_abdn_H, datatype = "Gaussian", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(abdn~  type_picture+  (1 | subj_id), grname = "subj_id", data = DATA_ALL_STROOP_abdn_L, datatype = "Gaussian", nboot =1000, npermut =  1000)


   #######


DATA_ALL_STROOP_abdn_H<- subset(DATA_ALL_STROOP_abdn, tolerance=="high")

Modelpoisstonk<- glmer(abdn~ age + sex +experience+type_picture + (1|subj_id), data = DATA_ALL_STROOP_abdn_H,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoisstonk)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_H,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoisstonk, Modelpoiss0)

##
t<-coef(summary(Modelpoissrhes))

write.table(t, file = "tonk.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

##
DATA_ALL_STROOP_abdn_M<- subset(DATA_ALL_STROOP_abdn, tolerance=="medium")

Modelpoissfasci<- glmer(abdn~ age + sex +type_picture + (1|subj_id), data = DATA_ALL_STROOP_abdn_M,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissfasci)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_M,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoissfasci, Modelpoiss0)
##
DATA_ALL_STROOP_abdn_L<- subset(DATA_ALL_STROOP_abdn, tolerance=="low")

Modelpoissrhes<- glmer(abdn~ age + sex +type_picture + (1|subj_id), data = DATA_ALL_STROOP_abdn_L,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissrhes)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_L,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoissrhes, Modelpoiss0)


postHocs<-glht(Modelpoissrhes, linfct=mcp(type_picture="Tukey"))
summary(postHocs)

psych::describeBy(DATA_ALL_STROOP_abdn_L$abdn, group = DATA_ALL_STROOP_abdn_L$type_picture) 
##
##############RANK EFFECT
#FASCI
DATA_ALL_STROOP_abdn_M_F<- subset(DATA_ALL_STROOP_abdn_M, sex=="         F")
Modelpoissfasci<- glmer(abdn~ RANK + (1|subj_id), data = DATA_ALL_STROOP_abdn_M_F,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissfasci)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_M_F,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoissfasci, Modelpoiss0)

#####RHESUS
DATA_ALL_STROOP_abdn_L_F<- subset(DATA_ALL_STROOP_abdn_L, sex=="         F")
Modelpoissrhes<- glmer(abdn~ RANK + (1|subj_id), data = DATA_ALL_STROOP_abdn_L_F,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissrhes)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_L_F,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoissrhes, Modelpoiss0)

##TONK
DATA_ALL_STROOP_abdn_H_M<- subset(DATA_ALL_STROOP_abdn_H, sex=="         M")
Modelpoisstonk<- glmer(abdn~ RANK + (1|subj_id), data = DATA_ALL_STROOP_abdn_H_M,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoisstonk)
Modelpoiss0<- glmer(abdn~ 1 + (1|subj_id), data = DATA_ALL_STROOP_abdn_H_M,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoiss)
anova(Modelpoisstonk, Modelpoiss0)



  ###BACKWARD SELECTION
Modelpoisssex<- glmer(abdn~ age + sex + tolerance+ type_picture +(1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoisssex)
Modelpoissage<- glmer(abdn~ age +  tolerance+ type_picture +(1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissage)
Modelpoissopti<- glmer(abdn~  tolerance+ type_picture +(1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(Modelpoissopti)

Modelpoisstol<- glmer(abdn~ type_picture +(1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelpoisspic<- glmer(abdn~  tolerance+ (1|subj_id), data = DATA_ALL_STROOP_abdn,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
anova (Modelpoissopti,Modelpoisstol)

anova (Modelpoissopti,Modelpoisspic)




########RANK EFFECT
DATA_ALL_STROOP_rank_tonk<-subset(DATA_ALL_STROOP_abdn, tolerance =="high" & sex=="         F")
DATA_ALL_STROOP_rank_fasci<-subset(DATA_ALL_STROOP_abdn, tolerance =="medium" & sex=="         F")
DATA_ALL_STROOP_rank_rhesus<-subset(DATA_ALL_STROOP_abdn, tolerance =="low" & sex=="         F")
 
#tonk
Modelpoiss<- glmer(abdn~ age + RANK + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_tonk,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelpoiss0<- glmer(abdn~ age + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_tonk,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

anova(Modelpoiss,Modelpoiss0)
#fasci
Modelpoiss<- glmer(abdn~ age + RANK + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_fasci,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelpoiss0<- glmer(abdn~ age + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_fasci,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

#rhesus
Modelpoiss<- glmer(abdn~ age + RANK + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_rhesus,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
Modelpoiss0<- glmer(abdn~ age + type_picture +(1|subj_id), data = DATA_ALL_STROOP_rank_rhesus,family=poisson,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(Modelpoiss)
#########



postHocs<-glht(Modelpoiss, linfct=mcp(tolerance="Tukey"))
summary(postHocs)

postHocs<-glht(Modelpoiss, linfct=mcp(type_picture="Tukey"))
summary(postHocs)
###########

t<-coef(summary(Modelpoiss))

write.table(t, file = "abd.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#################
   


#############


p <- DATA_ALL_STROOP_test%>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  mutate(type_picture = fct_relevel(type_picture, 
                                    "None", "Object","Neutral", "Threat"
  )) %>%
  
  ggplot(aes(x=tolerance, fill=type_picture)) + geom_bar()  + 
  labs(x = "\nDegree of Tolerance", y = "Number of time a response was witheld \n") +scale_x_discrete(labels=c("Low", "Medium", "High"))+labs(fill="Type of picture")+
  theme_classic(base_size = 22)

p +  scale_fill_brewer(palette="Set2")
##


##


################################################################################################################
##############################################################################################################
citation("ggplot2")

##AVERAGE NONE

DATAAVER<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/strooprhesus2.csv", header=TRUE)
summary(DATAAVER)
DATAAVER<-subset(DATAAVER,react_time>200 & react_time<35000)
DATAAVERNONE<-subset(DATAAVER,DATAAVER$stim=='None')
summary(DATAAVERNONE)
psych::describeBy(DATAAVERNONE$react_time, group = DATAAVERNONE$subj_id)



##REACTIONS

library(pscl)
library(emmeans)
data<- matrix(c(0,0,0,0,0,0,6,4,0,18,4,0),ncol=3,byrow=TRUE)
 colnames(data) <- c("Low","Medium","High")
 rownames(data) <- c("None","Object","Neutral", "Threat")


DATAREACTION<-read.csv("C:/Users/Marine/Desktop/DATA_STRAS/react.csv", header=TRUE)

DATAREACTION$tolerance <- factor(DATAREACTION$tolerance, levels = c("low", "medium","high"))
summary(DATAREACTION)

.none_vs_pic<-c(1, -3, 1, 1)
.object_vs_face<-c(1, 0, -2, 1)
.threat_vs_neutral<-c(1, 0, 0, -1)
m<-cbind(.none_vs_pic,.object_vs_face,.threat_vs_neutral)


contrasts(DATAREACTION$type_pict)<-m
##
##RFEPEATABILITY
rpt(Reaction~ (1 | subj_id), grname = "subj_id", data = DATAREACTION, datatype = "Gaussian", nboot =1000, npermut =  1000)

#ADJUSTED R
rpt(Reaction~  tolerance+  (1 | subj_id), grname = "subj_id", data = DATAREACTION, datatype = "Gaussian", nboot =1000, npermut =  1000)

########FOR TOLERANCE
model_tol<- zeroinfl(Reaction ~ type_pict + tolerance ,
                 data =DATAREACTION,
                 dist = "poisson")
summary(model_tol)

mnull <- update(model_tol, . ~ tolerance)
summary(mnull)
pchisq(2 * (logLik(model_tol) - logLik(mnull)), df = 3, lower.tail = FALSE)

library(rcompanion)
nagelkerke(model_tol)

Anova(model_tol,
      type="II",
      test="Chisq")
summary(model_tol)


marginal = emmeans(model_tol,
                   ~ tolerance)

pairs(marginal,
      adjust="tukey")

cld(marginal,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")

#########FOR TYPE OF PICTURES

model_pict<- zeroinfl(Reaction ~ type_pict,
                     data =DATAREACTION,
                     dist = "poisson")


summary(model_pict)
nagelkerke(model_pict)

Anova(model_pict,
      type="II",
      test="Chisq")
summary(model_pict)


marginal = emmeans(model_pict,
                   ~ type_pict)

pairs(marginal,
      adjust="tukey")

cld(marginal,
    alpha=0.05,
    Letters=letters,  ### Use lower-case letters for .group
    adjust="tukey")

#############
ici
summary(DATAREACTION)

p <- DATAREACTION %>%
  mutate(tolerance = fct_relevel(tolerance, 
                                 "low", "medium", "high"
  )) %>%
  
  mutate(type_pict = fct_relevel(type_pict, 
                                 "None", "Object", "Neutral","Threat"
  )) %>%
  
ggplot( mapping = aes(x=tolerance, y=Reaction, fill=type_pict)) +
  scale_fill_viridis(discrete = T, option="plasma")+
  geom_bar(stat="identity", position = position_stack(reverse=TRUE)) + 
  guides(fill = guide_legend(reverse=TRUE))+


  labs(x = "\nDegree of Tolerance", y = "Number of visible emotional reactions \n") + labs(fill="Type of picture")+
  theme_classic(base_size = 22)+
  scale_x_discrete(labels=c("Low", "Medium", "High")) +
  theme_classic(base_size = 22)
  
p +
  theme(axis.title.x = element_text(vjust=5),
        axis.title.y = element_text(vjust=0.01))+
theme(legend.position =  c(0.9, 0.7),  legend.title = element_text(size=12),   legend.text = element_text(size=11))
##


###
p <- ggplot(DATASTROOPS, aes(x=sex, y=control_Stroop, fill=sex)) + 
  geom_boxplot()+  
  
  geom_hline(yintercept=20, linetype="dashed", color = "black")+
  geom_jitter(color="black", size=0.4, alpha=0.3) +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=6,show_guide = FALSE)+
  
  labs(x = "\nSex of the subjects", y = "Distraction control score\n") +
  theme_classic(base_size = 17) 
p
###