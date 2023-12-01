# Lipshutz et al. Phylotranscriptomics reveals the convergent evolution of aggression is associated with 
# both shared and unique patterns of gene expression evolution in cavity-nesting songbirds

# Behavioral analyses and visualizations

# Load packages
library(ggplot2)

library(lattice)
library(MASS)
require(pscl) # alternatively can use package ZIM for zero-inflated models

library(lmtest)
library(lme4)
library(lmerTest)

require (officer)
require(rvg)
require(ggpubr)

# Load in behavioral data - you'll have to specify your own path
ten.species <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Analyses/Cavity.Agg.Collect.Allyears.csv", header = TRUE)

# Subset by sex
F.agg <- subset(ten.species, Sex == "F")
M.agg <- subset(ten.species, Sex == "M")

# Do cavity-nesting species have higher aggression than their close relatives that don't nest in cavities?

# Aggression - attack rate - both sexes

# Plot by species - both sexes
attack <- ggboxplot(ten.species, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1) + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                          "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
attack

# Export plot as ppt for Figure 1A
editable_graph <- dml(ggobj = attack)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "attack.pptx")

# Plot by nest strategy - both sexes
attack.nest <- ggboxplot(ten.species, x = "Nest.Type.Specific", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1) + scale_x_discrete(limits=c("Open", "Facultative Cavity", "Obligate.Cavity")) + theme_classic()
attack.nest
attack.nest_narrow <- attack.nest + theme(aspect.ratio = 2) 
attack.nest_narrow

# Export plot as ppt for Figure 1B
editable_graph <- dml(ggobj = attack.nest_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "attack.nest_narrow.pptx")


# Attack rate - females only
# Plot by species - females only
ggboxplot(F.agg, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1) + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                     "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()

# Plot by nest strategy - females only
ggboxplot(ten.species, x = "Nest.Type.Specific", y = "Attack.Rate", color = "Nest.Type.Specific", add = "jitter") + ylim(0,1) + theme_classic()

# Attack rate - males only
ggboxplot(M.agg, x = "Species", y = "Attack.Rate", color = "Nest.Type.Specific", add = "jitter") + ylim(0,1) + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                         "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
# Plot by nest strategy - females only
ggboxplot(M.agg, x = "Nest.Type.Specific", y = "Attack.Rate", color = "Nest.Type.Specific", add = "jitter") + ylim(0,1) + theme_classic()



# Distance from mount by species - both sexes
distance <- ggboxplot(ten.species, x = "Species", y = "Distance..raw.epochs.", color = "Sex", add = "jitter")
distance + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                               "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
distance

editable_graph <- dml(ggobj = distance)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "distance.pptx")


# Plot by distance from mount by nest strategy - both sexes
distance.nest <- ggboxplot(ten.species, x = "Nest.Type.Specific", y = "Distance..raw.epochs.", color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Open", "Facultative Cavity", "Obligate.Cavity")) + theme_classic()
distance.nest
distance.nest_narrow <- distance.nest + theme(aspect.ratio = 2) 
distance.nest_narrow

editable_graph <- dml(ggobj = distance.nest_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "distance.nest_narrow.pptx")


###################################################################################################################################

# Is aggression biased towards same-sex intruders? - Figure S1A

ggboxplot(ten.species, x = "Mount.Sex", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1) + theme_classic()
agg.sex <- ggboxplot(ten.species, x = "Mount.Sex", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1) + theme_classic()
agg.sex_narrow <- agg.sex + theme(aspect.ratio = 2) 
agg.sex_narrow

editable_graph <- dml(ggobj = agg.sex_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "agg.sex_narrow.pptx")


# Preliminary analyses prior to PGLM
aggression.sex.lmer <- lmer(Attack.Rate ~ Mount.Sex + Sex + Mount.Sex:Sex + (1|Family), data = ten.species)
anova(aggression.sex.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#Mount.Sex     0.12951 0.12951     1 282.32  3.1321 0.077842 . 
#Sex           0.23416 0.23416     1 282.27  5.6632 0.017989 * 
#Mount.Sex:Sex 0.38587 0.38587     1 282.31  9.3323 0.002467 **


# Is distance biased towards same-sex intruders? - Figure S1B
ggboxplot(ten.species, x = "Mount.Sex", y = "Distance..raw.epochs.", color = "Sex", add = "jitter") + theme_classic()
distance.sex <- ggboxplot(ten.species, x = "Mount.Sex", y = "Distance..raw.epochs.", color = "Sex", add = "jitter") + theme_classic()
distance.sex_narrow <- distance.sex + theme(aspect.ratio = 2) 
distance.sex_narrow

editable_graph <- dml(ggobj = distance.sex_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "distance.sex_narrow.pptx")


# Preliminary analyses prior to PGLM
distance.sex.lmer <- lmer(Distance..raw.epochs. ~ Mount.Sex + Sex + Mount.Sex:Sex + (1|Family), data = ten.species)
anova(distance.sex.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#Mount.Sex      0.470   0.470     1 283.35  0.0159 0.8997
#Sex           70.859  70.859     1 283.28  2.3994 0.1225
#Mount.Sex:Sex 18.718  18.718     1 283.22  0.6338 0.4266


######################################################################################################################

# Preliminary analyses of aggression prior to PGLM

# Linear mixed model for attack rate across nest strategies

Attack.lmer <- lmer(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex + (1|Family), data = ten.species)
anova(Attack.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Nest.Type.Specific     0.90854 0.45427     2 240.21 11.6729 1.452e-05 ***
#Sex                    0.32661 0.32661     1 290.27  8.3926  0.004054 ** 
#Nest.Type.Specific:Sex 0.14992 0.07496     2 290.23  1.9261  0.147568  

# Zero inflated Poisson
library(GLMMadaptive)
fm1 <- mixed_model(Physical.Contact ~ Nest.Type.Specific * Sex + Nest.Type.Specific:Sex, random = ~ 1 | Family, data = ten.species,
                   family = zi.poisson(), zi_fixed = ~ Nest.Type.Specific)
summary(fm1)
# Fixed effects:
#   Estimate Std.Err z-value   p-value
# (Intercept)                              0.9731  0.3243  3.0009 0.0026922
# Nest.Type.SpecificObligate.Cavity        1.4312  0.2918  4.9043   < 1e-04
# Nest.Type.SpecificOpen                  -1.9601  0.6633 -2.9549 0.0031282
# SexM                                     1.8401  0.2902  6.3413   < 1e-04
# Nest.Type.SpecificObligate.Cavity:SexM  -1.5870  0.2957 -5.3664   < 1e-04
# Nest.Type.SpecificOpen:SexM              0.7287  0.6711  1.0858 0.2775678

anova(attack.glmer,fm1)
AIC(attack.glmer,fm1)
# df      AIC
# attack.glmer  7 3943.119
# fm1          10 2485.070

library(performance)

r2_nakagawa(attack.glmer, by_group = FALSE, tolerance = 1e-05)
#Conditional R2: 0.929
#Marginal R2: 0.715

fm1 <- mixed_model(Physical.Contact ~ Nest.Type.Specific * Sex + Nest.Type.Specific:Sex, random = ~ 1 | Family, data = ten.species,
                   family = zi.poisson(), zi_fixed = ~ Nest.Type.Specific)
r2_nakagawa(fm1, by_group = FALSE, tolerance = 1e-05)
#Conditional R2: 0.945
#Marginal R2: 0.871

library(emmeans)
emmeans(attack.glmer, list(pairwise ~ Nest.Type.Specific), adjust = "tukey")
# $`pairwise differences of Nest.Type.Specific`
# 1                                    estimate    SE  df z.ratio p.value
# Facultative Cavity - Obligate.Cavity   -0.873 0.125 Inf -6.983  <.0001 
# Facultative Cavity - Open               1.968 0.286 Inf  6.887  <.0001 
# Obligate.Cavity - Open                  2.842 0.257 Inf 11.039  <.0001 

#####################################################################################################################

# Prelminary analyses of distance prior to PGLM

# Distance from mount
bp <- ggboxplot(ten.species, x = "Species", y = "Distance..raw.epochs.", color = "Sex", add = "jitter")
#bp <- ggboxplot(ten.species, x = "Species", y = "Distance..raw.epochs.", fill = "Sex", add = "dotplot")
bp + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                               "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()

Distance.lmer <- lmer(Distance..raw.epochs. ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex + (1|Family), data = ten.species)
anova(Distance.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# Nest.Type.Specific     87.255  43.627     2 220.35  1.5023 0.2249
# Sex                    34.205  34.205     1 291.40  1.1778 0.2787
# Nest.Type.Specific:Sex 26.077  13.038     2 291.34  0.4490 0.6387



##################################################################################################################################


# Subset by families
hirundinidae <- subset(ten.species, Family =="Hirundinidae")
parulidae <- subset(ten.species, Family =="Parulidae")
passeridae <- subset(ten.species, Family =="Passeridae")
turdidae <- subset(ten.species, Family =="Turdidae")
troglodytidae <- subset(ten.species, Family =="Troglodytidae")

ggboxplot(hirundinidae, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter")
agg.swallow.lmer <- lm(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex , data = hirundinidae)
anova (agg.swallow.lmer)
#                     Df  Sum Sq  Mean Sq F value   Pr(>F)   
#   Nest.Type.Specific  1 0.19664 0.196644  8.5413 0.005030 **
#   Sex                 1 0.17793 0.177934  7.7287 0.007426 **
#   Sex:Nest.Type       1 0.10512 0.105125  4.5662 0.037073 * 
#   Residuals          55 1.26624 0.023023             

ggboxplot(parulidae, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter")
agg.warbler.lmer <- lm(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex , data = parulidae)
anova (agg.warbler.lmer)
#                   Df Sum Sq Mean Sq F value    Pr(>F)    
# Nest.Type.Specific  1 0.9064 0.90638 14.1719 0.0003966 ***
# Sex                 1 0.1081 0.10810  1.6903 0.1987967    
# Sex:Nest.Type       1 0.0011 0.00109  0.0170 0.8967579    
# Residuals          57 3.6455 0.06396                      

ggboxplot(passeridae, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter")
agg.sparrow.lmer <- lm(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex , data = passeridae)
anova (agg.sparrow.lmer)
# Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Nest.Type.Specific      1 0.000782 0.0007816  0.1398 0.71013  
# Sex                     1 0.025762 0.0257619  4.6074 0.03692 *
#   Nest.Type.Specific:Sex  1 0.001059 0.0010594  0.1895 0.66531  
# Residuals              48 0.268385 0.0055914 

ggboxplot(turdidae, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter")
agg.thrush.lmer <- lm(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex , data = turdidae)
anova (agg.thrush.lmer)
# Df  Sum Sq  Mean Sq F value  Pr(>F)  
# Nest.Type.Specific      1 0.09967 0.099670  4.1056 0.04799 *
# Sex                     1 0.06246 0.062459  2.5728 0.11489  
# Nest.Type.Specific:Sex  1 0.03009 0.030087  1.2394 0.27081  
# Residuals              51 1.23809 0.024276 

ggboxplot(troglodytidae, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter")
agg.wren.lmer <- lm(Attack.Rate ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex , data = troglodytidae)
anova (agg.wren.lmer)
# Df Sum Sq Mean Sq F value    Pr(>F)    
# Nest.Type.Specific      1 0.0076 0.00758  0.1363 0.7130928    
# Sex                     1 0.7710 0.77104 13.8737 0.0003958 ***
# Nest.Type.Specific:Sex  1 0.0053 0.00527  0.0947 0.7591572    
# Residuals              69 3.8347 0.05558 


# Sparrows
Sparrows <-subset(ten.species, Family == "Passeridae")

ggboxplot(Sparrows, x = "Species", y = "Attack.Rate", color = "Sex", add = "jitter") + ylim(0,1)
Sparrow.Attack.lm <- lm(Attack.Rate ~ Species + Sex, data = Sparrows)
anova(Sparrow.Attack.lm) # Males attack more than females

ggboxplot(Sparrows, x = "Species", y = "Distance..raw.epochs.", color = "Sex", add = "jitter")
Sparrow.Attack.lm <- lm(Distance..raw.epochs. ~ Species + Sex, data = Sparrows)
anova(Sparrow.Attack.lm)


## Cavity vs facultative for HOSP and CARW
ten.species <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Analyses/Cavity.Agg.Collect.Allyears.csv", header = TRUE)

CARW <- subset(ten.species, Species == "Carolina Wren")

ggboxplot(CARW, x = "Nest.Type.Ind", y = "Attack.Rate", color = "Sex", add = "jitter")
CARW.facultative.Attack.lm <- lm(Attack.Rate ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = CARW)
anova(CARW.facultative.Attack.lm)
# Response: Attack.Rate
# Df  Sum Sq Mean Sq F value   Pr(>F)   
# Nest.Type.Ind      1 0.00003 0.00003  0.0006 0.981227   
# Sex                1 0.46394 0.46394  7.9269 0.008154 **
#   Nest.Type.Ind:Sex  1 0.00457 0.00457  0.0780 0.781760   
# Residuals         33 1.93139 0.05853                    

ggboxplot(CARW, x = "Nest.Type.Ind", y = "Distance..raw.epochs.", color = "Sex", add = "jitter")
CARW.facultative.Distance.lm <- lm(Distance..raw.epochs. ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = CARW)
anova(CARW.facultative.Distance.lm)


HOSP <- subset(ten.species, Species == "House Sparrow")

ggboxplot(HOSP, x = "Nest.Type.Ind", y = "Attack.Rate", color = "Sex", add = "jitter")
HOSP.facultative.Attack.lm <- lm(Attack.Rate ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = HOSP)
anova(HOSP.facultative.Attack.lm)
# Response: Attack.Rate
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Nest.Type.Ind      1 0.004762 0.0047623  0.6158 0.4394
# Sex                1 0.019762 0.0197625  2.5555 0.1216
# Nest.Type.Ind:Sex  1 0.001762 0.0017618  0.2278 0.6370
# Residuals         27 0.208803 0.0077334     

ggboxplot(HOSP, x = "Nest.Type.Ind", y = "Distance..raw.epochs.", color = "Sex", add = "jitter")
HOSP.facultative.Distance.lm <- lm(Distance..raw.epochs. ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = HOSP)
anova(HOSP.facultative.Distance.lm)



