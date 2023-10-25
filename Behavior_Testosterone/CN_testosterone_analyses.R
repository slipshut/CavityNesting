# Lipshutz et al. Phylotranscriptomics reveals the convergent evolution of aggression is associated with 
# both shared and unique patterns of gene expression evolution in cavity-nesting songbirds

# Testosterone analyses and visualizations

# Load packages
library(ggplot2)

library(lattice)
library(MASS)

library(lmtest)
library(emmeans)
library(lme4)
library(lmerTest)

require (officer)
require(rvg)
require(ggpubr)

# Load in testosterone data - you'll have to specify your own path
ten.species <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Analyses/Cavity.Agg.Collect.Allyears.csv", header = TRUE)

# Subset data by sex
F.agg <- subset(ten.species, Sex == "F")
M.agg <- subset(ten.species, Sex == "M")

# Plot log testosterone by species - both sexes
ggboxplot(ten.species, x = "Species", y = "logTestosterone", color = "Sex", size = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                                    "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
testosterone <- ggboxplot(ten.species, x = "Species", y = "logTestosterone", color = "Sex", size = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                                    "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
testosterone

# Export both sexes testosterone plot to powerpoint
editable_graph <- dml(ggobj = testosterone)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "testosterone.pptx")

# Plot log testosterone by nest strategy - both sexes
ggboxplot(ten.species, x = "Nest.Type.Specific", y = "logTestosterone", color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Open", "Facultative Cavity", "Obligate.Cavity")) + theme_classic()


# Plot log testosterone females
ggboxplot(F.agg, x = "Species", y = "logTestosterone", color = "Nest.Type.Specific", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                 "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
F.testosterone <- ggboxplot(F.agg, x = "Species", y = "logTestosterone", color = "Sex", size = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                                "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
F.testosterone

# Plot log testosterone males
ggboxplot(M.agg, x = "Species", y = "logTestosterone", color = "Nest.Type.Specific", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                 "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
M.testosterone <- ggboxplot(M.agg, x = "Species", y = "logTestosterone", color = "Nest.Type.Specific", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                 "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
M.testosterone


#########################################################################################################################################

### Sampling time and approaches - testosterone - females and males (prior to PGLM)

ggboxplot(ten.species, x = "Time.Point", y = "logTestosterone", color = "Sex", add = "jitter")

Time.Point.T.lmer <- lmer(logTestosterone ~ Time.Point + Sex + Time.Point:Sex + (1|Family), data = ten.species)
anova(Time.Point.T.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF  F value Pr(>F)    
#Time.Point      0.025   0.012     2 183.97   0.1055 0.8999    
#Sex            42.307  42.307     1 186.81 359.0413 <2e-16 ***
#Time.Point:Sex  0.269   0.135     2 187.16   1.1420 0.3214   

Immediate <- subset(ten.species, Time.Point == "Immediate")
cor.test(Immediate$logTestosterone, Immediate$Latency.Sample, method = "pearson") 
#t = 0.20486, df = 45, p-value = 0.8386, cor = 0.03052403 



#########################################################################################################################################

# Preliminary analyses of log testosterone (prior to PGLM)

ggboxplot(ten.species, x = "Nest.Type.Specific", y = "logTestosterone", color = "Sex", add = "jitter")

logT.lmer <- lmer(logTestosterone ~ Nest.Type.Specific + Sex + Nest.Type.Specific:Sex + (1|Family), data = ten.species)
anova(logT.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF  F value  Pr(>F)    
# Nest.Type.Specific      0.782   0.391     2 166.06   3.4173 0.03513 *  
# Sex                    48.727  48.727     1 186.39 426.0016 < 2e-16 ***
# Nest.Type.Specific:Sex  0.046   0.023     2 186.41   0.1990 0.81976     


emmeans(logT.lmer, list(pairwise ~ Nest.Type.Specific), adjust = "tukey")
# $`pairwise differences of Nest.Type.Specific`
# 1                                    estimate     SE  df t.ratio p.value
# Facultative Cavity - Obligate.Cavity    0.154 0.0850 171 1.809   0.1695 
# Facultative Cavity - Open               0.257 0.1039 144 2.473   0.0384 
# Obligate.Cavity - Open                  0.103 0.0598 189 1.721   0.2000 

FlogT.lmer <- lmer(logTestosterone ~ Nest.Type.Specific + (1|Family), data = F.agg)
anova(FlogT.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# Nest.Type.Specific 0.24935 0.12467     2 65.363  1.2911 0.2819

MlogT.lmer <- lmer(logTestosterone ~ Nest.Type.Specific + (1|Family), data = M.agg)
anova(MlogT.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# Nest.Type.Specific 0.52615 0.26308     2 88.119  2.0582 0.1338


# 10 species, both sexes
T.aov <- aov(logTestosterone ~ Nest.Type + Sex + Nest.Type:Sex, data = ten.species)
summary(T.aov)
#                Df Sum Sq Mean Sq F value Pr(>F)    
# Nest.Type       1   0.12    0.12   0.875  0.351    
# Sex             1  58.94   58.94 441.243 <2e-16 ***
# Nest.Type:Sex   1   0.02    0.02   0.185  0.668    
# Residuals     192  25.65    0.13    

# 10 species, females only
T.F.aov <- aov(logTestosterone ~ Nest.Type, data = F.agg)
summary(T.F.aov)
#             Df Sum Sq Mean Sq F value Pr(>F)
# Nest.Type    1  0.123  0.1233   1.146  0.287
# Residuals   90  9.690  0.1077 


# Does aggression correlate with T across species?
ggplot(ten.species, aes(x=logTestosterone, y=Attack.Rate, col=Species, shape = Sex)) + ylim(0,1) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cor.test(ten.species$logTestosterone,ten.species$Attack.Rate)
#t = 2.1098, df = 62, p-value = 0.03892, cor = 0.2588183

cor.test(F.agg$logTestosterone,F.agg$Attack.Rate) 
# t = 1.873, df = 25, p-value = 0.07281, cor = 0.3507881  

cor.test(M.agg$logTestosterone,M.agg$Attack.Rate)
# t = 1.1701, df = 35, p-value = 0.2499, cor = 0.1940174 


# Do songs correlate with T across species?
ggplot(ten.species, aes(x=logTestosterone, y=Songs, col=Species, shape = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cor.test(ten.species$logTestosterone,ten.species$Songs) 
#t = 2.61, df = 62, p-value = 0.01134, r = 0.3146373

cor.test(F.agg$logTestosterone,F.agg$Songs) 
# t = 0.6117, df = 25, p-value = 0.5463; r = 0.1214348 

cor.test(M.agg$logTestosterone,M.agg$Songs) 
# t = -0.55445, df = 35, p-value = 0.5828; r = -0.09331084 


# # Do flyovers correlate with T across species?
ggplot(ten.species, aes(x=logTestosterone, y=Flyovers, col=Species, shape = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# House wrens
HOWR <-subset(ten.species, Species == "House Wren")
HOWR.M <-subset(HOWR, Sex == "M")

ggplot(HOWR, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(HOWR.M$logTestosterone,HOWR.M$Attack.Rate) 
#t = -0.50158, df = 7, p-value = 0.6313


ggplot(HOWR, aes(x=logTestosterone, y=Songs, col = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(HOWR.M$logTestosterone,HOWR.M$Songs) 
#t = 0.65673, df = 7, p-value = 0.5323, cor = 0.24

# Carolina wren
CARW <-subset(ten.species, Species == "Carolina Wren")
CARW.M <-subset(CARW, Sex == "M")

ggplot(CARW, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(CARW.M$logTestosterone,CARW.M$Attack.Rate) 
#t = 1.0703, df = 5, p-value = 0.3334; r = 0.4317267


ggplot(CARW, aes(x=logTestosterone, y=Songs, col = Sex)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(CARW.M$logTestosterone,CARW.M$Songs) 
#t = 0.69029, df = 5, p-value = 0.5207, cor = 0.2949701


# Sparrows
Sparrows <-subset(ten.species, Family == "Passeridae")

f.passeridae <- subset(ten.species, Family =="Passeridae" & Sex == "F")
m.passeridae <- subset(ten.species, Family =="Passeridae" & Sex == "M")


ggboxplot(Sparrows, x = "Species", y = "logTestosterone", color = "Sex", add = "jitter")
Sparrow.T.lm <- lm(logTestosterone ~ Species*Sex, data = Sparrows)
anova(Sparrow.T.lm)
#              Df  Sum Sq Mean Sq  F value   Pr(>F)    
# Species      1  0.0227  0.0227   0.3494    0.559    
# Sex          1 14.4361 14.4361 221.7584 4.05e-15 ***
# Species:Sex  1  0.0169  0.0169   0.2600    0.614    
# Residuals   29  1.8878  0.0651  


HOSP <- subset(ten.species, Species == "House Sparrow")

ggboxplot(HOSP, x = "Nest.Type.Ind", y = "logTestosterone", color = "Sex", add = "jitter")
HOSP.facultative.T.lm <- lm(logTestosterone ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = HOSP)
anova(HOSP.facultative.T.lm)
#                  Df Sum Sq Mean Sq  F value    Pr(>F)    
# Nest.Type.Ind      2 0.6536  0.3268   5.4326   0.01928 *  
# Sex                1 7.7393  7.7393 128.6624 4.094e-08 ***
# Nest.Type.Ind:Sex  1 0.0042  0.0042   0.0702   0.79522    
# Residuals         13 0.7820  0.0602  

# Wrens
CARW <- subset(ten.species, Species == "Carolina Wren")

ggboxplot(CARW, x = "Nest.Type.Ind", y = "logTestosterone", color = "Sex", add = "jitter")
CARW.facultative.T.lm <- lm(logTestosterone ~ Nest.Type.Ind + Sex + Nest.Type.Ind:Sex, data = CARW)
anova(CARW.facultative.T.lm)
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# Nest.Type.Ind      1 0.17906 0.17906  5.3261   0.04367 *  
# Sex                1 2.86657 2.86657 85.2641 3.284e-06 ***
# Nest.Type.Ind:Sex  1 0.12675 0.12675  3.7700   0.08086 .  
# Residuals         10 0.33620 0.03362  

