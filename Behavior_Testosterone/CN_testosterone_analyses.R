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

T.time <- ggscatter(Immediate, x = "Latency.Sample", y = "logTestosterone", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "log Testosterone", xlab = "latency to sample")
T.time

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


###############################################################################################

# Does attack rate correlate with T across species?
ggplot(ten.species, aes(x=logTestosterone, y=Attack.Rate, col=Species, shape = Sex)) + ylim(0,1) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cor.test(ten.species$logTestosterone,ten.species$Attack.Rate)
#t = 2.1098, df = 62, p-value = 0.03892, r = 0.2588183 

cor.test(F.agg$logTestosterone,F.agg$Attack.Rate) 
#t = 1.873, df = 25, p-value = 0.07281, r = 0.3507881 

cor.test(M.agg$logTestosterone,M.agg$Attack.Rate)
#t = 1.1701, df = 35, p-value = 0.2499, r = 0.1940174 

# Does distance correlate with T across species?

ggplot(ten.species, aes(x=logTestosterone, y=Distance..raw.epochs., col=Species, shape = Sex))  + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

cor.test(ten.species$logTestosterone,ten.species$Distance..raw.epochs.)
# t = 0.049391, df = 63, p-value = 0.9608, r = 0.006222525 

cor.test(F.agg$logTestosterone,F.agg$Distance..raw.epochs.) 
# t = 0.12999, df = 25, p-value = 0.8976, r = 0.02599015 

cor.test(M.agg$logTestosterone,M.agg$Distance..raw.epochs.)
# t = -1.5146, df = 36, p-value = 0.1386, r = -0.244749 


# Looking into each species separately

# Tree swallow
TRES <-subset(ten.species, Species == "Tree Swallow")
ggplot(TRES, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(TRES$logTestosterone,TRES$Attack.Rate, method = 'spearman') 
#S = 420.37, p-value = 0.1233, rho = -0.4698253 

TRES.M <-subset(TRES, Sex == "M")
cor.test(TRES.M$logTestosterone,TRES.M$Attack.Rate, method = 'spearman') 
#S = 27.071, p-value = 0.5594, rho = -0.3535534

TRES.F <-subset(TRES, Sex == "F")
cor.test(TRES.F$logTestosterone,TRES.F$Attack.Rate, method = 'spearman') 
#S = 2, p-value = 0.002778, rho = 0.9642857 

# Barn swallow
BARS <-subset(ten.species, Species == "Barn Swallow")
ggplot(BARS, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(BARS$logTestosterone,BARS$Attack.Rate, method = 'spearman') 
# Not enough observations

# Eastern Bluebird
EABL <-subset(ten.species, Species == "Eastern Bluebird")
ggplot(EABL, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(EABL$logTestosterone,EABL$Attack.Rate, method = 'spearman') 
#S = 23.312, p-value = 0.5177, rho = 0.3339472 

EABL.M <-subset(EABL, Sex == "M")
cor.test(EABL.M$logTestosterone,EABL.M$Attack.Rate, method = 'spearman') 
#S = 6.8377, p-value = 0.6838, rho = 0.3162278 

EABL.F <-subset(EABL, Sex == "F")
cor.test(EABL.F$logTestosterone,EABL.F$Attack.Rate, method = 'spearman') 
#S = 2.2204e-16, p-value = 1, rho = 1

# American Robin
AMRO <-subset(ten.species, Species == "American Robin")
ggplot(AMRO, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(AMRO$logTestosterone,AMRO$Attack.Rate, method = 'spearman') 
#Not enough data


# House Sparrow
HOSP <-subset(ten.species, Species == "House Sparrow")
ggplot(HOSP, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(EABL$logTestosterone,HOSP$Attack.Rate, method = 'spearman') 
# Not enough data

# Eurasian Tree Sparrow
EUTS <-subset(ten.species, Species == "Eurasian Tree Sparrow")
ggplot(EUTS, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(EABL$logTestosterone,EABL$Attack.Rate, method = 'spearman') 
# Not enough data

# House wrens
HOWR <-subset(ten.species, Species == "House Wren")
ggplot(HOWR, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(HOWR$logTestosterone,HOWR$Attack.Rate, method = 'spearman') 
#S = 378.1, p-value = 0.9, rho = -0.03872845 

HOWR.M <-subset(HOWR, Sex == "M")
cor.test(HOWR.M$logTestosterone,HOWR.M$Attack.Rate, method = 'spearman') 
#S = 158.16, p-value = 0.4043, rho = -0.3179944

HOWR.F <-subset(HOWR, Sex == "F")
cor.test(HOWR.F$logTestosterone,HOWR.F$Attack.Rate, method = 'spearman') 
#S = 16, p-value = 0.4167, rho = -0.6 

# Carolina wren
CARW <-subset(ten.species, Species == "Carolina Wren")
ggplot(HOWR, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(CARW$logTestosterone,CARW$Attack.Rate, method = 'spearman') 
#S = 71.588, p-value = 0.02279, rho = 0.6746011 

CARW.M <-subset(CARW, Sex == "M")
cor.test(CARW.M$logTestosterone,CARW.M$Attack.Rate, method = 'spearman') 
#S = 22, p-value = 0.1667, rho = 0.6071429 

CARW.F <-subset(CARW, Sex == "F")
cor.test(CARW.F$logTestosterone,CARW.F$Attack.Rate, method = 'spearman') 
#S = 2.254, p-value = 0.2254, rho = 0.7745967 

# Prothonotary Warbler
PROW <-subset(ten.species, Species == "Prothonotary Warbler")
ggplot(PROW, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(PROW$logTestosterone,PROW$Attack.Rate, method = 'spearman') 
#S = 92, p-value = 0.5517, rho = 0.2333333 

PROW.M <-subset(PROW, Sex == "M")
cor.test(PROW.M$logTestosterone,PROW.M$Attack.Rate, method = 'spearman') 
#S = 18, p-value = 0.95, rho = 0.1 

PROW.F <-subset(PROW, Sex == "F")
cor.test(PROW.F$logTestosterone,PROW.F$Attack.Rate, method = 'spearman') 
#S = 14, p-value = 0.75, rho = -0.4

# Yellow Warbler
YEWA <-subset(ten.species, Species == "Yellow Warbler")
ggplot(YEWA, aes(x=logTestosterone, y=Attack.Rate, col = Sex)) + geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                                                                                           panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cor.test(YEWA$logTestosterone,YEWA$Attack.Rate, method = 'spearman') 
#S = 63.88, p-value = 0.5678, rho = 0.2395253 

YEWA.M <-subset(YEWA, Sex == "M")
cor.test(YEWA.M$logTestosterone,YEWA.M$Attack.Rate, method = 'spearman') 
#S = 36, p-value = 0.1333, rho = -0.8 

YEWA.F <-subset(YEWA, Sex == "F")
cor.test(YEWA.F$logTestosterone,YEWA.F$Attack.Rate, method = 'spearman') 
#S = 2, p-value = 1, rho = 0.5
