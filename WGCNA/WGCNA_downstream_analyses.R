# Lipshutz et al. Phylotranscriptomics reveals the convergent evolution of aggression is associated with 
# both shared and unique patterns of gene expression evolution in cavity-nesting songbirds

# WGCNA analyses and visualizations


# Load in eigengene data - you'll have to specify your own path
eigengenes <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA/WGCNA CN Signed Hybrid/Both sexes/eigengenes_all_ind.csv")

# WGCNA modules - both sexes

# Red module correlated with Sex
red<- ggboxplot(eigengenes, x = "Species", y = "MEred", color = "Sex", width = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",                                                                                                                   "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
red

red_narrow <- red + theme(aspect.ratio = 2) 
red_narrow

editable_graph <- dml(ggobj = red_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "red_narrow.pptx")

# Brown module correlated with Nesting Strategy, Attack Rate, Sex, and Interactions
brown<- ggboxplot(eigengenes, x = "Species", y = "MEbrown", width = 0.5, color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                                    "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
brown
brown_narrow <- brown + theme(aspect.ratio = 2) 
brown_narrow

# looks like this is just a house sparrow difference?
editable_graph <- dml(ggobj = brown_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "brown_narrow.pptx")

# What about brown module and nesting strategy?
brown.nesting <- ggboxplot(eigengenes, x = "Nesting.Strategy", y = "MEbrown", color = "Sex", add = "jitter") + theme_classic()
brown.nesting
editable_graph <- dml(ggobj = brown.nesting)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "brown.nesting.pptx")

# Dark green correlated with sex, attack interactions
darkgreen<- ggboxplot(eigengenes, x = "Species", y = "MEdarkgreen", weight = 0.5, color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",                                                                                                                       "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
darkgreen

darkgreen_narrow <- darkgreen + theme(aspect.ratio = 2) 
darkgreen_narrow

editable_graph <- dml(ggobj = darkgreen_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "darkgreen_narrow.pptx")


# Tan4 correlated with obligate CN
tan4 <- ggboxplot(eigengenes, x = "Species", y = "MEtan4", color = "Sex", width = 0.5, add = "jitter", ylim = c(-0.2, 0.2)) + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
tan4 
tan4_narrow <- tan4 + theme(aspect.ratio = 2) 
tan4_narrow

editable_graph <- dml(ggobj = tan4_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "tan4_narrow_0.5_axis.pptx")

# What about tan4 module and nesting strategy?
tan4.nesting <- ggboxplot(eigengenes, x = "Nesting.Strategy", y = "MEtan4", color = "Sex", add = "jitter") + theme_classic()
tan4.nesting

editable_graph <- dml(ggobj = tan4.nesting)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "tan4.nesting.pptx")

# Honeydew1 correlated with obligate CN
honeydew1<- ggboxplot(eigengenes, x = "Species", y = "MEhoneydew1", color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                               "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
honeydew1
editable_graph <- dml(ggobj = honeydew1)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "honeydew1.pptx")

honeydew1.nesting <- ggboxplot(eigengenes, x = "Nesting.Strategy", y = "MEhoneydew1", color = "Sex", add = "jitter") + theme_classic()
honeydew1.nesting
editable_graph <- dml(ggobj = honeydew1.nesting)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "honeydew1.nesting.pptx")

# Salmon1 correlated with obligate CN
salmon1<- ggboxplot(eigengenes, x = "Species", y = "MEsalmon1", color = "Sex", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                           "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
salmon1
editable_graph <- dml(ggobj = salmon1)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "salmon1.pptx")

salmon1.nesting <- ggboxplot(eigengenes, x = "Nesting.Strategy", y = "MEsalmon1", color = "Sex", add = "jitter") + theme_classic()
salmon1.nesting
editable_graph <- dml(ggobj = salmon1.nesting)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "salmon1.nesting.pptx")

#################################################################################

# WGCNA modules - females

F.eigengenes <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA/WGCNA CN Signed Hybrid/Females/eigengenes_female_plotting.csv")


# Females - Pink4
pink4 <- ggboxplot(F.eigengenes, x = "Species", y = "MEpink4", width = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
pink4 <- ggboxplot(F.eigengenes, x = "Species", y = "MEpink4", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
pink4
pink4_narrow <- pink4 + theme(aspect.ratio = 2) 
pink4_narrow

editable_graph <- dml(ggobj = pink4_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "pink4_narrow.pptx")


pink4.agg <- ggplot(F.eigengenes, aes(x=MEpink4, y=Aggression, col = Species)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pink4.agg
pink4.agg_narrow <- pink4.agg + theme(aspect.ratio = 2) 
pink4.agg_narrow


editable_graph <- dml(ggobj = pink4.agg_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "pink4.agg_narrow.pptx")


# Females - Sienna4
ggboxplot(F.eigengenes, x = "Species", y = "MEsienna4", width = 0.5, add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
sienna4 <- ggboxplot(F.eigengenes, x = "Species", y = "MEsienna4", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
sienna4
sienna4_narrow <- sienna4 + theme(aspect.ratio = 2) 
sienna4_narrow


sienna4.agg <- ggplot(F.eigengenes, aes(x=MEsienna4, y=Aggression, col = Species)) + 
  geom_point(size = 4) + theme_bw() + theme(panel.border = element_blank(), 
                                            panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sienna4.agg
sienna4.agg_narrow <- sienna4.agg + theme(aspect.ratio = 2) 
sienna4.agg_narrow

editable_graph <- dml(ggobj = sienna4.agg_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "sienna4.agg_narrow.pptx")


# Females - turquoise
ggboxplot(F.eigengenes, x = "Species", y = "MEturquoise", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
turquoise <- ggboxplot(F.eigengenes, x = "Species", y = "MEturquoise", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
turquoise_narrow <- turquoise + theme(aspect.ratio = 2) 
turquoise_narrow

editable_graph <- dml(ggobj = turquoise_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "turquoise_narrow.pptx")


# Females - coral3
ggboxplot(F.eigengenes, x = "Species", y = "MEcoral3", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                   "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
coral3 <- ggboxplot(F.eigengenes, x = "Species", y = "MEcoral3", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow", "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
coral3_narrow <- coral3 + theme(aspect.ratio = 2) 
coral3_narrow

editable_graph <- dml(ggobj = coral3_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "coral3_narrow.pptx")


# Females - royal blue
royalblue <- ggboxplot(F.eigengenes, x = "Species", y = "MEroyalblue", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                   "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
royalblue_narrow <- royalblue + theme(aspect.ratio = 2) 
royalblue_narrow

editable_graph <- dml(ggobj = royalblue_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "royalblue_narrow.pptx")


#################################################################################

# WGCNA module - males

setwd("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA/WGCNA CN Signed Hybrid/Males")
M.eigengenes <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA/WGCNA CN Signed Hybrid/Males/eigengenes_male.csv")

# Males - blue
blue <- ggboxplot(M.eigengenes, x = "Species", y = "MEblue", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                         "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
blue_narrow <- blue + theme(aspect.ratio = 2) 
blue_narrow

editable_graph <- dml(ggobj = blue_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "blue_narrow.pptx")


# Males - darkorange2
darkorange2 <- ggboxplot(M.eigengenes, x = "Species", y = "MEdarkorange2", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                       "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
darkorange2
darkorange2_narrow <- darkorange2 + theme(aspect.ratio = 2) 
darkorange2_narrow

editable_graph <- dml(ggobj = darkorange2_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "darkorange2_narrow.pptx")


# Males - purple
purple <- ggboxplot(M.eigengenes, x = "Species", y = "MEpurple", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                             "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
purple
purple_narrow <- purple + theme(aspect.ratio = 2) 
purple_narrow

editable_graph <- dml(ggobj = purple_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "purple_narrow.pptx")


# Males - yellowgreen
yellowgreen <- ggboxplot(M.eigengenes, x = "Species", y = "MEyellowgreen", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                                       "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
yellowgreen
yellowgreen_narrow <- yellowgreen + theme(aspect.ratio = 2) 
yellowgreen_narrow

editable_graph <- dml(ggobj = yellowgreen_narrow)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "yellowgreen_narrow.pptx")


# Males - darkolivegreen
ggboxplot(M.eigengenes, x = "Species", y = "MEdarkolivegreen", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                           "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()
# Males - antiquewhite1
ggboxplot(M.eigengenes, x = "Species", y = "MEantiquewhite1", add = "jitter") + scale_x_discrete(limits=c("Tree Swallow", "Barn Swallow", "Prothonotary Warbler", "Yellow Warbler", "Eurasian Tree Sparrow",
                                                                                                        "House Sparrow","Eastern Bluebird", "American Robin",  "House Wren", "Carolina Wren")) + theme_classic()

###################################################################################################


## Time diff in gene expression - WGCNA eigengenes 
WGCNA.latency <- read.csv("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA_latency.csv", header = TRUE)
Immediate <- subset(WGCNA.latency, Time.Point == "Immediate")


## Was there a difference in sampling types among nest strategies?
mosaicplot(data = WGCNA.latency, Time.Point ~ Nest.Type.Specific)
########## Not sure how to plot this??

# Red eigengene - both sexes
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEred", color = "Sex", add = "jitter")
Time.Point.Red.lmer <- lmer(MEred ~ Time.Point + Sex + Time.Point:Sex + (1|Family), data = WGCNA.latency)
anova(Time.Point.Red.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF  F value  Pr(>F)    
#Time.Point     0.00038 0.00019     2 112.82   0.2529 0.77701    
#Sex            0.72933 0.72933     1 111.82 962.0379 < 2e-16 ***
#Time.Point:Sex 0.00521 0.00260     2 111.69   3.4358 0.03564 *

cor.test(Immediate$MEred, Immediate$Latency.Sample, method = "pearson") 
#t = -1.4333, df = 36, p-value = 0.1604, cor = -0.2323434 


# Tan4 eigengene - both sexes
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEtan4", color = "Sex", add = "jitter")
Time.Point.Tan.lmer <- lmer(MEtan4 ~ Time.Point + Sex + Time.Point:Sex + (1|Family), data = WGCNA.latency)
anova(Time.Point.Tan.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
#Time.Point     0.154169 0.077085     2 105.52 11.5961 2.803e-05 ***
#Sex            0.000035 0.000035     1 112.44  0.0053    0.9420    
#Time.Point:Sex 0.011979 0.005990     2 112.30  0.9010    0.4091 

cor.test(Immediate$MEtan4, Immediate$Latency.Sample, method = "pearson") 
#t = 1.6361, df = 36, p-value = 0.1105, cor = 0.2630741


# MEpink4 eigengene - both sexes
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEpink4", color = "Sex", add = "jitter")
Time.Point.MEpink4.lmer <- lmer(MEpink4 ~ Time.Point + Sex + Time.Point:Sex + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEpink4.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                  Sum Sq   Mean Sq NumDF  DenDF F value Pr(>F)
#Time.Point     0.0138115 0.0069058     2 113.96  1.3381 0.2664
#Sex            0.0000329 0.0000329     1 111.57  0.0064 0.9365
#Time.Point:Sex 0.0164073 0.0082037     2 111.49  1.5896 0.2086

cor.test(Immediate$MEpink4, Immediate$Latency.Sample, method = "pearson") 
#t = 0.93034, df = 36, p-value = 0.3584, cor = 0.1532251


# MEsienna4 eigengene - both sexes
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEsienna4", color = "Sex", add = "jitter")
Time.Point.MEsienna4.lmer <- lmer(MEsienna4 ~ Time.Point + Sex + Time.Point:Sex + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEsienna4.lmer)
#Type III Analysis of Variance Table with Satterthwaite's method
#                  Sum Sq   Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point     0.0166762 0.0083381     2 112.07  3.5811 0.03108 *
#Sex            0.0012949 0.0012949     1 111.14  0.5561 0.45740  
#Time.Point:Sex 0.0029034 0.0014517     2 111.11  0.6235 0.53794  

cor.test(Immediate$MEsienna4, Immediate$Latency.Sample, method = "pearson") 
#t = 0.29925, df = 36, p-value = 0.7665, cor = 0.04981289


# MEroyalblue eigengene - females
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEroyalblue", add = "jitter")
Time.Point.MEroyalblue.lmer <- lmer(MEroyalblue ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEroyalblue.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point 0.057402 0.028701     2 54.379  2.9042 0.06333 

cor.test(Immediate$MEroyalblue, Immediate$Latency.Sample, method = "pearson") 
#t = -1.1368, df = 12, p-value = 0.2778, cor = -0.3117951


# MEturquoise eigengene - females
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEturquoise", add = "jitter")
Time.Point.MEturquoise.lmer <- lmer(MEturquoise ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEturquoise.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point 0.021507 0.010754     2 54.427  0.9133 0.4073

cor.test(Immediate$MEturquoise, Immediate$Latency.Sample, method = "pearson") 
#t = -0.71607, df = 12, p-value = 0.4877, p-value = 0.2778, cor = -0.2024314


# MEcoral3 eigengene - females
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEcoral3", add = "jitter")
Time.Point.MEcoral3.lmer <- lmer(MEcoral3 ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEcoral3.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point 0.014506 0.007253     2 54.456  0.6514 0.5253

cor.test(Immediate$MEcoral3, Immediate$Latency.Sample, method = "pearson") 
#t = -0.58636, df = 12, p-value = 0.5685, cor = -0.1668929


# MEblue eigengene - males
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEblue", add = "jitter")
Time.Point.MEblue.lmer <- lmer(MEblue ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEblue.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point 0.018821 0.0094106     2 56.146  0.9107 0.4081

cor.test(Immediate$MEblue, Immediate$Latency.Sample, method = "pearson") 
#t = -0.36453, df = 22, p-value = 0.7189, cor = -0.07748425 


# MEdarkorange2 eigengene - males
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEdarkorange2", add = "jitter")
Time.Point.MEdarkorange2.lmer <- lmer(MEdarkorange2 ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEdarkorange2.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#Time.Point 0.067817 0.033909     2 56.076  3.5255 0.03612 *

cor.test(Immediate$MEdarkorange2, Immediate$Latency.Sample, method = "pearson") 
#t = 1.5926, df = 22, p-value = 0.1255, cor = 0.3215172 


# MEpurple eigengene - males
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEpurple", add = "jitter")
Time.Point.MEpurple.lmer <- lmer(MEpurple ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEpurple.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F) 
# Time.Point 0.12144 0.060721     2 53.373  6.9501 0.002077 **

cor.test(Immediate$MEpurple, Immediate$Latency.Sample, method = "pearson") 
#t = -0.40433, df = 22, p-value = 0.6899, cor = -0.08588554 


# MEyellowgreen eigengene - males
ggboxplot(WGCNA.latency, x = "Time.Point", y = "MEyellowgreen", add = "jitter")
Time.Point.MEyellowgreen.lmer <- lmer(MEyellowgreen ~ Time.Point + (1|Family), data = WGCNA.latency)
anova(Time.Point.MEyellowgreen.lmer)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F) 
# Time.Point 0.00090925 0.00045463     2 56.134  0.0532 0.9482

cor.test(Immediate$MEyellowgreen, Immediate$Latency.Sample, method = "pearson") 
#t = -1.3625, df = 22, p-value = 0.1868, cor = -0.2789483 




