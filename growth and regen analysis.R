### Main body of code for manuscript: Time travelling seeds reveal that plant regeneration and growth traits are responding to climate change ###

#lines 64 - 202 factorial ANOVAs for each trait with p-adjustment
#lines 205 - 343 creating tidy data frame to perform remaining analyses and create figures - called "mydata"
#lines 345 - 378 analysing growth form and average trait change data and plotting figure 2
#lines 380 - 531 creating figure 3
#lines 535 - 539 determining collinearity between climate variables
#lines 540 - 713 running stepwise AICc meta-analytic models on each trait with all climate variables
#lines 716 - 851 creating figure 4

# opening libraries #

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)
library(metafor)
library(lme4)
library(tidyverse)
library(gridExtra)
library(broom)
library(tidyr)
library(purrr)
library(mvabund)
library(Hmisc)
library(phytools)
library(ape)
library(maps)
library(glmnet)
library(selectiveInference)
library(tidytext)
library(data.table)
library(rJava)
library(glmulti)
library(MuMIn)
library(ggthemes)

# reading data #

seedmassdata <- read.csv("seedmassdata.csv", stringsAsFactors = FALSE)
seeddimensiondata <- read.csv("seeddimensiondata.csv")
heightdata <- read.csv("height.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(heightlog = log10(Heightcm))
stemdata <- read.csv("stemdensity.csv") %>% 
  dplyr::rename("ModOld" = "Age") %>%
  dplyr::rename("stemdensity"="Density")
germdata <- read.csv("germinationdata.csv", header = T)
germdataraw <- read.csv("GerminationRaw.csv", header = T, na.strings = c(""))
biomassdata <- read.csv("Biomassdata.csv") %>%
  dplyr::mutate(totalbiomasslog = log10(Totalplant)) %>%
  dplyr::mutate(rootshootlog = log10(Roottoshoot))
tempdata <- read.csv("temperature_data_SuzEveringham.csv")
precipdata <- read.csv("precipitation_data_SuzEveringham.csv")
heatwavedata <- read.csv("heatwave_duration_data_SuzEveringham.csv")
aridityandvpddata <- read.csv("aridityandvpddata.csv", header=T)
growthformdata <- read.csv("growthform.csv")


###Factorial ANOVAs (linear models) for each trait to determine signficance between old and modern populations trait by trait for each species###

## for seed mass

#testing shape of data using hists
seedmasshists <- ggplot(seedmassdata,aes(x=Mass))+geom_histogram()+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

seedmasslm <- seedmassdata %>%
  group_by(Species) %>% 
  do(tidy(lm(Masslog ~ ModOld, data=.))) %>%
  dplyr::rename("Seedmass" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic") %>%
  filter(term=="ModOldOld") %>%
  select(-"term")

## for seed dimension variance

#testing shape of data using hists
seedvarhists <- ggplot(seeddimensiondata,aes(x=Variance1))+geom_histogram(binwidth = 0.001)+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

seeddimensionvarlm <- seeddimensiondata %>%
  group_by(Species) %>%
  do(tidy(lm(Variance1~ModOld, data=.))) %>%
  dplyr::rename("Seedshape" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic")%>%
  filter(term=="ModOldOld") %>%
  select(-"term")

## for seed viability 
## **NB for dormancy and seed viability I had to create a new dataframe with no NAs as model was not working with missing values**
viabanddormdata <- na.omit(germdataraw, cols=c("Viable", "Dormancy"))

#testing shape of data using hists
seedviabhists <- ggplot(viabanddormdata,aes(x=Viable))+geom_histogram(binwidth = 0.001)+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

seedviabilitylm <- viabanddormdata %>%
  group_by(Species) %>%
  do(tidy(glm(Viable~ModOld, family=binomial, data=.))) %>%
  dplyr::rename("Viability" ="p.value") %>%
  select(-"estimate", -"std.error", -"statistic") %>%
  filter(term=="ModOldOld") %>%
  select(-"term")

## for germination success

#testing shape of data using hists
seedviabhists <- ggplot(germdataraw,aes(x=Germ))+geom_histogram(binwidth = 0.001)+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

germsuccesslm <- germdataraw %>%
  group_by(Species) %>%
  do(tidy(glm(Germ~ModOld, family=binomial, data=.))) %>%
  dplyr::rename("GerminationSuccess" ="p.value") %>%
  select(-"estimate", -"std.error", -"statistic") %>%
  filter(term=="ModOldOld") %>%
  select(-"term")

#testing shape of data using hists
seedviabhists <- ggplot(viabanddormdata,aes(x=Dormancy))+geom_histogram(binwidth = 0.001)+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

dormancylm <- viabanddormdata %>%
  group_by(Species) %>%
  do(tidy(glm(Dormancy~ModOld, family=binomial, na.action = "na.exclude", data=.))) %>%
  dplyr::rename("Dormancy" ="p.value") %>%
  select(-"estimate", -"std.error", -"statistic") %>%
  filter(term=="ModOldOld") %>%
  select(-"term")

# for plant height

#testing shape of data using hists
heighthists <- ggplot(heightdata,aes(x=heightlog))+geom_histogram()+facet_wrap(~Species, ncol=5, scales='free')+theme_bw()

plantheightlm <- heightdata %>%
  group_by(Species) %>%
  do(tidy(lm(heightlog~ModOld, data=.))) %>%
  dplyr::rename("Plantheight" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic")%>%
  filter(term=="ModOldOld") %>%
  select(-"term")

# for plant biomass (total)

#testing shape of data using hists
biomasshists <- ggplot(biomassdata,aes(x=totalbiomasslog))+geom_histogram()+facet_wrap(~Species, scales='free')+theme_bw()

totalbiomasslm <- biomassdata %>%
  group_by(Species) %>%
  do(tidy(lm(totalbiomasslog~ModOld, data=.))) %>%
  dplyr::rename("Totalbiomass" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic")%>%
  filter(term=="ModOldOld") %>%
  select(-"term")

# for root to shoot ratio

#testing shape of data using hists
rootshoothists <- ggplot(biomassdata,aes(x=rootshootlog))+geom_histogram()+facet_wrap(~Species, scales='free')+theme_bw()

# step 1= create a new df because the zeros are creating issues for the logged values- remove these and then do rest of same steps
rootshootdata <- biomassdata[-grep("-Inf",biomassdata$rootshootlog),]

# then do all the same usual steps as other regressions
roottoshootlm <- rootshootdata %>%
  group_by(Species) %>%
  do(tidy(lm(rootshootlog~ModOld, data=., na.action=na.exclude))) %>%
  dplyr::rename("Roottoshoot" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic")%>%
  filter(term=="ModOldOld") %>%
  select(-"term")

##for stem density

#testing shape of data using hists
stemdensityhists <- ggplot(stemdata,aes(x=Density))+geom_histogram()+facet_wrap(~Species, scales='free')+theme_bw()

# for stem density
stemdensitylm <- stemdata %>%
  group_by(Species) %>%
  do(tidy(lm(stemdensity~ModOld, data=.))) %>%
  dplyr::rename("Stemdensity" = "p.value") %>%
  select(-"estimate", -"std.error", -"statistic")%>%
  filter(term=="ModOldOld") %>%
  select(-"term")

#now stitching
traitmodels <- list(seedmasslm, germsuccesslm, seedviabilitylm, dormancylm, seeddimensionvarlm, totalbiomasslm, roottoshootlm, plantheightlm, stemdensitylm) %>%
  purrr::reduce(left_join, by= c("Species")) %>%
  as.data.frame()

write.csv(traitmodels, file = "./traitmodels.csv")

## converting from wide dataframe to long to get it ready for p adjustment

pvalues <- data.matrix(traitmodels) #convert to matrix
print(pvalues)

#performing p.adjustment
adjustedpvalues <- p.adjust(pvalues, method="holm")
print(adjustedpvalues)

#creating dataframe from adjusted p value matrix

traitmodelsadjusted <- data.frame(matrix(adjustedpvalues, ncol=10, nrow=43))
write.csv(traitmodelsadjusted, file = "./traitmodelsadjusted.csv")


### Tidying data for regeneration and growth manuscript and creating new data frame with all necessary variables - need this to answer hypotheses 2-4 ###
### For each data frame I am calculating the lnRR and vi using the escalc function, this all doens't work in loops it only really works in long form renaming of data sets

#for seed mass
seedmass1 <- ddply(seedmassdata, c("Species", "ModOld"), summarise,  
      avg.wt = mean(Mass), massSD = sd(Mass), massnumber = n()) ## sometimes the n() function doesn't work if you accidentally reload plyr after dplyr
seedmass2 <- dcast(setDT(seedmass1), Species ~ ModOld, value.var = c("avg.wt", "massSD", "massnumber")) %>%
  mutate(seedmassLRR = log(avg.wt_Modern/avg.wt_Old)) ##the reason why I've kept this calculation here is to double check that the lnRR in metafor is working, which it is
seedmass3 <- escalc(measure="ROM", n1i = massnumber_Modern, n2i = massnumber_Old, m1i = avg.wt_Modern, m2i = avg.wt_Old, sd1i = massSD_Modern, sd2i = massSD_Old, data = seedmass2)
seedmass <- list(seedmass3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(seedmass.p.value = cut(Seedmass,
                       breaks = c(-Inf, 0.05, Inf),
                       labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for seed dimension variance

seeddimensionvariance1 <- ddply(seeddimensiondata, c("Species", "ModOld"), summarise,
                          avg.shape = mean(Variance1), shapeSD = sd(Variance1), shapenumber= n())
seeddimensionvariance2 <- dcast(setDT(seeddimensionvariance1), Species ~ ModOld, value.var = c("avg.shape", "shapeSD", "shapenumber")) %>%
  mutate(seedmassLRR = log(avg.shape_Modern/avg.shape_Old))
seeddimensionvariance3 <- escalc(measure="ROM", n1i = shapenumber_Modern, n2i = shapenumber_Old, m1i = avg.shape_Modern, m2i = avg.shape_Old, sd1i = shapeSD_Modern, sd2i = shapeSD_Old, data=seeddimensionvariance2)
seeddimensionvariance <- list(seeddimensionvariance3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species")%>%
  mutate(seeddimension.p.value = cut(Seedshape,
                       breaks = c(-Inf, 0.05, Inf),
                       labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2
  
#for germination success
# NB for seed viability, germination success and seed dormancy I 
#needed to calculate proportions rather than log odds ratio because it is binary data ###

seedgerminationsuccess1 <- germdataraw %>%
  mutate(Species=as.character(Species), ModOld=as.character(ModOld)) %>%
  group_by(Species, ModOld) %>%
  summarise(numbergerm= sum(Germ == 'yes'), germnumber = n())
seedgerminationsuccess2 <- dcast(setDT(seedgerminationsuccess1), Species ~ ModOld, value.var = c("numbergerm", "germnumber"))
seedgerminationsuccess3 <- escalc(measure="OR", ai=numbergerm_Mod, ci=numbergerm_Old, n1i=germnumber_Mod, n2i=germnumber_Old, data=seedgerminationsuccess2)
seedgerminationsuccess <- list(seedgerminationsuccess3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(germination.p.value = cut(GerminationSuccess,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

# for seed viability

seedviability1 <- germdataraw %>%
  mutate(Species=as.character(Species), ModOld=as.character(ModOld)) %>%
  group_by(Species, ModOld) %>%
  summarise(numberviab= sum(Viable == 'yes'), viablenumber = n())
seedviability2 <- dcast(setDT(seedviability1), Species ~ ModOld, value.var = c("numberviab", "viablenumber"))
seedviability3 <- escalc(measure="OR", ai=numberviab_Mod, ci=numberviab_Old, n1i=viablenumber_Mod, n2i=viablenumber_Old, data=seedviability2)
seedviability <- list(seedviability3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(viability.p.value = cut(Viability,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

# for seed dormancy
dormancy1 <- germdataraw %>%
  mutate(Species=as.character(Species), ModOld=as.character(ModOld)) %>%
  group_by(Species, ModOld) %>%
  summarise(numberdormant= sum(Dormancy == 'yes'), dormancynumber= n())
dormancy2 <- dcast(setDT(dormancy1), Species ~ ModOld, value.var = c("numberdormant", "dormancynumber"))
dormancy3 <- escalc(measure="OR", ai=numberdormant_Mod, ci=numberdormant_Old, n1i=dormancynumber_Mod, n2i=dormancynumber_Old, data=dormancy2)
dormancy <- list(dormancy3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(dormancy.p.value = cut(Dormancy,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for plant height

plantheight1 <- ddply(heightdata, c("Species", "ModOld"), summarise,
                                avg.height = mean(Heightcm), 
                                heightSD = sd(Heightcm), 
                                heightnumber= n())
plantheight2 <- dcast(setDT(plantheight1), Species ~ ModOld, value.var = c("avg.height", "heightSD", "heightnumber")) %>%
  mutate(heightLRR = log(avg.height_Mod/avg.height_Old))
plantheight3 <- escalc(measure="ROM", n1i = heightnumber_Mod, n2i = heightnumber_Old, m1i = avg.height_Mod, m2i = avg.height_Old, sd1i = heightSD_Mod, sd2i = heightSD_Old, data=plantheight2)
plantheight <- list(plantheight3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(height.p.value = cut(Plantheight,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for stem density

stemdensity1 <- ddply(stemdata, c("Species", "ModOld"), summarise, 
                                  avg.stem = mean(stemdensity), 
                                  stemSD = sd(stemdensity), 
                                  stemnumber= n())
stemdensity2 <- dcast(setDT(stemdensity1), Species ~ ModOld, value.var = c("avg.stem", "stemSD", "stemnumber")) %>%
  mutate(stemLRR = log(avg.stem_Mod/avg.stem_Old))
stemdensity3 <- escalc(measure = "ROM", n1i = stemnumber_Mod, n2i = stemnumber_Old, m1i = avg.stem_Mod, m2i = avg.stem_Old, sd1i = stemSD_Mod, sd2i = stemSD_Old, data = stemdensity2)
stemdensity <- list(stemdensity3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by = "Species") %>%
  mutate(stemdensity.p.value = cut(Stemdensity,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

# for total biomass

biomass1 <- ddply(biomassdata, c("Species", "ModOld"), summarise,
                  avg.biomass= mean(Totalplant),
                  biomassSD = sd(Totalplant),
                  biomassnumber = n())
biomass2 <- dcast(setDT(biomass1), Species ~ ModOld, value.var = c("avg.biomass", "biomassSD", "biomassnumber")) %>%
  mutate(biomassLRR = log(avg.biomass_Mod/avg.biomass_Old))
biomass3 <- escalc(measure="ROM", n1i = biomassnumber_Mod, n2i = biomassnumber_Old, m1i = avg.biomass_Mod, m2i = avg.biomass_Old, sd1i = biomassSD_Mod, sd2i = biomassSD_Old, data = biomass2)
biomass <- list(biomass3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(biomass.p.value = cut(Totalbiomass,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for root to shoot ratio

roottoshoot1 <-ddply(biomassdata, c("Species", "ModOld"), summarise,
                     avg.roottoshoot = mean(Roottoshoot),
                     roottoshootSD = sd(Roottoshoot),
                     roottoshootnumber = n())
roottoshoot2 <- dcast(setDT(roottoshoot1), Species ~ ModOld, value.var = c("avg.roottoshoot", "roottoshootSD", "roottoshootnumber")) %>%
  mutate(roottoshootLRR = log(avg.roottoshoot_Mod/avg.roottoshoot_Old))
roottoshoot3 <- escalc(measure = "ROM", n1i = roottoshootnumber_Mod, n2i = roottoshootnumber_Old, m1i = avg.roottoshoot_Mod, m2i = avg.roottoshoot_Old, sd1i= roottoshootSD_Mod, sd2i = roottoshootSD_Old, data=roottoshoot2)
roottoshoot <- list(roottoshoot3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(roottoshoot.p.value = cut(Roottoshoot,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2


## now want to join these dataframes with all the weather data frames to create data file that I will work with for rest of code

mydata <- list(seedmass, seeddimensionvariance, seedviability, dormancy, seedgerminationsuccess, plantheight, biomass, roottoshoot, stemdensity, tempdata, precipdata, heatwavedata, aridityandvpddata, growthformdata) %>%
  purrr::reduce(left_join, by="Species") %>%
  as.data.frame()

write.csv(mydata, file = "./mydata.csv") 


### Analysing growth form against average trait change 

## what I've now done is taken this data frame and in excel
## I have made a column for average trait LRR, this file is now called "average trait "averagetraitchange.csv" in the folder and 
## can be used for growth form vs. average trait change analyses
#read in data
averagetraitchange <- read.csv("averagetraitchange.csv") %>%
  left_join(growthformdata, by="Species") %>% ## join with growth form data
  mutate(absolutetraitchange = abs(Averagetraitchange)) ## add a column for absolute change

## linear model - NB must log because it is not normally distributed
growthformlm <- lm(log(absolutetraitchange) ~ growthform, data=averagetraitchange)
summary(growthformlm)

r.squaredGLMM(growthformlm)

plot(growthformlm)

## plot 
growthformplot <- ggplot(aes(x=growthform, y=log(absolutetraitchange)), data = averagetraitchange) +
  geom_violin(fill="sandybrown") +
  stat_summary(fun.y=mean, geom="point", size=2, color="salmon") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black") +
  theme_classic(base_size = 15) +
  xlab("Growth Form") +
  ylab("Absolute trait change (ln)") +
  annotate("text", x=1, y=1.06, label="a") +
  annotate("text", x=2, y=0.9, label="a") +
  annotate("text", x=3, y=0.41, label="a") +
  annotate("text", x=4, y=0.7, label="a")
tiff("figure2.tiff")
plot(growthformplot)
dev.off()

### Plotting trait change for figure 3

mydatafigure3 <- {mydata[is.na(mydata)] <- 0; mydata} ## had to make NAs zero because that makes the species line up properly on figure 3 
  
mydatafigure3 <- reshape::rename(mydatafigure3, c(yi= "stemyi", yi.x = "seedmassyi", yi.y= "seedshapeyi", yi.x.x = "viableyi", yi.y.y = "dormancyyi", yi.x.x.x = "germinationyi", yi.y.y.y = "heightyi", yi.x.x.x.x = "biomassyi", yi.y.y.y.y = "roottoshootyi",
                                                  vi= "stemvi", vi.x = "seedmassvi", vi.y= "seedshapevi", vi.x.x = "viablevi", vi.y.y = "dormancyvi", vi.x.x.x = "germinationvi", vi.y.y.y = "heightvi", vi.x.x.x.x = "biomassvi", vi.y.y.y.y = "roottoshootvi"))

#create colour palette for different species for each graph
mycolours <- c("sandybrown", "darkolivegreen") ## creating my own colour vector

seedmassplot <- ggplot(mydatafigure3, aes(x=Species, y=seedmassyi, color=seedmass.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=seedmassyi-sqrt(seedmassvi), ymax=seedmassyi+sqrt(seedmassvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Seed Mass (g)") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank()) +
  coord_flip()
print(seedmassplot)

seeddimensionplot <- ggplot(mydatafigure3, aes(x=Species, y=seedshapeyi, color=seeddimension.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) + 
  geom_errorbar(aes(ymin=seedshapeyi-sqrt(seedshapevi), ymax=seedshapeyi+sqrt(seedshapevi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Seed Shape") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
print(seeddimensionplot)

viabilityplot <- ggplot(mydatafigure3, aes(x=Species, y=viableyi, color=viability.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=viableyi-sqrt(viablevi), ymax=viableyi+sqrt(viablevi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Seed Viability") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
plot(viabilityplot)

germsuccessplot <- ggplot(mydatafigure3, aes(x=Species, y=germinationyi, color=germination.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=germinationyi-sqrt(germinationvi), ymax=germinationyi+sqrt(germinationvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Germination Success") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
print(germsuccessplot)

dormancyplot <- ggplot(mydatafigure3, aes(x=Species, y=dormancyyi, col=dormancy.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=dormancyyi-sqrt(dormancyvi), ymax=dormancyyi+sqrt(dormancyvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Seed Dormancy") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
print(dormancyplot)

plantheightplot <- ggplot(mydatafigure3, aes(x=Species, y=heightyi, col=height.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=heightyi-sqrt(heightvi), ymax=heightyi+sqrt(heightvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Plant Height (cm)") +
  theme_classic(base_size=12) + 
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
print(plantheightplot)

biomassplot <- ggplot(mydatafigure3, aes(x=Species, y=biomassyi, col=biomass.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=biomassyi-sqrt(biomassvi), ymax=biomassyi+sqrt(biomassvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Total Biomass (g)") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()

plot(biomassplot)

roottoshootplot <- ggplot(mydatafigure3, aes(x=Species, y=roottoshootyi, col=roottoshoot.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=roottoshootyi-sqrt(roottoshootvi), ymax=roottoshootyi+sqrt(roottoshootvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Root:Shoot (cm)") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  coord_flip()
plot(roottoshootplot)

stemdensityplot <- ggplot(mydatafigure3, aes(x=Species, y=stemyi, col=stemdensity.p.value)) +
  geom_point() +
  scale_color_manual(values=mycolours) +
  geom_errorbar(aes(ymin=stemyi-sqrt(stemvi), ymax=stemyi+sqrt(stemvi)), width=.2) +
  geom_hline(yintercept=0) +
  ylab("Stem Density (g/cm3)") +
  labs(color= "             Linear Model 
       Significance (<0.05)") +
  theme_classic(base_size=12) +
  theme(axis.line.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "top") +
  coord_flip()
plot(stemdensityplot)

# now stitching all the individual trait plots together
traitchangeplot <- grid.arrange(seedmassplot, seeddimensionplot, germsuccessplot, viabilityplot, dormancyplot,
                                plantheightplot, biomassplot, roottoshootplot, 
                                stemdensityplot, ncol=9, widths= c(2,1,1,1,1,1,1,1,1))
tiff("figure3.tiff")
print(traitchangeplot)
dev.off()


### Determining collinearity between climate variables ###

climatedata <- dplyr::select(seedmass, Changeavtemp, Changetemprange, Changetempvar, Changeav, Changerange, Changevar, Changemaxseasonal, Changeminseasonal, Changedrought, Changeheatwave, Changevpd, Changemaxdryspell)
climate.cor <- cor(climatedata)
print(climate.cor)

## create function for rma.glmulti- rma with yi and vi
rma.glmulti <- function(formula, data, ...) {
  rma(as.formula(paste(deparse(formula))), vi, data=data, method = "ML", ...)
}

##glmulti rma with seed mass
seedmassmodelselection <- glmulti(seedmassLRR ~ Changeavtemp + Changetemprange + Changetempvar + 
                                    Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                    Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=seedmass,
                                  level=1, fitfunction=rma.glmulti, crit="aicc")

print(seedmassmodelselection) ## look at best model
plot(seedmassmodelselection) ## plot model selection  
summary(seedmassmodelselection@objects[[1]]) ##summary of best model selected from model selection
plot(seedmassmodelselection, type="s") #plot of relative importance of various model terms


mmi <- as.data.frame(coef(seedmassmodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)


## with seed dimension

seedshapemodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                     Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                     Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform + CO2change, data=seeddimensionvariance,
                                   level=1, fitfunction=rma.glmulti, crit="aicc")
summary(seedshapemodelselection@objects[[1]])
print(seedshapemodelselection)
plot(seedshapemodelselection)

mmi <- as.data.frame(coef(seedshapemodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## with seed germination success

germinationmodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                       Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                       Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=seedgerminationsuccess,
                                     level=1, fitfunction=rma.glmulti, crit="aicc")

print(germinationmodelselection)
summary(germinationmodelselection@objects[[1]])

mmi <- as.data.frame(coef(germinationmodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## with seed viability
viabilitymodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                     Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                     Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=seedviability,
                                   level=1, fitfunction=rma.glmulti, crit="aicc")

print(viabilitymodelselection)
summary(viabilitymodelselection@objects[[1]])
mmi <- as.data.frame(coef(viabilitymodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## for dormancy model selection

dormancymodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                    Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                    Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=dormancy,
                                  level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize = 8192)

print(dormancymodelselection)

summary(dormancymodelselection@objects[[1]])
mmi <- as.data.frame(coef(dormancymodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## for plant height model selection

plantheightmodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                       Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                       Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=plantheight,
                                     level=1, fitfunction=rma.glmulti, crit="aicc")
print(plantheightmodelselection)
summary(plantheightmodelselection@objects[[1]])
mmi <- as.data.frame(coef(plantheightmodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

## for biomass model selection

biomassmodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                   Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                   Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=biomass,
                                 level=1, fitfunction=rma.glmulti, crit="aicc")
print(biomassmodelselection)
summary(biomassmodelselection@objects[[1]])
mmi <- as.data.frame(coef(biomassmodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)
## stem density model selection

stemdensitymodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                       Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                       Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=stemdensity,
                                     level=1, fitfunction=rma.glmulti, crit="aicc")
print(stemdensitymodelselection)
mmi <- as.data.frame(coef(stemdensitymodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)
## for root to shoot model selection

roottoshootmodelselection <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                       Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                       Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought + growthform, data=roottoshoot,
                                     level=1, fitfunction=rma.glmulti, crit="aicc")
print(roottoshootmodelselection)
summary(roottoshootmodelselection@objects[[1]])
mmi <- as.data.frame(coef(roottoshootmodelselection))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)


### plots for figure 4 (draft)

#creating point size for seed shape
sizeshape <- 1/sqrt(seeddimensionvariance$vi) #create a size for points based on the sampling variance
sizeshape <- sizeshape/max(sizeshape) #proportional to the max sampling variance

## seed shape vs dry spell duration
seedshapedryspell <- ggplot(seeddimensionvariance, aes(Changemaxdryspell, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizeshape*5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(cex=8, x="Change in maximum dry spell (days)", y="Change in seed shape") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(seedshapedryspell)

#creating point size for seed viability
sizeviable <- 1/sqrt(seedviability$vi) #create a size for points based on the sampling variance #proportional to the max sampling variance

#seed viability vs temperature variability
seedviabilitytempvar <- ggplot(seedviability, aes(Changetempvar, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizeviable*1.5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(x="Change temperature variability", y="Change in seed viability") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(panel.border= element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(seedviabilitytempvar)


# seed viability vs max dry spell
seedviabilitydryspell <- ggplot(seedviability, aes(Changemaxdryspell, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizeviable*1.5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(cex=8, x="Change maximum dry spell duration (days)", y="Change in seed viability") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(seedviabilitydryspell)




#creating point size for seed germination graph
sizegerm <- 1/sqrt(seedgerminationsuccess$vi)

# seed germination vs temperature variability

germinationtempvar <- ggplot(seedgerminationsuccess, aes(Changetempvar, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizegerm*1.5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(cex=8, x="Change in temperature variability", y="Change in germination success") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(germinationtempvar)


#creating point size for plant height graph
sizeheight <- 1/sqrt(plantheight$vi) #create a size for points based on the sampling variance
sizeheight <- sizeheight/max(sizeheight) #proportional to the max sampling variance

# plant height vs max heatwave duration
plantheightheatwave <- ggplot(plantheight, aes(Changeheatwave, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizeheight*5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(cex=8, x="Change maximum heatwave duration (days)") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(plantheightheatwave)

#creating size for biomass graph

sizebiomass <-  1/sqrt(biomass$vi) #create a size for points based on the sampling variance
sizebiomass <- sizebiomass/max(sizebiomass) #proportional to the max sampling variance

# Biomass vs max heatwave duration

biomassmaxheatwave <- ggplot(biomass, aes(Changeheatwave, yi)) + 
  geom_point(colour="darkolivegreen4", size=sizebiomass*5) + 
  geom_smooth(method="lm", se=FALSE, colour="sandybrown") +
  labs(cex=8, x="Change max heatwave duration (days)", y="Change in plant biomass (g)") +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  theme(panel.border= element_blank())+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  theme(text = element_text(size=8), axis.text.x = element_text(size=6))
plot(biomassmaxheatwave)


## creating size of points for root to shoot ratios
sizeroottoshoot <- 1/sqrt(roottoshoot$vi)
sizeroottoshoot <- sizeroottoshoot/max(sizeroottoshoot)

figure4 <- grid.arrange(seedshapedryspell, seedviabilitytempvar, seedviabilitydryspell,
                        germinationtempvar, plantheightheatwave, biomassmaxheatwave, ncol=3)
tiff("figure4.tiff")
plot(figure4)
dev.off()