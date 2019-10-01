## analysing CO2 and trait change for manuscript:
## Time travelling seeds reveal that plant regeneration and growth traits are responding to climate change

#read data
seedmassdata <- read.csv("seedmassdata.csv", stringsAsFactors = FALSE)
seeddimensiondata <- read.csv("seeddimensiondata.csv")
heightdata <- read.csv("height.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(heightlog = log10(Heightcm))
stemdata <- read.csv("stemdensity.csv") %>%
  dplyr::rename("ModOld" = "Age") %>%
  dplyr::rename("stemdensity"="Density")
germdata <- read.csv("germinationdata.csv", header = T)
biomassdata <- read.csv("Biomassdata.csv") %>%
  dplyr::mutate(totalbiomasslog = log10(Totalplant)) %>%
  dplyr::mutate(rootshootlog = log10(Roottoshoot))
CO2data <- read.csv("CO2.csv")

#cleaning data
#for seed mass
seedmass1 <- ddply(seedmassdata, c("Species", "ModOld"), summarise,  
                   avg.wt = mean(Mass), massSD = sd(Mass), massnumber = dplyr::n()) ## sometimes the n() function doesn't work if you accidentally reload plyr after dplyr
seedmass2 <- dcast(setDT(seedmass1), Species ~ ModOld, value.var = c("avg.wt", "massSD", "massnumber")) %>%
  mutate(seedmassLRR = log(avg.wt_Modern/avg.wt_Old)) ##the reason why I've kept this calculation here is to double check that the lnRR in metafor is working, which it is
seedmass3 <- escalc(measure="ROM", n1i = massnumber_Modern, n2i = massnumber_Old, m1i = avg.wt_Modern, m2i = avg.wt_Old, sd1i = massSD_Modern, sd2i = massSD_Old, data = seedmass2)
seedmass <- list(seedmass3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(seedmass.p.value = cut(Seedmass,
                                breaks = c(-Inf, 0.05, Inf),
                                labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for seed dimension variance

seeddimensionvariance1 <- ddply(seeddimensiondata, c("Species", "ModOld"), summarise,
                                avg.shape = mean(Variance1), shapeSD = sd(Variance1), shapenumber= n())
seeddimensionvariance2 <- dcast(setDT(seeddimensionvariance1), Species ~ ModOld, value.var = c("avg.shape", "shapeSD", "shapenumber")) %>%
  mutate(seeddimLRR = log(avg.shape_Modern/avg.shape_Old))
seeddimensionvariance3 <- escalc(measure="ROM", n1i = shapenumber_Modern, n2i = shapenumber_Old, m1i = avg.shape_Modern, m2i = avg.shape_Old, sd1i = shapeSD_Modern, sd2i = shapeSD_Old, data=seeddimensionvariance2)
seeddimensionvariance <- list(seeddimensionvariance3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
  purrr::reduce(left_join, by="Species")%>%
  mutate(seeddimension.p.value = cut(Seedshape,
                                     breaks = c(-Inf, 0.05, Inf),
                                     labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2

#for germination success
# NB for seed viability, germination success and seed dormancy I 
#needed to calculate proportions so these look a lot different###

seedgerminationsuccess1 <- germdataraw %>%
  mutate(Species=as.character(Species), ModOld=as.character(ModOld)) %>%
  group_by(Species, ModOld) %>%
  summarise(numbergerm= sum(Germ == 'yes'), germnumber = n())
seedgerminationsuccess2 <- dcast(setDT(seedgerminationsuccess1), Species ~ ModOld, value.var = c("numbergerm", "germnumber"))
seedgerminationsuccess3 <- escalc(measure="OR", ai=numbergerm_Mod, ci=numbergerm_Old, n1i=germnumber_Mod, n2i=germnumber_Old, data=seedgerminationsuccess2)
seedgerminationsuccess <- list(seedgerminationsuccess3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
seedviability <- list(seedviability3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
dormancy <- list(dormancy3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
plantheight <- list(plantheight3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
stemdensity <- list(stemdensity3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
biomass <- list(biomass3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
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
roottoshoot <- list(roottoshoot3, tempdata, precipdata, heatwavedata, growthformdata, aridityandvpddata, traitmodels, CO2data) %>%
  purrr::reduce(left_join, by="Species") %>%
  mutate(roottoshoot.p.value = cut(Roottoshoot,
                                   breaks = c(-Inf, 0.05, Inf),
                                   labels= c("Significant", "Non-signficant"))) ## adding this bit in, in order to colour the graphs in figure 2


CO2seedmass <- lm(seedmass$seedmassLRR ~ seedmass$CO2change)
summary(CO2seedmass)
plot(CO2seedmass)

CO2seedshape <- lm(seeddimensionvariance$seeddimLRR ~ seeddimensionvariance$CO2change)
summary(CO2seedshape)
plot(CO2seedshape)

CO2seedviab <- lm(seedviability$yi ~ seedviability$CO2change)
summary(CO2seedviab)
plot(CO2seedviab)

CO2seeddorm <- lm(dormancy$yi ~ dormancy$CO2change)
summary(CO2seeddorm)
plot(CO2seeddorm)

CO2seedgerm <- lm(seedgerminationsuccess$yi ~ seedgerminationsuccess$CO2change)
summary(CO2seedgerm)
plot(CO2seedgerm)

CO2plantheight <- lm(plantheight$heightLRR ~ plantheight$CO2change)
summary(CO2plantheight)
plot(CO2plantheight)

CO2biomass <- lm(biomass$biomassLRR ~ biomass$CO2change)
summary(CO2biomass)
plot(CO2biomass)

CO2roottoshoot <- lm(roottoshoot$roottoshootLRR ~ roottoshoot$CO2change)
summary(CO2roottoshoot)
plot(CO2roottoshoot)

Co2stemdensity <- lm(stemdensity$stemLRR ~ stemdensity$CO2change)
summary(Co2stemdensity)
plot(Co2stemdensity)

##running model selection meta-analytic methods with CO2 data

seedmass <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                                       Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                                       Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                                      + growthform + CO2change, data=seedmass,
                                      level=1, fitfunction=rma.glmulti, crit="aicc")
print(seedmass)
summary(seedmass@objects[[1]])
mmi <- as.data.frame(coef(seedmass))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

seedshape <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=seedimensionvariance,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(seedshape)
summary(seedshape@objects[[1]])
mmi <- as.data.frame(coef(seedshape))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

seedviab <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=seedviability,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(seedviab)
summary(seedviab@objects[[1]])
mmi <- as.data.frame(coef(seedviab))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

germination <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=seedgerminationsuccess,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(germination)
summary(germination@objects[[1]])
mmi <- as.data.frame(coef(germination))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

dormancy <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=dormancy,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(dormancy)
summary(dormancy@objects[[1]])
mmi <- as.data.frame(coef(dormancy))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

height <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=plantheight,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(height)
summary(height@objects[[1]])
mmi <- as.data.frame(coef(height))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

biomass <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=biomass,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(biomass)
summary(biomass@objects[[1]])
mmi <- as.data.frame(coef(biomass))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

roottoshoot <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=roottoshoot,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(roottoshoot)
summary(roottoshoot@objects[[1]])
mmi <- as.data.frame(coef(roottoshoot))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)

stemdensity <- glmulti(yi ~Changeavtemp + Changetemprange + Changetempvar + 
                      Changeav + Changerange + Changevar + Changemaxseasonal + Changeminseasonal + 
                      Changedrought + Changeheatwave + Changevpd + Changemaxdryspell + Changedrought 
                    + growthform + CO2change, data=stemdensity,
                    level=1, fitfunction=rma.glmulti, crit="aicc")
print(stemdensity)
summary(stemdensity@objects[[1]])
mmi <- as.data.frame(coef(stemdensity))
mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
mmi$z <- mmi$Estimate / mmi$SE
mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]
round(mmi, 4)