# life history analysis - Sadeq, Mills, Beckerman

# libraries
library(tidyverse)
library(readxl)
library(car)
library(gridExtra)
library(broom)

# get the data ----
lh_dat <- read_xlsx("Life history data.RLH_SS.xlsx")
glimpse(lh_dat)

# preliminary graphing ----
ggplot(lh_dat, aes(x = Cu, y = Size_at_Maturity, colour = JuJu))+
  geom_point()+
  facet_grid(JuJu~Antibiotic)

# Data Management and Summarisation ----

# summarise all of the traits to get mean ± se
sumDat <- lh_dat %>% 
  group_by(Cu, JuJu, Antibiotic) %>% 
  summarise(meanSize = mean(Size_at_Maturity, na.rm = TRUE),
            lwr_S = meanSize - sd(Size_at_Maturity, na.rm = TRUE)/sqrt(sum(!is.na(Size_at_Maturity))),
            upr_S = meanSize + sd(Size_at_Maturity, na.rm = TRUE)/sqrt(sum(!is.na(Size_at_Maturity))),
            meanAge = mean(Age_at_Maturity, na.rm = TRUE),
            lwr_A = meanAge - sd(Age_at_Maturity, na.rm = TRUE)/sqrt(sum(!is.na(Age_at_Maturity))),
            upr_A = meanAge + sd(Age_at_Maturity, na.rm = TRUE)/sqrt(sum(!is.na(Age_at_Maturity))),
            meanGR = mean(Growth_Rate, na.rm = TRUE),
            lwr_GR = meanGR - sd(Growth_Rate, na.rm = TRUE)/sqrt(sum(!is.na(Growth_Rate))),
            upr_GR = meanGR + sd(Growth_Rate, na.rm = TRUE)/sqrt(sum(!is.na(Growth_Rate))),
            meanRepro = mean(`Clutch Size`, na.rm = TRUE),
            lwr_Repro = meanRepro - sd(`Clutch Size`, na.rm = TRUE)/sqrt(sum(!is.na(`Clutch Size`))),
            upr_Repro = meanRepro + sd(`Clutch Size`, na.rm = TRUE)/sqrt(sum(!is.na(`Clutch Size`))),
            meanInd = mean(Neckteeth_Induction, na.rm = TRUE),
            lwr_ind = meanInd - sd(Neckteeth_Induction, na.rm = TRUE)/sqrt(sum(!is.na(Neckteeth_Induction))),
            upr_ind = meanInd + sd(Neckteeth_Induction, na.rm = TRUE)/sqrt(sum(!is.na(Neckteeth_Induction))),
            meanLipid = mean(Lipid_Index, na.rm = TRUE),
            lwr_lip = meanLipid - sd(Lipid_Index, na.rm = TRUE)/sqrt(sum(!is.na(Lipid_Index))),
            upr_lip = meanLipid + sd(Lipid_Index, na.rm = TRUE)/sqrt(sum(!is.na(Lipid_Index))))

# adjust factor levels for plotting            
sumDat <- sumDat %>% data.frame() %>% 
  mutate(JuJu = case_when(
    JuJu == "No" ~ "Control",
    JuJu == "JuJu" ~ "Kairomone")) %>% 
  mutate(JuJu = factor(JuJu, levels = c("Control", "Kairomone"))) %>% 
  mutate(Antibiotic = factor(Antibiotic, 
                             levels = c("Control", "Antibiotic"))) %>% 
  rename(Kairomone = JuJu)

sumDat

# make the six plot figure - factorial design ----
pSize <- ggplot(sumDat, aes(x = Cu, y = meanSize, 
                   group = Kairomone,
                   colour = Kairomone,
                   ymin = lwr_S,
                   ymax = upr_S))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Size at Maturity")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  xlab("Copper (µg/l)")+
  scale_x_continuous(breaks = c(0,5))+
  theme_bw()+ggtitle("A")+
  theme(legend.title=element_blank())

pAge <- ggplot(sumDat, aes(x = Cu, y = meanAge, 
                         group = Kairomone,
                         colour = Kairomone,
                         ymin = lwr_A,
                         ymax = upr_A))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Age at Maturity")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  xlab("Copper (µg/l)")+
  scale_x_continuous(breaks = c(0,5))+
  theme_bw()+ggtitle("B")+
  theme(legend.title=element_blank())

pGR <- ggplot(sumDat, aes(x = Cu, y = meanGR, 
                         group = Kairomone,
                         colour = Kairomone,
                         ymin = lwr_GR,
                         ymax = upr_GR))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Growth Rate")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  xlab("Copper (µg/l)")+
  scale_x_continuous(breaks = c(0,5))+
  theme_bw()+ggtitle("C")+
  theme(legend.title=element_blank())

pRepro <- ggplot(sumDat, aes(x = Cu, y = meanRepro, 
                         group = Kairomone,
                         colour = Kairomone,
                         ymin = lwr_Repro,
                         ymax = upr_Repro))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Reproduction")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  xlab("Copper (µg/l)")+
  scale_x_continuous(breaks = c(0,5))+
  theme_bw()+ggtitle("D")+
  theme(legend.title=element_blank())

pInd <- ggplot(sumDat, aes(x = Cu, y = meanInd, 
                         group = Kairomone,
                         colour = Kairomone,
                         ymin = lwr_ind,
                         ymax = upr_ind))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Induction")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  xlab("Copper (µg/l)")+
  scale_x_continuous(breaks = c(0,5))+
  theme_bw()+ggtitle("E")+
  theme(legend.title=element_blank())

pLipid <- ggplot(sumDat, aes(x = Cu, y = meanLipid, 
                         group = Kairomone,
                         colour = Kairomone,
                         ymin = lwr_lip,
                         ymax = upr_lip))+
  geom_point(position = position_dodge(width = 0.9))+
  geom_line(position = position_dodge(width = 0.9))+
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.9))+
  facet_wrap(~Antibiotic)+ylab("Mean Condition")+
  scale_colour_manual(values = c(Control = "black", Kairomone = "Red"))+
  theme_bw()+ggtitle("F")+
  xlab("Copper (µg/ml)")+
  scale_x_continuous(breaks = c(0,5))+
  theme(legend.title=element_blank())

# arrange the plots
grid.arrange(pSize,pAge,pGR,pRepro,pInd,pLipid, nrow = 2)

# Statistical Models -----

# overall MANOVA - three way interaction on phenotype
mv_LH <- lm(cbind(Size_at_Maturity, Age_at_Maturity, Growth_Rate,
                  `Clutch Size`, Neckteeth_Induction, Lipid_Index) ~ Cu*JuJu*Antibiotic, data = lh_dat)
Man_out <- Manova(mv_LH, univariate = TRUE)
Man_out
summary(Man_out, univariate = TRUE, multivariate = FALSE)



# univariate models
mod_size <- lm(Size_at_Maturity ~ Cu*JuJu*Antibiotic, data = lh_dat)
mod_age <- lm(Age_at_Maturity ~ Cu*JuJu*Antibiotic, data = lh_dat)
mod_GR <- lm(Growth_Rate ~ Cu*JuJu*Antibiotic, data = lh_dat)
mod_repro <- lm(`Clutch Size`~ Cu*JuJu*Antibiotic, data = lh_dat)
mod_ind <- lm(Neckteeth_Induction ~ Cu*JuJu*Antibiotic, data = lh_dat)
mod_lipid <- lm(Lipid_Index ~ Cu*JuJu*Antibiotic, data = lh_dat)

a1 <- Anova(mod_size)
a2 <- Anova(mod_age)
a3 <- Anova(mod_GR)
a4 <- Anova(mod_repro)
a5 <- Anova(mod_ind)
a6 <- Anova(mod_lipid)

##  collect anova tables
# out_table <- tidy(bind_rows(a1, a2, a3, a4, a5, a6)) %>% 
#   mutate(term = rep(c("Cu","Predation","Antibiotic", 
#                       "Cu:Predation", "Predation:Antibiotic", "Cu:Antibiotic",
#                       "Cu:Predation:Antibiotic","Residuals"), 6),
#          model = rep(c("size","age","gr","repr","ind","lipid"), each = 8)) %>% 
#   select(model, term, `Sum Sq`, Df, `F value`, `Pr(>F)`) %>% 
#   mutate(sig = case_when(
#     `Pr(>F)` < 0.05 & `Pr(>F)` >= 0.01 ~ "*",
#     `Pr(>F)` < 0.01 & `Pr(>F)` >= 0.001 ~ "**",
#     `Pr(>F)` < 0.001 ~ "***",
#     `Pr(>F)` >= 0.05 ~ "")) %>% 
#   filter(., !grepl("Residuals", term))

out_table <- tidy(bind_rows(a1, a2, a3, a4, a5, a6)) %>% 
  mutate(term = rep(c("Cu","Predation","Antibiotic", 
                      "Cu:Predation", "Predation:Antibiotic", "Cu:Antibiotic",
                      "Cu:Predation:Antibiotic","Residuals"), 6),
         model = rep(c("size","age","gr","repr","ind","lipid"), each = 8)) %>%
  select(model, term, sumsq, df, statistic, p.value) %>% 
  mutate(sig = case_when(
        p.value < 0.05 & p.value >= 0.01 ~ "*",
        p.value < 0.01 & p.value >= 0.001 ~ "**",
        p.value < 0.001 ~ "***",
        p.value >= 0.05 ~ "")) %>% 
  filter(., !grepl("Residuals", term))

out_table
write_csv(out_table, file = "allTrait_Model_univariateANOVA.csv")


# CONTROL/NO ANTIBIOTIC ONLY MODELS ----

mv_LH_II <- lm(cbind(Size_at_Maturity, Age_at_Maturity, Growth_Rate,
                  `Clutch Size`, Neckteeth_Induction, Lipid_Index) ~ Cu*JuJu, 
               data = filter(lh_dat, Antibiotic == "Control"))

Man_out_controls <- Manova(mv_LH_II, univariate = TRUE)
Man_out_controls

# univariate ANOVAs - no antibiotics
mod_size_II <- lm(Size_at_Maturity ~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_size_II)

mod_age_II <- lm(Age_at_Maturity ~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_age_II)

mod_GR_II <- lm(Growth_Rate ~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_GR_II)

mod_repro_II <- lm(`Clutch Size`~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_repro_II)

mod_ind_II <- lm(Neckteeth_Induction ~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_ind_II)

mod_lipid_II <- lm(Lipid_Index ~ Cu*JuJu, data = filter(lh_dat, Antibiotic == "Control"))
Anova(mod_lipid_II)

a1_c <- Anova(mod_size_II)
a2_c <- Anova(mod_age_II)
a3_c <- Anova(mod_GR_II)
a4_c <- Anova(mod_repro_II)
a5_c <- Anova(mod_ind_II)
a6_c <- Anova(mod_lipid_II)

out_table_control <- tidy(bind_rows(a1_c, a2_c, a3_c, a4_c, a5_c, a6_c)) %>% 
  mutate(term = rep(c("Cu","Predation","Cu:Predation","Residuals"), 6),
         model = rep(c("size","age","gr","repr","ind","lipid"), each = 4)) %>%
  select(model, term, sumsq, df, statistic, p.value) %>% 
  mutate(sig = case_when(
    p.value < 0.05 & p.value >= 0.01 ~ "*",
    p.value < 0.01 & p.value >= 0.001 ~ "**",
    p.value < 0.001 ~ "***",
    p.value >= 0.05 ~ "")) %>% 
  filter(., !grepl("Residuals", term))

out_table_control
write_csv(out_table_control, file = "allTrait_Model_univariateANOVA_NoAntibiotic.csv")
