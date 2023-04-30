# R code for multispecies occupancy analysis of lynx, wildcat, and wolf published in Dyck et al., 2022
# Citation: Dyck, M. A., Iosif, R., Promberger–Fürpass, B., & Popescu, V. D. (2022). Dracula’s ménagerie: A multispecies occupancy analysis of lynx, wildcat, and wolf in the Romanian Carpathians. *Ecology and evolution*, 12(5), e8921.
# Original code provided by Marissa Dyck (https://github.com/marissadyck/Multi_occ_LWW) and amended by Olivier Gimenez

# Summary -----------------------------------------------------------------

# This code is a simplified version of the code the authors used for Dyck et al., 2022. 
# Current script contains analyses for the winter season only. 

# Libraries ---------------------------------------------------------------

library(unmarked)
library(tidyverse)

# Source ------------------------------------------------------------------

# source code provided from Ken Kellner to fix package bugs related to the predict function
source('scripts/om_predict_fix.R')

# Data (winter) -------------------------------------------------------------

# read in data

# detection history
winter_spp <- 
  read.csv('data/species_matrix_winter.csv')

# set Trapcode as factor
winter_spp$TrapCode <- 
  as.factor(winter_spp$TrapCode)

# trap effort/observation covaraites
winter_traps <- 
  read.csv('data/trap_effort_winter.csv')

# set Trapcode as factor
winter_traps$TrapCode <-
  as.factor(winter_traps$TrapCode)

# site covariates
winter_sites <- 
  read.csv('data/cams_data_winter.csv')
  
# alter variable structure
winter_sites$TrapCode <- 
  as.factor(winter_sites$TrapCode)

winter_sites$Z <- 
  as.numeric(winter_sites$Z)

winter_sites$Impact <- 
  as.factor(winter_sites$Impact)

winter_sites$CLC2018 <- 
  as.factor(winter_sites$CLC2018)

# Format data (winter) ----------------------------------------------------

#creating and unmarkedFrameOccuMulti
#unmarkedFrameOccuMulti(y, siteCovs=NULL, obsCovs=NULL, mapInfo)

#creating y: A list (optionally a named list) of length S where each element is an MxJ matrix of the detection, non-detection data for one species, where M is the number of sites, J is the maximum number of sampling periods per site, and S is the number of species in the analysis.
y_multi_winter <- 
  list(
    matrix(unlist(winter_spp[2:9]),   ncol = 8,  byrow = F), # Lynx
    matrix(unlist(winter_spp[10:17]), ncol = 8,  byrow = F), # Wildcat
    matrix(unlist(winter_spp[18:25]), ncol = 8,  byrow = F)# Wolf
    ) 

# add species names
names(y_multi_winter) <- c("Lynx", "Wildcat", "Wolf")

# observation covariates/trap effort
obs_covs_winter <- 
  list(winter_traps[, 2:9])

# site covariates
# combining impact at site for detection covariate and calling it 'some'. 
# Not enough samples per each category to keep them separate.
# Combining also simplifies variable reducing potential bias given that some rangers may report "low" while others more specifically report "logging" etc.
some <- c("Isolated buildings", "Logging", "Village", "Low")

# Combine specific CORINE Land Cover types into broader categories (CLC_ag, CLC_open, CLC_forest)
winter_sites <-
  winter_sites %>% 
  mutate(CLC_ag = (CLC211 + CLC222 + CLC242 + CLC243),# proportion agriculture
         CLC_open = (CLC231 + CLC321 + CLC322 + CLC324), # proportion open habitat
         CLC_forest = (CLC311 + CLC312 + CLC313),# proportion CLC
         Impact.2 = case_when(Impact == "None" ~ 0,
                              Impact %in% some ~ 1))

# change Impact.2 to factoor
winter_sites$Impact.2 <- 
  as.factor(winter_sites$Impact.2)

# create new object winter_sites.scaled for scaled variables 
winter_sites.scaled <- 
  winter_sites

# scale variables
winter_sites.scaled$denslocalr <- scale(winter_sites.scaled$denslocalr)
winter_sites.scaled$distlocalr <- scale(winter_sites.scaled$distlocalr)
winter_sites.scaled$distnatlro <- scale(winter_sites.scaled$distnatlro)
winter_sites.scaled$distsettle <- scale (winter_sites.scaled$distsettle)
winter_sites.scaled$diststream <- scale(winter_sites.scaled$diststream)
winter_sites.scaled$TRI5 <- scale(winter_sites.scaled$TRI5)
winter_sites.scaled$Z <- scale(winter_sites.scaled$Z)
winter_sites.scaled$CLC311 <- scale(winter_sites.scaled$CLC311)
winter_sites.scaled$CLC312 <- scale(winter_sites.scaled$CLC312)
winter_sites.scaled$CLC313 <- scale(winter_sites.scaled$CLC313)
winter_sites.scaled$CLC_forest <- scale(winter_sites.scaled$CLC_forest)
winter_sites.scaled$CLC_ag <- scale(winter_sites.scaled$CLC_ag)
winter_sites.scaled$CLC_open <- scale(winter_sites.scaled$CLC_open)

# create site covarioates data for the unmarkedFrameOccuMulti 
winter_site_covs <- 
  data.frame(winter_sites.scaled)

# adding lynx presence/absence as potential covariate for detection of other species
Lynx <- data.frame(winter_spp[2:9]) # pull lynx detection history out of full detection history

# adding wildcat presence/absence as potential covariate for detection of other species
Wildcat <- data.frame(winter_spp[10:17])

# adding wolf presence/absence as potential covariate for detection of other species
Wolves <- data.frame(winter_spp[18:25]) 

# combine observation covariates plus species detection histories
winter_obs_covs_species <- 
  c(list(effort = data.frame(winter_traps[,2:9])),
    list(Lynx = Lynx),
    list(Wildcat = Wildcat),
    list(Wolf = Wolves))

# create unmarkedFrameOccuMulti object for analysis
winter_occ_data <- 
  unmarkedFrameOccuMulti(
    y_multi_winter, 
    siteCovs = winter_site_covs, 
    obsCovs = winter_obs_covs_species)


# Explore data (winter) ---------------------------------------------------

summary(winter_occ_data)
# naive occupancy for each species (sites with at least 1 detection / total sites)
# lynx 0.67
# wildcat 0.34
# wolf 0.51

plot(winter_occ_data)


# Detection function (winter) ---------------------------------------------

winter_LWW_det <- c(
  '~Wolf + Wildcat + diststream', 
  '~Lynx + Wolf+ diststream', 
  '~Lynx +Wildcat+ diststream')

# Final model (winter) ----------------------------------------------------

#lynx:wildcat movement, wolf habitat
winter_LWW_OF_205 <-  c('~denslocalr', '~Z', '~denslocalr',       # marginal occupancy lynx, wildcat, wolf
                        '~TRI5', '~CLC_forest', '~CLC_forest',    # co-occupancy lynx:wildcat, lynx:wolf, wildcat:wolf
                        '~0')                                     # three species occupancy (not tested, set up to 0 here)

# fit model
winter_LWW_205 <- occuMulti(winter_LWW_det, winter_LWW_OF_205, winter_occ_data)
summary(winter_LWW_205)


# Marginal occupancy results (winter) -------------------------------------

# calculate marginal occupancy using the predict function

# lynx
winter_L_mo<- 
  (predict(winter_LWW_205, #the model we want to use to make predictions
           'state', #specifies we want to predict the occupancy state
           species ='Lynx')) #specifies which species

# saving predicted values for each species to plot
winter_L_pred <- 
  mean(winter_L_mo$Predicted)

winter_L_low <- 
  mean(winter_L_mo$lower)

winter_L_up <- 
  mean(winter_L_mo$upper)

# wildcat
winter_Wc_mo <-
  (predict(winter_LWW_205,
           'state',
           species ='Wildcat'))

winter_Wc_pred <- 
  mean(winter_Wc_mo$Predicted)

winter_Wc_low <-
  mean(winter_Wc_mo$lower)

winter_Wc_up <- 
  mean(winter_Wc_mo$upper)

# wolf
winter_W_mo <- 
  (predict(winter_LWW_205,
           'state',
           species='Wolf'))

winter_W_pred <- 
  mean(winter_W_mo$Predicted)

winter_W_low <- 
  mean(winter_W_mo$lower)

winter_W_up <- 
  mean(winter_W_mo$upper)

# save predicted values for all three species and assign variable names to plot together
winter_pred <- 
  c(winter_L_pred , winter_Wc_pred, winter_W_pred)

winter_upper <- 
  c(winter_L_up, winter_Wc_up, winter_W_up)

winter_lower <- 
  c(winter_L_low, winter_Wc_low, winter_W_low)

# create variables names for graph 
winter_species <- 
  c("Lynx", "Wildcat", "Wolf")

winter_season <- 
  c("Winter", "Winter", "Winter")

# combine into data frame
winter_pred.df<-
  data.frame(cbind(winter_species, winter_season, winter_pred, winter_upper, winter_lower))

# change variable structure
winter_pred.df$winter_pred <- 
  as.numeric(winter_pred.df$winter_pred)

winter_pred.df$winter_upper <- 
  as.numeric(winter_pred.df$winter_upper)

winter_pred.df$winter_lower <- 
  as.numeric(winter_pred.df$winter_lower)

ggplot(
  data = winter_pred.df, 
  aes(x = winter_species, 
      y = winter_pred))+ 
  geom_point()+
  geom_errorbar(
    data = winter_pred.df, 
    aes(x = winter_species, 
        ymin = winter_lower, 
        ymax = winter_upper), 
    width = 0.25,
    color = '#003E83')+
  scale_y_continuous(
    breaks = seq(0,1,0.1),
    limits = c(0,1),
    expand = c(0,0))+
  labs(
    title = 'Predicted marginal occupancy for lynx, wildcat, and wolf in Winter',
    x = 'Species',
    y = 'Occupancy')+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5,
                              size = 18),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15))


# Marginal detection results (winter) -------------------------------------

# calculate marginal detection using the predict function

# lynx
winter_L_md <- 
  (predict(winter_LWW_205, # the model we want to use to make predictions
           'det', # specifies we want the predicted detection
           species='Lynx')) # specifies the species

mean(winter_L_md$Predicted)

# wildcat
winter_Wc_md <- 
  (predict(winter_LWW_205,
           'det',
           species ='Wildcat'))

mean(winter_Wc_md $Predicted)

# wolf
winter_W_md <- 
  (predict(winter_LWW_205,
           'det',
           species ='Wolf'))

mean(winter_W_md$Predicted)



# Conditional occupancy results (winter) ----------------------------------

#calculate conditional occupancy results (predicted occupancy of one species conditional on the presence/absence of another species) using the predict function

# Lynx

# lynx | wildcat present
winter_L_co_Wc <- 
  (predict(winter_LWW_205,
           'state',
           species ='Lynx',
           cond ='Wildcat')) # add the conditional species

mean(winter_L_co_Wc$Predicted)

# lynx | wildcat absent
winter_L_noWc <- 
  (predict(winter_LWW_205,
           'state',
           species ='Lynx',
           cond ='-Wildcat'))

mean(winter_L_noWc$Predicted)

# lynx | wolf present
winter_L_co_W <- 
  (predict(winter_LWW_205,
           'state',
           species='Lynx',
           cond='Wolf')) 

mean(winter_L_co_W$Predicted)

# lynx | wolf absent
winter_L_noW <- 
  (predict(winter_LWW_205,
           'state',
           species ='Lynx',
           cond ='-Wolf'))

mean(winter_L_noW$Predicted)

# Wildcat

# wildcat | lynx present
winter_Wc_co_L <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wildcat',
           cond ='Lynx')) 

mean(winter_Wc_co_L$Predicted)

# wildcat | lynx absent
winter_Wc_noL <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wildcat',
           cond ='-Lynx')) 

mean(winter_Wc_noL$Predicted)

# wildcat | wolf present
winter_Wc_co_W <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wildcat',
           cond ='Wolf'))

mean(winter_Wc_co_W$Predicted)

# Wildcat | Wolf absent
winter_Wc_noW <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wildcat',
           cond ='-Wolf')) 

mean(winter_Wc_noW$Predicted)

# Wolf

# wolf | lynx present
winter_W_co_L<- 
  (predict(winter_LWW_205,
           'state',
           species ='Wolf',
           cond ='Lynx'))

mean(winter_W_co_L$Predicted)

# wolf | lynx absent
winter_W_noL <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wolf',
           cond ='-Lynx')) 

mean(winter_W_noL$Predicted)

# wolf | wildcat present
winter_W_co_Wc <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wolf',
           cond ='Wildcat')) 

mean(winter_W_co_Wc$Predicted)

# wolf | wildcat absent
winter_W_noWc <- 
  (predict(winter_LWW_205,
           'state',
           species ='Wolf',
           cond ='-Wildcat'))

mean(winter_W_noWc$Predicted)




# Conditional occupancy graph ---------------------------------------------

# reproduce Figure 4 in paper

# create new data frame holding all variables constant except the one we want to plot using expand.grid function. 
# This new data frame must include all variables in the model including both detection on occupancy variables for all species and species combinations
winter_205_CLC.df <- 
  data.frame(
    expand.grid(
      CLC_forest = seq(
        min(winter_sites.scaled$CLC_forest), 
        max(winter_sites.scaled$CLC_forest), 
        0.01),
      Z = mean(winter_sites.scaled$Z),
      denslocalr = mean(winter_sites.scaled$denslocalr),
      TRI5 = mean(winter_sites.scaled$TRI5),
      diststream = mean(winter_sites.scaled$diststream)))

CLC_forest <-
  matrix(
  seq(from = 0.08505, to = 1.00, length.out = 509),
  ncol = 1)

# create variable for predicted occupancy of each species conditional on every other species presence/absence

# lynx | wolf
winter_Epsi_L_W <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Lynx', 
          cond = 'Wolf', 
          newdata = winter_205_CLC.df)

# lynx | wolf absent
winter_Epsi_L_noW <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Lynx', 
          cond = '-Wolf', 
          newdata= winter_205_CLC.df)

#lynx | wildcat
winter_Epsi_L_Wc <-
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Lynx', 
          cond = 'Wildcat', 
          newdata = winter_205_CLC.df)

# lynx | wildcat absent
winter_Epsi_L_noWc <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Lynx', 
          cond = '-Wildcat', 
          newdata = winter_205_CLC.df)

# wildcat | lynx
winter_Epsi_Wc_L <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wildcat', 
          cond = 'Lynx', 
          newdata = winter_205_CLC.df)

# wildcat | lynx absent
winter_Epsi_Wc_noL <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wildcat', 
          cond = '-Lynx', 
          newdata = winter_205_CLC.df)

# wildcat | wolf
winter_Epsi_Wc_W <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wildcat', 
          cond = 'Wolf', 
          newdata = winter_205_CLC.df)

# wildcat | wolf absent
winter_Epsi_Wc_noW <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wildcat', 
          cond = '-Wolf', 
          newdata = winter_205_CLC.df)

# wolf | lynx
winter_Epsi_W_L <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wolf', 
          cond = 'Lynx', 
          newdata = winter_205_CLC.df)

# wold | lynx absent
winter_Epsi_W_noL <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wolf', 
          cond = '-Lynx', 
          newdata = winter_205_CLC.df)

# wolf | wildcat
winter_Epsi_W_Wc <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wolf', 
          cond = 'Wildcat', 
          newdata = winter_205_CLC.df)

# wolf | wildcat absent
winter_Epsi_W_noWc <- 
  predict(winter_LWW_205, 
          type ="state", 
          species = 'Wolf', 
          cond = '-Wildcat', 
          newdata= winter_205_CLC.df)


# add conditional column for each of the conditional occupancy data frames that specifies present/absent and another for species, as well as column for 1 SE, and column for the unscaled variable CLC_forest

# lynx | wolf
winter_Epsi_L_W <- 
  winter_Epsi_L_W %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Present",
         c.spp = "Wolf")

winter_Epsi_L_W <- 
  cbind(winter_Epsi_L_W, CLC_forest)

# lynx | wolf absent
winter_Epsi_L_noW <- 
  winter_Epsi_L_noW %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Absent",
         c.spp = "Wolf")

winter_Epsi_L_noW <- 
  cbind(winter_Epsi_L_noW, CLC_forest)

# lynx | wildcat
winter_Epsi_L_Wc <- 
  winter_Epsi_L_Wc %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Present",
         c.spp = "Wildcat")

winter_Epsi_L_Wc <- 
  cbind(winter_Epsi_L_Wc, CLC_forest)

# lynx | wildcat absent
winter_Epsi_L_noWc <- 
  winter_Epsi_L_noWc %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Lynx",
         conditional = "Absent",
         c.spp = "Wildcat")

winter_Epsi_L_noWc <- 
  cbind(winter_Epsi_L_noWc, CLC_forest)

# wildcat | lynx
winter_Epsi_Wc_L <- 
  winter_Epsi_Wc_L %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wildcat",
         conditional = "Present",
         c.spp = "Lynx")

winter_Epsi_Wc_L <- 
  cbind(winter_Epsi_Wc_L, CLC_forest)

# wildcat | lynx absent
winter_Epsi_Wc_noL <- 
  winter_Epsi_Wc_noL %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wildcat",
         conditional = "Absent",
         c.spp = "Lynx")

winter_Epsi_Wc_noL <- 
  cbind(winter_Epsi_Wc_noL, CLC_forest)

# wildcat | wolf
winter_Epsi_Wc_W <- 
  winter_Epsi_Wc_W %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wildcat",
         conditional = "Present",
         c.spp = "Wolf")

winter_Epsi_Wc_W <- 
  cbind(winter_Epsi_Wc_W, CLC_forest)

# wildcat | wolf absent
winter_Epsi_Wc_noW <- 
  winter_Epsi_Wc_noW %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wildcat",
         conditional = "Absent",
         c.spp = "Wolf")

winter_Epsi_Wc_noW <- 
  cbind(winter_Epsi_Wc_noW, CLC_forest)

# wolf | lynx
winter_Epsi_W_L <- 
  winter_Epsi_W_L %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Present",
         c.spp = "Lynx")

winter_Epsi_W_L <- 
  cbind(winter_Epsi_W_L, CLC_forest)

# wolf | lynx absent
winter_Epsi_W_noL <- 
  winter_Epsi_W_noL %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Absent",
         c.spp = "Lynx")

winter_Epsi_W_noL <- 
  cbind(winter_Epsi_W_noL, CLC_forest)

# wolf | wildcat
winter_Epsi_W_Wc <- 
  winter_Epsi_W_Wc %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Present",
         c.spp = "Wildcat")

winter_Epsi_W_Wc <- 
  cbind(winter_Epsi_W_Wc, CLC_forest)

# wolf | wildcat absent
winter_Epsi_W_noWc <- 
  winter_Epsi_W_noWc %>%
  mutate(lwr = Predicted-SE,
         upr = Predicted + SE,
         spp = "Wolf",
         conditional = "Absent",
         c.spp = "Wildcat")

winter_Epsi_W_noWc <- 
  cbind(winter_Epsi_W_noWc, CLC_forest)

winter_Epsi_condOcc_all <- 
  rbind(winter_Epsi_L_W, 
        winter_Epsi_L_noW, 
        winter_Epsi_L_Wc, 
        winter_Epsi_L_noWc,
        winter_Epsi_W_L, 
        winter_Epsi_W_noL, 
        winter_Epsi_W_Wc, 
        winter_Epsi_W_noWc,
        winter_Epsi_Wc_L, 
        winter_Epsi_Wc_noL, 
        winter_Epsi_Wc_W, 
        winter_Epsi_Wc_noW)

#View(winter_Epsi_condOcc_all)


# graph

# labels and color vectors
colors <- 
  c("plum4", "lightseagreen")

label <- 
  c("Absent", "Present")

lines <- 
  c("dotted", "solid")

occ.labels <- 
  c("Probability of lynx", 
     "Probability of wildcat", 
     "Probability of wolf")

names(occ.labels) <- 
  c("Lynx", 
    "Wildcat", 
    "Wolf")

cond.labels <- 
  c("Conditional on lynx", 
    "Conditional on wildcat", 
    "Conditional on wolf")

names(cond.labels) <- 
  c("Lynx", 
    "Wildcat", 
    "Wolf")

winter_condOcc.plot <- 
  ggplot(
    data = winter_Epsi_condOcc_all, 
    aes(x = CLC_forest, 
        y = Predicted, 
        group = conditional))+
  geom_ribbon(
    aes(ymin = lwr, 
        ymax = upr, 
        fill = conditional), 
    alpha = 0.7)+
  geom_line(
    aes(x = CLC_forest, 
        y = Predicted, 
        linetype = conditional), 
    linewidth = 0.5, 
    color = "black")+
  labs(x = expression ("Proportion forest "~km^2),
       y = "Occupancy Probability")+
  facet_grid(c.spp ~ spp, 
             labeller = labeller(c.spp = cond.labels, spp = occ.labels))+
  scale_linetype_manual(values = lines, labels=label)+
  scale_fill_manual(values = colors, labels = label)+
  coord_cartesian(ylim = c(0,1.1))+
  scale_y_continuous(breaks = seq(from = 0, to =1, by = 0.25))+
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 15)); winter_condOcc.plot

# columns are spp it is predicting occupancy for, and rows are the species the predictions are conditional on

