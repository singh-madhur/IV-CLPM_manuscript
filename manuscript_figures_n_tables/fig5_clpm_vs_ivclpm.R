## Figure: LRT in CLPM vs. IV-CLPM
## Madhur Singh
## Aug 30, 2022

## Contrast the multiple LRTs in IV-CLPM and the lagged effect LRT in CLPM.


# clear workspace
rm(list=ls(all=TRUE))

# SET THE WORKING DIRECTORY
# getwd()

# Libraries
library(tidyverse)
library(RColorBrewer)

# Data --------------------------------------------------------------------------

## IV-CLPM
## From 1_IV-CLPM_varying_time_interval.R
keep <- read.table("resfile_varyingLag_bidir_IV-CLPM.dta")
colnames(keep) <-  c(
  ## Data-generating parameters
  'bxy','byx','bxx','byy','rxy','Tdis',
  ## LRT Statistics
  ## Unidirectional Effect
  'NCP_distal_xy','NCP_prox_xy1','NCP_prox_xy2',
  'NCP_prox12_xy','NCP_dist_prox1_xy','NCP_dist_prox2_xy','NCP_all_xy',
  ## LRTs for the bidirectional effects
  'NCP_distal_bd','NCP_prox1_bd','NCP_prox2_bd',
  'NCP_prox12_bd','NCP_dist_prox1_bd','NCP_dist_prox2_bd','NCP_all_bd',
  ## R-sq - to check that the simulated variables have sensible relationship
  'IVx->X1','IVx->X2','IVy->Y1','IVy->Y2',
  'X1->Y1','X2->Y2','Y1->X1','Y2->X2','X1->Y2','Y1->X2',
  # ## Recovered model estimates
  'b_distal_xy','b_prox1_xy','b_prox2_xy','b_distal_yx','b_prox1_yx','b_prox2_yx',
  ## Standardized path coefficients
  'std_distal_xy','std_prox1_xy','std_prox2_xy','std_distal_yx','std_prox1_yx','std_prox2_yx'
)

glimpse(keep)
psych::describe(keep)

## CLPM
## From 2_CLPM_varying_time_interval.R
keep2 <- read.table("resfile_varyingLag_bidir_CLPM.dta")
colnames(keep2) <-  c(
  ## Data-generating parameters
  'bxy','byx','bxx','byy','rxy','Tdis',
  ## R-sq - to check that the simulated variables have sensible relationship
  'X1->Y1','X2->Y2','Y1->X1','Y2->X2','X1->Y2','Y1->X2',
  ## LRT Statistics
  ## Unidirectional Effect
  'NCP_distal_xy',
  ## LRTs for the bidirectional effects
  'NCP_distal_bd',
  # ## Recovered model estimates
  'b_distal_xy','b_distal_yx',
  ## Standardized path coefficients
  'std_distal_xy','std_distal_yx'
)

glimpse(keep2)
psych::describe(keep2)

# ensure the rows are in the same order in both datasets

keep <- keep |> arrange(bxy, byx, bxx, byy, rxy, Tdis)
keep2 <- keep2 |> arrange(bxy, byx, bxx, byy, rxy, Tdis)

keep$NCP_clpm_lag_xy <- keep2$NCP_distal_xy
keep$NCP_clpm_lag_bd <- keep2$NCP_distal_bd

# Plot -----------------------------------------------------------------------------

clpm_v_ivclpm <- keep |> 
  filter(bxx==.7, byy==.7, bxy==0.2, byx==0.2, rxy==0.3) |>  
  pivot_longer(cols = c('NCP_clpm_lag_xy','NCP_dist_prox2_xy','NCP_prox_xy2','NCP_distal_xy'), 
               names_to = "parameter", names_prefix = "NCP_",
               values_to = "ncp") |> 
  mutate( ref_par = ifelse(parameter %in% c("distal_xy","prox_xy2"),0,1),
          parameter = as_factor(parameter),
          parameter = fct_recode(parameter, "CLPM: Distal Effect (1df)" = "clpm_lag_xy",
                                 "IV-CLPM: Distal + Proximal_W2 (2df)" = "dist_prox2_xy",
                                 "IV-CLPM: Proximal_W2 (1df)" = "prox_xy2",
                                 "IV-CLPM: Distal (1df)" = "distal_xy") ) |> 
  ggplot(aes(x = Tdis, y = ncp, color = parameter)) +
  geom_line(stat = "summary") +
  geom_point(size=1, stat = "summary") +
  coord_cartesian(ylim = c(0,100)) + 
  scale_color_manual(values = c(brewer.pal(11,"RdYlBu")[c(4,1,11,8)])) +
  labs(title = "IV-CLPM vs. CLPM",#
       y = "LRT Statistic", x = "Time Interval", 
       color = "Unidirectional Effects Tested") +
  theme_bw(12) + 
  theme(legend.position = c(0.75,0.87),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))


pdf("fig5_clpm_vs_ivclpm_ncp.pdf", width = 6, height = 6)
clpm_v_ivclpm
dev.off()

