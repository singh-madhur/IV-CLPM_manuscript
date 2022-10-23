## Figure: CLPM vs. IV-CLPM

## Aim: Examine how changes in the data-generating parameters impact
## the effect estimates in IV-CLPM and CLPM.


# clear workspace
rm(list=ls(all=TRUE))

# SET THE WORKING DIRECTORY
# getwd()

# Libraries
library(tidyverse)
library(patchwork)
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

keep2$mod <- factor(rep("CLPM: Distal", nrow(keep2)))


## Plots ------------------------------------------------------------------------


### 1. IV-CLPM by Causal Effect Size --------------------------------------------

ivclpm_bxy <-  keep %>% 
  filter(byx==.2, bxx==.7, byy==.7, rxy==.3) |> 
  pivot_longer(cols = c("b_prox1_xy","b_distal_xy","b_prox2_xy"), 
               names_to = "parameter", names_prefix = "b_",
               values_to = "estimate") %>% 
  mutate(parameter = as_factor(parameter),
         parameter = fct_recode(parameter, "Proximal, Wave 1" = "prox1_xy",
                                "Distal" = "distal_xy",
                                "Proximal, Wave 2" = "prox2_xy"),
         bxy = fct_rev(as_factor(bxy))) %>% 
  ggplot(aes(Tdis, estimate, color = bxy)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~ parameter) +
  scale_color_brewer(palette = "Set2") +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  labs(title = "By Causal Effect Size",
       y = "Estimate", x = "Time Interval", color = "Lagged\nEffect") +
  theme_bw(12) + 
  theme(legend.position = "left")


### 2. IV-CLPM by Predictor AR --------------------------------------------------

ivclpm_bxx <- keep %>% 
  filter(bxy==.2, byx==.2, byy==.7, rxy==.3) |> 
  pivot_longer(cols = c("b_prox1_xy","b_distal_xy","b_prox2_xy"), 
               names_to = "parameter", names_prefix = "b_",
               values_to = "estimate") %>% 
  mutate(parameter = as_factor(parameter),
         parameter = fct_recode(parameter, "Proximal, Wave 1" = "prox1_xy",
                                "Distal" = "distal_xy",
                                "Proximal, Wave 2" = "prox2_xy"),
         bxx = fct_rev(as_factor(bxx))) %>% 
  ggplot(aes(Tdis, estimate, color = bxx)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~ parameter) +
  scale_color_brewer(palette = "Set2") +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  labs(title = "By Predictor Autoregression",
       y = "Estimate", x = "Time Interval", color = "Predictor\nAR1") +
  theme_bw(12) + 
  theme(legend.position = "left")



### 3. IV-CLPM by Outcome AR ----------------------------------------------------

ivclpm_byy <- keep %>% 
  filter(bxy==.2, byx==.2, bxx==.7, rxy==.3) |> 
  pivot_longer(cols = c("b_prox1_xy","b_distal_xy","b_prox2_xy"), 
               names_to = "parameter", names_prefix = "b_",
               values_to = "estimate") %>% 
  mutate(parameter = as_factor(parameter),
         parameter = fct_recode(parameter, "Proximal, Wave 1" = "prox1_xy",
                                "Distal" = "distal_xy",
                                "Proximal, Wave 2" = "prox2_xy"),
         byy = fct_rev(as_factor(byy))) %>% 
  ggplot(aes(Tdis, estimate, color = byy)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~ parameter) +
  scale_color_brewer(palette = "Set2") +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  labs(title = "By Outcome Autoregression",
       y = "Estimate", x = "Time Interval", color = "Outcome\nAR1") +
  theme_bw(12) + 
  theme(legend.position = "left")



### 4. CLPM by Causal Effect Size -----------------------------------------------s

clpm_bxy <- keep2 %>%
  filter(byx==.2, bxx==.7, byy==.7, rxy==.3) |> 
  mutate(bxy = as_factor(bxy),
         bxy = fct_rev(bxy)) %>% 
  ggplot(aes(Tdis, b_distal_xy, color = bxy)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~mod) +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  scale_color_brewer(palette = "Set2") +
  labs(y = NULL, x = "Time Interval") +
  theme_bw(12) + 
  theme(legend.position = "none")

### 5. CLPM by Predictor AR -----------------------------------------------------

clpm_bxx <-  keep2 %>%
  filter(bxy==.2, byx==.2, byy==.7, rxy==.3) |> 
  mutate(bxx = as_factor(bxx),
         bxx = fct_rev(bxx)) %>% 
  ggplot(aes(Tdis, b_distal_xy, color = bxx)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~mod) +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  scale_color_brewer(palette = "Set2") +
  labs(y = NULL, x = "Time Interval") +
  theme_bw(12) + 
  theme(legend.position = "none")

### 6. CLPM by Outcome AR -------------------------------------------------------

clpm_byy <- keep2 %>%
  filter(bxy==.2, byx==.2, bxx==.7, rxy==.3) |> 
  mutate(byy = as_factor(byy),
         byy = fct_rev(byy)) %>% 
  ggplot(aes(Tdis, b_distal_xy, color = byy)) +
  geom_line(stat = "summary") +
  geom_point(size=0.5, stat = "summary") +
  facet_wrap(~mod) +
  coord_cartesian(xlim = c(0,40), ylim = c(0,.75)) +
  scale_color_brewer(palette = "Set2") +
  labs(y = NULL, x = "Time Interval") +
  theme_bw(12) + 
  theme(legend.position = "none")



## Compile ----------------------------------------------------------------------

row1 <- ivclpm_bxy + clpm_bxy + plot_layout(nrow = 1,widths = c(3,1))

row2 <- ivclpm_bxx + clpm_bxx + plot_layout(nrow = 1,widths = c(3,1))

row3 <- ivclpm_byy + clpm_byy + plot_layout(nrow = 1,widths = c(3,1))

pdf("fig3_effect_by_data_parameters.pdf", width = 12, height = 9)
row1 / row2 / row3 +
  plot_annotation(tag_levels = 'A', tag_suffix = ".")
dev.off()


