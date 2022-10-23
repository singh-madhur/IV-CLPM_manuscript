## Figure: CLPM vs. IV-CLPM

## Contrast the lagged effect with the three causal paths in IV-CLPM
## and how they change with varying time intervals


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
## 2_CLPM_varying_time_interval.R
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


## Plots ------------------------------------------------------------------------

# clpm
clpm <-keep2 |> 
  filter(bxx==.7, byy==.7, bxy==0.2, byx==0.2, rxy==0.3) |> 
  ggplot(aes(Tdis, b_distal_xy)) +
  geom_line(stat = "summary", color = brewer.pal(3,"Dark2")[2]) +
  geom_point(size=1, stat = "summary", color = brewer.pal(3,"Dark2")[2]) +
  coord_cartesian(ylim = c(0,max(keep$b_prox1_xy)) 
  ) +
  labs(title = "CLPM: Distal Effect",
       y = "Causal Estimate", x = "Time Interval between Study Waves") +
  theme_bw(12)

# iv-clpm
ivclpm <- keep %>% 
  filter(bxx==.7, byy==.7, bxy==0.2, byx==0.2, rxy==0.3) |>
  pivot_longer(cols = c('b_prox1_xy','b_distal_xy','b_prox2_xy'), 
               names_to = "parameter", names_prefix = "b_",
               values_to = "estimate") %>% 
  mutate(parameter = as_factor(parameter),
         parameter = fct_recode(parameter, "Proximal, Wave 1" = "prox1_xy",
                                "Distal" = "distal_xy",
                                "Proximal, Wave 2" = "prox2_xy"),
         parameter = fct_relevel(parameter, "Proximal, Wave 1")
  ) %>% 
  ggplot(aes(Tdis, estimate, color = parameter)) +
  geom_line(stat = "summary") +
  geom_point(size = 1, stat = "summary") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0,max(keep$b_prox1_xy))
  ) +
  labs(title = "IV-CLPM: Proximal and Distal Effects",
       y = "Causal Estimate", x = "Time Interval between Study Waves", color = "Effect Type") +
  theme_bw(12) + 
  theme(legend.position = c(0.8,0.5),
        legend.text = element_text(size = 8))


## Combine with {patchwork} ----------------------------------------------------

pdf("fig2_clpm_vs_ivclpm.pdf", width = 10, height = 5)
clpm + ivclpm +
  plot_annotation(tag_levels = 'A', tag_suffix = ".") 
dev.off()

