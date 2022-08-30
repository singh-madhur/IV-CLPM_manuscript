## Figure: CLPM vs. IV-CLPM
## Madhur Singh
## Aug 30, 2022

## Aim: Examine how changes in the residual correlation in simulated data impact
## the effect estimates in IV-CLPM and CLPM.
## the model-estimated residual correlation in either model.


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
## From 3_corr_residuals_ivclpm.R
keep <- read.table("resfile_rescor_IV-CLPM.dta")
colnames(keep) <-  c(
  ## Data-generating parameters
  'rxy','Tdis',
  # ## Recovered model estimates
  'b_distal_xy','b_prox1_xy','b_prox2_xy',
  'r_xy1','r_xy2'
)

psych::describe(keep)

## CLPM
## From 4_corr_residuals_clpm.R
keep2 <- read.table("resfile_rescor_CLPM.dta")
colnames(keep2) <-  c(
  ## Data-generating parameters
  'rxy','Tdis',
  # ## Recovered model estimates
  'b_distal_xy','r_xy1','r_xy2'
)

psych::describe(keep2)



## Plots ------------------------------------------------------------------------

# facet panel labels
rxy_labs = c("r_exy = 0.1","r_exy = 0.3")
names(rxy_labs) = c("0.1","0.3")

### 1. CLPM ---------------------------------------------------------------------

clpm <-
  keep2 %>% 
  pivot_longer(cols = c('b_distal_xy','r_xy1','r_xy2'), 
               names_to = "parameter",
               values_to = "estimate") %>% 
  mutate(`Parameter Type` = as_factor(ifelse(parameter == "b_distal_xy",
                                             "Regression", "Correlation")),
         Parameter = as_factor(parameter),
         Parameter = fct_recode(Parameter, "Distal" = "b_distal_xy",
                                "r_exy1" = "r_xy1",
                                "r_exy2" = "r_xy2")) %>% 
  ggplot(aes(Tdis, estimate, color = Parameter)) +
  geom_line(aes(linetype = `Parameter Type`), stat = "summary") +
  geom_point(aes(shape = `Parameter Type`), size = 1, stat = "summary") +
  facet_wrap(~ rxy, labeller = labeller(rxy = rxy_labs)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = brewer.pal(n=7,name="Dark2")[c(2,6,8)]) +
  coord_cartesian(xlim = c(0,40), 
                  ylim = c(min(keep$r_xy2), max(keep2$r_xy2))) +
  scale_y_continuous(breaks = seq(-.6,.6,0.3)) +
  labs(title = "CLPM", y = "Estimate", x = "Time Interval between Study Waves") +
  theme_bw(14) + 
  theme(legend.position = "right")


### 2. IV-CLPM ------------------------------------------------------------------

ivclpm <-
  keep %>% 
  pivot_longer(cols = c('b_prox1_xy','b_distal_xy','b_prox2_xy','r_xy1','r_xy2'), 
               names_to = "parameter",
               values_to = "estimate") %>% 
  mutate(`Parameter Type` = as_factor(ifelse(parameter %in% c("b_prox1_xy","b_distal_xy","b_prox2_xy"),
                                             "Regression", "Correlation")),
         Parameter = as_factor(parameter),
         Parameter = fct_recode(Parameter, "Proximal_W1" = "b_prox1_xy",
                                "Distal" = "b_distal_xy",
                                "Proximal_W2" = "b_prox2_xy",
                                "r_exy1" = "r_xy1",
                                "r_exy2" = "r_xy2"),
         Parameter = fct_relevel(Parameter, "Proximal_W1")) %>% 
  ggplot(aes(Tdis, estimate, color = Parameter)) +
  geom_line(aes(linetype = `Parameter Type`), stat = "summary") +
  geom_point(aes(shape = `Parameter Type`), size = 1, stat = "summary") +
  facet_wrap(~ rxy, labeller = labeller(rxy = rxy_labs)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = brewer.pal(n=7,name="Dark2")[c(1:3,6,8)]) +
  coord_cartesian(xlim = c(0,40), 
                  ylim = c(min(keep$r_xy2), max(keep2$r_xy2))) +
  scale_y_continuous(breaks = seq(-.6,.6,0.3)) +
  labs(title = "IV-CLPM", y = "Estimate", x = "Time Interval between Study Waves") +
  theme_bw(14) + 
  theme(legend.position = "right")



## Compile ----------------------------------------------------------------------

pdf("figs2_effect_by_residual_corr.pdf", width = 10, height = 10)
clpm / ivclpm +
  plot_annotation(tag_levels = 'A', tag_suffix = ".")
dev.off()


