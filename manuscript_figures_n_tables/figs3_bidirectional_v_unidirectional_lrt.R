## Figure: Bidirectional vs Unidirectional IV-CLPM

## Contrast the LRTs in bidirectional and unidirectional IV-CLPM, given unidirectional causal effect in the true model.


# clear workspace
rm(list=ls(all=TRUE))

# SET THE WORKING DIRECTORY
# getwd()

# Libraries
library(tidyverse)
library(RColorBrewer)

# Data --------------------------------------------------------------------------

## From 5_bidir_v_unidir_IV-CLPM.R
keep <- read.table("resfile_bidir_v_unidir_IV-CLPM.dta")
colnames(keep) <-  c(
  ## Data-generating parameters
  'model','Tdis',
  ## LRT Statistics for the Effect of X on Y
  'NCP_distal','NCP_prox1','NCP_prox2',
  ## LRT Statistics for Bidirectional
  'NCP_distal_bd','NCP_prox1_bd','NCP_prox2_bd'
)

psych::describe(keep)

keep <- as_tibble(keep)

# Plots -------------------------------------------------------------------------


bidir_v_unidir <- keep |> 
  pivot_longer(cols = c("NCP_prox1","NCP_distal","NCP_prox2"), 
               names_to = "parameter", names_prefix = "NCP_",
               values_to = "ncp") %>% 
  mutate(parameter = fct_relevel(as_factor(parameter), "prox1"),
         parameter = fct_recode(parameter,"Proximal_W1" = "prox1",
                                "Distal" = "distal",
                                "Proximal_W2" = "prox2"),
         model = fct_recode(as_factor(model), "Unidirectional IV-CLPM" = "unidir", 
                            "Bidirectional IV-CLPM" = "bidir")) %>% 
  ggplot(aes(x = Tdis, y = ncp, color = parameter)) +
  geom_line() +
  geom_point(size=1) +
  facet_wrap(~ model) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0,max(keep$NCP_prox1_bd))) +
  labs(title = "1-DF LRTs of Unidirectional Effects",
       y = "LRT Statistic", x = "Time Interval", color = "Effect Tested") +
  theme_bw(12) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))



pdf("figs3_bidir_v_unidir_ncp.pdf", width = 10, height = 6)
bidir_v_unidir
dev.off()
