#' ---
#' title: "R code for the analysis conducted in 'Characterizing brain age the Alzheimer\\'s disease connectome project using a deep neural network pre-trained on the UK biobank'"
#' author: "Nagesh Adluru"
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#' ---
#'
#' # Initialization
# Loading the libraries =========
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)
library(GGally)
library(forcats)
library(stringr)
library(latex2exp)
library(modelr)
library(tidytext)
library(lemon)
library(scales)
library(purrr)

# Initializing variables ======
rm(list = ls(all = T))
csvroot = 'CSVs/'
figroot = 'Figures/'

# ggplot theme ====
dodge = position_dodge(width = 0.9)
txtSize = 12
gtheme = theme(
  legend.key = element_rect(colour = "black"),
  legend.title = element_text(size = txtSize),
  legend.text = element_text(size = txtSize),
  legend.background = element_blank(),
  legend.position = "top", 
  strip.text.x = element_text(size = txtSize),
  strip.text.y = element_text(size = txtSize),
  strip.background = element_blank(),
  axis.text = element_text(colour = "black", size = txtSize),
  axis.title.x = element_text(size = txtSize),
  axis.title.y = element_text(size = txtSize),
  plot.title = element_text(size = txtSize),
  axis.title = element_text(size = txtSize),
  axis.line = element_line(),
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "gray"),
  panel.grid.minor = element_line(size = 0.2, linetype = 'solid', colour = "gray"),
  # ticks
  axis.ticks.length = unit(0.25, "cm")
)

#' # figure 1 (ca vs. ba)
set.seed(1)
# Loading the predicted brain age values
df = read.csv(paste0(csvroot, 'obadcp_binrange2.csv')) %>%
  mutate(BrainAGE = EstimatedAge - CalendarAge) %>% 
  inner_join(read.csv(paste0(csvroot, 'NIHToolBoxReshaped_01062021.csv')), 
             by = c('ID' = 'SubjectID'))

# Performing the bias correction
df %<>% crossv_loo %>% 
  mutate(model = map(train, ~ lm(BrainAGE ~ CalendarAge, data = .)),
         predicted = map2(model, test, ~ augment(.x, newdata = .y))) %>% 
  unnest(predicted) %>% 
  mutate(BiasCorrectedEstimatedAge = EstimatedAge - .fitted,
         BiasCorrectedBrainAGE = BiasCorrectedEstimatedAge - CalendarAge) %>% 
  select(-train, -test)

dfsample = df %>% 
  select(CalendarAge, Sex, ConsensusDiagnosis) %>% 
  group_by(Sex, ConsensusDiagnosis) %>% 
  summarise(n = n(), 
            ma = mean(CalendarAge), 
            sa = sd(CalendarAge), .groups = 'drop')

#+ fig.width=7.85, fig.height=5.05, warning=F
p = df %>% select(ID, CalendarAge, EstimatedAge, BiasCorrectedEstimatedAge, ConsensusDiagnosis, Sex) %>% 
  pivot_longer(c(BiasCorrectedEstimatedAge, EstimatedAge), 
               names_to = 'BA',
               values_to = 'BAValue') %>% 
  mutate(BA = plyr::mapvalues(BA, 
                              from = c('EstimatedAge',
                                       'BiasCorrectedEstimatedAge'), 
                              to = c('Uncorrected', 
                                     'Bias corrected'))) %>% 
  ggplot(aes(x = CalendarAge, 
             y = BAValue, 
             shape = BA)) + 
  geom_vline(data = dfsample, aes(xintercept = ma), linetype = 'longdash') +
  geom_vline(data = dfsample, aes(xintercept = ma - sa), alpha = 0.75) +
  geom_vline(data = dfsample, aes(xintercept = ma + sa), alpha = 0.75) +
  geom_rect(data = dfsample, aes(xmin = ma - sa, xmax = ma + sa, 
                                 ymin = -Inf, ymax = Inf), 
            alpha = 0.1, fill = 'gray64', inherit.aes = F) + 
  geom_point() + 
  geom_line(aes(group = ID), 
            alpha = 0.5, size = 0.05) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text(data = dfsample, aes(x = Inf, y = -Inf, 
                                 label = paste0('n = ', n)), 
            size = txtSize/2, hjust = 1, vjust = -1, inherit.aes = F) +
  geom_text(data = df %>% 
              select(BiasCorrectedBrainAGE, ConsensusDiagnosis, Sex) %>%
              group_by(ConsensusDiagnosis, Sex) %>% 
              summarise(mae = mean(abs(BiasCorrectedBrainAGE)), 
                        .groups = 'drop'), 
            aes(x = Inf, y = -Inf, 
                label = paste0('MAE = ', round(mae, 1), ' [y]')), 
            size = txtSize / 2.5, hjust = 1, vjust = -3.0, inherit.aes = F) +
  facet_rep_grid(Sex ~ fct_relevel(ConsensusDiagnosis, 'AD', 
                                   after = 2)) + 
  gtheme + scale_shape_manual(values = c(0, 1)) + 
  labs(x = 'Calendar age [y]', y = 'FCNN brain age [y]') + 
  theme(legend.title = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(5, -10, -10, -10), 
        strip.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        strip.text.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16), 
        legend.key = element_blank(), 
        legend.background = element_blank()) + 
  scale_y_continuous(breaks = pretty_breaks(n = 5))
p
ggsave(paste0(figroot, 'BA_ADCP_AAIC2021.pdf'), 
       p, 
       width = 7.85, 
       height = 5.05)

#' # figure 2 (box plots)
#+ fig.width=4.05, fig.height=3.80, warning=F
p = df %>% 
  ggplot(aes(y = BiasCorrectedBrainAGE, 
             x = fct_relevel(ConsensusDiagnosis, 'AD', after = 2), 
             color = Sex)) + 
  geom_boxplot(alpha = 1, outlier.shape = 1, outlier.size = 4) + 
  geom_hline(yintercept = 0, size = 1.0, linetype = 'longdash') + 
  gtheme + 
  theme(legend.title = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(5, -10, -20, -10), 
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size = 16), 
        axis.text.x = element_text(size = 16), 
        legend.key = element_blank(), 
        legend.background = element_blank()) + 
  labs(x = '', y = 'Excess aging of the brain [y]') + 
  scale_color_brewer(palette = 'Set1') + 
  scale_y_continuous(breaks = pretty_breaks(10))
p
ggsave(paste0(figroot, 'BoxPlot_ADCP_AAIC2021.pdf'), 
       p, 
       width = 4.05, 
       height = 3.80)