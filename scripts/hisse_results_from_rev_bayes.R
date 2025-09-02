library(tidyverse)
library(devtools)
library(ggplot2)
library(RevGadgets)
library(stableGR)

HiSSE_file <- paste0("output/reptiles_HiSSE_plasticity1_50000_time.log")
pdata <- processSSE(HiSSE_file)
pdata

speciation <- pdata[pdata['rate'] == 'speciation', ]
extinction <- pdata[pdata['rate'] == 'extinction', ]
diversification <- pdata[pdata['rate'] == 'net-diversification', ]

speciation2 <- speciation |>
  pivot_wider(names_from = label, values_from = "value", 
              id_cols= Iteration, names_prefix = "speci")

extinction2 <- extinction |>
  pivot_wider(names_from = label, values_from = "value", 
              id_cols= Iteration, names_prefix = "extinc")

diversification2 <- diversification |>
  pivot_wider(names_from = label, values_from = "value", 
              id_cols= Iteration, names_prefix = "divers")

general_rates <- inner_join(speciation2, extinction2, by="Iteration")
general_rates <- inner_join(general_rates, diversification2, by="Iteration")

round(colMeans(general_rates), 5)

mcmc_result_pivoted <- general_rates %>% 
  pivot_longer(starts_with('spec'),
               names_to="Speciation",
               values_to = 'LambdaPosterior')

temporary <- general_rates %>% 
  pivot_longer(starts_with('extin'),
               names_to="Extinction",
               values_to = 'ExtinctionPosterior')

temporary2 <- general_rates %>% 
  pivot_longer(starts_with('divers'),
               names_to="Diversification",
               values_to = 'DiversificationPosterior')

mcmc_result_pivoted$Extinction <- temporary$Extinction
mcmc_result_pivoted$ExtinctionPosterior <- temporary$ExtinctionPosterior

mcmc_result_pivoted$Diversification <- temporary2$Diversification
mcmc_result_pivoted$DiversificationPosterior <- temporary2$DiversificationPosterior

### bf diversification ###
bf_diversif <- mcmc_result_pivoted |>
  select(Diversification, DiversificationPosterior) |>
  group_by(Diversification) |>
  mutate(row = row_number()) |>
  pivot_wider(names_from=Diversification,
              values_from=DiversificationPosterior) |>
  select(-row)
diversification_means <- colMeans(bf_diversif)

#################### plots ############################
# plot the rates
# plot1 <- ggplot2::ggplot(pdata, ggplot2::aes(x = value, fill = observed_state)) + 
#   ggplot2::geom_density(alpha = 0.8) + ggplot2::scale_fill_manual(values = colFun(length(unique(pdata$observed_state))), 
#                                                                   name = "Observed state") + ggplot2::facet_wrap(rate ~ hidden_state, 
#                                                                                                                  scales = "free", ncol = 2) + 
#   ggplot2::xlab("Rate") + ggplot2::ylab("Posterior density") + 
#   ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
#                                        panel.grid.minor = ggplot2::element_blank(), 
#                                        strip.background = ggplot2::element_blank())
# ggsave(paste0("HiSSE_density_50000.png"),plot1, width=10, height=5)

#plot2 <- ggplot2::ggplot(pdata, ggplot2::aes(x = observed_state, y = value, fill = observed_state)) + 
#  ggplot2::geom_boxplot(alpha = 0.8) + 
#  ggplot2::scale_fill_manual(values = colFun(length(unique(pdata$observed_state))), 
#                             name = "Observed state") + 
#  ggplot2::facet_wrap(rate ~ hidden_state, scales = "free", ncol = 2) + 
#  ggplot2::xlab("Observed state") + 
#  ggplot2::ylab("Posterior density") + 
#  ggplot2::theme_bw() + 
#  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
#                 panel.grid.minor = ggplot2::element_blank(), 
#                 strip.background = ggplot2::element_blank())


fill_colors <- colFun(length(unique(pdata$observed_state)))  # ex: vermelho e azul

# cores um pouco mais escuras usando alpha
line_colors <- alpha(fill_colors, 0.8)  # deixa ligeiramente mais escuro

plot2 <- ggplot(pdata, aes(x = observed_state, y = value, fill = observed_state)) + 
  geom_boxplot(aes(color = observed_state), alpha = 0.8, size = 0.7) + 
  scale_fill_manual(values = fill_colors, name = "Observed state") +
  scale_color_manual(values = line_colors, guide = "none") +  # linhas e pontos
  facet_wrap(rate ~ hidden_state, scales = "free", ncol = 2) + 
  xlab("Observed state") + 
  ylab("Posterior density") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

ggsave(paste0("HiSSE_boxplot_50000.png"),plot2, width=8, height=6)
write.csv(pdata, 
          file = "hisse_data.csv", 
          row.names = FALSE)
