library(devtools)
library(ggplot2)
library(RevGadgets)

HiSSE_file <- paste0("output/reptiles_HiSSE_plasticity12.log")
pdata <- processSSE(HiSSE_file)

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

#mcmc_max <- nrow(general_rates)
#mcmc_out_burn_in <- round(nrow(general_rates) * 0.2) + 1
#general_rates <- general_rates[mcmc_out_burn_in:mcmc_max, ]
#speciation2 <- speciation2[mcmc_out_burn_in:mcmc_max, ]
#extinction2 <- extinction2[mcmc_out_burn_in:mcmc_max, ]
#diversification2 <- diversification2[mcmc_out_burn_in:mcmc_max, ]

n.eff(as.matrix(general_rates[, -1]))
n.eff(as.matrix(speciation2[, -1]))
n.eff(as.matrix(extinction2[, -1]))
n.eff(as.matrix(diversification2[, -1]))

###### BAYES FACTOR CALCULATION ######
bf_mean <- function(x, y) x / y
bf_timestep <- function(x, y) mean(x / y)
mean_posteriors <- colMeans(general_rates)[2:ncol(general_rates)]
round(mean_posteriors, 5)
# lamb 1 x lamb 2 #
bf_mean(mean_posteriors[1], mean_posteriors[2])
# lamb 1 x lamb 3 #
bf_mean(mean_posteriors[1], mean_posteriors[3])
# lamb 2 x lamb 3 #
bf_mean(mean_posteriors[2], mean_posteriors[3])

# mu 1 x mu 2 #
bf_mean(mean_posteriors[4], mean_posteriors[5])
# mu 1 x mu 3 #
bf_mean(mean_posteriors[4], mean_posteriors[6])
# mu 2 x mu 3 #
bf_mean(mean_posteriors[5], mean_posteriors[6])

# transition #
# q12 x q13 #
bf_mean(mean_posteriors[7], mean_posteriors[8])
# q12 x q21 #
bf_mean(mean_posteriors[7], mean_posteriors[9])
# q12 x q23 #
bf_mean(mean_posteriors[7], mean_posteriors[10])
# q12 x q31 #
bf_mean(mean_posteriors[7], mean_posteriors[11])
# q12 x q32 #
bf_mean(mean_posteriors[7], mean_posteriors[12])
# q13 x q21 #
bf_mean(mean_posteriors[8], mean_posteriors[9])
# q13 x q23 #
bf_mean(mean_posteriors[8], mean_posteriors[10])
# q13 x q31 #
bf_mean(mean_posteriors[8], mean_posteriors[11])
# q13 x q32 #
bf_mean(mean_posteriors[8], mean_posteriors[12])
# q21 x q23 #
bf_mean(mean_posteriors[9], mean_posteriors[10])
# q21 x q31 #
bf_mean(mean_posteriors[9], mean_posteriors[11])
# q21 x q32 #
bf_mean(mean_posteriors[9], mean_posteriors[12])
# q23 x q31 #
bf_mean(mean_posteriors[10], mean_posteriors[11])
# q23 x q32 #
bf_mean(mean_posteriors[10], mean_posteriors[12])
# q31 x q32 #
bf_mean(mean_posteriors[11], mean_posteriors[12])

mean_posteriors[7:12]

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
# diversification #
# diversi 1 x diversi 2 #
bf_mean(diversification_means[1], diversification_means[2])
# diversi 1 x diversi 3 #
bf_mean(diversification_means[1], diversification_means[3])
# diversi 2 x diversi 3 #
bf_mean(diversification_means[2], diversification_means[3])

#################### plots ############################
# temporal series #
save_result <- TRUE
markov <- ggplot(mcmc_result_pivoted, aes(Iteration, DiversificationPosterior)) +
  geom_line() + 
  xlab('Time') + ylab('Log(L)') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('markov_result', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(markov)
  dev.off()
}

markov_lambda <- ggplot(speciation2, aes(x=Iteration)) +
  geom_line(aes(y=speci0A), color='blue') +
  geom_line(aes(y=speci0B), color='black') +
  geom_line(aes(y=speci1A), color='yellow') +
  geom_line(aes(y=speci1B), color='red') +
  xlab('Time') + ylab('Lambda') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('markov_result_lambda', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(markov_lambda)
  dev.off()
}

# speciation #
speciation <- mcmc_result_pivoted %>% 
  ggplot(aes(LambdaPosterior, fill = Speciation)) + 
  geom_density(alpha=0.7) +
  scale_x_log10() +
  scale_x_sqrt() +
  theme_bw() +
  scale_fill_hue(labels = c("0A", 
                            "0B", 
                            "1A",
                            "1B"),
                 name="States") +
  theme(legend.position = c(0.8, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
  xlab("Speciation") + ylab('Posterior Density') 


if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('speciation', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(speciation)
  dev.off()
}

# extinction #
extinction <- mcmc_result_pivoted %>% 
  ggplot(aes(ExtinctionPosterior, fill = Extinction)) + 
  geom_density(alpha=0.7) +
  scale_x_log10() +
  scale_x_sqrt() +
  theme_bw() +
  scale_fill_hue(labels = c("0A", 
                            "0B", 
                            "1A",
                            "1B"),
                 name="States") +
  theme(legend.position = c(0.8, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
  xlab("Extinction") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('extinction', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(extinction)
  dev.off()
}

# diversification #
diversification <- mcmc_result_pivoted %>% 
  ggplot(aes(DiversificationPosterior, fill = Diversification)) + 
  geom_density(alpha=0.7) +
  scale_x_log10() +
  scale_x_sqrt() +
  theme_bw() +
  scale_fill_hue(labels = c(labels = c("0A", 
                                       "0B", 
                                       "1A",
                                       "1B")),
                 name="States") +
  theme(legend.position = c(0.8, 0.75),
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) +
  xlab("Diversification") + ylab('Posterior Density') 

if(save_result == TRUE){
  tiff(filename = file.path(here('output'), paste0('diversification', ".tif")),
       #width = 1000,
       width = 800, #without abundance and occupancy
       height = 600,
       #height = 400, #abundance and occupancy
       units = "px",
       res = 100)
  print(diversification)
  dev.off()
}

# plot the rates

ggplot2::ggplot(pdata, ggplot2::aes(x = value, fill = observed_state)) + 
  ggplot2::geom_density(alpha = 0.8) + ggplot2::scale_fill_manual(values = colFun(length(unique(pdata$observed_state))), 
                                                                  name = "Observed state") + ggplot2::facet_wrap(rate ~ hidden_state, 
                                                                                                                 scales = "free", ncol = 2) + 
  ggplot2::xlab("Rate") + ggplot2::ylab("Posterior density") + 
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                                       panel.grid.minor = ggplot2::element_blank(), 
                                       strip.background = ggplot2::element_blank())

ggsave(paste0("HiSSE_div_rates_activity_period.png"),plot, width=10, height=5)
