library(devtools)
#install_github("GuangchuangYu/ggtree")
#install_github("revbayes/RevGadgets")
library(ggplot2)
library(RevGadgets)

HiSSE_file <- paste0("output/reptiles_HiSSE_anc_states_results.tree")
p_anc <- processAncStates(HiSSE_file)
plot <- plotAncStatesMAP(p_anc,
                         tree_layout = "rect",
                         tip_labels_size = 1) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.1,0.85),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.title = element_text(size=6), #change legend title font size
        legend.text = element_text(size=4))
plot
#ggsave(paste0("HiSSE_anc_states_activity_period.png"),plot, width=8, height=8)

HiSSE_file <- paste0("output/reptiles_HiSSE_plasticity.log")
pdata <- processSSE(HiSSE_file)

# plot the rates
plot <- plotMuSSE(pdata) +
  theme(legend.position = c(0.875,0.915),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))
plot

ggsave(paste0("HiSSE_div_rates_activity_period.png"),plot, width=5, height=5)
