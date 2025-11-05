# Combine satured plots
# 11/05/2025
# RCS

#setup working directories

repo <- "C:/Users/rcscott/MetaDisease"
graphs <- paste(repo,"/Graphs",sep="")

# Load library
library(ggplot2)
library(patchwork)

#load figs

setwd(graphs)
TwoPatches <- readRDS("alpha_sat_plot2.RDS")
TwoPatches

FivePatches <- readRDS("alpha_sat_plot5.RDS")
FivePatches

Sat_plots <- TwoPatches + FivePatches + plot_annotation(tag_levels = "A")
Sat_plots
ggsave(Sat_plots, filename = "Sat_plots.png")
