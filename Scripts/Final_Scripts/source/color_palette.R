##Create color palettes for final figures 

library(tidyverse)
library(dplyr)
library(RColorBrewer)

brewer.pal(n = 7, name = "Set2")
display.brewer.pal(n = 7, name = "Set2")

## Genotype colors (8)

genotype_rgb_colors <- c("XYF" = "rgb(166,206,227)",
                     "XYM" = "rgb(31,120,180)",
                     "XXF" = "rgb(251,154,153)",
                     "XXM" = "rgb(227,26,28)",
                     "XXYF" = "rgb(253,191,111)",
                     "XXYM" = "rgb(255,127,0)",
                     "XYYF" = "rgb(178,223,138)",
                     "XYYM" = "rgb(51,160,44)")

genotype_hex_colors <- c("XYF" = "#A6CEE3",
                     "XYM" = "#1F78B4",
                     "XXF" = "#FB9A99",
                     "XXM" = "#E31A1C",
                     "XXYF" = "#FDBF6F",
                     "XXYM" = "#FF7F00",
                     "XYYF" = "#B2DF8A",
                     "XYYM" = "#33A02C")

## Gonad colors (2)

gonad_hex_colors <- c("Male" = "#007FFF",
                      "Female" = "#FF0000")

## SexChrom colors (4) 

sexchrom_hex_colors <- c("XY" = "#4A90E2",
                         "XX" = "#D0021B",
                         "XXY" = "#F5A623",
                         "XYY" = "#7ED321")

## Trisomy colors (2) 

trisomy_hex_colors <- c("Non-Trisomy" = "#1f77b4",
                      "Trisomy" = "#ff7f0e")

## Cell type colors (11)

cell_rgb_colors <- c("Astrocyte" = "rgb(141,211,199)",
                 "Ependymal" = "rgb(255,255,179)",
                 "GABA" = "rgb(190,186,218)",
                 "Chol" = "rgb(251,128,114)",
                 "Glut" = "rgb(128,177,211)",
                 "Inh-IMN" = "rgb(253,180,98)",
                 "Microglia" = "rgb((179,222,105)",
                 "Misc-NT" = "rgb(217,217,217)",
                 "Oligo" = "rgb(252,205,229)",
                 "OPC" = "rgb(188,128,189)",
                 "Vascular" = "rgb(204,235,197)")

cell_hex_colors <- c("Astrocyte" = "#8dd3c7",
                 "Ependymal" = "#ffeda0",
                 "GABA" = "#bebada",
                 "Chol" = "#fb8072",
                 "Glut" = "#80b1d3",
                 "Inh-IMN" = "#fdb462",
                 "Microglia" = "#b3de69",
                 "Misc-NT" = "#d9d9d9",
                 "Oligo" = "#fccde5",
                 "OPC" = "#bc80bd",
                 "Vascular" = "#ccebc5")

## Sex factor colors (7)

# effect_rgb_colors <- c("Normative" = "rgb(231,111,81)",
#                    "Gonad" = "rgb(244,162,97)",
#                    "SexChrom" = "rgb(138,177,125)",
#                    "Gonad:SexChrom" = "rgb(233,196,106)",
#                    "addX" = "rgb(132,155,203)",
#                    "addY1" = "rgb(42,157,143)",
#                    "addY2" = "rgb(59,172,182)")

# effect_rgb_colors <- c(
#   "Normative" = "rgb(252,141,98)",
#   "Gonad" = "rgb(231,138,195)",
#   "SexChrom" = "rgb(166,216,84)",
#   "Gonad:SexChrom" = "rgb(255,217,47)",
#   "addX" = "rgb(120,180,255)",
#   "addY1" = "rgb(102,194,165)",
#   "addY2" = "rgb(180,180,255)"
# )
# 
# effect_hex_colors <- c(
#   "Normative" = "#FC8D62",    # Light coral
#   "Gonad" = "#E78AC3",        # Light orange
#   "SexChrom" = "#A6D854",     # Light teal green
#   "Gonad:SexChrom" = "#FFD92F", # Light golden yellow
#   "addX" = "#78b4ff",         # Light sky blue
#   "addY1" = "#66C2A5",        # Light aqua
#   "addY2" = "#b4b4ff"         # Light indigo-blue
# )

effect_rgb_colors <- c(
  "Normative" = "rgb(252,141,98)",
  "GE" = "rgb(231,138,195)",
  "SCE" = "rgb(166,216,84)",
  "GExSCE" = "rgb(255,217,47)",
  "XCD" = "rgb(120,180,255)",
  "YCD1" = "rgb(102,194,165)",
  "YCD2" = "rgb(180,180,255)"
)

effect_hex_colors <- c(
  "Typical" = "#FC8D62",    # Light coral
  "GE" = "#E78AC3",        # Light orange
  "SCE" = "#A6D854",     # Light teal green
  "GExSCE" = "#FFD92F", # Light golden yellow
  "XCD" = "#78b4ff",         # Light sky blue
  "YCD1" = "#66C2A5",        # Light aqua
  "YCD2" = "#b4b4ff"         # Light indigo-blue
)
