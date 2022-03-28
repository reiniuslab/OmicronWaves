library(data.table)

# amplification curve
bxf.amp <- fread("data_decoded/BXF_PLATE_K_N1_AmplificationCurve.csv")

## PLOT
library(ggplot2)
library(cowplot)

p.bxf.amp <-
ggplot(bxf.amp[`Well Position` %in% c("A2", "E6")], aes(x=Cycle, y=log10(`Delta Rn`), col=`Target Name`)) +
  geom_line() +
  facet_wrap(~`Well Position`, ncol = 1) +
  scale_color_manual(values=c("#2B4B9B", "#0A652E", "#999999")) +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/bxf_ampcurve.pdf", width = 5,height = 4, p.bxf.amp)
