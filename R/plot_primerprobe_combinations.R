library(data.table)

# primer-probe combinations
fls.pp <- list.files("primer_probe_combinations", "^S.*.txt", full.names = T)
names(fls.pp) <- gsub("\\..*", "", basename(fls.pp))

dat.pp <- as.data.table(melt(lapply(fls.pp, read.delim), id.vars="X"))
colnames(dat.pp)[1:2] <- c("rv", "fw")
dat.pp[, c("target", "variant", "probe") := tstrsplit(L1, "_")]

## PLOT
library(ggplot2)
library(cowplot)

p.pp.heat <-
ggplot(dat.pp, aes(x=rv, y=fw, fill=as.numeric(value))) +
  geom_tile() +
  facet_grid(variant~probe, scales="free", space="free") +
  labs(x="Reverse primer", y="Forward primer") +
  scale_y_discrete(limits=rev) +
  scale_fill_viridis_c(name="Ct", direction = -1) +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.border = element_rect(color="black"))

ggsave2("plots/primerprobe_combinations_heatmap.pdf", width = 8, height = 4, p.pp.heat)