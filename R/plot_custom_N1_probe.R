library(data.table)

##

as.qpcr <- function(x, lim=50){
  out <- as.numeric(x)
  out[x == "Undetermined"] <- lim
  out
}

##

# C28311T

# amp curves
bxf.amp <- fread("data_decoded/BXF_PLATE_K_N1_AmplificationCurve.csv")
bxf.amp.custom <- fread("data_decoded/BXF_PLATE_L_N1CUSTOM_AmplificationCurve.csv")

bxf.amp.merge <- rbindlist(list("N1" = bxf.amp, "N1custom" = bxf.amp.custom), idcol = T)
  
# ct values
bxf.ct <- fread("data_decoded/BXF_PLATE_K_N1_Results.csv")
bxf.ct.custom <- fread("data_decoded/BXF_PLATE_L_N1CUSTOM_Results.csv")

bxf.ct.merge <- rbindlist(list("N1" = bxf.ct, "N1custom" = bxf.ct.custom), idcol = T)
bxf.ct.merge[, qc_pass := CT[`Target Name` == "RNAse P"] != "Undetermined", by=c("Well", ".id")]
bxf.ct.merge[, covid := CT[`Target Name` == "N Gene"] != "Undetermined", by=c("Well", ".id")]
bxf.ct.merge[, BA.1 := CT[`Target Name` == "OmiS"] != "Undetermined", by=c("Well", ".id")]

## PLOT
library(ggplot2)
library(cowplot)

# amp curves for modified N1 probe
p.n1custom.amp <-
ggplot(bxf.amp.merge[Well %in% bxf.ct[`Target Name` == "OmiS" & CT != "Undetermined", Well] & `Target Name` == "N Gene" & Cycle >= 15], aes(x=Cycle, y=`Delta Rn`, group=factor(.id), col=.id)) +
  geom_path() +
  geom_vline(data = bxf.ct.merge[Well %in% bxf.ct[`Target Name` == "OmiS" & CT != "Undetermined", Well] & `Target Name` == "N Gene"], show.legend = F, lty=2, aes(xintercept = as.numeric(CT), col=.id)) +
  facet_wrap(~`Well Position`, ncol = 7) +
  scale_color_brewer(palette="Set1", name="Probe") +
  theme_cowplot() +
  theme(strip.background = element_blank(), strip.text = element_blank())

# ct difference
p.n1custom.ct <-
ggplot(bxf.ct.merge[qc_pass == T & covid == T & `Target Name` == "N Gene"], aes(x=.id, y=as.qpcr(CT), col=.id)) +
  geom_line(col="black", aes(group=factor(Well))) +
  geom_point(show.legend = F) +
  stat_summary(fun.y="median", geom="text", col="black", aes(label=round(..y..,2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=38)), geom="text", show.legend = F, col="black" ) +
  facet_grid(~BA.1) +
  labs(x="Probe", y="Ct value") +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave2("plots/customN1_amplification.pdf", width = 6, height = 4, p.n1custom.amp)
ggsave2("plots/customN1_ct.pdf", width = 3, height = 4, p.n1custom.ct)