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
ggplot(bxf.amp.merge[Well %in% bxf.ct[`Target Name` == "OmiS" & CT != "Undetermined", Well] & `Target Name` == "N Gene" & Cycle >= 15], aes(x=Cycle, y=log10(`Delta Rn`), group=factor(.id), col=.id)) +
  geom_path() +
  #geom_vline(data = bxf.ct.merge[Well %in% bxf.ct[`Target Name` == "OmiS" & CT != "Undetermined", Well] & `Target Name` == "N Gene"], show.legend = F, lty=2, aes(xintercept = as.numeric(CT), col=.id)) +
  labs(y="Delta Rn (log10)") +
  facet_wrap(~`Well Position`, ncol = 7) +
  scale_color_brewer(palette="Set1", name="Probe") +
  theme_cowplot() +
  theme(strip.background = element_blank(), strip.text = element_blank())

# ct difference
p.n1custom.ct <-
ggplot(dcast(bxf.ct.merge[qc_pass == T & covid == T & `Target Name` == "N Gene" & BA.1 == T], Well~.id, value.var = "CT"), aes(x=as.qpcr(N1), y=as.qpcr(N1custom) )) + 
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  annotate("text", x=23,y=35, label=dcast(bxf.ct.merge[qc_pass == T & covid == T & `Target Name` == "N Gene" & BA.1 == T], Well~.id, value.var = "CT")[,paste("R^2 =",round(cor(as.numeric(N1), as.numeric(N1custom))^2, 3))] ) +
  labs(x= "Ct (CDC N1 probe)", y="Ct (Omicron N1 probe)") +
  coord_cartesian(xlim=c(20,40), ylim=c(20,40)) +
  scale_color_brewer(palette="Set1") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/customN1_amplification.pdf", width = 8, height = 4, p.n1custom.amp)
ggsave2("plots/customN1_ct.pdf", width = 3, height = 3, p.n1custom.ct)