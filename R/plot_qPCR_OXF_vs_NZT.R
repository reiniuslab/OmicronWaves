library(data.table)

##

as.qpcr <- function(x, lim=50){
  out <- as.numeric(x)
  out[x == "Undetermined"] <- lim
  out
}

##

# pangolin results
pangolin <- fread("data_decoded/Pangolin_results_denovo_DECODED.csv")
pangolin[, Sample := gsub("/.*", "", `Sequence name`)]

# qPCR data
oxf.nzt <- fread("data_decoded/OXF_NZT_20220308_DECODED.csv")
oxf.nzt[pangolin, c("Lineage","Lineage_group") := list(Lineage, `Scorpio call`), on = "Sample_Name == Sample"]
oxf.nzt[, N := .N, by=Lineage_group]

## PLOT
library(ggplot2)
library(cowplot)

# OXF vs NZT comparisons
p.pcr.oxf.nzt.comp1 <-
ggplot(oxf.nzt[grepl("negative|positive", OXF_Outcome, ignore.case = T) & grepl("negative|positive", NZT_Outcome, ignore.case = T)], aes(x=as.qpcr(OXF_N1_Ct), y=as.qpcr(NZT_N1_Ct) )) +
  geom_hline(yintercept=50, col="grey", lty=2) +
  geom_vline(xintercept=50, col="grey", lty=2) +
  geom_point() +
  annotate("text", x = 40,y = 15, label = oxf.nzt[,paste0("R^2 = ", round(cor(as.numeric(OXF_N1_Ct), as.numeric(NZT_N1_Ct), use = "complete.obs")^2, 4))]) +
  geom_smooth(method="lm", aes(x=as.numeric(OXF_N1_Ct), y=as.numeric(NZT_N1_Ct))) +
  labs(x="OXF N1 Ct", y="NZT N1 Ct") +
  coord_cartesian(xlim=c(10,50), ylim=c(10,50)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.pcr.oxf.nzt.comp2 <-
ggplot(oxf.nzt[grepl("negative|positive", OXF_Outcome, ignore.case = T) & grepl("negative|positive", NZT_Outcome, ignore.case = T)], aes(x=as.qpcr(OXF_N1_Ct), y=as.qpcr(NZT_RdRp_Ct) )) +
  geom_hline(yintercept=50, col="grey", lty=2) +
  geom_vline(xintercept=50, col="grey", lty=2) +
  geom_point() +
  annotate("text", x = 40,y = 15, label = oxf.nzt[,paste0("R^2 = ", round(cor(as.numeric(OXF_N1_Ct), as.numeric(NZT_RdRp_Ct), use = "complete.obs")^2, 4))]) +
  geom_smooth(method="lm", aes(x=as.numeric(OXF_N1_Ct), y=as.numeric(NZT_RdRp_Ct))) +
  labs(x="OXF N1 Ct", y="NZT RdRp Ct") +
  coord_cartesian(xlim=c(10,50), ylim=c(10,50)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.pcr.oxf.nzt.comp3 <-
ggplot(oxf.nzt[grepl("negative|positive", OXF_Outcome, ignore.case = T) & grepl("negative|positive", NZT_Outcome, ignore.case = T)], aes(x=as.qpcr(NZT_N1_Ct), y=as.qpcr(NZT_RdRp_Ct) )) +
  geom_hline(yintercept=50, col="grey", lty=2) +
  geom_vline(xintercept=50, col="grey", lty=2) +
  geom_point() +
  annotate("text", x = 40,y = 15, label = oxf.nzt[,paste0("R^2 = ", round(cor(as.numeric(NZT_N1_Ct), as.numeric(NZT_RdRp_Ct), use = "complete.obs")^2, 4))]) +
  geom_smooth(method="lm", aes(x=as.numeric(NZT_N1_Ct), y=as.numeric(NZT_RdRp_Ct))) +
  labs(x="NZT N1 Ct", y="NZT RdRp Ct") +
  coord_cartesian(xlim=c(10,50), ylim=c(10,50)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.pcr.oxf.nzt.comp4 <-
ggplot(oxf.nzt[grepl("negative|positive", OXF_Outcome, ignore.case = T) & grepl("negative|positive", NZT_Outcome, ignore.case = T)], aes(x=as.qpcr(OXF_RNaseP_Ct), y=as.qpcr(NZT_RNaseP_Ct) )) +
  geom_hline(yintercept=50, col="grey", lty=2) +
  geom_vline(xintercept=50, col="grey", lty=2) +
  geom_point() +
  annotate("text", x = 40,y = 15, label = oxf.nzt[,paste0("R^2 = ", round(cor(as.numeric(OXF_RNaseP_Ct), as.numeric(NZT_RNaseP_Ct), use = "complete.obs")^2, 4))]) +
  geom_smooth(method="lm", aes(x=as.numeric(OXF_RNaseP_Ct), y=as.numeric(NZT_RNaseP_Ct))) +
  labs(x="OXF RNaseP Ct", y="NZT RNaseP Ct") +
  coord_cartesian(xlim=c(10,50), ylim=c(10,50)) +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/oxf_nzt_comparison.pdf", width = 6, height = 6,
plot_grid(
  p.pcr.oxf.nzt.comp1, p.pcr.oxf.nzt.comp2, 
  p.pcr.oxf.nzt.comp3, p.pcr.oxf.nzt.comp4
))


##
library(ggpubr)
library(ggbeeswarm)
oxf.nzt.melt <- melt(oxf.nzt, id.vars = c("Sample_Name", "OXF_Outcome", "NZT_Outcome", "Lineage", "Lineage_group", "N"))

p.pcr.oxf.nzt.omi <-
ggplot(oxf.nzt.melt[grepl("positive", OXF_Outcome, ignore.case = T) & !grepl("RNaseP|OmiS", variable)], aes(x=variable, y=as.numeric(value), group=paste(variable, OXF_Outcome) )) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.25, position = position_dodge(0.9), aes(fill=OXF_Outcome)) +
  stat_summary(fun = "median", geom="text", position = position_dodge(0.9), angle=90, aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-14)), geom="text", position = position_dodge(0.9) ) +
  stat_compare_means(tip.length = 0, step.increase = 0.025, method="wilcox.test", method.args = list(alternative = "less"), aes(label = paste0("P = ", ..p.format..) ) ) +
  labs(y = "Ct value", x = "Assay") +
  scale_y_reverse() +
  scale_fill_brewer(palette="Paired", name="Omicron (BA.1)", labels= c("PCR Negative", "PCR Positive") ) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))

ggsave2("plots/oxf_nzt_omi.pdf", width = 5, height = 6, p.pcr.oxf.nzt.omi)

p.pcr.oxf.nzt.lin <-
ggplot(oxf.nzt.melt[grepl("positive", OXF_Outcome, ignore.case = T) & !grepl("RNaseP|OmiS", variable) & !is.na(Lineage)], aes(x=variable, y=as.numeric(value), group=paste(variable, Lineage_group) )) +
  geom_quasirandom(stroke=NA, alpha=0.25, dodge.width = 0.9, aes(col=Lineage_group)) +
  geom_boxplot(width=0.25, outlier.shape = NA, position = position_dodge(0.9), aes(fill=Lineage_group)) +
  stat_summary(fun = "median", geom="text", position = position_dodge(0.9), angle=90, aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-15)), geom="text", position = position_dodge(0.9) ) +
  #stat_compare_means(tip.length = 0, step.increase = 0.025, aes(label = paste0("P = ", ..p.format..) ) ) +
  labs(y = "Ct value", x = "Assay") +
  scale_y_reverse() +
  scale_fill_brewer(palette="Pastel1", name="Lineage", aesthetics = c("col", "fill") ) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))

ggsave2("plots/oxf_nzt_lineage.pdf", width = 5, height = 6, p.pcr.oxf.nzt.lin)
