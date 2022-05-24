library(data.table)

##

median_cl_boot <- function(x, conf = 0.95) {
  # function from http://rstudio-pubs-static.s3.amazonaws.com/28101_41a7995107d94c8dbb07bbf7cd7e8291.html
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}

##

# pangolin results
pangolin <- fread("data_decoded/Pangolin_results_denovo_DECODED.csv")
pangolin[, Sample := gsub("/.*", "", `Sequence name`)]

# qPCR data
oxf.clin <- fread("data_decoded/OXF_20220308_DECODED.csv")
oxf.clin[, Analysis_Date := as.Date(Analysis_Date, tryFormats = c("%d/%m/%Y"))]
oxf.cont <- fread("data_decoded/Control_samples_BXF_OXF_20220308_DECODED.csv")

# merge CT data
oxf.merge <- rbind(oxf.cont, oxf.clin, fill=T)
oxf.merge[pangolin, c("Lineage","Lineage_group") := list(Lineage, `Scorpio call`), on = "Sample_Name == Sample"]
oxf.merge[, N := .N, by=Lineage_group]

## PLOT
library(ggplot2)
library(cowplot)

# pcr timeseries
p.pcr.timeseries.swe <- 
ggplot(oxf.clin[PCR_Outcome == "POSITIVE"], aes(x= Analysis_Date, fill = OXF_Outcome)) + 
  geom_bar(position="fill") +
  geom_text(data = oxf.clin[PCR_Outcome == "POSITIVE", .N, by = Analysis_Date], inherit.aes = F, angle=90, size=3, hjust=1, vjust=0.5, aes(x=Analysis_Date, y=1, label=N)) +
  scale_fill_brewer(palette="Paired") +
  theme_cowplot() +
  theme(legend.position = "top")

p.pcr.timeseries.reg <-
ggplot(oxf.clin[PCR_Outcome == "POSITIVE"], aes(x=Analysis_Date, fill = OXF_Outcome)) + 
  geom_bar(position="fill") +
  geom_text(data = oxf.clin[PCR_Outcome == "POSITIVE", .N, by = c("Analysis_Date", "Region")], inherit.aes = F, angle=90, size=3, hjust=1, vjust=0.5, aes(x=Analysis_Date, y=1, label=N)) +
  scale_fill_brewer(palette="Paired") +
  facet_grid(Region~.) +
  theme_cowplot() +
  theme(strip.background = element_blank() )

ggsave2("plots/pcr_timeseries_sweden.pdf", width = 8, height = 3, p.pcr.timeseries.swe)
ggsave2("plots/pcr_timeseries_regions.pdf", width = 8, height = 6, p.pcr.timeseries.reg)

# agreement pcr/seq control experiments
p.pcr.seq.con.conf <-
ggplot(oxf.merge[grepl("CON", Sample_Name) & OXF_N1_Ct != "Undetermined" & !is.na(Lineage)], aes(x=OXF_OmiS_Ct != "Undetermined", y=Lineage_group == "Omicron (BA.1-like)")) + 
  geom_jitter(col="grey", stroke=NA, width = 0.3, height = 0.3) +
  geom_text(data = oxf.merge[grepl("CON", Sample_Name) & OXF_N1_Ct != "Undetermined" & !is.na(Lineage), melt(table(list(x = OXF_OmiS_Ct != "Undetermined", y = Lineage_group == "Omicron (BA.1-like)")))], aes(x=x, y=y, label=value)) +
  labs(x="BA.1 Positive (qPCR)", y="BA.1 Positive (sequencing)") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave2("plots/pcr_seq_control_confusionmatrix.pdf", width = 3, height = 3, p.pcr.seq.con.conf)

# agreement pcr/seq ALL samples
p.pcr.seq.agreement <-
ggplot(oxf.merge[OXF_N1_Ct != "Undetermined" & !is.na(Lineage) , list(Lineage_group == "Omicron (BA.1-like)" & OXF_OmiS_Ct != "Undetermined" | Lineage_group != "Omicron (BA.1-like)" & OXF_OmiS_Ct == "Undetermined")],aes(y="", fill=V1)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  coord_polar() +
  scale_fill_brewer(palette="Set1", name = NULL, labels = c("Disagreement", "Agreement")) +
  theme_void() +
  theme(aspect.ratio = 1)

ggsave2("plots/pcr_seq_agreement.pdf", width = 5, height = 4, p.pcr.seq.agreement)

pangolin[Sample %in% oxf.merge[grepl("CLI", Sample_Name) & OXF_N1_Ct != "Undetermined" & !is.na(Lineage) & !(Lineage_group == "Omicron (BA.1-like)" & OXF_OmiS_Ct != "Undetermined" | Lineage_group != "Omicron (BA.1-like)" & OXF_OmiS_Ct == "Undetermined"), Sample_Name]] # Disagreeing samples

# classified lineage clinical samples only
p.pcr.lineage <-
ggplot(oxf.merge[grepl("CLI", Sample_Name) & OXF_N1_Ct != "Undetermined" & !is.na(Lineage)], aes(x = OXF_OmiS_Ct != "Undetermined", fill = Lineage )) + 
  geom_bar() +
  stat_count(geom = "text", aes(label=..count..)) +
  labs(y="Count", x="Omicron (BA.1)") +
  scale_fill_brewer(palette="Pastel2") +
  theme_cowplot()

ggsave2("plots/pcr_lineage_assignment_clinical.pdf", width = 3, height = 3, p.pcr.lineage)

# ct differences
library(ggbeeswarm)
library(ggpubr)
oxf.merge[OXF_N1_Ct != "Undetermined", wilcox.test(as.numeric(OXF_N1_Ct)~OXF_OmiS_Ct != "Undetermined")]$p.value # p-value
oxf.merge[OXF_N1_Ct != "Undetermined", median_cl_boot(as.numeric(OXF_N1_Ct)), by = OXF_OmiS_Ct != "Undetermined"] # c.i

p.pcr.n.all <-
ggplot(oxf.merge[OXF_N1_Ct != "Undetermined"], aes(x=OXF_OmiS_Ct != "Undetermined", y = as.numeric(OXF_N1_Ct))) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.25, show.legend = F, aes(fill=OXF_OmiS_Ct != "Undetermined")) +
  stat_summary(fun = "median", geom="text", aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-15)), geom="text" ) +
  labs(y = "N gene Ct", x = "Omicron (BA.1)") +
  coord_flip(ylim=c(30,20)) +
  scale_y_reverse() +
  scale_fill_brewer(palette="Paired") +
  theme_cowplot()

p.pcr.rnasep.all <-
ggplot(oxf.merge[OXF_N1_Ct != "Undetermined"], aes(x=OXF_OmiS_Ct != "Undetermined", y = as.numeric(OXF_RNaseP_Ct))) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.25, show.legend = F, aes(fill=OXF_OmiS_Ct != "Undetermined")) +
  stat_summary(fun = "median", geom="text", aes(label = round(-..y.., 2) )) +
  labs(y = "RNaseP Ct", x = "Omicron (BA.1)") +
  coord_flip() +
  scale_y_reverse() +
  scale_fill_brewer(palette="Paired") +
  theme_cowplot()

p.pcr.n.lin <-
ggplot(oxf.merge[OXF_N1_Ct != "Undetermined" & N > 10 & !is.na(Lineage)], aes(x=Lineage_group, y = as.numeric(OXF_N1_Ct))) +
  geom_quasirandom(stroke=NA, alpha=0.25) +
  geom_boxplot(width=0.25, show.legend = F, outlier.shape = NA, aes(fill=Lineage_group)) +
  stat_summary(fun = "median", geom="text", aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-16)), geom="text" ) +
  stat_compare_means(comparisons = list(c("Delta (B.1.617.2-like)", "Omicron (BA.1-like)"), c("Delta (B.1.617.2-like)", "Omicron (BA.2-like)"), c("Omicron (BA.2-like)", "Omicron (BA.1-like)") ), tip.length = 0, step.increase = 0.025 ) +
  scale_fill_brewer(palette="Pastel1") +
  scale_y_reverse() +
  labs(y = "N gene Ct", x = "Assigned lineage") +
  theme_cowplot()

ggsave2("plots/pcr_n_fulldata.pdf", width = 5, height = 2, p.pcr.n.all)
ggsave2("plots/pcr_rnasep_fulldata.pdf", width = 5, height = 2, p.pcr.rnasep.all)
ggsave2("plots/pcr_n_lineage.pdf", width = 3, height = 5, p.pcr.n.lin)

# ct over time
library(quantreg)
oxf.clin[PCR_Outcome == "POSITIVE", summary(rq(as.numeric(OXF_N1_Ct)~Analysis_Date+OXF_Outcome ))] # quantile regression

p.pcr.timeseries.n <-
ggplot(oxf.clin[PCR_Outcome == "POSITIVE"], aes(y = as.numeric(OXF_N1_Ct), x = as.factor(Analysis_Date), fill=OXF_Outcome)) +
  geom_boxplot(show.legend = F) +
  stat_summary(fun = "median", geom="hline", show.legend = F, lty=2, aes(x=1, yintercept=-..y.., col=OXF_Outcome)) +
  labs(y = "N gene Ct") +
  scale_fill_brewer(palette="Paired", aesthetics = c("fill", "col")) +
  scale_y_reverse() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), strip.background = element_blank() )

p.pcr.timeseries.dens <-
ggplot(oxf.clin[PCR_Outcome == "POSITIVE"], aes(x = as.numeric(OXF_N1_Ct), col=OXF_Outcome )) +
  geom_density() +
  labs(x = "N gene Ct") +
  coord_flip() +
  scale_color_brewer(palette="Paired") +
  scale_x_reverse() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), strip.background = element_blank() )

ggsave2("plots/pcr_timeseries_control.pdf", width = 12, height = 3,
plot_grid(align = "hv", axis = "trbl", rel_widths = c(1,0.5),
  p.pcr.timeseries.n, p.pcr.timeseries.dens
))

# ct per region
oxf.clin[PCR_Outcome == "POSITIVE" & !grepl("^AAA_", Region), median(as.numeric(OXF_N1_Ct)[OXF_Outcome == "Omicron (BA.1) POSITIVE"], na.rm = T) - median(as.numeric(OXF_N1_Ct)[OXF_Outcome == "Omicron (BA.1) NEGATIVE"], na.rm = T), by="Region"][,c("avg" = mean(V1), "sd" = sd(V1), "range" = range(V1) )] # difference per region

p.pcr.region.n <-
ggplot(oxf.clin[PCR_Outcome == "POSITIVE" & !grepl("^AAA_", Region)], aes(y = as.numeric(OXF_N1_Ct), x="", fill = OXF_Outcome)) +
  geom_boxplot() +
  stat_summary(fun = "median", geom="text", show.legend = F, position = position_dodge(0.75), aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-16)), geom="text", aes(col=OXF_Outcome) ) +
  stat_compare_means(tip.length = 0, aes(label = signif(..p.., 2) ) ) +
  facet_wrap(~Region, nrow=2) +
  labs(y = "N gene Ct") +
  #coord_cartesian(ylim =c(32,18)) +
  scale_fill_brewer(palette="Paired", aesthetics = c("fill", "col")) +
  scale_y_reverse() +
  theme_cowplot() +
  theme(strip.background = element_blank(), panel.background = element_rect(colour="black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave2("plots/pcr_region_control.pdf", width = 6, height = 3, p.pcr.region.n)

ggplot(oxf.clin[PCR_Outcome == "POSITIVE" & !grepl("^AAA_", Region)], aes(y = as.numeric(OXF_N1_Ct), x=Analysis_Date > as.IDate("2022-02-09"), col = OXF_Outcome)) +
  geom_boxplot(show.legend = T) +
  stat_summary(fun = "median", geom="text", show.legend = F, position = position_dodge(0.75), aes(label = round(-..y.., 2) )) +
  #stat_compare_means(tip.length = 0, step.increase = 0.025 ) +
  facet_grid(~Region, scales="free_x") +
  labs(y = "N gene Ct") +
  scale_fill_brewer(palette="Paired", aesthetics = c("fill", "col")) +
  scale_x_discrete(drop=FALSE) +
  scale_y_reverse() +
  theme_cowplot() +
  theme(strip.background = element_blank())