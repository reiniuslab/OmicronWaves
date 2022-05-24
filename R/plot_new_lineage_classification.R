library(data.table)

# data
pangolin.old <- fread("data_decoded/Pangolin_results_denovo_DECODED.csv")
pangolin.new <- fread("data_decoded/Pangolin_results_denovo_20220517_DECODED.csv")

pangolin.merge <- pangolin.new[, list("Name" = `Sequence name`, "Lineage.New" = Lineage, "Scorpio.New" = `Scorpio call`)][pangolin.old[, list("Name" = gsub("/.*","",`Sequence name`), "Lineage.Old" = Lineage, "Scorpio.Old" = `Scorpio call`)], on = "Name == Name"]
pangolin.merge[, N := .N, by="Lineage.Old"]

# qPCR data
oxf.clin <- fread("data_decoded/OXF_20220308_DECODED.csv")
oxf.clin[, Analysis_Date := as.Date(Analysis_Date, tryFormats = c("%d/%m/%Y"))]
oxf.clin[pangolin.new, c("Lineage_new","Lineage_group") := list(Lineage, `Scorpio call`), on = "Sample_Name == `Sequence name`"]
oxf.clin[, N_lineage := .N, by=Lineage_new]

## PLOT
library(ggplot2)
library(cowplot)
p.pang.clade <-
ggplot(pangolin.merge, aes(x=gsub(" ","\n",Scorpio.Old), y=gsub(" ","\n",Scorpio.New) )) + 
  geom_jitter(width = 0.3, height = 0.3) +
  geom_hline(yintercept = seq(1.5,7.5,1)) +
  geom_vline(xintercept = seq(1.5,6.5,1)) +
  labs(x="Pangolin 3.1.20 (2022-03-08)", y="Pangolin 4.0.6 (2022-05-17)") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

p.pang.lineage <-
ggplot(pangolin.merge[N>10], aes(x=Lineage.Old, y=factor(Lineage.New, levels=paste0("BA.",c("1", "1.13", "1.14", "1.15", "1.16", "1.17", "1.17.2", "1.18", "1.20", "1.21", "1.1", "1.1.1", "1.1.4", "1.1.8", "1.1.11", "1.1.18", "2", "2.6", "2.9", "2.25"))) )) + 
  geom_jitter(width = 0.3, height = 0.3) +
  geom_hline(yintercept = seq(1.5,20.5,1)) +
  geom_vline(xintercept = seq(1.5,3.5,1)) +
  labs(x="Pangolin 3.1.20 (2022-03-08)", y="Pangolin 4.0.6 (2022-05-17)") +
  theme_cowplot() +
  theme(aspect.ratio = 2.5)

ggsave2("plots/lineage_update_confusion_matrix.pdf",width = 9, height = 4,
plot_grid(rel_widths = c(1,0.5),
  p.pang.clade, p.pang.lineage
))

##
library(ggbeeswarm)
library(ggpubr)
p.pang.n <-
ggplot(oxf.clin[grepl("BA.2", Lineage_new) & N_lineage > 10], aes(x=Lineage_new, y=as.numeric(OXF_N1_Ct))) +
  geom_quasirandom(groupOnX = T, width = 0.25) +
  geom_boxplot() +
  #stat_summary(fun.data="median_cl_boot", col="red", position = position_nudge(0.4)) +
  stat_summary(fun = "median", geom="text", position=position_nudge(-0.45), angle=90, aes(label = round(-..y.., 2) )) +
  stat_summary(fun.data = function(x) return(data.frame(label=paste0("n = ",length(x)), y=-32)), geom="text" ) +
  stat_compare_means(aes(label=paste("P (MWU) =",..p.format..) )) +
  scale_y_reverse() +
  labs(x="Pangolin 4.0.6 (2022-05-17)", y="N1 Ct") +
  theme_cowplot()

ggsave2("plots/lineage_update_ba2_n_ct.pdf", width = 3, height = 4, p.pang.n)

##

p.pang.region <-
ggplot(oxf.clin[!is.na(Lineage_new) & !grepl("^AAA", Region), list("BA.2.9" = mean(Lineage_new == "BA.2.9"), "N" = paste0(sum(Lineage_new == "BA.2.9"),"/",.N) ), by="Region"], aes(y=Region, x="", fill=BA.2.9, label=N)) +
  geom_tile() +
  geom_text() +
  scale_fill_viridis_c() +
  expand_limits(fill=0) +
  theme_void() +
  theme(axis.text.y = element_text())

ggsave2("plots/lineage_update_region.pdf", width = 3, height = 3, p.pang.region)