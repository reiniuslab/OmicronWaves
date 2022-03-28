library(data.table)

# qPCR data
oxf.clin <- fread("data_decoded/OXF_20220308_DECODED.csv")
oxf.clin[, Analysis_Date := as.Date(Analysis_Date, tryFormats = c("%d/%m/%Y"))]

# nextclade mutations
nextclade <- fread("data_decoded/Nextclade_results_denovo_DECODED.csv")
nextclade[, sample := gsub("/.*","",seqName)]
nextclade[, H78Y := grepl("ORF3a:H78Y", aaSubstitutions)]

# pangolin lineage classification
pangolin <- fread("data_decoded/Pangolin_results_denovo_DECODED.csv")
pangolin[, Sample := gsub("/.*","",`Sequence name`)]

# add data
oxf.clin[nextclade, c("Clade", "H78Y") := list(clade, H78Y), on = "Sample_Name == sample"]
oxf.clin[pangolin, c("Lineage", "Lineage_group") := list(Lineage, `Scorpio call`), on = "Sample_Name == Sample"]

## PLOT
library(ggplot2)
library(cowplot)

# ORF3a:H78Y status
p.h78y.cnt <-
ggplot(oxf.clin[!is.na(H78Y) & !is.na(Lineage) & Lineage != "None"], aes(x=Lineage_group, fill=H78Y)) +
  geom_bar() +
  stat_count(geom="text", aes(label=..count..)) +
  labs(y="Count", x="Lineage") +
  scale_fill_brewer(palette="Set1") +
  theme_cowplot()

ggsave2("plots/h78y_count.pdf", width = 3, height = 3, p.h78y.cnt)

# by region
#oxf.clin[, Month := month(Analysis_Date)]

p.h78y.region <- 
ggplot(oxf.clin[!is.na(H78Y) & !is.na(Lineage) & Lineage == "BA.2", list("H78Y" = mean(H78Y), N = .N ), by=c("Region")], aes(y=Region, x="", fill=H78Y)) + 
  geom_tile() +
  geom_text(aes(label=N)) +
  scale_fill_viridis_c() +
  theme_cowplot()

ggsave2("plots/h78y_region.pdf", width = 4, height = 4, p.h78y.region)