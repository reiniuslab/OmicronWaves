library(data.table)

# qPCR results
oxf.clin <- fread("data_decoded/OXF_20220308_DECODED.csv")
oxf.cont <- fread("data_decoded/Control_samples_BXF_OXF_DECODED.csv")
oxf.merge <- rbind(oxf.cont, oxf.clin, fill=T)

# pangolin results
pangolin <- fread("data_decoded/Pangolin_results_denovo_DECODED.csv")
pangolin[, Sample := gsub("/.*", "", `Sequence name`)]

# n50
dat.n50 <- fread("denovo/n50.txt", sep = ":", header = F)
dat.n50[, c("Sample_Name", "N50") := list(gsub("\\..*", "", basename(V1) ), as.numeric(sapply(strsplit(V2, " "), "[[", 2)) )]

# coverage
fls.cov <- list.files("alignment_stats", "coverage_.*.txt", full.names = T)
names(fls.cov) <- gsub("coverage_|\\.txt", "", basename(fls.cov) )

dat.cov <- rbindlist(lapply(fls.cov, fread), idcol = T)

# minimap2 qc results (flagstat)
fls.stats <- list.files("alignment_stats", "flagstat_.*.txt", full.names = T)
names(fls.stats) <- gsub("flagstat_|\\.txt", "", basename(fls.stats) )

dat.stats <- rbindlist(lapply(fls.stats, fread, sep = "", header = F), idcol = T)
dat.stats[, type := gsub("[[:digit:]]* \\+ [[:digit:]]* | \\(.*", "", V1)]

dat.stats <- dat.stats[type %in% c("in total", "mapped"), 
         list("sample" = .id[type == "in total"],
              "total_reads" = gsub(" .*", "", V1[type == "in total"]),
              "mapped_reads" = gsub(" .*", "", V1[type == "mapped"]),
              "percentage_mapped" = gsub(".*\\(|%.*", "", V1[type == "mapped"]) 
              )]

# merge data
dat.stats[dat.cov, c("coverage", "depth") := list(coverage, meandepth), on = "sample == .id"]
dat.stats[dat.n50, "n50" := N50, on = "sample == Sample_Name"]
dat.stats[oxf.merge, OXF_N1_Ct := as.numeric(OXF_N1_Ct), on = "sample == Sample_Name"]
dat.stats[, PCR := !is.na(OXF_N1_Ct)]
dat.stats[pangolin, "lineage" := Lineage, on = "sample == Sample"]
dat.stats[, sequencing := !is.na(lineage)]
dat.stats[grepl("NEGATIVE_CONTROL", sample), c("sequencing", "PCR") := list(FALSE, FALSE)]

dat.stats.melt <- melt(dat.stats, id.vars = c("sample", "sequencing", "PCR", "lineage"))
dat.stats.melt[, value := as.numeric(value)]
dat.stats.melt[is.na(value) & grepl("percentage", variable), value := 0]

## PLOT
library(ggplot2)
library(cowplot)

#
library(ggbeeswarm)
p.seq.stats <-
ggplot(dat.stats.melt, aes(y=value, x= paste(PCR, sequencing) )) +
  geom_quasirandom(stroke=NA, alpha=0.25, show.legend = F, aes(col = paste(PCR, sequencing) )) +
  geom_boxplot(outlier.shape = NA, width=0.25, aes(fill = paste(PCR, sequencing) )) +
  facet_grid(variable~., scales = "free", switch = "y") +
  labs(x=NULL, y=NULL) +
  scale_color_brewer(palette="Pastel1", name = "PCR Seq", aesthetics = c("fill", "colour")) +
  theme_cowplot() +
  theme(strip.placement = "outside", strip.background = element_blank(), panel.border = element_rect(colour="black") )

ggsave2("plots/seq_stats_boxplot.pdf", width = 4, height = 12, p.seq.stats)

# pca
library(ggfortify)

pca.stats <- prcomp(data.frame(dcast(dat.stats.melt, sample~variable, value.var = "value"), row.names = 1), scale. = T)

autoplot(pca.stats, data=dat.stats[, list("SequencingPCR" = paste(sequencing, PCR))], colour="SequencingPCR")