library(data.table)

##

splitgrep <- function(x, split, pattern, ... ){
  y <- unlist(strsplit(x, split, ... ))
  idx <- grepl(pattern, y)
  y[idx]
}

##

## nextclade mutations
nextclade <- fread("data_decoded/Nextclade_results_denovo_DECODED.csv")
nextclade[, sample := gsub("/.*","",seqName)]

## table of mutations from Colson et al. Table1 (doi:10.1002/jmv.27789)
colson <- fread("Colson_Table1.txt")

## all mutations
# remove samples with erroneous frameshifts, format as TRUE/FALSE table and remove mutations only found in 10 sample
dat.all <- nextclade[qc.frameShifts.status != "bad" & clade != "20B", strsplit(paste(substitutions, deletions, insertions, sep = ","), ","), by=c("sample", "clade", "Nextclade_pango")]
dat.all.cast <- dcast(dat.all, sample+clade+Nextclade_pango ~ V1, value.var = "V1")
dat.all.cast <- dat.all.cast[,lapply(.SD, function(x) !is.na(x) ), by=c("sample", "clade", "Nextclade_pango")]

dat.all.filt <- dat.all.cast[, c(T,T,T,F, apply(dat.all.cast[,-1:-4],2,sum)>10), with=F]

# melt data and add positions (for indel ranges, just selects the first position) 
dat.all.melt <- melt(dat.all.filt, id.vars = c("sample", "clade", "Nextclade_pango"))
dat.all.melt[, position := as.numeric(gsub("[[:alpha:]]|[[:punct:]]|-.*","", variable))]
dat.all.melt[value == T, call := variable]

# PCA
pca.all <- prcomp(dat.all.filt[,-1:-3], scale. = T)

## Same but for S mutations only, threshold = mutation found in > 1 samples
dat.s <- nextclade[qc.frameShifts.status != "bad", splitgrep(paste(aaSubstitutions, aaDeletions, aaInsertions, sep = ","), ",", "^S"), by=c("sample", "clade", "Nextclade_pango")]
dat.s.cast <- dcast(dat.s, sample+clade+Nextclade_pango ~ V1, value.var = "V1")
dat.s.cast <- dat.s.cast[,lapply(.SD, function(x) !is.na(x) ), by=c("sample", "clade", "Nextclade_pango")]

dat.s.filt <- dat.s.cast[, c(T,T,T, apply(dat.s.cast[,-1:-3],2,sum)>1), with=F]

# melt and add pos
dat.s.melt <- melt(dat.s.filt, id.vars = c("sample", "clade", "Nextclade_pango"))
dat.s.melt[, position := as.numeric(gsub("[[:alpha:]]|[[:punct:]]","", variable))]
dat.s.melt[value == T, call := variable]

# pca
pca.s <- prcomp(dat.s.filt[,-1:-3], scale. = T)


## PLOT
library(ggplot2)
library(cowplot)

# mutation heatmap
p.mut.heat.all <-
ggplot(dat.all.melt, aes(x=position, y=sample, fill=call)) +
  geom_tile(show.legend = F, width=150) +
  facet_grid(clade~., scales="free", space="free") +
  ggtitle("All mutations") +
  scale_fill_hue(na.value = NA) +
  theme_void() +
  theme(panel.border = element_rect(color="black", fill=NA), panel.spacing = unit(0, "cm"))

p.mut.heat.s <-
ggplot(dat.s.melt, aes(x=position, y=sample, fill=call)) +
  geom_tile(show.legend = F, width=10) +
  facet_grid(clade~., scales="free", space="free") +
  ggtitle("S mutations") +
  scale_fill_hue(na.value = NA) +
  theme_void() +
  theme(panel.border = element_rect(color="black", fill=NA), panel.spacing = unit(0, "cm"))

ggsave2("plots/mutation_heatmap_all.pdf", width = 12, height = 6, p.mut.heat.all)
ggsave2("plots/mutation_heatmap_s.pdf", width = 12, height = 6, p.mut.heat.s)

# pca
library(ggfortify)
p.mut.pca.all <- 
  autoplot(pca.all, data=dat.all.filt[,1:3], colour="clade" ) + ggtitle("All mutations") + theme_cowplot() + theme(aspect.ratio = 1)
p.mut.pca.s <- 
  autoplot(pca.s, data=dat.s.filt[,1:3], colour="clade" ) + ggtitle("S mutations") + theme_cowplot() + theme(aspect.ratio = 1)

ggsave2("plots/mutation_pca_all.pdf", width = 7, height = 6, p.mut.pca.all)
ggsave2("plots/mutation_pca_s.pdf", width = 7, height = 6, p.mut.pca.s)

# specific mutation patterns from Colson et al.
p.mut.colson <-
ggplot(melt(dat.all.filt[grepl("21I|21J|21K|21L", clade),c(1:3, which(colnames(dat.all.filt) %in% colson$Nucleotide_changes)), with=F], id.vars = c("sample", "clade", "Nextclade_pango")), aes(y=reorder(variable, -as.numeric(gsub("[[:alpha:]]","",variable))), x=sample, fill=value)) + 
  geom_tile() + 
  facet_grid(~clade, scales="free") +
  labs(y="Mutation spectrum from Colson et al.") +
  scale_fill_grey(start = 0.8, end = 0.2, name="Mutation") +
  theme_void() +
  theme(axis.text.y = element_text(), axis.title.y = element_text(angle=90), panel.border = element_rect(fill=NA))

p.mut.colson.exp <- 
ggplot(melt(colson[Nucleotide_changes %in% colnames(dat.all.filt)], id.vars=c("Region", "Nucleotide_changes", "AA_changes")), aes(x=variable, y=reorder(Nucleotide_changes , -as.numeric(gsub("[[:alpha:]]","",Nucleotide_changes ))), fill=value == "X" ))+
  geom_tile(show.legend = F) +
  facet_grid(~variable, scales="free") +
  scale_fill_grey(start = 0.8, end = 0.2, name="Mutation") +
  theme_void() +
  theme(panel.border = element_rect(fill=NA))

ggsave2("plots/mutation_heatmap_colson.pdf", width = 5, height = 3,
plot_grid(rel_widths = c(1,0.25),
  p.mut.colson, p.mut.colson.exp
))