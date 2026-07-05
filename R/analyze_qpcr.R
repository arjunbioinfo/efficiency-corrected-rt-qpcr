# Efficiency-corrected RT-qPCR: normalized target/reference ratio
#
# SYNTHETIC DATA. Every Cq value in data/raw is randomly generated and does not
# correspond to any real experiment. This demonstrates the calculation only.
#
# Pfaffl-style efficiency correction:
#   Efficiency (%)   = (-1 + 10^(-1 / slope)) * 100
#   E                = 1 + Efficiency(%) / 100
#   Ct(100%)         = Cq * log2(E)
#   Quantity         = 100 * 2^( mean(std100 Ct(100%)) - sample Ct(100%) )
#   Normalized ratio = Target quantity / Reference quantity   (per technical replicate)
#
# Technical ratios are averaged per biological replicate, then summarized
# (mean, SD) across replicates within each sample group.

suppressWarnings(suppressMessages(library(readxl)))

root     <- getwd()
raw_path <- file.path(root, "data", "raw", "qpcr_raw_ct.xlsx")
out_dir  <- file.path(root, "results", "generated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- as.data.frame(read_excel(raw_path))

# efficiency per assay (from the recorded standard-curve fit)
eff   <- tapply(raw$EfficiencyPercent, raw$Assay, function(x) x[1])
log2E <- function(assay) log2(1 + eff[[assay]] / 100)

mean_std100 <- sapply(names(eff), function(a) {
  s <- raw[raw$Assay == a & raw$Role == "standard" & raw$StandardQuantity == 100, ]
  mean(s$Cq)
})
quantity_of <- function(assay, cq) {
  l2 <- log2E(assay); 100 * 2 ^ ((mean_std100[[assay]] * l2) - (cq * l2))
}

samp <- raw[raw$Role == "sample", ]
samp$Ct100    <- mapply(function(a, cq) cq * log2E(a), samp$Assay, samp$Cq)
samp$quantity <- mapply(quantity_of, samp$Assay, samp$Cq)

sq <- samp[, c("Assay", "BiologicalReplicate", "Group", "TechnicalReplicate", "Cq", "quantity")]
sq$Cq <- round(sq$Cq, 4); sq$quantity <- round(sq$quantity, 3)
write.csv(sq, file.path(out_dir, "sample_quantities.csv"), row.names = FALSE)

trm <- aggregate(Cq ~ Assay + BiologicalReplicate + Group, data = samp, FUN = mean)
names(trm)[names(trm) == "Cq"] <- "mean_Cq"; trm$mean_Cq <- round(trm$mean_Cq, 4)
write.csv(trm, file.path(out_dir, "technical_replicate_means.csv"), row.names = FALSE)

tgt <- samp[samp$Assay == "target", ]
ref <- samp[samp$Assay == "reference", ]
m <- merge(tgt, ref, by = c("BiologicalReplicate", "Group", "TechnicalReplicate"),
           suffixes = c("_target", "_reference"))
m$target_reference_ratio <- m$quantity_target / m$quantity_reference
trr <- m[order(m$Group, m$BiologicalReplicate, m$TechnicalReplicate),
         c("Group", "BiologicalReplicate", "TechnicalReplicate", "target_reference_ratio")]
trr$target_reference_ratio <- round(trr$target_reference_ratio, 6)
write.csv(trr, file.path(out_dir, "technical_replicate_ratios.csv"), row.names = FALSE)

nr <- aggregate(target_reference_ratio ~ BiologicalReplicate + Group, data = m, FUN = mean)
names(nr)[names(nr) == "target_reference_ratio"] <- "mean_ratio"
nr <- nr[order(nr$Group, nr$BiologicalReplicate), ]; nr$mean_ratio <- round(nr$mean_ratio, 6)
write.csv(nr, file.path(out_dir, "normalized_ratios.csv"), row.names = FALSE)

groups <- sort(unique(nr$Group))
group <- do.call(rbind, lapply(groups, function(g) {
  v <- nr$mean_ratio[nr$Group == g]
  data.frame(group = g,
             rep1_ratio = round(v[1], 6), rep2_ratio = round(v[2], 6), rep3_ratio = round(v[3], 6),
             mean_ratio = round(mean(v), 6), sd_ratio = round(sd(v), 6))
}))
write.csv(group, file.path(out_dir, "group_summary.csv"), row.names = FALSE)

ref_path <- file.path(root, "results", "efficiency_corrected_ratio_summary.csv")
if (file.exists(ref_path)) {
  r <- read.csv(ref_path)
  cat(sprintf("Group summary matches committed snapshot: %s\n",
              isTRUE(all.equal(group$mean_ratio, r$mean_ratio, tolerance = 1e-4))))
}
print(group)
cat("\nWrote CSVs to", out_dir, "\n")
