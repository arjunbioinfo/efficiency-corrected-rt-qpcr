# Efficiency-corrected RT-qPCR analysis: CCAT1 intron1 normalized to MYC intron1
#
# DUMMY / SYNTHETIC DATA. Every Cq value in data/raw is randomly generated and
# does not correspond to any real experiment. The script exists to demonstrate
# the calculation, not to report a biological result.
#
# PPT / Pfaffl-style formulas implemented:
#   Efficiency (%)   = (-1 + 10^(-1 / slope)) * 100
#   E                = 1 + Efficiency(%) / 100
#   Ct(100%)         = Cq * log2(E)
#   Quantity         = 100 * 2^( mean(std100 Ct(100%)) - sample Ct(100%) )
#   Normalized ratio = CCAT1 quantity / MYC quantity   (per matched technical replicate)
#
# The two technical ratios are averaged per biological replicate, then the group
# mean and SD are taken across the three biological replicates per confluency.

suppressWarnings(suppressMessages(library(readxl)))

root       <- getwd()
raw_path   <- file.path(root, "data", "raw", "ccat1_myc_intron1_raw_ct.xlsx")
out_dir    <- file.path(root, "results", "generated")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- as.data.frame(read_excel(raw_path))

# ---- efficiency per target (taken from the recorded standard-curve fit) ----
eff <- tapply(raw$EfficiencyPercent, raw$Target, function(x) x[1])
log2E <- function(target) log2(1 + eff[[target]] / 100)

# ---- Ct(100%) and quantity relative to the 100-copy standard ----
mean_std100 <- sapply(names(eff), function(tg) {
  s <- raw[raw$Target == tg & raw$Role == "standard" & raw$StandardQuantity == 100, ]
  mean(s$Cq)
})

quantity_of <- function(target, cq) {
  l2 <- log2E(target)
  100 * 2 ^ ((mean_std100[[target]] * l2) - (cq * l2))
}

samp <- raw[raw$Role == "sample", ]
samp$Ct100    <- mapply(function(tg, cq) cq * log2E(tg), samp$Target, samp$Cq)
samp$quantity <- mapply(quantity_of, samp$Target, samp$Cq)

# ---- sample quantities ----
sq <- samp[, c("Target", "BiologicalReplicate", "Confluency",
               "TechnicalReplicate", "Cq", "quantity")]
sq$Cq <- round(sq$Cq, 4); sq$quantity <- round(sq$quantity, 3)
write.csv(sq, file.path(out_dir, "sample_quantities.csv"), row.names = FALSE)

# ---- technical-replicate mean Cq per target/rep/confluency ----
trm <- aggregate(Cq ~ Target + BiologicalReplicate + Confluency, data = samp, FUN = mean)
names(trm)[names(trm) == "Cq"] <- "mean_Cq"
trm$mean_Cq <- round(trm$mean_Cq, 4)
write.csv(trm, file.path(out_dir, "technical_replicate_means.csv"), row.names = FALSE)

# ---- CCAT1/MYC ratio per matched technical replicate ----
ccat <- samp[samp$Target == "ccat1 intron1", ]
myc  <- samp[samp$Target == "MYC intron1", ]
key  <- function(d) paste(d$BiologicalReplicate, d$Confluency, d$TechnicalReplicate)
m    <- merge(ccat, myc, by = c("BiologicalReplicate", "Confluency", "TechnicalReplicate"),
             suffixes = c("_ccat1", "_myc"))
m$ccat1_myc_ratio <- m$quantity_ccat1 / m$quantity_myc
trr <- m[, c("BiologicalReplicate", "Confluency", "TechnicalReplicate", "ccat1_myc_ratio")]
trr <- trr[order(trr$BiologicalReplicate, trr$Confluency, trr$TechnicalReplicate), ]
trr$ccat1_myc_ratio <- round(trr$ccat1_myc_ratio, 6)
write.csv(trr, file.path(out_dir, "technical_replicate_ratios.csv"), row.names = FALSE)

# ---- normalized ratio per biological replicate (mean of technical ratios) ----
nr <- aggregate(ccat1_myc_ratio ~ BiologicalReplicate + Confluency, data = m, FUN = mean)
names(nr)[names(nr) == "ccat1_myc_ratio"] <- "mean_ccat1_myc_ratio"
nr <- nr[order(nr$BiologicalReplicate, nr$Confluency), ]
nr$mean_ccat1_myc_ratio <- round(nr$mean_ccat1_myc_ratio, 6)
write.csv(nr, file.path(out_dir, "normalized_ratios.csv"), row.names = FALSE)

# ---- group summary per confluency across biological replicates ----
confs <- sort(unique(nr$Confluency))
group <- do.call(rbind, lapply(confs, function(cf) {
  v <- nr$mean_ccat1_myc_ratio[nr$Confluency == cf]
  v <- v[order(unique(nr$BiologicalReplicate))]
  data.frame(confluency = cf,
             R1_ratio = round(v[1], 6), R2_ratio = round(v[2], 6), R3_ratio = round(v[3], 6),
             mean_ratio = round(mean(v), 6), sd_ratio = round(sd(v), 6))
}))
write.csv(group, file.path(out_dir, "group_summary.csv"), row.names = FALSE)

# ---- verify against the committed snapshot ----
ref_path <- file.path(root, "results", "ppt_formula_corrected_workbook_summary.csv")
if (file.exists(ref_path)) {
  ref <- read.csv(ref_path)
  ok  <- isTRUE(all.equal(group$mean_ratio, ref$mean_ratio, tolerance = 1e-4))
  cat(sprintf("Group summary matches committed snapshot: %s\n", ok))
}
print(group)
cat("\nWrote CSVs to", out_dir, "\n")
