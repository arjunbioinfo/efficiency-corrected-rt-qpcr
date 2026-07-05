# Efficiency correction in RT-qPCR — method note

A short, neutral reference for the calculation used in this repository. It
contains no experiment-, sample-, or lab-specific information.

## 1. Why efficiency matters

Relative quantification with the popular 2^(−ΔΔCt) method assumes every assay
amplifies with **100% efficiency** — the amount of product doubles each cycle.
Real primer/probe assays rarely reach exactly 100%, and a target assay and a
reference assay usually differ from each other. If those efficiencies are not
equal, an uncorrected ratio is systematically biased: the assay with the higher
efficiency is over-weighted. Efficiency correction (the Pfaffl approach) removes
that bias by rescaling every Cq onto a common 100%-efficiency basis before
quantities are compared.

## 2. Getting efficiency from a standard curve

Run a dilution series of a standard for each assay (e.g. 10-fold steps) in
duplicate, and fit Cq against log10(input):

```text
slope          = slope of Cq vs log10(quantity)
Efficiency (%) = (-1 + 10^(-1 / slope)) * 100
E              = 1 + Efficiency(%) / 100
```

Interpretation:

- A slope of **−3.32** corresponds to **100%** efficiency (E = 2.0).
- Efficiencies roughly within **90–110%** are generally considered acceptable.
- **R²** of the linear fit should be high (commonly ≥ 0.98); a poor fit means the
  slope — and therefore the efficiency — is unreliable.
- Exclude standard points that are undetermined or clearly off the line (often
  the most dilute point).

## 3. Efficiency-corrected quantification

For each measured Cq, convert to a 100%-efficiency-equivalent value and then to a
relative quantity against the highest standard:

```text
Ct(100%) = Cq * log2(E)
Quantity = 100 * 2^( mean(std100 Ct(100%)) - sample Ct(100%) )
```

The normalized readout is the ratio of the two assays’ quantities:

```text
Normalized ratio = Target quantity / Reference quantity
```

Compute the ratio per matched technical replicate, average technical replicates
within a biological replicate, then summarize (mean, SD) across biological
replicates within each group.

## 4. Controls worth including

- **No-RT control** — the reverse-transcriptase step omitted. Any signal here
  indicates amplification from genomic DNA rather than cDNA.
- **No-template control (NTC)** — water in place of template. Signal indicates
  contamination or primer-dimer.
- **Inter-plate calibrator** — a shared sample run on every plate if results are
  combined across plates.

## 5. Practical notes

- Reverse transcription is a **conversion**, not an amplification; a 1:1 RNA→cDNA
  conversion is assumed downstream and is rarely exact, so keep RT conditions
  consistent across samples.
- Assess RNA integrity (e.g. an RIN-type metric) before quantification;
  degradation is non-random and biases results.
- High-GC or structured templates can under-convert; a brief denaturation step
  and a strand-displacing reverse transcriptase help.

## Reference

Pfaffl, M.W. (2001). *A new mathematical model for relative quantification in
real-time RT-PCR.* Nucleic Acids Research 29(9):e45.
