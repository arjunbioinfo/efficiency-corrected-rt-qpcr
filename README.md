# Efficiency-corrected RT-qPCR analysis (CCAT1 intron1 / MYC intron1)

> **Synthetic demo repository.** Every number in this repository is randomly
> generated dummy data. It contains **no real experimental values** and reports
> **no biological result**. The purpose is to demonstrate an efficiency-corrected
> RT-qPCR calculation and a reproducible analysis workflow.

This repository reproduces an efficiency-corrected RT-qPCR analysis for
CCAT1 intron1 normalized to MYC intron1 in HCT116 WT samples across three cell
confluencies (0.1, 0.5, 1.0), with three biological replicates and two technical
replicates each.

## Result on the dummy data

| Cell confluency | R1 ratio | R2 ratio | R3 ratio | Mean CCAT1/MYC | SD |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0.446925 | 0.495239 | 0.436762 | 0.459642 | 0.031244 |
| 0.5 | 0.821351 | 0.893508 | 0.580211 | 0.765024 | 0.164068 |
| 1.0 | 0.221490 | 0.299057 | 0.675234 | 0.398594 | 0.242697 |

These values are meaningless (synthetic) and are shown only to confirm that the
workbook formulas, the R script, and the committed summary all agree.

## Calculation

The formulas follow `references/reverse_transcription_presentation.pptx`:

```text
Efficiency (%)   = (-1 + 10^(-1 / slope)) * 100
E                = 1 + Efficiency(%) / 100
Ct(100%)         = Cq * log2(E)
Quantity         = 100 * 2^( mean(std100 Ct(100%)) - sample Ct(100%) )
Normalized ratio = CCAT1 quantity / MYC quantity
```

CCAT1/MYC is computed for each matched technical-replicate pair. The two technical
ratios are averaged per biological replicate, then the group mean and SD are taken
across the three biological replicates for each confluency.

## Standard-curve quality (dummy)

| Target | R-squared | Efficiency | Assessment |
|---|---:|---:|---|
| CCAT1 intron1 | 0.999754 | 94.92% | Within the usual 90-110% range |
| MYC intron1 | 0.999960 | 98.84% | Within the usual 90-110% range |

The lowest standard (0.1) is undetermined for both assays and is not used in the fit.

## Reproduce the analysis

Requirements: R and the `readxl` package.

```r
install.packages("readxl")
```

Run from the repository root:

```bash
Rscript R/analyze_qpcr.R
```

The script reads `data/raw/ccat1_myc_intron1_raw_ct.xlsx`, writes CSV files to
`results/generated/`, and verifies the generated group summary against
`results/ppt_formula_corrected_workbook_summary.csv`.

## Regenerate the dummy dataset

The synthetic data is produced by a seeded generator, so the whole repository is
reproducible end to end:

```bash
python3 tools/generate_dummy_data.py
```

## Repository layout

```text
R/             Reproducible R analysis
tools/         Seeded dummy-data generator
data/raw/      Synthetic raw Cq workbook (input to the analysis)
references/    PPT method and Pfaffl template
results/       Live-formula workbook, static values snapshot, summary CSV
docs/          Static HTML report (optional GitHub Pages source)
```

The `*_live.xlsx` workbook keeps the original formulas so you can inspect and edit
the calculation in Excel. The `*_values.xlsx` workbook is the same sheet with every
formula evaluated to a static number. The summary CSV is the validated snapshot the
R script reproduces.

## GitHub Pages (optional)

A static report lives in `docs/`. Because it is dummy data, it is safe to publish.
To enable it, push the repository and open **Settings → Pages → Deploy from a branch**,
then choose `main` and the `/docs` folder. Preview locally with:

```bash
python3 -m http.server 8000   # then open http://localhost:8000/docs/
```
