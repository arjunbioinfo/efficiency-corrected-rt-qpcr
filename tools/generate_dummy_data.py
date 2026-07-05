#!/usr/bin/env python3
"""Seeded generator for the synthetic RT-qPCR dataset.

Self-contained and path-independent: writes the raw Cq workbook, the committed
summary CSV, and the generated result CSVs relative to the repository root.
Requires: numpy, openpyxl.

Run from the repository root:
    python3 tools/generate_dummy_data.py
"""
import math, os, csv
import numpy as np
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 20260705
rng = np.random.default_rng(SEED)

# ---- standard-curve parameters (realistic-looking efficiencies) ----
SLOPE_C, SLOPE_M = -3.45, -3.35
eff = lambda s: (10 ** (-1.0 / s) - 1.0) * 100.0
EFF_C, EFF_M = eff(SLOPE_C), eff(SLOPE_M)
L2C, L2M = math.log2(1 + EFF_C / 100), math.log2(1 + EFF_M / 100)

CONF = {"10": 0.1, "50": 0.5, "100": 1.0}
CLS = ["10", "50", "100"]
BRS = [1, 2, 3]

def std_curve(intercept, slope):
    return {q: [round(intercept + slope * math.log10(q) + rng.normal(0, 0.03), 4)
                for _ in range(2)] for q in (100, 10, 1)}

STD_C, STD_M = std_curve(30.4, SLOPE_C), std_curve(33.1, SLOPE_M)
MSC, MSM = float(np.mean(STD_C[100])), float(np.mean(STD_M[100]))

def sample_cts(center):
    d = {}
    for b in BRS:
        for c in CLS:
            gm = rng.normal(center, 0.45)
            d[(b, c)] = [round(gm + rng.normal(0, 0.03), 4) for _ in range(2)]
    return d

CT_C, CT_M = sample_cts(24.6), sample_cts(26.3)

def r2(std, slope, intercept):
    xs = [math.log10(q) for q in std for _ in std[q]]
    ys = [v for q in std for v in std[q]]
    xs, ys = np.array(xs), np.array(ys)
    pred = intercept + slope * xs
    return 1 - np.sum((ys - pred) ** 2) / np.sum((ys - np.mean(ys)) ** 2)

R2_C, R2_M = round(float(r2(STD_C, SLOPE_C, 30.4)), 6), round(float(r2(STD_M, SLOPE_M, 33.1)), 6)
qty = lambda ct, mstd, l2: 100.0 * 2 ** ((mstd * l2) - (ct * l2))

# ---- raw workbook ----
os.makedirs(os.path.join(ROOT, "data", "raw"), exist_ok=True)
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "raw_ct"
ws.append(["Well", "Target", "Sample", "Confluency", "BiologicalReplicate",
           "TechnicalReplicate", "Role", "Cq", "StandardQuantity",
           "SlopeOrNA", "EfficiencyPercent", "R2"])
n = 0
def add(*row):
    global n; n += 1
    ws.append([f"W{n}"] + list(row))
for q in (100, 10, 1):
    for i, ct in enumerate(STD_C[q]):
        add("ccat1 intron1", f"std {q}", "", "", i + 1, "standard", ct, q, SLOPE_C, EFF_C, R2_C)
for q in (100, 10, 1):
    for i, ct in enumerate(STD_M[q]):
        add("MYC intron1", f"std {q}", "", "", i + 1, "standard", ct, q, SLOPE_M, EFF_M, R2_M)
for (b, c), cts in sorted(CT_C.items()):
    for i, ct in enumerate(cts):
        add("ccat1 intron1", f"BR{b} {c}", CONF[c], b, i + 1, "sample", ct, "", SLOPE_C, EFF_C, R2_C)
for (b, c), cts in sorted(CT_M.items()):
    for i, ct in enumerate(cts):
        add("MYC intron1", f"R{b} {c} WT", CONF[c], b, i + 1, "sample", ct, "", SLOPE_M, EFF_M, R2_M)
wb.save(os.path.join(ROOT, "data", "raw", "ccat1_myc_intron1_raw_ct.xlsx"))

# ---- summary + generated CSVs ----
order = [(b, c) for b in BRS for c in CLS]
ratios = {}
for (b, c) in order:
    cq = [qty(CT_C[(b, c)][t], MSC, L2C) for t in range(2)]
    mq = [qty(CT_M[(b, c)][t], MSM, L2M) for t in range(2)]
    ratios[(b, c)] = float(np.mean([cq[t] / mq[t] for t in range(2)]))

summ = []
for c in CLS:
    v = [ratios[(b, c)] for b in BRS]
    summ.append([CONF[c], round(v[0], 6), round(v[1], 6), round(v[2], 6),
                 round(float(np.mean(v)), 6), round(float(np.std(v, ddof=1)), 6)])

def wcsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

wcsv(os.path.join(ROOT, "results", "ppt_formula_corrected_workbook_summary.csv"),
     ["confluency", "R1_ratio", "R2_ratio", "R3_ratio", "mean_ratio", "sd_ratio"], summ)
g = os.path.join(ROOT, "results", "generated")
os.makedirs(g, exist_ok=True); open(os.path.join(g, ".gitkeep"), "w").close()
wcsv(os.path.join(g, "group_summary.csv"),
     ["confluency", "R1_ratio", "R2_ratio", "R3_ratio", "mean_ratio", "sd_ratio"], summ)

print(f"Efficiency: CCAT1={EFF_C:.2f}%  MYC={EFF_M:.2f}%   R2: {R2_C} / {R2_M}")
for s in summ:
    print("  conf %.1f  mean=%.6f  sd=%.6f" % (s[0], s[4], s[5]))
print("Wrote raw workbook + summary CSVs under", ROOT)
