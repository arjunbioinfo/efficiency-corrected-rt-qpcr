#!/usr/bin/env python3
"""Seeded generator for the synthetic RT-qPCR dataset (method demo).

Self-contained and path-independent: writes the raw Cq workbook, the committed
summary CSV, and the generated result CSVs relative to the repo root.
Requires: numpy, openpyxl.  Run:  python3 tools/generate_dummy_data.py
"""
import math, os, csv
import numpy as np
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(20260705)

SLOPE_T, SLOPE_R = -3.45, -3.35
eff = lambda s: (10 ** (-1.0 / s) - 1.0) * 100.0
EFF_T, EFF_R = eff(SLOPE_T), eff(SLOPE_R)
L2T, L2R = math.log2(1 + EFF_T / 100), math.log2(1 + EFF_R / 100)
GROUPS, REPS = ["A", "B", "C"], [1, 2, 3]

def std_curve(intercept, slope):
    return {q: [round(intercept + slope * math.log10(q) + rng.normal(0, 0.03), 4) for _ in range(2)]
            for q in (100, 10, 1)}
STD_T, STD_R = std_curve(30.4, SLOPE_T), std_curve(33.1, SLOPE_R)
MST, MSR = float(np.mean(STD_T[100])), float(np.mean(STD_R[100]))

def samp(center):
    d = {}
    for rep in REPS:
        for g in GROUPS:
            gm = rng.normal(center, 0.45)
            d[(rep, g)] = [round(gm + rng.normal(0, 0.03), 4) for _ in range(2)]
    return d
CT_T, CT_R = samp(24.6), samp(26.3)

def r2(std, slope, intercept):
    xs = [math.log10(q) for q in std for _ in std[q]]; ys = [v for q in std for v in std[q]]
    xs, ys = np.array(xs), np.array(ys); pred = intercept + slope * xs
    return round(float(1 - np.sum((ys - pred) ** 2) / np.sum((ys - np.mean(ys)) ** 2)), 6)
R2T, R2R = r2(STD_T, SLOPE_T, 30.4), r2(STD_R, SLOPE_R, 33.1)
qty = lambda ct, mstd, l2: 100.0 * 2 ** ((mstd * l2) - (ct * l2))

os.makedirs(os.path.join(ROOT, "data", "raw"), exist_ok=True)
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "raw_ct"
ws.append(["Well", "Assay", "Sample", "Group", "BiologicalReplicate",
           "TechnicalReplicate", "Role", "Cq", "StandardQuantity",
           "SlopeOrNA", "EfficiencyPercent", "R2"])
n = 0
def add(*row):
    global n; n += 1; ws.append([f"W{n}"] + list(row))
for q in (100, 10, 1):
    for i, ct in enumerate(STD_T[q]):
        add("target", f"std {q}", "", "", i + 1, "standard", ct, q, SLOPE_T, EFF_T, R2T)
for q in (100, 10, 1):
    for i, ct in enumerate(STD_R[q]):
        add("reference", f"std {q}", "", "", i + 1, "standard", ct, q, SLOPE_R, EFF_R, R2R)
for (rep, g), cts in sorted(CT_T.items()):
    for i, ct in enumerate(cts):
        add("target", f"Grp {g} Rep{rep}", g, rep, i + 1, "sample", ct, "", SLOPE_T, EFF_T, R2T)
for (rep, g), cts in sorted(CT_R.items()):
    for i, ct in enumerate(cts):
        add("reference", f"Grp {g} Rep{rep}", g, rep, i + 1, "sample", ct, "", SLOPE_R, EFF_R, R2R)
wb.save(os.path.join(ROOT, "data", "raw", "qpcr_raw_ct.xlsx"))

order = [(rep, g) for rep in REPS for g in GROUPS]
ratios = {}
for (rep, g) in order:
    cq = [qty(CT_T[(rep, g)][t], MST, L2T) for t in range(2)]
    rq = [qty(CT_R[(rep, g)][t], MSR, L2R) for t in range(2)]
    ratios[(rep, g)] = float(np.mean([cq[t] / rq[t] for t in range(2)]))
summ = []
for g in GROUPS:
    v = [ratios[(rep, g)] for rep in REPS]
    summ.append([g, round(v[0], 6), round(v[1], 6), round(v[2], 6),
                 round(float(np.mean(v)), 6), round(float(np.std(v, ddof=1)), 6)])
def wcsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
HED = ["group", "rep1_ratio", "rep2_ratio", "rep3_ratio", "mean_ratio", "sd_ratio"]
wcsv(os.path.join(ROOT, "results", "efficiency_corrected_ratio_summary.csv"), HED, summ)
g = os.path.join(ROOT, "results", "generated")
os.makedirs(g, exist_ok=True); open(os.path.join(g, ".gitkeep"), "w").close()
wcsv(os.path.join(g, "group_summary.csv"), HED, summ)

print(f"Efficiency: target={EFF_T:.2f}%  reference={EFF_R:.2f}%   R2: {R2T} / {R2R}")
for s in summ:
    print("  group %s  mean=%.6f  sd=%.6f" % (s[0], s[4], s[5]))
