#!/usr/bin/env python3
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from math import sqrt
from scipy.stats import norm

# ============================================================
# DeLong helpers
# DeLong Implementation (from Sun & Xu 2014 adaptation)
# ============================================================

def compute_midrank(x):
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T + 1
    return T2

def fastDeLong(preds_sorted, label_1_count):
    m = label_1_count
    n = preds_sorted.shape[1] - m

    Tx = np.apply_along_axis(compute_midrank, 1, preds_sorted)

    Ty = Tx[:, :m]
    Tz = Tx[:, m:]

    V10 = (Ty - Ty.mean(axis=1, keepdims=True)) / m
    V01 = (Tz - Tz.mean(axis=1, keepdims=True)) / n

    sx = np.var(V10.sum(axis=0), ddof=1)
    sy = np.var(V01.sum(axis=0), ddof=1)

    auc = (Tx[:, :m].mean() - (m+1)/2) / n
    return auc, sx/m + sy/n

def delong_auc_and_var(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    # correct ordering for DeLong: positives first
    order = np.argsort(-y_true)
    y_sorted = y_true[order]
    preds_sorted = y_pred[order]

    pos_count = np.sum(y_sorted == 1)
    preds_matrix = preds_sorted.reshape(1, -1)

    auc = roc_auc_score(y_true, y_pred)
    auc_delong, auc_var = fastDeLong(preds_matrix, pos_count)
    # ensure no division by zero
    if auc_var == 0:
        auc_var = 1e-16
    return auc, auc_var


# ============================================================
# Random-effects meta-analysis (DerSimonian–Laird)
# ============================================================

def random_effects_meta(effects, variances):
    effects = np.asarray(effects)
    variances = np.asarray(variances)
    weights = 1 / variances

    fixed = np.sum(weights * effects) / np.sum(weights)
    Q = np.sum(weights * (effects - fixed)**2)
    df = len(effects) - 1
    tau2 = max(0, (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights)))

    weights_re = 1 / (variances + tau2)
    pooled = np.sum(weights_re * effects) / np.sum(weights_re)
    pooled_var = 1 / np.sum(weights_re)

    return pooled, pooled_var, tau2


# ============================================================
# Load replicate files
# ============================================================

def load_file(path):
    df = pd.read_csv(path,sep="\t")
    # ground_truth    pred_prob
    if not {"ground_truth", "pred_prob"}.issubset(df.columns):
        raise ValueError(f"Input {path} must contain ground_truth and pred columns.")
    return df["ground_truth"].astype(int).values, df["pred_prob"].astype(float).values



# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputA", nargs="+", required=True)
    parser.add_argument("--inputB", nargs="+", required=True)
    args = parser.parse_args()

    aucA_list, varA_list = [], []
    aucB_list, varB_list = [], []

    print("\n=== GROUP A REPLICATES ===")
    for path in args.inputA:
        y, p = load_file(path)
        auc, var = delong_auc_and_var(y, p)
        aucA_list.append(auc)
        varA_list.append(var)
        print(f"{path}: AUROC={auc:.4f}, Var={var:.6f}")

    print("\n=== GROUP B REPLICATES ===")
    for path in args.inputB:
        y, p = load_file(path)
        auc, var = delong_auc_and_var(y, p)
        aucB_list.append(auc)
        varB_list.append(var)
        print(f"{path}: AUROC={auc:.4f}, Var={var:.6f}")

    # Meta-analysis within groups
    aucA_pool, varA_pool, tauA = random_effects_meta(aucA_list, varA_list)
    aucB_pool, varB_pool, tauB = random_effects_meta(aucB_list, varB_list)

    # Compare pooled AUROCs
    diff = aucA_pool - aucB_pool
    se_diff = sqrt(varA_pool + varB_pool)
    z = diff / se_diff
    p = 2 * (1 - norm.cdf(abs(z)))

    print("\n=== META-ANALYSIS SUMMARY ===")
    print(f"Pooled AUROC A = {aucA_pool:.4f}")
    print(f"Pooled AUROC B = {aucB_pool:.4f}")
    print(f"Difference A-B = {diff:.6f}")
    print(f"P-value = {p:.2e}")
    print(f"tau²_A = {tauA:.2e}, tau²_B = {tauB:.2e}")
    print("Lower bound estimate for 2-tailed T tests: 2e-33 (report 0.00+00 as <2.2e-16)")

if __name__ == "__main__":
    main()

