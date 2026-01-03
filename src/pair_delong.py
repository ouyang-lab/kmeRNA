#!/usr/bin/env python3
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from math import sqrt
import sys

# ---------------------------------------------------------
# DeLong Implementation (from Sun & Xu 2014 adaptation)
# ---------------------------------------------------------
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
    T2[J] = T + 1  # positions starting at 1
    return T2

def fastDeLong(preds_sorted, label_1_count):
    m = label_1_count
    n = preds_sorted.shape[1] - m

    pos_preds = preds_sorted[:, :m]
    neg_preds = preds_sorted[:, m:]

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
    """
    Return:
      auc: standard AUROC (sklearn)
      auc_var: variance of AUROC from DeLong
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)

    # sort by label: positives (1) first, then negatives (0)
    # this is what fastDeLong expects
    order = np.argsort(-y_true)  # y_true=1 come first
    y_sorted = y_true[order]
    preds_sorted = y_pred[order]

    pos_count = np.sum(y_sorted == 1)
    if pos_count == 0 or pos_count == len(y_sorted):
        raise ValueError("Need both positive and negative samples for AUC.")

    preds_matrix = preds_sorted.reshape(1, -1)  # shape (n_classifiers=1, n_samples)

    auc = roc_auc_score(y_true, y_pred)  # actual AUROC
    auc_est_delong, auc_var = fastDeLong(preds_matrix, pos_count)

    # auc_est_delong should be ~equal to auc; we return auc for clarity
    return auc, auc_var


# ---------------------------------------------------------
# Meta-analysis using random-effects model (DerSimonian-Laird)
# ---------------------------------------------------------
def random_effects_meta(diffs, ses):
    weights_fixed = 1 / (ses**2)
    fixed_effect = np.sum(weights_fixed * diffs) / np.sum(weights_fixed)
    Q = np.sum(weights_fixed * (diffs - fixed_effect)**2)
    df = len(diffs) - 1
    tau2 = max(0, (Q - df) / (np.sum(weights_fixed) - np.sum(weights_fixed**2) / np.sum(weights_fixed)))
    weights_random = 1 / (ses**2 + tau2)
    pooled = np.sum(weights_random * diffs) / np.sum(weights_random)
    se_pooled = sqrt(1 / np.sum(weights_random))

    z = pooled / se_pooled
    from scipy.stats import norm
    pval = 2 * (1 - norm.cdf(abs(z)))
    return pooled, se_pooled, pval, tau2

# ---------------------------------------------------------
# Load replicates
# ---------------------------------------------------------
def load_file(path):
    df = pd.read_csv(path,sep="\t")
    # ground_truth    pred_prob
    if not {"ground_truth", "pred_prob"}.issubset(df.columns):
        raise ValueError(f"Input {path} must contain ground_truth and pred columns.")
    return df["ground_truth"].values, df["pred_prob"].values

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputA", nargs="+", required=True, help="Replicate files for group A")
    parser.add_argument("--inputB", nargs="+", required=True, help="Replicate files for group B")
    args = parser.parse_args()

    if len(args.inputA) != len(args.inputB):
        print("ERROR: groups must have same number of biological replicates.", file=sys.stderr)
        sys.exit(1)

    n_rep = len(args.inputA)
    aucA, aucB = [], []
    diffs, ses = [], []

    for i in range(n_rep):
        yA, pA = load_file(args.inputA[i])
        yB, pB = load_file(args.inputB[i])

        aucA_i, varA = delong_auc_and_var(yA, pA)
        aucB_i, varB = delong_auc_and_var(yB, pB)

        diff = aucA_i - aucB_i
        se = sqrt(varA + varB)

        aucA.append(aucA_i)
        aucB.append(aucB_i)
        diffs.append(diff)
        ses.append(se)

        print(f"Replicate {i+1}:")
        print(f"  AUROC A = {aucA_i:.4f}")
        print(f"  AUROC B = {aucB_i:.4f}")
        print(f"  Difference = {diff:.4f}, SE = {se:.4f}")
        print("")

    diffs = np.array(diffs)
    ses = np.array(ses)

    pooled, se_pooled, pval, tau2 = random_effects_meta(diffs, ses)

    print("===== META-ANALYSIS RESULTS =====")
    print(f"Pooled AUROC difference (A - B): {pooled:.4f}")
    print(f"SE: {se_pooled:.4f}")
    print(f"p-value: {pval:.2e}")
    print(f"Between-replicate variance tau²: {tau2:.2e}")

if __name__ == "__main__":
    main()

