#!/usr/bin/env python3
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.
import argparse
import numpy as np
import pandas as pd
import shap
import pickle
import joblib
import keras
from itertools import product


# --------------------------------------------------------
# Generate k-mer feature names
# --------------------------------------------------------
def initialize_kmers(min_k, max_k):
    """Initialize all possible kmers."""
    bases = ['A', 'C', 'G', 'U']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return kmers


# --------------------------------------------------------
# Load embedded data
# --------------------------------------------------------
def load_embedded_data(path):
    with open(path, "rb") as f:
        raw = pickle.load(f)
    ids = np.array(raw[::2])
    data = np.array(raw[1::2])
    return ids, data


# --------------------------------------------------------
# Compute SHAP using KernelExplainer
# --------------------------------------------------------
def compute_kernel_shap(model, data, ids, sample_size=200, bg_size=100):
    background = data[np.random.choice(len(data), bg_size, replace=False)]

    # prediction wrapper
    def predict_fn(X):
        X = X.astype(np.float32)
        out = model.predict(X)
        return out.reshape(-1)

    print(f"[INFO] Initializing KernelExplainer with background size = {bg_size}")
    explainer = shap.KernelExplainer(predict_fn, background) #silent=True)

    X = data[:sample_size]

    print("[INFO] Computing SHAP values...")
    shap_values = explainer.shap_values(X)
    shap_values = np.array(shap_values)

    print(f"[INFO] SHAP matrix shape = {shap_values.shape}")

    # positive vs negative groups
    neg_mask = np.char.endswith(ids[:sample_size], "_neg")
    pos_mask = ~neg_mask

    pos_shap = shap_values[pos_mask]
    neg_shap = shap_values[neg_mask]

    mean_diff = pos_shap.mean(axis=0) - neg_shap.mean(axis=0)

    return shap_values, mean_diff


# --------------------------------------------------------
# Main
# --------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Compute SHAP values using KernelExplainer with unified k-mer naming.")
    parser.add_argument("--model", required=True, help="Path to .keras or sklearn model (.joblib/.pkl)")
    parser.add_argument("--data", required=True, help="Path to embedded input data (.pkl)")
    parser.add_argument("--output", required=True, help="Output prefix for saving SHAP results")
    parser.add_argument("--sample_size", type=int, default=100, help="Number of samples for SHAP computation")
    parser.add_argument("--bg_size", type=int, default=50, help="Background size for Kernel SHAP")
    parser.add_argument("--min_k", type=int, default=1, help="Minimum k-mer length")
    parser.add_argument("--max_k", type=int, default=4, help="Maximum k-mer length")

    args = parser.parse_args()

    # ---------------- Load Model ----------------
    print("[INFO] Loading model:", args.model)
    try:
        if args.model.endswith(".keras"):
            model = keras.models.load_model(args.model)
        else:
            model = joblib.load(args.model)
    except Exception as e:
        print("[ERROR] Failed to load model:", e)
        return

    # ---------------- Load Data ----------------
    print("[INFO] Loading embedded data:", args.data)
    ids, data = load_embedded_data(args.data)

    # ---------------- Build Feature Names ----------------
    print("[INFO] Generating k-mer feature names...")
    kmers = initialize_kmers(args.min_k, args.max_k)

    feature_names = (
        [f"{k}-query" for k in kmers] +
        [f"{k}-target" for k in kmers] +
        ["RNAduplex", "RNAcofold"]
    )

    # ---------------- Assign Feature Names to ALL Models ----------------
    model.feature_names_in_ = np.array(feature_names, dtype=object)

    # ---------------- Compute SHAP Values ----------------
    shap_values, mean_diff = compute_kernel_shap(
        model,
        data,
        ids,
        sample_size=args.sample_size,
        bg_size=args.bg_size
    )

    # ---------------- Validate Feature Dimensions ----------------
    n_features = shap_values.shape[1]

    if n_features != len(feature_names):
        print(f"[WARN] SHAP output has {n_features} features but expected {len(feature_names)}.")
        print("[WARN] Truncating or padding feature names is NOT recommended.")
        feature_names = feature_names[:n_features]

    # ---------------- Build SHAP Output DataFrame ----------------
    df = pd.DataFrame(
        shap_values.T,
        index=feature_names,
        columns=[f"sample_{i}" for i in range(shap_values.shape[0])]
    )

    df["mean_diff"] = mean_diff

    out_path = args.output + ".tsv"
    print("[INFO] Saving SHAP results:", out_path)
    df.to_csv(out_path, sep="\t")

    print("[INFO] All done.")


if __name__ == "__main__":
    main()

