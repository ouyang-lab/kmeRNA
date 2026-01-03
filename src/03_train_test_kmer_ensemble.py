#!/usr/bin/env python3
# MIT License
# Copyright (c) 2026 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import os
import argparse
import pickle
import sys
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier, AdaBoostClassifier
import multiprocessing as mp
from multiprocessing import Queue
import numpy as np
from scipy import stats
import random
from datetime import datetime
import joblib
import csv
import time
from itertools import product, permutations, combinations_with_replacement
import traceback

########################
# Utility / safety helpers
########################

def sanitize_X(X):
    """
    Ensure X is a float64 numpy array with no NaN/inf.
    """
    X = np.asarray(X, dtype=np.float64)
    X = np.nan_to_num(X, nan=0.0, posinf=1e9, neginf=-1e9)
    return X


def safe_predict_proba(model, X, model_name="MODEL"):
    """
    Safe wrapper around predict_proba that:
      - converts X → float64 and removes NaN/inf
      - handles 1-class models
      - ensures returned probs length matches n_samples
      - reports which rows were excluded by the estimator
    """

    X = sanitize_X(X)
    n_samples = X.shape[0]

    try:
        probs = model.predict_proba(X)
    except Exception as e:
        print(f"[{model_name}] predict_proba FAILED. Returning zeros. Error:\n{e}")
        traceback.print_exc()
        return np.zeros(n_samples, dtype=np.float64), list(range(n_samples))

    # probs should be (n_samples, n_classes). If shorter, identify dropped rows.
    dropped_rows = []

    if probs.shape[0] != n_samples:
        print(f"[{model_name}] WARNING: predict_proba returned {probs.shape[0]} rows for {n_samples} inputs")

        # Build a map of rows that survived
        surviving_indices = np.arange(probs.shape[0])

        # Missing rows = the difference
        dropped_rows = list(sorted(set(range(n_samples)) - set(surviving_indices)))

        print(f"[{model_name}] DROPPED ROWS: {dropped_rows}")

        # Fix length mismatch (pad zeros)
        if probs.shape[0] < n_samples:
            pad_len = n_samples - probs.shape[0]
            probs = np.pad(probs, ((0, pad_len), (0, 0)), mode="constant", constant_values=0.0)

        # If it's longer (rare), truncate
        elif probs.shape[0] > n_samples:
            probs = probs[:n_samples, :]

    # Now choose the correct probability column
    if probs.ndim == 1:
        # flatten case
        p1 = probs.astype(np.float64)

    elif probs.shape[1] == 1:
        # One-class model
        if hasattr(model, "classes_") and model.classes_[0] == 1:
            p1 = np.ones(n_samples, dtype=np.float64)
        else:
            p1 = np.zeros(n_samples, dtype=np.float64)

    else:
        # Multi-class, choose column for class 1 when available
        if hasattr(model, "classes_") and 1 in model.classes_:
            col = list(model.classes_).index(1)
        else:
            col = probs.shape[1] - 1
        p1 = probs[:, col].astype(np.float64)

    return p1, dropped_rows

def mp_train_wrapper(args):
    """
    Worker wrapper for multiprocessing.
    Each worker trains ONE model and returns:
    (model_name, model_object, train_probabilities, logs)
    """
    model_type, X, y, seed = args

    if model_type == "rf":
        model, probs, logs = run_random_forest_classifier(X, y, seed)
    elif model_type == "et":
        model, probs, logs = run_extra_trees_classifier(X, y, seed)
    elif model_type == "gb":
        model, probs, logs = run_gradient_boosting_classifier(X, y, seed)
    else:
        return (model_type, None, None, [f"Unknown model type {model_type}"])

    return (model_type, model, probs, logs)


########################
# I/O helpers
########################

# Function to load embedded matrices from pickle file
def load_embedded_data(file_path):
    with open(file_path, 'rb') as f:
        embedded_in = pickle.load(f)
    ids = embedded_in[::2]
    embedded_out = embedded_in[1::2]
    return ids, embedded_out

def save_model_data(file_path,model):
    model.save(file_path)

def load_model_data(file_path):
    model = keras.models.load_model(file_path)
    return model

def prepare_data(embedded_data, labels):
    X = sanitize_X(embedded_data)
    y = np.asarray(labels)
    return X, y

def initialize_kmers(min_k,max_k):
    """Initialize all possible kmers with count 0"""
    bases = ['A', 'C', 'G', 'U']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return {kmer: 0 for kmer in kmers}

def get_feature_importance(model, feature_names):
    importance_scores = model.feature_importances_
    feature_importance_dict = dict(zip(feature_names, importance_scores))
    sorted_importance = sorted(feature_importance_dict.items(), key=lambda x: x[1], reverse=True)
    return sorted_importance

########################
# Model definitions
########################

def random_forest_classifier(train_data, train_labels, seed_val):
    train_data = sanitize_X(train_data)
    train_labels = np.asarray(train_labels)

    start = time.time()
    logs = []
    rf_classifier = RandomForestClassifier(
        n_estimators=512,
        criterion='gini',
        max_features=0.8,
        max_depth=100,
        min_impurity_decrease=0,
        random_state=seed_val,
        n_jobs=32
    )
    logs.append(f"RF Classifier initialized {datetime.now()}")
    logs.append("Training the RF model...")
    rf_classifier.fit(train_data, train_labels)
    logs.append("RF Training complete")
    logs.append("Predicting class 1 probabilities using RF model...")
    train_predictions = safe_predict_proba(rf_classifier, train_data)
    stop = time.time()
    logs.append(f"RF Predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return rf_classifier, train_predictions, logs

def extra_trees_classifier(train_data, train_labels, seed_val):
    train_data = sanitize_X(train_data)
    train_labels = np.asarray(train_labels)

    start = time.time()
    logs = []
    et_classifier = ExtraTreesClassifier(
        n_estimators=1024,
        criterion='gini',
        max_features=0.8,
        max_depth=64,
        min_impurity_decrease=0,
        random_state=seed_val,
        n_jobs=32
    )
    logs.append(f"ExtraTrees Classifier initialized {datetime.now()}")
    logs.append("Training the ET model...")
    et_classifier.fit(train_data, train_labels)
    logs.append("ET Training complete")
    logs.append("Predicting class 1 probabilities using ET model...")
    train_predictions = safe_predict_proba(et_classifier, train_data)
    stop = time.time()
    logs.append(f"ET predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return et_classifier, train_predictions, logs

def gradient_boosting_classifier(train_data, train_labels, seed_val):
    train_data = sanitize_X(train_data)
    train_labels = np.asarray(train_labels)

    logs = []
    start = time.time()
    gb_classifier = GradientBoostingClassifier(
        n_estimators=768,
        max_features=0.8,
        subsample=0.8,
        max_depth=5,
        learning_rate=0.1,
        min_impurity_decrease=0,
        n_iter_no_change=5,
        random_state=seed_val
    )
    logs.append(f"GradientBoosting Classifier initialized {datetime.now()}")
    logs.append("Training the GB model...")
    gb_classifier.fit(train_data, train_labels)
    print("GB Training complete")
    print("Predicting class 1 probabilities using GB model...")
    logs.append("GB Training complete")
    logs.append("Predicting class 1 probabilities using GB model...")
    train_predictions = safe_predict_proba(gb_classifier, train_data)
    stop = time.time()
    logs.append(f"GB predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return gb_classifier, train_predictions, logs

def histogram_gradient_boosting_classifier(train_data, train_labels, seed_val):
    train_data = sanitize_X(train_data)
    train_labels = np.asarray(train_labels)

    logs = []
    hgb_classifier = HistGradientBoostingClassifier(
        max_features=None,
        random_state=seed_val
    )
    print(f"HistGradientBoosting Classifier initialized {datetime.now()}")
    print("Training the HGB model...")
    hgb_classifier.fit(train_data, train_labels)
    print("HGB Training complete")
    print("Predicting class 1 probabilities using HGB model...")
    train_predictions = safe_predict_proba(hgb_classifier, train_data)
    print(f"HGB predictions complete {datetime.now()}")
    logs.append("HGB Training complete")
    logs.append("HGB predictions complete")
    return hgb_classifier, train_predictions, logs

def adaboost_classifier(train_data, train_labels, seed_val):
    train_data = sanitize_X(train_data)
    train_labels = np.asarray(train_labels)

    logs = []
    ab_estimator = GradientBoostingClassifier(
        max_features=0.67,
        max_depth=12,
        min_impurity_decrease=0
    )
    ab_classifier = AdaBoostClassifier(
        estimator=ab_estimator,
        n_estimators=20,
        learning_rate=1.0,
        random_state=seed_val
    )
    logs.append(f"AdaBoost Classifier initialized {datetime.now()}")
    logs.append("Training the AB model...")
    ab_classifier.fit(train_data, train_labels)
    logs.append("AB Training complete")
    logs.append("Predicting class 1 probabilities using AB model...")
    train_predictions = safe_predict_proba(ab_classifier, train_data)
    logs.append(f"AB predictions complete {datetime.now()}")
    return ab_classifier, train_predictions, logs

########################
# Testing wrappers (if used)
########################

def random_forest_testing(rf_classifier, test_data):
    return safe_predict_proba(rf_classifier, test_data)

def extra_trees_testing(et_classifier, test_data):
    return safe_predict_proba(et_classifier, test_data)

def gradient_boosting_testing(gb_classifier, test_data):
    return safe_predict_proba(gb_classifier, test_data)

def adaboost_testing(ab_classifier, test_data):
    return safe_predict_proba(ab_classifier, test_data)

########################
# Run wrappers with crash protection
########################

def run_random_forest_classifier(train_data, train_labels, seed_val):
    logs = []
    train_labels = np.asarray(train_labels)
    rf_out = None
    train_probabilities_rf = np.zeros(len(train_labels), dtype=np.float64)
    try:
        rf_out, train_probabilities_rf, logs = random_forest_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in RandomForest: {e}")
        traceback.print_exc()
        logs.append(f"Error in RandomForest: {e}")
    return (rf_out, train_probabilities_rf, logs)

def run_extra_trees_classifier(train_data, train_labels, seed_val):
    logs = []
    train_labels = np.asarray(train_labels)
    et_out = None
    train_probabilities_et = np.zeros(len(train_labels), dtype=np.float64)
    try:
        et_out, train_probabilities_et, logs = extra_trees_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in ExtraTrees: {e}")
        traceback.print_exc()
        logs.append(f"Error in ExtraTrees: {e}")
    return (et_out, train_probabilities_et, logs)

def run_gradient_boosting_classifier(train_data, train_labels, seed_val):
    logs = []
    train_labels = np.asarray(train_labels)
    gb_out = None
    train_probabilities_gb = np.zeros(len(train_labels), dtype=np.float64)
    try:
        gb_out, train_probabilities_gb, logs = gradient_boosting_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in GradientBoosting: {e}")
        traceback.print_exc()
        logs.append(f"Error in GradientBoosting: {e}")
    return (gb_out, train_probabilities_gb, logs)

def run_histogram_gradient_boosting_classifier(train_data, train_labels, seed_val):
    logs = []
    train_labels = np.asarray(train_labels)
    hgb_out = None
    train_probabilities_hgb = np.zeros(len(train_labels), dtype=np.float64)
    try:
        hgb_out, train_probabilities_hgb, logs = histogram_gradient_boosting_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in HistGradientBoosting: {e}")
        traceback.print_exc()
        logs.append(f"Error in HistGradientBoosting: {e}")
    return "hgb", (hgb_out, train_probabilities_hgb, logs)

def run_ada_gradient_boosting_classifier(train_data, train_labels, seed_val):
    logs = []
    train_labels = np.asarray(train_labels)
    ab_out = None
    train_probabilities_ab = np.zeros(len(train_labels), dtype=np.float64)
    try:
        ab_out, train_probabilities_ab, logs = adaboost_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in AdaBoosting: {e}")
        traceback.print_exc()
        logs.append(f"Error in AdaBoosting: {e}")
    return "ab", (ab_out, train_probabilities_ab, logs)

def wrapper_function(wrap_args):
    """Wraps the worker function and stores the result in the shared dict."""
    results = wrap_args[0]
    out_args = wrap_args[1:4]
    name_arg = wrap_args[-1]
    if name_arg == "rf":
        key, value = run_random_forest_classifier(*out_args)
    elif name_arg == "et":
        key, value = run_extra_trees_classifier(*out_args)
    elif name_arg == "gb":
        key, value = run_gradient_boosting_classifier(*out_args)
    elif name_arg == "hgb":
        key, value = run_histogram_gradient_boosting_classifier(*out_args)
    elif name_arg == "ab":
        key, value = run_ada_gradient_boosting_classifier(*out_args)
    else:
        print(f"{name_arg} is not recognized")
        return
    results[key] = value

########################
# Main
########################

def main():
    start = time.time()
    parser = argparse.ArgumentParser(description="Train CNN on embedded data")
    parser.add_argument("--train_data", help="Path to embedded training data pickle file", required=False)
    parser.add_argument("--train_labels", help="Path to embedded training labels text file", required=False)
    parser.add_argument("--test_data", help="Path to embedded test data pickle file", required=False)
    parser.add_argument("--test_labels", help="Path to embedded test labels text file", required=False)
    parser.add_argument("--model_out", help="Path to trained model information", required=False)
    parser.add_argument("--model_in", help="Path to pre-trained model information", required=False)
    parser.add_argument("--seed_val", type=int, help="input seed value", required=False)
    parser.add_argument("--output_train_prob", help="Path to training output", required=False)
    parser.add_argument("--output_test_prob", help="Path to testing output", required=False)
    parser.add_argument("--min_k", type=int, help="Min k-mer length in bases", default=1)
    parser.add_argument("--max_k", type=int, help="Max k-mer length in bases", default=4)
    args = parser.parse_args()

    print(f'Random seed: {args.seed_val}')
    os.environ['PYTHONHASHSEED'] = str(args.seed_val)
    random.seed(args.seed_val)
    np.random.seed(args.seed_val)
    seed_val = int(args.seed_val)

    # Load training and test embedded data
    if args.train_data:
        train_ids, train_data = load_embedded_data(args.train_data)
        train_labels = np.loadtxt(args.train_labels)
        train_data = sanitize_X(train_data)

    if args.test_data:
        test_ids, test_data = load_embedded_data(args.test_data)
        test_labels = np.loadtxt(args.test_labels)
        test_data = sanitize_X(test_data)

    kmer_init = initialize_kmers(args.min_k, args.max_k)
    feature_names = [kmer1 + "-query" for kmer1 in kmer_init] + \
                    [kmer2 + "-target" for kmer2 in kmer_init] + \
                    ["RNAduplex_energy", "RNAcofold_energy"]

    if args.model_in:
        print("Loading RF models from file...")
        rf_in = joblib.load(args.model_in + ".RF.joblib")
        et_in = joblib.load(args.model_in + ".ET.joblib")
        gb_in = joblib.load(args.model_in + ".GB.joblib")

        if args.output_train_prob and args.train_data:
            print("Beginning training predictions...")
            train_probabilities_rf = safe_predict_proba(rf_in, train_data)
            print("RF model training predictions complete")
            train_probabilities_et = safe_predict_proba(et_in, train_data)
            print("ET model training predictions complete")
            train_probabilities_gb = safe_predict_proba(gb_in, train_data)
            print("GB model training predictions complete")

        if args.output_test_prob and args.test_data:
            print("Beginning testing predictions...")
            #test_probabilities_rf = safe_predict_proba(rf_in, test_data)
            test_probabilities_rf, dropped_rf = safe_predict_proba(rf_in, test_data, model_name="RF")
            if dropped_rf:
                print(f"RF model dropped rows: {dropped_rf}")
            print("RF model testing predictions complete")
            #test_probabilities_et = safe_predict_proba(et_in, test_data)
            test_probabilities_et, dropped_et = safe_predict_proba(et_in, test_data, model_name="ET")
            if dropped_et:
                print(f"RF model dropped rows: {dropped_et}")
            print("ET model testing predictions complete")
            #test_probabilities_gb = safe_predict_proba(gb_in, test_data)
            test_probabilities_gb, dropped_gb = safe_predict_proba(gb_in, test_data, model_name="GB")
            if dropped_rf:
                print(f"RF model dropped rows: {dropped_gb}")
            print("GB model testing predictions complete")

    elif args.model_out and args.train_data:

        print("Launching parallel training for RF, ET, GB...")

        # Prepare job list for workers
        job_list = [
            ("rf", train_data, train_labels, seed_val),
            ("et", train_data, train_labels, seed_val),
            ("gb", train_data, train_labels, seed_val),
        ]

        # Train 3 models in parallel
        #import multiprocessing as mp
        with mp.Pool(processes=3) as pool:
            results = pool.map(mp_train_wrapper, job_list)

        # Collect results from workers
        model_dict = {}
        for model_type, model_obj, train_probs, logs in results:
            model_dict[model_type] = {
                "model": model_obj,
                "probs": train_probs,
                "logs": logs
            }

        # Unpack final results
        rf_out = model_dict["rf"]["model"]
        train_probabilities_rf = model_dict["rf"]["probs"]
        logs_rf = model_dict["rf"]["logs"]

        et_out = model_dict["et"]["model"]
        train_probabilities_et = model_dict["et"]["probs"]
        logs_et = model_dict["et"]["logs"]

        gb_out = model_dict["gb"]["model"]
        train_probabilities_gb = model_dict["gb"]["probs"]
        logs_gb = model_dict["gb"]["logs"]

        # Print logs
        for elem in logs_rf: print("[RF]", elem)
        for elem in logs_et: print("[ET]", elem)
        for elem in logs_gb: print("[GB]", elem)

        print("Parallel training complete.")

        # Save trained models
        print("Saving models...")
        joblib.dump(rf_out, args.model_out + ".RF.joblib")
        joblib.dump(et_out, args.model_out + ".ET.joblib")
        joblib.dump(gb_out, args.model_out + ".GB.joblib")
        print("Models saved.")

        # Save importance scores
        print("Saving model importances...")
        sorted_importance_rf = get_feature_importance(rf_out, feature_names)
        sorted_importance_et = get_feature_importance(et_out, feature_names)
        sorted_importance_gb = get_feature_importance(gb_out, feature_names)

        with open(args.model_out + ".RF.feature_importance.tsv", 'w') as f:
            for feature, importance in sorted_importance_rf:
                f.write(f"{feature}\t{importance}\n")

        with open(args.model_out + ".ET.feature_importance.tsv", 'w') as f:
            for feature, importance in sorted_importance_et:
                f.write(f"{feature}\t{importance}\n")

        with open(args.model_out + ".GB.feature_importance.tsv", 'w') as f:
            for feature, importance in sorted_importance_gb:
                f.write(f"{feature}\t{importance}\n")

        print("Feature importances saved.")
 
    ###########################
    # Writing out probabilities
    ###########################
    if args.train_data and args.output_train_prob and 'train_probabilities_rf' in locals():
        N = min(
            len(train_labels),
            len(train_probabilities_rf),
            len(train_probabilities_et),
            len(train_probabilities_gb)
        )
        if N < len(train_labels):
            print(f"Warning: training length mismatch, truncating to {N} samples.")
        scores_train_probs = []
        for pos in range(N):
            scores_train_probs.append((
                str(int(train_labels[pos])),
                str(train_probabilities_rf[pos]),
                str(train_probabilities_et[pos]),
                str(train_probabilities_gb[pos])
            ))
        print("Writing training data predictions to file...")
        with open(args.output_train_prob, 'w') as f:
            f.write('ground_truth\tpred_prob_RF\tpred_prob_ET\tpred_prob_GB\n')
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(scores_train_probs)
        print("RF training results written.")

    if args.test_data and args.output_test_prob and 'test_probabilities_rf' in locals():
        N = min(
            len(test_labels),
            len(test_probabilities_rf),
            len(test_probabilities_et),
            len(test_probabilities_gb)
        )
        if N < len(test_labels):
            print(f"Warning: testing length mismatch, truncating to {N} samples.")
        scores_test_probs = []
        for pos in range(N):
            scores_test_probs.append((
                str(int(test_labels[pos])),
                str(test_probabilities_rf[pos]),
                str(test_probabilities_et[pos]),
                str(test_probabilities_gb[pos])
            ))
        print("Writing testing data predictions to file...")
        with open(args.output_test_prob, 'w') as f:
            f.write('ground_truth\tpred_prob_RF\tpred_prob_ET\tpred_prob_GB\n')
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(scores_test_probs)
        print("Ensemble model testing results written")
        end = time.time()
        print("Ensemble model time elapsed: ", round((end-start), 1), "sec", round((end-start)/60, 2), "min")

if __name__ == "__main__":
    main()

