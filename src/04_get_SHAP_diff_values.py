#!/usr/bin/env python3
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import argparse
import shap
import numpy as np
import tensorflow as tf
import pandas as pd
import keras
import pickle
from sklearn.ensemble import ExtraTreesClassifier 
from sklearn.inspection import permutation_importance
from scipy import stats
#from tensorflow.keras.models import load_model
import matplotlib.pyplot as plt
from itertools import product, permutations, combinations_with_replacement
import joblib

def load_embedded_data(file_path):
    with open(file_path, 'rb') as f:
        embedded_in = pickle.load(f)
    ids = embedded_in[::2]
    embedded_out = np.array(embedded_in[1::2])
    return ids, embedded_out

def save_model_data(file_path,model):
    model.save(file_path)

def load_model_data(file_path):
    model = keras.models.load_model(file_path)
    return model
def load_tree_model_data(file_path):
    model = joblib.load(file_path)
    return model

def prepare_data(embedded_data, labels):
    return np.array(embedded_data), np.array(labels)

def reshape_input_for_model(model, data):
    expected_input_shape = model.input_shape  # (None, 1457) or (None, 1457, 1)

    if len(expected_input_shape) == 2:
        # Model expects 2D input: (batch, features)
        if len(data.shape) == 3:
            # Strip the trailing dimension
            data = data.squeeze(-1)
    elif len(expected_input_shape) == 3:
        # Model expects 3D input: (batch, features, channels)
        if len(data.shape) == 2:
            # Add a trailing dimension
            data = np.expand_dims(data, axis=-1)

    return data

def reshape_for_model(model, input_array):
    """Ensure the input array shape matches the model's expected input."""
    model_input_shape = model.input_shape

    if len(model_input_shape) == 2:
        # Model expects (batch_size, features)
        if len(input_array.shape) == 3:
            input_array = input_array.squeeze(-1)
    elif len(model_input_shape) == 3:
        # Model expects (batch_size, features, channels)
        if len(input_array.shape) == 2:
            input_array = np.expand_dims(input_array, -1)
    return input_array

def reshape_for_tree_model(mode, input_array):
    """
    Ensure the input array shape is compatible with scikit-learn tree-based models.
    These models expect input of shape (n_samples, n_features).
    """
    if input_array.ndim == 1:
        # Single sample with 1D input, reshape to (1, n_features)
        input_array = input_array.reshape(1, -1)
    elif input_array.ndim == 3 and input_array.shape[2] == 1:
        # Convert (n_samples, n_features, 1) to (n_samples, n_features)
        input_array = input_array.squeeze(-1)
    elif input_array.ndim > 2:
        raise ValueError(f"Unexpected input shape {input_array.shape} for scikit-learn model.")

    return input_array

def initialize_kmers(min_k,max_k):
    """Initialize all possible kmers with count 0"""
    bases = ['A', 'C', 'G', 'T']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return {kmer: 0 for kmer in kmers}



def compute_shap_values_diff(model, ids, data, sample_size, bg_sample_size=100):
    # Updated to ensure background and data do not overlap
    # Shuffle data and ids together
    indices = np.random.permutation(len(data))
    data = data[indices]
    ids = np.array(ids)[indices]

    background = data[np.random.choice(data.shape[0], bg_sample_size, replace=False)]

    # Exclude background samples from data and ids
    mask = ~np.any(np.all(data[:, None] == background, axis=2), axis=1)
    remaining_data = data[mask]
    remaining_ids = ids[mask]

    # Restrict remaining data and ids to sample size
    data = remaining_data[:sample_size]
    ids = remaining_ids[:sample_size]

    # Reshape background and data for SHAP
    background = reshape_for_model(model, background)
    data = reshape_for_model(model, data)

    print("Model expects:", model.input_shape)
    print("Background shape:", background.shape)
    print("Data shape:", data.shape)

    # Compute SHAP values
    explainer = shap.GradientExplainer(model, background)
    print("GradientExplainer finished")
    try:
        shap_values = explainer.shap_values(data)
    except Exception as e:
        print(f"Failed: {e}")
        quit()

    # Convert to array in case it's a multi-class list of arrays
    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    # Split SHAP values into pos/neg based on ID suffix (containined "_neg")
    ids = np.array(ids)
    neg_mask = np.char.endswith(ids, "_neg")
    pos_mask = ~neg_mask

    pos_shap = shap_values[pos_mask]
    neg_shap = shap_values[neg_mask]

    # Compute mean difference
    #mean_diff = reshape_for_model(np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0))
    mean_diff_raw = np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0)
    mean_diff = mean_diff_raw.reshape(-1)
    mean_diff = mean_diff.reshape(1, shap_values.shape[1]) 

    return shap_values, data, mean_diff


def compute_shap_values_diff_old(model, ids, data, sample_size, bg_sample_size=100):
    # Shuffle data and ids together
    indices = np.random.permutation(len(data))
    data = data[indices]
    ids = np.array(ids)[indices]

    # Limit to sample size
    data = data[:sample_size]
    ids = ids[:sample_size]

    # Prepare background for SHAP
    background = data[np.random.choice(data.shape[0], bg_sample_size, replace=False)]
    background = reshape_for_model(model, background)
    data = reshape_for_model(model, data)

    print("Model expects:", model.input_shape)
    print("Background shape:", background.shape)
    print("Data shape:", data.shape)

    # Compute SHAP values
    explainer = shap.GradientExplainer(model, background)
    print("GradientExplainer finished")
    try:
        shap_values = explainer.shap_values(data)
    except Exception as e:
        print(f"Failed: {e}")
        quit()

    # Convert to array in case it's a multi-class list of arrays
    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    # Split SHAP values into pos/neg based on ID suffix
    ids = np.array(ids)
    neg_mask = np.char.endswith(ids, "_neg")
    pos_mask = ~neg_mask

    pos_shap = shap_values[pos_mask]
    neg_shap = shap_values[neg_mask]

    # Compute mean difference
    #mean_diff = reshape_for_model(np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0))
    mean_diff_raw = np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0)
    mean_diff = mean_diff_raw.reshape(-1)
    mean_diff = mean_diff.reshape(1, shap_values.shape[1]) 

    return shap_values, data, mean_diff

def compute_perm_values_diff(model, ids, data, labels, sample_size, n_repeats=10):
    # Shuffle data, ids, and labels in sync
    indices = np.random.permutation(len(data))
    data = data[indices][:sample_size]
    ids = np.array(ids)[indices][:sample_size]
    labels = np.array(labels)[indices][:sample_size]

    # Split data by class using "_neg" suffix
    neg_mask = np.char.endswith(ids, "_neg")
    pos_mask = ~neg_mask

    X_pos = data[pos_mask]
    y_pos = labels[pos_mask]
    X_neg = data[neg_mask]
    y_neg = labels[neg_mask]

    print(f"Positive samples: {X_pos.shape[0]}, Negative samples: {X_neg.shape[0]}")

    # Compute permutation importance on each group
    pos_result = permutation_importance(model, X_pos, y_pos, n_repeats=n_repeats, random_state=42)
    neg_result = permutation_importance(model, X_neg, y_neg, n_repeats=n_repeats, random_state=42)

    # Compute the mean difference in importances
    perm_diff = pos_result.importances_mean - neg_result.importances_mean

    return pos_result, neg_result, perm_diff


def compute_shap_values(model, data, sample_size, bg_sample_size=100):
    # Use a small background set for SHAP (important for efficiency)
    background = data[np.random.choice(data.shape[0], bg_sample_size, replace=False)]
    background = reshape_for_model(model, background) 
    mask = ~np.any(np.all(data[:, None] == background, axis=2), axis=1)
    remaining_data = data[mask]
    data_out = remaining_data[np.random.permutation(len(remaining_data))]
    data_out = data_out[:sample_size]
 
    data_out = reshape_for_model(model, data_out)
    
    print("Model expects:", model.input_shape)
    print("Background shape:", background.shape)
    print("Data shape:", data_out.shape)
    # Use Explainer (works fine with Keras models)
    explainer = shap.GradientExplainer(model, background)
    #explainer = shap.DeepExplainer(model, background)
    print("GradientExplainer finished")
    try:
        shap_values = explainer.shap_values(data_out)
    except Exception as e:
        print(f"Failed: {e}")
        quit()

    return shap_values, data_out

def compute_tree_shap_values(model, data, sample_size, bg_sample_size=100):
    print("Input Shape",data.shape)
    # Use a small background set for SHAP (important for efficiency)
    background = data[np.random.choice(data.shape[0], bg_sample_size, replace=False)]
    mask = ~np.any(np.all(data[:, None] == background, axis=2), axis=1)
    remaining_data = data[mask]
    data_out = remaining_data[np.random.permutation(len(remaining_data))]
    data_out = data_out[:sample_size]
    
    print("Background shape:", background.shape)
    print("Data shape:", data_out.shape)
    #quit()
    background = np.array(background).astype(np.float32)
    data_out = np.array(data_out).astype(np.float32)
    # Use Explainer (works well with Keras models)
    explainer = shap.TreeExplainer(model, background) # model_output="probability")
    #explainer = shap.DeepExplainer(model, background)
    #quit()
    print("TreeExplainer finished")
    try:
        shap_values = explainer.shap_values(data_out)
    except:
        print(f"Failed")
        quit()

    return shap_values, data_out

def compute_tree_shap_values_diff(model, ids, data, kmer_init, sample_size, bg_sample_size=100):
    # Updated to ensure background and data do not overlap
    # Shuffle data and ids together
    indices = np.random.permutation(len(data))
    data = data[indices]
    ids = np.array(ids)[indices]

    background = data[np.random.choice(data.shape[0], bg_sample_size, replace=False)]

    # Exclude background samples from data and ids
    mask = ~np.any(np.all(data[:, None] == background, axis=2), axis=1)
    remaining_data = data[mask]
    remaining_ids = ids[mask]

    # Restrict remaining data and ids to sample size
    data = remaining_data[:sample_size]
    ids = remaining_ids[:sample_size]

    # Reshape background and data for SHAP
    background = reshape_for_tree_model(model, background)
    data = reshape_for_tree_model(model, data)

    #print("Model expects:", model.input_shape)
#    print("Background shape:", background.shape)
#    print("Data shape:", data.shape)
    background = np.array(background).astype(np.float32)
    data_out = np.array(data).astype(np.float32)
    
    model.feature_names_in_ = [kmer1+"-query" for kmer1 in kmer_init]+[kmer2+"-target" for kmer2 in kmer_init]

    # Compute SHAP values
    #explainer = shap.TreeExplainer(model, background)
    explainer = shap.TreeExplainer(model, model_output="raw")
    #explainer = shap.Explainer(model)
    print("TreeExplainer finished")
    try:
        shap_values = explainer.shap_values(data_out)
    except Exception as e:
        print(f"Failed: {e}")
        quit()

    # Convert to array in case it's a multi-class list of arrays
    if isinstance(shap_values, list):
        shap_values = shap_values[0]

    # Split SHAP values into pos/neg based on ID suffix (containined "_neg")
    ids = np.array(ids)
    neg_mask = np.char.endswith(ids, "_neg")
    pos_mask = ~neg_mask

    pos_shap = shap_values[pos_mask]
    neg_shap = shap_values[neg_mask]

    # Compute mean difference
    #mean_diff = reshape_for_model(np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0))
    mean_diff_raw = np.mean(pos_shap, axis=0) - np.mean(neg_shap, axis=0)
    mean_diff = mean_diff_raw.reshape(-1)
    mean_diff = mean_diff.reshape(1, shap_values.shape[1],2) 

    return shap_values, data_out, mean_diff



def main():
    parser = argparse.ArgumentParser(description="Compute SHAP values for a Keras model or scikit-learn Tree model.")
    parser.add_argument("--model", required=True, help="Path to the Keras model (.keras) or scikit-learn Tree model")
    parser.add_argument("--data", required=True, help="Path to the input data (.pkl)")
    parser.add_argument("--output", help="shap_summary[.tsv,.png] "+"| Output path for SHAP summary plot")
    parser.add_argument("--min_k",default=1,type=int)
    parser.add_argument("--max_k",default=4,type=int)
    parser.add_argument("--sample_size", default=1000,type=int)
    parser.add_argument("--bg_sample_size", default=100,type=int)
    parser.add_argument("--tree",action='store_true',default=False)
    args = parser.parse_args()
   
    sample_size=args.sample_size
    bg_sample_size=args.bg_sample_size
    #sample_size=1000 #2500
    #bg_sample_size=100 #250

    # Load model and data
    print("Loading model...")
    if args.tree:
        model = load_tree_model_data(args.model)
    else:
        model = load_model_data(args.model)

    print("Loading data...")
    ids, data = load_embedded_data(args.data)
    #prep_data, prep_ids = prepare_data(ids,data)
    #data = reshape_input_for_model(model, data)
    kmer_init = initialize_kmers(args.min_k,args.max_k)
    if args.tree:
        print("Computing SHAP values...")
        shap_values, data_out, mean_diff = compute_tree_shap_values_diff(model, ids, data,kmer_init, sample_size, bg_sample_size)
    else:
        # since we use dropout
        model.trainable = False
        model.training = False
        print("Computing SHAP values...")
        shap_values, data_out, mean_diff = compute_shap_values_diff(model, ids, data, sample_size, bg_sample_size)
    
    # get feature names
    #kmer_init = initialize_kmers(args.min_k,args.max_k)
    #kmer_init_pairs = initialize_kmer_pairs(args.max_k)
    
    feature_names = [kmer1+"-query" for kmer1 in kmer_init]+[kmer2+"-target" for kmer2 in kmer_init]
    #print("Saving SHAP values to file ()")
    #np.savetxt(f"{args.output}.tsv", shap_values, delimiter="\t")
    #df_out = pd.DataFrame(shap_values.T, index=feature_names)
    #df_out.to_csv(f"{args.output}.tsv", sep="\t", index=False)

    # Squeeze if necessary (e.g., shape (1, 2000, 1457) -> (2000, 1457))
    #shap_array = np.squeeze(shap_values)
    #shap_array = shap_values.squeeze(-1)
    shap_array = shap_values.reshape((sample_size, len(feature_names))) # 1457))
    # For classification: use shap_array[class_index] or aggregate all classes
    #if shap_array.ndim == 3:
        # Aggregate across all classes (2 classes)
    #    shap_array = np.mean(np.abs(shap_array), axis=0)
    #elif shap_array.ndim == 4:
    #    shap_array = np.mean(np.abs(shap_array), axis=0)
    #elif shap_array.ndim == 2:
    #    shap_array = shap_array[1]
    
    print(f"shap_values type: {type(shap_values)}")
    print(f"shap_values shape: {np.array(shap_values).shape}")
    print(f"shap_array type: {type(shap_array)}")
    print(f"shap_array shape: {np.array(shap_array).shape}")
    
    # Compute mean absolute SHAP value per feature
    #global_importance = np.mean(np.abs(shap_array), axis=0)

    # Create DataFrame with feature names and importance scores
    #importance_df = pd.DataFrame({
    #    "feature": feature_names,
    #    "importance": global_importance
    #})

    # Sort by importance (descending)
    #importance_df = importance_df.sort_values(by="importance", ascending=False)
    #sv = shap_values.squeeze()  # removes singleton dimensions
    #importance_df = pd.DataFrame(shap_array.T, columns=[f"sample_{i}" for i in range(shap_array.shape[-1])], index=feature_names)
    importance_df = pd.DataFrame(shap_array.T, 
                             columns=[f"sample_{i}" for i in range(shap_array.shape[0])],  # 100 samples
                             index=feature_names)
    importance_df["sum"] = importance_df.sum(axis=1)
    importance_df["mean"] = importance_df.mean(axis=1)
    importance_df["stderr"] = importance_df.sem(axis=1)
    importance_df["mean_diff"] = mean_diff.T
    print(f"importance_df type: {type(importance_df)}")
    print(f"importance_df shape: {np.array(importance_df).shape}")

    print("Saving SHAP values to file")
    # Save to TSV
    try:
        importance_df.to_csv(f"{args.output}.shap_values.tsv", sep="\t")
        #importance_df.to_csv(f"{args.output}.global.csv", sep=",", index=False)
    except:
        print("Failed to save SHAP values to file")
    #print("type:", type(shap_values))
    #print("shape:", np.array(shap_values).shape)
    #print("Plotting SHAP summary...")
    
    # Use first class if shap_values is a list
    #if isinstance(shap_values, list):
    #    shap_to_plot = shap_values[0]
    #else:
    #    shap_to_plot = shap_values
    #print(f"type: {type(shap_values)}")
    #print(f"shape: {np.array(shap_values).shape}")
    # Squeeze and plot
    #shap_to_plot = np.squeeze(shap_to_plot)
    #shap.summary_plot(shap_to_plot, data, show=False, feature_names=feature_names)
    #try:
    
    #shap.summary_plot(shap_values[], data_out, show=False, feature_names=feature_names)
    #plt.savefig(f"{args.output}.png")
    #print(f"SHAP summary saved to {args.output}.tsv and {args.output}.png")
    #except:
        #print("SHAP summary failed to save")

if __name__ == "__main__":
    main()

