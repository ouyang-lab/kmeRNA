#!/usr/bin/env python3
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import os
import argparse
import pickle
import sys
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
import numpy as np
from scipy import stats
import random
from datetime import datetime
import joblib
import csv
import time
from itertools import product, permutations, combinations_with_replacement
import traceback


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
    return np.array(embedded_data), np.array(labels)

def initialize_kmers(min_k,max_k):
    """Initialize all possible kmers with count 0"""
    bases = ['A', 'C', 'G', 'T']
    kmers = []
    for length in range(min_k, max_k + 1):
        kmers.extend([''.join(p) for p in product(bases, repeat=length)])
    return {kmer: 0 for kmer in kmers}


def get_feature_importance(model, feature_names):
    """
    Save feature importance scores from a RandomForestClassifier model into a text file.
    
    Parameters:
        - model: The trained RandomForestClassifier model.
        - feature_names: List containing names of features in the same order as they were used in training.
    """
    # Get feature importances from the model
    importance_scores = model.feature_importances_
    
    # Create a dictionary to map feature names to their importance scores
    feature_importance_dict = dict(zip(feature_names, importance_scores))
    
    # Sort the dictionary by importance scores (descending order)
    sorted_importance = sorted(feature_importance_dict.items(), key=lambda x: x[1], reverse=True)
    return sorted_importance

def random_forest_classifier(train_data, train_labels, seed_val):
    start=time.time()
    logs = []
    # Initialize the Random Forest Classifier
    rf_classifier = RandomForestClassifier(
            n_estimators=512,  #512+256
            criterion='gini',
            max_features=0.8,
            max_depth=100,
            min_impurity_decrease=0, #1e-7,
            random_state=seed_val,
            n_jobs=64)
    logs.append(f"RF Classifier intitialized {datetime.now()}")
    logs.append("Training the RF model...")
    # Train the classifier
    rf_classifier.fit(train_data, train_labels)
    logs.append("RF Training complete")
    logs.append("Predicting class 1 probabilities using RF model...")
    # Make predictions on the test data
    train_predictions = rf_classifier.predict_proba(train_data)[:, 1]
    
    stop=time.time()
    logs.append(f"RF Predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return rf_classifier, train_predictions, logs

def extra_trees_classifier(train_data, train_labels, seed_val):
    start=time.time()
    logs = []
    # Initialize the Extra Trees Classifier
    et_classifier = ExtraTreesClassifier(
            n_estimators=1024, 
            criterion='gini',
            max_features=0.8,
            max_depth=None,
            min_impurity_decrease=0,
            random_state=seed_val,
            n_jobs=64)
    
    logs.append(f"ExtraTrees Classifier intitialized {datetime.now()}")
    logs.append("Training the ET model...")
    # Train the classifier
    et_classifier.fit(train_data, train_labels)
    logs.append("ET Training complete")
    logs.append("Predicting class 1 probabilities using ET model...")
    # Make predictions on the test data
    train_predictions = et_classifier.predict_proba(train_data)[:, 1]
    
    stop=time.time()
    logs.append(f"ET predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return et_classifier, train_predictions, logs

def gradient_boosting_classifier(train_data, train_labels, seed_val):
    logs = []
    start=time.time()
    # Initialize the Gradient Boosting Classifier
    gb_classifier = GradientBoostingClassifier(
            n_estimators=768,
            max_features=0.8,
            subsample=0.8,
            max_depth=5,
            learning_rate=0.1,
            min_impurity_decrease=0,
            n_iter_no_change=5,
            random_state=seed_val)
    logs.append(f"GradientBoosting Classifier intitialized {datetime.now()}")
    logs.append("Training the GB model...")
    # Train the classifier
    gb_classifier.fit(train_data, train_labels)
    print("GB Training complete")
    print("Predicting class 1 probabilities using GB model...")
    logs.append("GB Training complete")
    logs.append("Predicting class 1 probabilities using GB model...")
    # Make predictions on the test data
    train_predictions = gb_classifier.predict_proba(train_data)[:, 1]
    stop=time.time()
    logs.append(f"GB predictions complete {datetime.now()}\t{(stop-start)/60} min elapsed")
    return gb_classifier, train_predictions, logs


def random_forest_testing(rf_classifier, test_data):
    test_predictions = rf_classifier.predict_proba(test_data)[:, 1]
    return test_predictions

def extra_trees_testing(et_classifier, test_data):
    test_predictions = et_classifier.predict_proba(test_data)[:, 1]
    return test_predictions

def gradient_boosting_testing(gb_classifier, test_data):
    test_predictions = gb_classifier.predict_proba(test_data)[:, 1]
    return test_predictions


def run_random_forest_classifier(train_data, train_labels, seed_val):
    try:
        rf_out, train_probabilities_rf, logs = random_forest_classifier(train_data, train_labels, seed_val) 
        # Model training and prediction logic
    except Exception as e:
        print(f"Error in RandomForest: {e}")
        traceback.print_exc()
    #    results["rf"] = random_forest_classifier(train_data, train_labels, test_data)
    return (rf_out, train_probabilities_rf, logs)

def run_extra_trees_classifier(train_data, train_labels, seed_val):
    try:
        et_out, train_probabilities_et, logs = extra_trees_classifier(train_data, train_labels, seed_val)
        # Model training and prediction logic
    except Exception as e:
        print(f"Error in ExtraTrees: {e}")
        traceback.print_exc()
    return (et_out, train_probabilities_et, logs)
 
def run_gradient_boosting_classifier(train_data, train_labels, seed_val):
    try:
        gb_out, train_probabilities_gb, logs = gradient_boosting_classifier(train_data, train_labels, seed_val)
    except Exception as e:
        print(f"Error in GradientBoosting: {e}")
        traceback.print_exc()
    return (gb_out, train_probabilities_gb, logs)



def wrapper_function(wrap_args):
    """Wraps the worker function and stores the result in the shared dict."""
    results = wrap_args[0]
    out_args = wrap_args[1:4]
    name_arg = wrap_args[-1]
    #results = {}
    if name_arg == "rf":
        key, value = run_random_forest_classifier(*out_args)
    elif name_arg == "et":
        key, value = run_extra_trees_classifier(*out_args)
    elif name_arg == "gb":
        key, value = run_gradient_boosting_classifier(*out_args)
    else:
        print(f"{name_arg} is not recognized")
    results[key] = value
    #return results    

def main():
    start = time.time()
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Train ensemble models on embedded data")
    parser.add_argument("--train_data", help="Path to embedded training data pickle file", required=False)
    parser.add_argument("--train_labels", help="Path to embedded training labels text file", required=False)
    parser.add_argument("--test_data", help="Path to embedded test data pickle file", required=False)
    parser.add_argument("--test_labels", help="Path to embedded test labels text file", required=False)
    parser.add_argument("--model_out", help="Path to trained model information",required=False)
    parser.add_argument("--model_in",help="Path to pre-trained model information",required=False)
    parser.add_argument("--seed_val",type=int,help="input seed value",required=False)
    parser.add_argument("--output_train_prob", help="Path to training output",required=False)
    parser.add_argument("--output_test_prob", help="Path to testing output",required=False)
    parser.add_argument("--min_k",type=int,help="Min k-mer length in bases",default=1)
    parser.add_argument("--max_k",type=int,help="Max k-mer length in bases",default=4)
    #parser.add_argument("--output_mini_prob",help="Path to mini_output",required=False)
    args = parser.parse_args()
    print(f'Random seed: {args.seed_val}')
    os.environ['PYTHONHASHSEED']=str(args.seed_val)
    random.seed(args.seed_val)
    np.random.seed(args.seed_val)
    seed_val = int(args.seed_val)

    # Load training and test embedded data
    if args.train_data:
        train_ids, train_data = load_embedded_data(args.train_data)
        #train_data = train_tmp_data.reshape(-1, train_tmp_data.shape[-1])
        train_labels = np.loadtxt(args.train_labels)
    if args.test_data:
        test_ids, test_data = load_embedded_data(args.test_data)
        #test_data = test_tmp_data.reshape(-1, train_tmp_data.shape[-1])
        test_labels = np.loadtxt(args.test_labels)
    kmer_init = initialize_kmers(args.min_k,args.max_k)
    kmer_init_pairs = initialize_kmer_pairs(args.min_k,args.max_k)
    feature_names = [kmer1+"-query" for kmer1 in kmer_init]+[kmer2+"-target" for kmer2 in kmer_init]
    if args.model_in:
        print("Loading RF models from file...") 
        rf_in = joblib.load(args.model_in+".RF.joblib")
        et_in = joblib.load(args.model_in+".ET.joblib")
        gb_in = joblib.load(args.model_in+".GB.joblib")
        if args.output_train_prob:
            print(f"Beginning training predictions...")
            train_probabilities_rf = rf_in.predict_proba(train_data)[:, 1]
            print(f"RF model predictions complete")
            train_probabilities_et = et_in.predict_proba(train_data)[:, 1]
            print(f"ET model predictions complete")
            train_probabilities_gb = gb_in.predict_proba(train_data)[:, 1]
            print(f"GB model predictions complete")
        if args.output_test_prob:
            print(f"Beginning testing predictions...")
            test_probabilities_rf = rf_in.predict_proba(test_data)[:, 1]
            print(f"RF model predictions complete")
            test_probabilities_et = et_in.predict_proba(test_data)[:, 1]
            print(f"ET model predictions complete")
            test_probabilities_gb = gb_in.predict_proba(test_data)[:, 1]
            print(f"GB model predictions complete")
    elif args.model_out:
        rf_out, train_probabilities_rf, logs_rf = run_random_forest_classifier(train_data, train_labels, seed_val)
        for elem in logs_rf:
            print(elem)
        
        et_out, train_probabilities_et, logs_et = run_extra_trees_classifier(train_data, train_labels, seed_val)
        for elem in logs_et:
            print(elem)
        
        gb_out, train_probabilities_gb, logs_gb = run_gradient_boosting_classifier(train_data, train_labels, seed_val)
        for elem in logs_gb:
            print(elem)

        print("Saving RF models...")
        joblib.dump(rf_out, args.model_out+".RF.joblib")
        joblib.dump(et_out, args.model_out+".ET.joblib")
        joblib.dump(gb_out, args.model_out+".GB.joblib")
        print("RF models Saved")
            
        print("Saving RF model importances...")
        sorted_importance_rf = get_feature_importance(rf_out, feature_names)
        importance_file_rf = args.model_out+'.RF.feature_importance.tsv'
        sorted_importance_et = get_feature_importance(et_out, feature_names)
        importance_file_et = args.model_out+'.ET.feature_importance.tsv'
        sorted_importance_gb = get_feature_importance(gb_out, feature_names)
        importance_file_gb = args.model_out+'.GB.feature_importance.tsv'
        with open(importance_file_rf,'w') as file:
            for feature, importance in sorted_importance_rf:
                file.write(f"{feature}\t{importance}\n")
        print("RF model importances written to file.")
        with open(importance_file_et,'w') as file:
            for feature, importance in sorted_importance_et:
                file.write(f"{feature}\t{importance}\n")
        print("ET model importances written to file.")
        with open(importance_file_gb,'w') as file:
            for feature, importance in sorted_importance_gb:
                file.write(f"{feature}\t{importance}\n")
        print("GB model importances written to file.")
 
#    # Create and compile the model
    
    if args.train_data and args.output_train_prob:
        scores_train_probs = []
        for pos in range(len(train_labels)):
            scores_train_probs += [(str(int(train_labels[pos])),str(train_probabilities_rf[pos]),str(train_probabilities_et[pos]),str(train_probabilities_gb[pos]))]
        print("Writing training data predictions to file...")
        with open(args.output_train_prob, 'w') as f:
            f.write(str('ground_truth\tpred_prob_RF\tpred_prob_ET\tpred_prob_GB\n'))
            writer = csv.writer(f, delimiter='\t')  # Using tab delimiter for a 4-column text file
            writer.writerows(scores_train_probs)
        print("RF training results written.")

    if args.test_data and args.output_test_prob:
        scores_test_probs = []
        for pos in range(len(test_labels)):
            scores_test_probs += [(str(int(test_labels[pos])),str(test_probabilities_rf[pos]),str(test_probabilities_et[pos]),str(test_probabilities_gb[pos]))]
        print("Writing testing data predictions to file...")
        with open(args.output_test_prob, 'w') as f:
            f.write(str('ground_truth\tpred_prob_RF\tpred_prob_ET\tpred_prob_GB\n'))
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(scores_test_probs)
        print("Ensemble model testing results written")
        end = time.time()
        print("Ensemble model time elapsed: ",round((end-start),1),"sec",round((end-start)/60,2),"min")
    
if __name__ == "__main__":
    main()

