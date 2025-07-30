#!/usr/bin/env python3
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.
import argparse
import numpy as np
import scipy
from sklearn.metrics import classification_report, roc_auc_score, roc_curve, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d 


def interpolate_roc(thresholds, tprs, fprs, target_threshold):
    """
    Interpolates TPR and 1-FPR for a given threshold using linear interpolation.
    
    Parameters:
    thresholds (list or array-like): The list of threshold values.
    tprs (list or array-like): The list of True Positive Rates corresponding to the thresholds.
    fprs (list or array-like): The list of False Positive Rates corresponding to the thresholds.
    target_threshold (float): The threshold value at which to interpolate TPR and 1-FPR.
    
    Returns:
    tuple: Interpolated TPR and 1-FPR at the given threshold.
    """
    
    if not (min(thresholds) <= target_threshold <= max(thresholds)):
        raise ValueError("The target threshold is out of the range of provided thresholds.")
    
    # Create interpolation functions
    interp_tpr = interp1d(thresholds, tprs, kind='linear')
    interp_fpr = interp1d(thresholds, fprs, kind='linear')
    
    # Interpolate TPR and FPR
    interpolated_tpr = interp_tpr(target_threshold)
    interpolated_fpr = interp_fpr(target_threshold)
    
    # Calculate 1 - FPR
    #interpolated_1_fpr = 1 - interpolated_fpr
    
    return interpolated_tpr, interpolated_fpr


def print_roc_plot(tpr,fpr,roc_auc,optimal_threshold,plot_filepath):
    plot_filepath_out = plot_filepath+"_roc.pdf"
    #idx_half = np.where(thresholds==0.5)[0][0] 
    plt.figure()
    lw = 0.5
    plt.plot(1-fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.4f)' % roc_auc)
    plt.plot([1, 0], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.plot([0, 0.5], [0.5,0.5], color='gray', lw=lw, linestyle='--')
    plt.plot([0, 1-fpr[optimal_threshold]], [0.5,0.5], color='gray', lw=lw, linestyle='--')
    plt.plot([1,1-fpr[optimal_threshold]],[1,tpr[optimal_threshold]],color='yellow', lw=lw, linestyle='--')
    plt.scatter(1-fpr[optimal_threshold], tpr[optimal_threshold], marker='o', color='red', label='Optimal Threshold')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('1 - False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    #plt.show()
    plt.savefig(plot_filepath_out,format="pdf")

def print_prc_plot(precision,recall,auprc,plot_filepath):
    #idx_half = np.where(thresholds==0.5)[0][0] 
    plot_filepath_out = plot_filepath+"_prc.pdf"
    plt.figure()
    lw = 0.5
    plt.plot(recall, precision, color='darkorange', lw=lw, label='PR curve (area = %0.4f)' % auprc)
    plt.plot([1, 0], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.plot([0, 1], [0.5,0.5], color='gray', lw=lw, linestyle='--')
    #plt.plot([0, 1-fpr[optimal_threshold]], [0.5,0.5], color='gray', lw=lw, linestyle='--')
    #plt.plot([1,1-fpr[optimal_threshold]],[1,tpr[optimal_threshold]],color='yellow', lw=lw, linestyle='--')
    #plt.scatter(1-fpr[optimal_threshold], tpr[optimal_threshold], marker='o', color='red', label='Optimal Threshold')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision Recall Curve')
    plt.legend(loc="lower right")
    #plt.show()
    plt.savefig(plot_filepath_out,format="pdf")


def calculate_optimal_threshold(thresholds,tpr,fpr):
    """
    Calculate optimal threshold given list of thresholds, tpr, fpr from roc_curve.
    Optimal threshold maximizes TPR and minimizes FPR. This is approximated by
    maximizing the difference between TPR and FPR.
    
    Parameters: thresholds, tpr, fpr

    Returns:
    float: optimal_threshold
    """
    return thresholds[np.argmin(abs((1-fpr) - tpr))]
    #return thresholds[np.argmin(np.sqrt((1-tpr)**2+(1-(1-fpr)**2)))]

def calculate_constrained_threshold(args,thresholds,tpr,fpr):
    # Calculate the Euclidean distance from (TPR, 1-FPR) to (1, 1) 
    distances = np.sqrt((1 - tpr)**2 + (1-(1-fpr)**2))

    # Define a penalty term for thresholds far from 0.5
    penalty = np.abs(thresholds - float(args.threshold_constr))

    # Combine the distance and penalty to create a modified cost function
    cost = distances + penalty

    # Find the threshold that minimizes the modified cost function
    optimal_threshold = thresholds[np.argmin(cost)]

    return optimal_threshold

def calculate_roc_auc(args,ground_truth, predicted_probabilities):
    """
    Calculate ROC AUC score given ground truth labels and predicted probabilities.

    Parameters:
    ground_truth (list): List of true labels (0 or 1).
    predicted_probabilities (list): List of predicted probabilities.

    Returns:
    float: ROC AUC score.
    """
    auc_out = roc_auc_score(ground_truth, predicted_probabilities)
    fpr, tpr, thresholds = roc_curve(ground_truth, predicted_probabilities)
    if args.threshold_constr:
        optimal_threshold = calculate_constrained_threshold(args,thresholds,tpr, fpr)
    else:
        optimal_threshold = calculate_optimal_threshold(thresholds,tpr, fpr)
    #class_report = classification_report(ground_truth, predicted_probabilities)
    #return roc_auc_score(ground_truth, predicted_probabilities)
    #ppv = precision_score(ground_truth, predicted_probabilities)
    return auc_out, fpr, tpr, thresholds, optimal_threshold

def calculate_prc_auc(ground_truth, predicted_probabilities):
    """
    Calculate ROC AUC score given ground truth labels and predicted probabilities.

    Parameters:
    ground_truth (list): List of true labels (0 or 1).
    predicted_probabilities (list): List of predicted probabilities.

    Returns:
    float: ROC AUC score.
    """
    auprc_out = average_precision_score(ground_truth, predicted_probabilities)
    precision, recall, thresholds = precision_recall_curve(ground_truth, predicted_probabilities)
    #optimal_threshold = calculate_optimal_threshold(thresholds,tpr, fpr)
    #class_report = classification_report(ground_truth, predicted_probabilities)
    #return roc_auc_score(ground_truth, predicted_probabilities)
    #ppv = precision_score(ground_truth, predicted_probabilities)
    return auprc_out, precision, recall, thresholds #, optimal_threshold

def find_index_within_tolerance(lst, value, tolerance):
    """
    Finds the index of the first element in the list that is within the given tolerance of the specified value.

    Parameters:
    lst (list of float): The list of elements to search through.
    value (float): The value to compare against.
    tolerance (float): The tolerance within which an element is considered a match.

    Returns:
    int: The index of the first matching element, or -1 if no such element is found.
    """
    for index, element in enumerate(lst):
        if abs(element - value) <= tolerance:
            return index
    return -1

def load_data(file_path):
    """
    Load data from a text file where each line contains two columns:
    ground truth label (0 or 1) and predicted probability.

    Parameters:
    file_path (str): Path to the input file.

    Returns:
    tuple: Tuple containing lists of ground truth labels and predicted probabilities.
    """
    ground_truth = []
    predicted_probabilities = []
    
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            label = columns[0]
            prob = columns[1]
            
            if label[0] != 'g':  # Skip rows with 'g' as the first character in the label
                ground_truth.append(int(label))
                predicted_probabilities.append(float(prob))
    return ground_truth, predicted_probabilities

def load_data_ens(file_path):
    """
    Load data from a text file where each line contains five columns:
    ground truth label (0 or 1) and 4 predicted probabilities.

    Parameters:
    file_path (str): Path to the input file.

    Returns:
    tuple: Tuple containing lists of ground truth labels and predicted probabilities.
    """
    ground_truth = []
    rf_predicted_probabilities = []
    et_predicted_probabilities = []
    gb_predicted_probabilities = []
    ab_predicted_probabilities = []
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            label = columns[0]
            rf_prob = columns[1]
            et_prob = columns[2]
            gb_prob = columns[3]
            #ab_prob = columns[4]
            
            if label[0] != 'g':  # Skip rows with 'g' as the first character in the label
                ground_truth.append(int(label))
                rf_predicted_probabilities.append(float(rf_prob))
                et_predicted_probabilities.append(float(et_prob))
                gb_predicted_probabilities.append(float(gb_prob))
                #ab_predicted_probabilities.append(float(ab_prob))
    #return ground_truth, rf_predicted_probabilities, et_predicted_probabilities, gb_predicted_probabilities, ab_predicted_probabilities
    return ground_truth, rf_predicted_probabilities, et_predicted_probabilities, gb_predicted_probabilities

#    with open(file_path, 'r') as file:
#        for line in file:
#            label, prob = line.strip().split()
#            if label[0] != 'g':
#                ground_truth.append(int(label))
#                predicted_probabilities.append(float(prob))
#    return ground_truth, predicted_probabilities
def calc_RF_perf(ground_truth,rf_predicted_probabilities,args):
    # Calculate RF ROC AUC score
    roc_auc, fpr, tpr, thresholds, optimal_threshold = calculate_roc_auc(args,ground_truth, rf_predicted_probabilities)
    opt_idx = np.where(thresholds==optimal_threshold)[0][0]
    auprc_out, precision, recall, prc_thresholds = calculate_prc_auc(ground_truth, rf_predicted_probabilities)
    if args.threshold:
        tpr_out, fpr_out = interpolate_roc(thresholds, tpr, fpr, args.threshold)

    auprc = average_precision_score(ground_truth, rf_predicted_probabilities)
    if args.plot_roc:
        print_roc_plot(tpr,fpr,roc_auc,opt_idx,args.plot_roc+'.rf')
        print_prc_plot(precision, recall,auprc_out,args.plot_roc+'.rf')
    N = len(ground_truth)
    # Output ROC statistics
    print(N,"predictions")
    print("RF PRC AUC:", auprc)
    print("RF ROC AUC:", roc_auc)
    if args.threshold:
        print("RF ROC selected threshold:", args.threshold)
        print("RF ROC optimal TPR:",tpr_out)
        print("RF ROC optimal 1-FPR:", 1-fpr_out)
    else:    
        print("RF ROC optimal threshold:", optimal_threshold)
        print("RF ROC optimal TPR:",tpr[opt_idx])
        print("RF ROC optimal 1-FPR:", 1-fpr[opt_idx])

def calc_ET_perf(ground_truth,et_predicted_probabilities,args):
    # Calculate ET ROC AUC score
    roc_auc, fpr, tpr, thresholds, optimal_threshold = calculate_roc_auc(args,ground_truth, et_predicted_probabilities)
    opt_idx = np.where(thresholds==optimal_threshold)[0][0]
    auprc_out, precision, recall, prc_thresholds = calculate_prc_auc(ground_truth, et_predicted_probabilities)
    if args.threshold:
        tpr_out, fpr_out = interpolate_roc(thresholds, tpr, fpr, args.threshold)

    auprc = average_precision_score(ground_truth, et_predicted_probabilities)
    if args.plot_roc:
        print_roc_plot(tpr,fpr,roc_auc,opt_idx,args.plot_roc+'.et')
        print_prc_plot(precision, recall,auprc_out,args.plot_roc+'.et')
    N = len(ground_truth)
    # Output ROC statistics
    print(N,"predictions")
    print("ET PRC AUC:", auprc)
    print("ET ROC AUC:", roc_auc)
    if args.threshold:
        print("ET ROC selected threshold:", args.threshold)
        print("ET ROC optimal TPR:",tpr_out)
        print("ET ROC optimal 1-FPR:", 1-fpr_out)
    else:    
        print("ET ROC optimal threshold:", optimal_threshold)
        print("ET ROC optimal TPR:",tpr[opt_idx])
        print("ET ROC optimal 1-FPR:", 1-fpr[opt_idx])

def calc_GB_perf(ground_truth,gb_predicted_probabilities,args):
    # Calculate GB ROC AUC score
    roc_auc, fpr, tpr, thresholds, optimal_threshold = calculate_roc_auc(args,ground_truth, gb_predicted_probabilities)
    opt_idx = np.where(thresholds==optimal_threshold)[0][0]
    auprc_out, precision, recall, prc_thresholds = calculate_prc_auc(ground_truth, gb_predicted_probabilities)
    if args.threshold:
        tpr_out, fpr_out = interpolate_roc(thresholds, tpr, fpr, args.threshold)

    auprc = average_precision_score(ground_truth, gb_predicted_probabilities)
    if args.plot_roc:
        print_roc_plot(tpr,fpr,roc_auc,opt_idx,args.plot_roc+'.gb')
        print_prc_plot(precision, recall,auprc_out,args.plot_roc+'.gb')
    N = len(ground_truth)
    # Output ROC statistics
    print(N,"predictions")
    print("GB PRC AUC:", auprc)
    print("GB ROC AUC:", roc_auc)
    if args.threshold:
        print("GB ROC selected threshold:", args.threshold)
        print("GB ROC optimal TPR:",tpr_out)
        print("GB ROC optimal 1-FPR:", 1-fpr_out)
    else:    
        print("GB ROC optimal threshold:", optimal_threshold)
        print("GB ROC optimal TPR:",tpr[opt_idx])
        print("GB ROC optimal 1-FPR:", 1-fpr[opt_idx])

def calc_AB_perf(ground_truth,ab_predicted_probabilities,args):
    # Calculate AB ROC AUC score
    roc_auc, fpr, tpr, thresholds, optimal_threshold = calculate_roc_auc(args,ground_truth, ab_predicted_probabilities)
    opt_idx = np.where(thresholds==optimal_threshold)[0][0]
    auprc_out, precision, recall, prc_thresholds = calculate_prc_auc(ground_truth, ab_predicted_probabilities)
    if args.threshold:
        tpr_out, fpr_out = interpolate_roc(thresholds, tpr, fpr, args.threshold)

    auprc = average_precision_score(ground_truth, ab_predicted_probabilities)
    if args.plot_roc:
        print_roc_plot(tpr,fpr,roc_auc,opt_idx,args.plot_roc+'.ab')
        print_prc_plot(precision, recall,auprc_out,args.plot_roc+'.ab')
    N = len(ground_truth)
    # Output ROC statistics
    print(N,"predictions")
    print("AB PRC AUC:", auprc)
    print("AB ROC AUC:", roc_auc)
    if args.threshold:
        print("AB ROC selected threshold:", args.threshold)
        print("AB ROC optimal TPR:",tpr_out)
        print("AB ROC optimal 1-FPR:", 1-fpr_out)
    else:    
        print("AB ROC optimal threshold:", optimal_threshold)
        print("AB ROC optimal TPR:",tpr[opt_idx])
        print("AB ROC optimal 1-FPR:", 1-fpr[opt_idx])




def main():
    parser = argparse.ArgumentParser(description="Calculate ROC AUC score from a 2-column input file.")
    parser.add_argument("--input", help="Path to the input file containing ground truth labels and predicted probabilities.")
    parser.add_argument("--plot_roc")
    parser.add_argument("--threshold",type=float)
    parser.add_argument("--gt",type=int)
    parser.add_argument("--pred",type=int)
    parser.add_argument("--threshold_constr")
    #parser.add_argument("-tc","--threshold_constr",action='store_true',default=False)
    args = parser.parse_args()

    # Load data from input file
    #ground_truth, predicted_probabilities = load_data(args.input)
    #ground_truth, rf_predicted_probabilities, et_predicted_probabilities, gb_predicted_probabilities, ab_predicted_probabilities = load_data_ens(args.input)
    ground_truth, rf_predicted_probabilities, et_predicted_probabilities, gb_predicted_probabilities = load_data_ens(args.input)
    calc_RF_perf(ground_truth,rf_predicted_probabilities,args)
    calc_ET_perf(ground_truth,et_predicted_probabilities,args)
    calc_GB_perf(ground_truth,gb_predicted_probabilities,args)
    #calc_AB_perf(ground_truth,ab_predicted_probabilities,args)

if __name__ == "__main__":
    main()


