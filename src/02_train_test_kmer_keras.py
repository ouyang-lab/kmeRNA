#!/usr/bin/env python3
# MIT License
# Copyright (c) 2025 Eric Nels Pederson, University of Massachusetts Amherst
# Licensed under the MIT License. See LICENSE file in the project root for full license information.

import os
import argparse
import pickle
import tensorflow as tf
import keras
from keras import initializers, layers, regularizers, optimizers, losses, models
import numpy as np
from scipy import stats
import random
from datetime import datetime
import joblib
import csv
import time
from matplotlib import pyplot as plt

# declare classes for Keras callbacks:
# for storing all of the losses during training
class BatchMetrics(keras.callbacks.Callback):
    def __init__(self):
        super().__init__()
        self.batch_losses = []
        self.epoch_steps = []
        self.current_epoch_steps = 0
        self.pred_batch_losses = []
    #  custom callback functions available:
    #def on_train_begin(self, logs=None) | fit()
    #def on_train_end(self, logs=None) | fit()
    #def on_test_begin(self, logs=None) | evaluate()
    #def on_test_end(self, logs=None) | evaluate()
    #def on_predict_begin(self, logs=None) | predict()
    #def on_predict_end(self, logs=None) | predict()
    #def on_epoch_begin(self, epoch, logs=None)
    #def on_epoch_end(self, epoch, logs=None)
    #def on_batch_begin(self, batch, logs=None)
    #def on_batch_end(self, batch, logs=None)
    def on_train_batch_end(self, batch, logs=None):
        self.batch_losses.append(logs['loss'])
        self.current_epoch_steps += 1 # increase step counter

    def on_epoch_begin(self, epoch, logs=None):
        self.current_epoch_steps = 0  # reset step counter for the new epoch

    def on_epoch_end(self, epoch, logs=None):
        self.epoch_steps.append(self.current_epoch_steps) 

    def get_steps_per_epoch(self):
        return self.epoch_steps
    
    def on_predict_batch_end(self, batch, logs=None):
        self.pred_batch_losses.append(logs['loss'])

    def get_batch_losses(self, logs=None):
        return self.batch_losses

    def get_batch_pred_losses(self, logs=None):
        return self.pred_batch_losses

class SaveModelCallback(keras.callbacks.Callback):
    def __init__(self, model_path="saved_models", model_name="RNA-RNA"):
        super().__init__()
        self.model_path = model_path
        self.model_name = model_name
        self.min_loss = float("inf")
        #os.makedirs(self.save_dir, exist_ok=True) # update later with path check

    def on_epoch_end(self, epoch, logs=None):
        if logs is None:
            return
        current_loss = logs.get("loss")
        if current_loss is not None and current_loss < self.min_loss:
            self.min_loss = current_loss
            save_path = f"{self.model_path}.{self.model_name}_epoch_{epoch+1}.keras"
            self.model.save(save_path)
            print(f"\nModel saved at: {save_path}")

class LinearDecayLR(keras.callbacks.Callback):
    def __init__(self, initial_lr=0.001, min_lr=1e-5, k=1.2):
        super(LinearDecayLR, self).__init__()
        self.min_lr = min_lr
        self.initial_lr = initial_lr
        self.k = k  # Scaling factor for controlling the rate of decay. Higher values indicate a steeper drop in LR 
        self.initial_loss = None
        self.old_lr = None
        self.time_start=time.time()
    
    def on_train_batch_end(self, batch, logs=None):
        logs = logs or {}
        loss = logs.get("val_loss", logs.get("loss"))
        if self.initial_loss is None and loss is not None:
            self.initial_loss = loss
            print(f"Initial loss set to {self.initial_loss:.7f}")
    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        current_loss = logs.get("val_loss", logs.get("loss"))
        if current_loss is None or self.initial_loss is None:
            return

        
        # Compute new learning rate
        loss_reduction = max(0, self.initial_loss - current_loss)
        new_lr = max(0, self.initial_lr - self.k * loss_reduction * self.initial_lr)
        if new_lr <= self.min_lr:
            new_lr = self.min_lr
        # Apply new learning rate
        self.model.optimizer.learning_rate.assign(new_lr)
        if self.old_lr != None:
            print(f"End of Epoch {epoch + 1}: Adjusting learning rate from {self.old_lr:.7f} to {new_lr:.7f}")
        else:
            print(f"End of Epoch {epoch + 1}: Adjusting learning rate from {self.initial_lr:.7f} to {new_lr:.7f}")
        print(f"End of Epoch {epoch + 1}: {round((time.time()-self.time_start)/60,3)} min elapsed training")
        self.old_lr = new_lr
        #self.time_end = time.time()

def plot_AUPRC(history, file_path):
    """
    Plot training and validation loss from the history of a model.fit() call and save to a PDF.
    
    Parameters:
    - history: The history object returned by model.fit().
    - file_path: The path (including file name) where the plot will be saved as a PDF.
    """
    # Extract loss and validation loss
    dat_auprc = history.history['AUPRC']
    val_dat_auprc = history.history.get('val_AUPRC')  # Use .get() to handle cases without validation
    
    # Plot the loss
    plt.figure(figsize=(8, 5))
    plt.plot(dat_auprc, label='Training AUPRC', marker='o',color='k')
    if val_dat_auprc:  # Only plot if validation loss is available
        plt.plot(val_dat_auprc, label='Validation AUPRC', marker='_',color="r")
    plt.title('AUPRC vs. Epochs')
    plt.xlabel('Epochs')
    plt.ylabel('AUPRC')
    plt.legend()
    plt.grid(True)
    
    # Save the plot to a PDF
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.close()  # Close the figure to free memory

def plot_AUROC(history, file_path):
    """
    Plot training and validation loss from the history of a model.fit() call and save to a PDF.
    
    Parameters:
    - history: The history object returned by model.fit().
    - file_path: The path (including file name) where the plot will be saved as a PDF.
    """
    # Extract loss and validation loss
    dat_auprc = history.history['AUROC']
    val_dat_auprc = history.history.get('val_AUROC')  # Use .get() to handle cases without validation
    
    # Plot the loss
    plt.figure(figsize=(8, 5))
    plt.plot(dat_auprc, label='Training AUROC', marker='o',color='k')
    if val_dat_auprc:  # Only plot if validation loss is available
        plt.plot(val_dat_auprc, label='Validation AUROC', marker='_',color="r")
    plt.title('AUROC vs. Epochs')
    plt.xlabel('Epochs')
    plt.ylabel('AUROC')
    plt.legend()
    plt.grid(True)
    # Save the plot to a PDF
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.close()  # Close the figure to free memory


def plot_loss(history, file_path):
    """
    Plot training (and validation) epoch loss from the history of a model.fit() call and save to a PDF.
    
    Parameters:
    - history: The history object returned by model.fit().
    - file_path: The path (including file name) where the plot will be saved as a PDF.
    """
    # Extract loss and validation loss
    loss = history.history['loss']
    val_loss = history.history.get('val_loss')  # Use .get() to handle cases without validation
    
    # Plot the loss
    plt.figure(figsize=(8, 5))
    plt.plot(loss, label='Training Loss', marker='o',color='k')
    if val_loss:  # Only plot if validation loss is available
        plt.plot(val_loss, label='Validation Loss', marker='_',color="r")
    plt.title('Loss vs. Epochs')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    
    # Save the plot to a PDF
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.close()  # Close the figure to free memory

def plot_batch_loss(history, file_path, pred=False):
    """
    Plot training batch loss from the custom callback of a model.fit() call and save to a PDF.
    
    Parameters:
    - history: The history object returned by BatchMetrics class.
    - file_path: The path (including file name) where the plot will be saved as a PDF.
    """
    if pred:
         # Extract loss and validation loss
        #loss = history.pred_batch_losses
        loss = history.get_batch_losses()
        # Plot the loss
        plt.figure(figsize=(8, 5))
        plt.plot(loss, label='Prediction Loss', marker='o',color='k')
        
    else:
        # Extract loss and validation loss
        #loss = history.pred_batch_losses
        loss = history.get_batch_losses()
        epochs = history.get_steps_per_epoch()
        # Plot the loss
        plt.figure(figsize=(8, 5))
        for epoch in epochs:
            plt.axvline(x=epoch, color='r', linestyle='--', alpha=0.7, label='Epoch' if epoch == epochs[0] else "")
        plt.plot(loss, label='Training Loss', marker='o',color='k')
   
    plt.title('Loss vs. Batches')
    plt.xlabel('Batches')
    plt.ylabel('Loss')
    plt.legend()
    plt.grid(True)
    
    # Save the plot to a PDF
    plt.savefig(file_path, format='pdf', bbox_inches='tight')
    plt.close()  # Close the figure to free memory


# FILE I/O
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


def create_model_optimal(input_shape,args):    
    #depth_dim = max_depth_dim(input_shape)
    #''' initialize '''    
    model = models.Sequential()
    #''' input layer. kernel and bias intialization. 1024 seems to work well. tanh is the best activation func '''
    # input_shape: (batch,2234,1)
    #if args.max_k <= 5:
    nodes1=1024
    nodes2=512
    nodes3=256
    nodes4=64
    nodes5=32
    nodes6=32
    l1_val=5e-5
    l2_val=5e-5

    model.add(layers.Input(input_shape))
    model.add(layers.Dense(nodes1, activation='tanh', kernel_initializer='ones', bias_initializer='zeros'))
    model.add(layers.Dropout(float(args.drop1/(nodes2*(nodes1+1))),seed=args.seed_val+0))
    model.add(layers.Dense(nodes2, activation='gelu', kernel_constraint=tf.keras.constraints.MaxNorm(5)))
    model.add(layers.Dropout(float(args.drop2/(nodes3*(nodes2+1))),seed=args.seed_val+1))    
    model.add(layers.Dense(nodes3, activation='gelu', kernel_constraint=tf.keras.constraints.MaxNorm(5)))
    model.add(layers.Dropout(float(args.drop3/(nodes4*(nodes3+1))),seed=args.seed_val+2))
    model.add(layers.Dense(nodes4, activation='gelu', kernel_constraint=tf.keras.constraints.MaxNorm(5)))
    model.add(layers.Dropout(float(args.drop4/(nodes5*(nodes4+1))),seed=args.seed_val+3))
    model.add(layers.Dense(nodes5, activation='gelu', kernel_regularizer=tf.keras.regularizers.L1L2(l1=l1_val,l2=l2_val)))
    model.add(layers.Flatten())
    model.add(layers.Dense(nodes6, activation='gelu'))
    model.add(layers.Dense(1, activation='sigmoid'))
    model.compile(
        optimizer=tf.keras.optimizers.Adam(
            learning_rate=0.001, 
            beta_1=0.9, 
            beta_2=0.999), 
        loss='binary_crossentropy',
        metrics=[
            tf.keras.metrics.TruePositives(thresholds=0.5,name="TP"),
            tf.keras.metrics.FalseNegatives(thresholds=0.5,name="FN"),
            tf.keras.metrics.TrueNegatives(thresholds=0.5,name="TN"),
            tf.keras.metrics.FalsePositives(thresholds=0.5,name="FP"),
            tf.keras.metrics.BinaryAccuracy(threshold=0.5,name="Acc"),
            tf.keras.metrics.AUC(curve='ROC',name="AUROC"),
            tf.keras.metrics.AUC(curve='PR',name="AUPRC")])
            #tf.keras.metrics.AUC(from_logits=False,curve='ROC',name="AUROC"),
            #tf.keras.metrics.AUC(from_logits=False,curve='PR',name="AUPRC")])
    return model

def main():
    start = time.time()
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Train CNN on embedded data")
    parser.add_argument("--train_data", help="Path to embedded training data pickle file", required=False)
    parser.add_argument("--train_labels", help="Path to embedded training labels text file", required=False)
    parser.add_argument("--test_data", help="Path to embedded test data pickle file", required=False)
    parser.add_argument("--test_labels", help="Path to embedded test labels text file", required=False)
    parser.add_argument("--model_out", help="Path to trained model information",required=False)
    parser.add_argument("--model_in",help="Path to pre-trained model information",required=False)
    parser.add_argument("--seed_val",type=int,help="input seed value",required=False)
    parser.add_argument("--output_train_prob", help="Path to training output",required=False)
    parser.add_argument("--output_test_prob", help="Path to testing output",required=False)
    parser.add_argument("--plot_path", default=False, required=False)
    parser.add_argument("--batch_size", type=int,required=False)
    parser.add_argument("--max_epochs",type=int,default=8,required=False)
    parser.add_argument("--min_k",type=int,default=1,required=False)
    parser.add_argument("--max_k",type=int,default=4,required=False)
    parser.add_argument("--cv", action='store_true')
    parser.add_argument("--drop1",type=float,default=1,help="Number of nodes (float) to drop between 1st and 2nd dense layer")
    parser.add_argument("--drop2",type=float,default=1,help="Number of nodes (float) to drop between 2nd and 3rd dense layer")
    parser.add_argument("--drop3",type=float,default=1,help="Number of nodes (float) to drop between 3rd and 4th dense layer")
    parser.add_argument("--drop4",type=float,default=1,help="Number of nodes (float) to drop between 4th and 5th dense layer")
    parser.add_argument("--task",type=str,default="RNA-RNA",choices=["RNA-RNA","miRNA-RNA","sRNA-mRNA"],required=False)
    parser.add_argument("--save_best_epochs",action='store_true',help="Save the Keras models at each epoch if the loss or validation loss is lowest so far.")
    #parser.add_argument("--output_mini_prob",help="Path to mini_output",required=False)
    args = parser.parse_args()
    
    print(''.join('NN seed: ' + str(args.seed_val)))
    os.environ['PYTHONHASHSEED']=str(args.seed_val)
    random.seed(args.seed_val)
    np.random.seed(args.seed_val)

    if args.cv and args.test_data and args.test_labels:
        VALIDATION_STEP_FLAG=True
    else:
        VALIDATION_STEP_FLAG=False
    # Load training and test embedded data
    

    # default to 64 if not given during testing
    if args.model_in and not args.batch_size:
        batch_size_out = 64
    elif args.model_in and args.batch_size:
        batch_size_out = args.batch_size

    if args.train_data:
        print("Loading training data from file...")
        train_pair_ids, train_data_in = load_embedded_data(args.train_data)
        train_labels = np.loadtxt(args.train_labels)
        train_data, train_ids = prepare_data(train_data_in, train_labels)
        print("Training data loaded")
    if (args.test_data and VALIDATION_STEP_FLAG) or (args.test_data and not args.train_data):
        print("Loading testing data from file...")
        test_pair_ids, test_data_in = load_embedded_data(args.test_data)

        test_labels = np.loadtxt(args.test_labels)
        test_data, test_ids = prepare_data(test_data_in, test_labels)
        print("Testing data loaded")
        # Create and compile the model on padded data
    if args.train_data:
        input_shape = (train_data[0].shape[0],1)
    elif args.test_data and not args.train_data:
        input_shape= (test_data[0].shape[0],1)
    if args.model_in:
        print("Loading NN model from file...")
        model = load_model_data(args.model_in+'.keras')
        print("NN model loaded")
        print("Model summary:")
        model.summary()
    elif args.model_out:
        if args.task=="RNA-RNA":
            if not args.batch_size:
                batch_size_out=48
            else:
                batch_size_out=args.batch_size
            model = create_model_optimal(input_shape,args)
            print("RNA-RNA NN model initialized")
        elif args.task=="miRNA-RNA":
            if not args.batch_size:
                batch_size_out=64
                print(f"Default batch size: {batch_size_out}")
                print("Consider adjusting this parameter using the --batch_size argument")
            else:
                batch_size_out=args.batch_size
            model = create_model_optimal(input_shape,args)
            print("miRNA-RNA NN model initialized")
        elif args.task=="sRNA-RNA":
            if not args.batch_size:
                batch_size_out=80
                print(f"Default batch size: {batch_size_out}")
                print("Consider adjusting this parameter using the --batch_size argument")
            else:
                batch_size_out=args.batch_size
            model = create_model_optimal(input_shape,args)
            print("sRNA-RNA NN model initialized")

    # Create and compile the model
    SHUFFLE_FLAG=True # 'batch', True or False
    if args.cv and args.test_data and args.test_labels:
        early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss',
                restore_best_weights=True,
                patience=2,
                verbose=1)
    else:
        early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='loss',
                restore_best_weights=True,
                patience=2,
                verbose=1)
    if args.save_best_epochs:
        callback_list = [
                early_stop_callback, # built in early stopping
                BatchMetrics(), # saves and plots batch loss during training and predicting
                LinearDecayLR(), # slowly decreases LR over time for stability and convergence
                SaveModelCallback(model_path=args.model_out, model_name=args.task)] # saves model to file if loss if loss is less than minimum
    else:
                callback_list = [
                    early_stop_callback, # built in early stopping
                    BatchMetrics(), # saves and plots batch loss during training and predicting
                    LinearDecayLR()]
    if not args.model_in:
        # Train the model
        print("Training the NN model...")
        if VALIDATION_STEP_FLAG:
            history = model.fit(train_data, train_ids, epochs=args.max_epochs, batch_size=batch_size_out,callbacks=callback_list, shuffle=SHUFFLE_FLAG,validation_data=(test_data, test_ids),verbose=2)
        else:
    
            history = model.fit(train_data, train_ids, epochs=args.max_epochs, batch_size=batch_size_out,callbacks=callback_list, shuffle=SHUFFLE_FLAG,verbose=2)
        #print("NN model training complete")
        
    #model.summary()
        if args.plot_path:
            plot_loss(history, args.plot_path+'.loss.pdf')
            plot_batch_loss(BatchMetrics(), args.plot_path+'.batch_loss.pdf')
            plot_AUROC(history,args.plot_path+'.auroc.pdf')
            plot_AUPRC(history,args.plot_path+'.auprc.pdf')
            print("NN model training history saved")
        print("NN model training complete")
 
    if args.model_out and args.train_data and not args.model_in:
        print("Saving NN model...")
        save_model_data(args.model_out+'.keras',model)
        print("NN model saved")
 
#    # Create and compile the model
    
    # Predict probabilities on test data

    if args.output_train_prob:
        print("Predicting NN training data probabilities")
        train_probabilities = model.predict(train_data,batch_size=batch_size_out,verbose=2)
        
        scores_train_probs = []
        prob = str()
        for pos in range(len(train_labels)):
            prob = str(float(train_probabilities[pos][0]))
            scores_train_probs.append((str(int(train_labels[pos])),prob))
        print("Writing NN training data predictions to file...")
        with open(args.output_train_prob, 'w') as f:
            f.write(str('ground_truth\tpred_prob\n'))
            writer = csv.writer(f, delimiter='\t')  # Using tab delimiter for a 2-column text file
            writer.writerows(scores_train_probs)
        print("NN training results written.")
    if args.test_data and args.train_data:
        print("Loading testing data from file...")
        test_pair_ids, test_data_in = load_embedded_data(args.test_data)

        test_labels = np.loadtxt(args.test_labels)
        test_data, test_ids = prepare_data(test_data_in, test_labels)
        print("Testing data loaded")

    if args.output_test_prob:
        print("Predicting NN testing data probabilities")
        test_probabilities = model.predict(test_data,batch_size=batch_size_out,verbose=2)
        scores_test_probs = []
        prob = str()
        for pos in range(len(test_labels)):
            prob = str(float(test_probabilities[pos][0]))
            scores_test_probs.append((str(int(test_labels[pos])),prob))
        print("Writing NN testing data predictions to file...")
        with open(args.output_test_prob, 'w') as f:
            f.write(str('ground_truth\tpred_prob\n'))
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(scores_test_probs)
        print("NN testing results written")
    end = time.time()
    print("NN time elapsed: ",round((end-start),1),"sec",round((end-start)/60,2),"min")

if __name__ == "__main__":
    main()

