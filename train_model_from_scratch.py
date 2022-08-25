#!/usr/bin/env python
# coding: utf-8

# import basic modules
import os
from datetime import datetime
from joblib import dump
from glob import glob
import numpy as np
import pandas as pd
from datetime import datetime
import seaborn as sn
from sklearn.utils.validation import FLOAT_DTYPES
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import keras
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
import argparse
from aggregate_results import aggregate_results_genomic_dimension
from __init__ import __version__

###--- utils functions ---###
# function to plot a complete confusion matrix with common metrics.
def make_confusion_matrix(cf,
                          group_names=None,
                          categories='auto',
                          count=True,
                          percent=True,
                          cbar=True,
                          xyticks=True,
                          xyplotlabels=True,
                          sum_stats=True,
                          figsize=None,
                          cmap='Blues',
                          title=None,
                          path=None):
    '''
    ###############################################################################################
    CITATION: taken from: https://github.com/DTrimarchi10/confusion_matrix/blob/master/cf_matrix.py
    ###############################################################################################
    
    This function will make a pretty plot of an sklearn Confusion Matrix cm using a Seaborn heatmap visualization.
    Arguments
    ---------
    cf:            confusion matrix to be passed in
    group_names:   List of strings that represent the labels row by row to be shown in each square.
    categories:    List of strings containing the categories to be displayed on the x,y axis. Default is 'auto'
    count:         If True, show the raw number in the confusion matrix. Default is True.
    normalize:     If True, show the proportions for each category. Default is True.
    cbar:          If True, show the color bar. The cbar values are based off the values in the confusion matrix.
                   Default is True.
    xyticks:       If True, show x and y ticks. Default is True.
    xyplotlabels:  If True, show 'True Label' and 'Predicted Label' on the figure. Default is True.
    sum_stats:     If True, display summary statistics below the figure. Default is True.
    figsize:       Tuple representing the figure size. Default will be the matplotlib rcParams value.
    cmap:          Colormap of the values displayed from matplotlib.pyplot.cm. Default is 'Blues'
                   See http://matplotlib.org/examples/color/colormaps_reference.html
                   
    title:         Title for the heatmap. Default is None.
    '''


    # CODE TO GENERATE TEXT INSIDE EACH SQUARE
    blanks = ['' for i in range(cf.size)]

    if group_names and len(group_names)==cf.size:
        group_labels = ["{}\n".format(value) for value in group_names]
    else:
        group_labels = blanks

    if count:
        group_counts = ["{0:0.0f}\n".format(value) for value in cf.flatten()]
    else:
        group_counts = blanks

    if percent:
        group_percentages = ["{0:.2%}".format(value) for value in cf.flatten()/np.sum(cf)]
    else:
        group_percentages = blanks

    box_labels = [f"{v1}{v2}{v3}".strip() for v1, v2, v3 in zip(group_labels,group_counts,group_percentages)]
    box_labels = np.asarray(box_labels).reshape(cf.shape[0],cf.shape[1])


    # CODE TO GENERATE SUMMARY STATISTICS & TEXT FOR SUMMARY STATS
    if sum_stats:
        #Accuracy is sum of diagonal divided by total observations
        accuracy  = np.trace(cf) / float(np.sum(cf))

        #if it is a binary confusion matrix, show some more stats
        if len(cf)==2:
            #Metrics for Binary Confusion Matrices
            precision = cf[1,1] / sum(cf[:,1])
            recall    = cf[1,1] / sum(cf[1,:])
            f1_score  = 2*precision*recall / (precision + recall)
            stats_text = "\n\nAccuracy={:0.3f}\nPrecision={:0.3f}\nRecall={:0.3f}\nF1 Score={:0.3f}".format(
                accuracy,precision,recall,f1_score)
        else:
            stats_text = "\n\nAccuracy={:0.3f}".format(accuracy)
    else:
        stats_text = ""


    # SET FIGURE PARAMETERS ACCORDING TO OTHER ARGUMENTS
    if figsize==None:
        #Get default figure size if not set
        figsize = plt.rcParams.get('figure.figsize')

    if xyticks==False:
        #Do not show categories if xyticks is False
        categories=False


    # MAKE THE HEATMAP VISUALIZATION
    plt.figure(figsize=figsize)
    sn.heatmap(cf,annot=box_labels,fmt="",cmap=cmap,cbar=cbar,xticklabels=categories,yticklabels=categories)

    if xyplotlabels:
        plt.ylabel('True label')
        plt.xlabel('Predicted label' + stats_text)
    else:
        plt.xlabel(stats_text)
    
    if title:
        plt.title(title)
    
    if path:
        plt.tight_layout()
        plt.savefig(path)


def retrieve_dataset(currents_features_folderpath):
    # retrieve and concat TT and CC training currents dataset 
    df_CC = pd.read_table(glob(os.path.join(currents_features_folderpath, "*CC*"))[0], header=None)
    df_CC["Label"] = ["CCcontext" for i in range(df_CC.shape[0])]
    df_CC.columns = ["region", "pos_1_based", "read_name", "strand", "cur-3", "cur-2", "cur-1", "cur0", "cur+1", "cur+2", "cur+3", "dw-3", "dw-2", "dw-1", "dw0", "dw+1", "dw+2", "dw+3", "label"]
    df_TT = pd.read_table(glob(os.path.join(currents_features_folderpath, "*TT*"))[0], header=None)
    df_TT["label"] = ["TTcontext" for i in range(df_TT.shape[0])]
    df_TT.columns = ["region", "pos_1_based", "read_name", "strand", "cur-3", "cur-2", "cur-1", "cur0", "cur+1", "cur+2", "cur+3", "dw-3", "dw-2", "dw-1", "dw0", "dw+1", "dw+2", "dw+3", "label"]
    df = pd.concat([df_CC, df_TT])
    df.fillna(np.float64(0), inplace=True)
    df.drop_duplicates(inplace=True, ignore_index=True)
    X = df.iloc[:,4:-1].copy().values
    X = np.asarray(X).astype(np.float64)
    y_info = df.loc[:,["region", "pos_1_based", "read_name", "strand", "label"]].copy().values

    # retrieve CTcontext features
    df_CT = pd.read_table(glob(os.path.join(currents_features_folderpath, "*CT*"))[0], header=None)
    df_CT["label"] = ["CTcontext" for i in range(df_CT.shape[0])]
    df_CT.columns = ["region", "pos_1_based", "read_name", "strand", "cur-3", "cur-2", "cur-1", "cur0", "cur+1", "cur+2", "cur+3", "dw-3", "dw-2", "dw-1", "dw0", "dw+1", "dw+2", "dw+3", "label"]
    df_CT.fillna(np.float64(0), inplace=True)
    df_CT.drop_duplicates(inplace=True, ignore_index=True)
    X_CT = df_CT.iloc[:,4:-1].copy().values
    X_CT = np.asarray(X_CT).astype(np.float64)
    y_info_CT = df_CT.loc[:,["region", "pos_1_based", "read_name", "strand", "label"]].copy().values

    return X, y_info, X_CT, y_info_CT



def train_model_from_scratch(model_name, features_folderpath, bam_filepath, Tproba_thr=0.8, ref_filepath=None, organism=None):
    '''
    Main function to train a model from scratch.
    Inputs:
        - model_name: a string indicating the name to use when saving the model.
        - features_folderpath: a string indicating the full path for the currents and dwells times extracted by the currents_dwell_retriever_multiprocessing.py script.
    '''
    start_time = datetime.now()
    ###--- define inputs ---###
    model_n = model_name

    features_folderpath = features_folderpath

    # retrieve sample name
    sample_name = os.path.basename(features_folderpath).split(".")[0]

    ###-- define outputs ---###
    model_path = os.path.join(os.path.dirname(features_folderpath), f"{sample_name}.model_CNN_{sample_name}_{model_n}")
    if not os.path.exists(model_path):
        os.mkdir(model_path)

    # define model name
    model_name = f'model_CNN_{sample_name}_{model_n}_{datetime.now().strftime("%d_%m_%Y_%H_%M_%S")}_from_scratch.h5'


    # ### sample TRAINING
    ###--- START LOADING OF DATA ---###

    print(f"Retrieving datasets from {features_folderpath}.", flush=True)
    X, y_info, X_CT, y_info_CT = retrieve_dataset(features_folderpath)
    print("\nX shape:", X.shape, flush=True)
    print("y (with info) shape:", y_info.shape, flush=True)
    print("\nX:\n", X, flush=True)
    print("\ny (with info):\n", y_info, flush=True)

    print("\nX CT shape:", X_CT.shape, flush=True)
    print("y CT (with info) shape:", y_info_CT.shape, flush=True)
    print("\nX CT:\n", X_CT, flush=True)
    print("\ny CT (with info):\n", y_info_CT, flush=True)

    print("\n\n", flush=True)

    print("\n### CC-TT reads Dataset (training dataset) Stats:\n", flush=True)
    print("Region counts:", np.unique(y_info[:,0], return_counts=True), flush=True)
    print("Region-Position couple counts:", pd.DataFrame(y_info[:,:2]).value_counts().shape[0], flush=True)
    print("Unique Reads counts:", pd.DataFrame(y_info[:,2]).value_counts().shape[0], flush=True)
    print("Cotenxt counts dataset:", np.unique(y_info[:,4], return_counts=True), flush=True)

    print("\n### CT reads Dataset Stats:\n", flush=True)
    print("Region counts:", np.unique(y_info_CT[:,0], return_counts=True), flush=True)
    print("Region-Position couple counts:", pd.DataFrame(y_info_CT[:,:2]).value_counts().shape[0], flush=True)
    print("Unique Reads counts:", pd.DataFrame(y_info_CT[:,2]).value_counts().shape[0], flush=True)
    print("Cotenxt counts dataset:", np.unique(y_info_CT[:,4], return_counts=True))

    # Standardizing data
    sc = StandardScaler()
    X = sc.fit_transform(X)
    # transforming into tensor (1st feature currents, 2nd feature dwell times)
    X = np.stack([X[:,:7], X[:,7:]],axis=2)
    print("\nStandardizing X dataset and transforming it into 3D tensor (1st feature currents, 2nd feature dwell times).", flush=True)
    print("X tensor shape:", X.shape, flush=True)



    # transforming y_info into a vector of categorically encoded label
    print("\nEncoding label as OHE version and appending to y_info array (0 class CCcontext, 1 class TTcontext)", flush=True)
    y_cat = tf.keras.utils.to_categorical([1 if (i=="TTcontext") else 0 for i in y_info[:,4]])
    y_info = np.concatenate([y_info, y_cat], axis=1)
    print("y (with info and OHE Label):\n", y_info, "\n", flush=True)


    # splitting into train, validation and test sets
    # produce test set
    print("\nProducing Training, Validation and Test dataset (for CC and TT reads).", flush=True)
    X_train, X_test, y_train_info, y_test_info = train_test_split(X, y_info, random_state=42, test_size=0.4, shuffle=True)
    # produce train and validation sets
    X_train, X_val, y_train_info, y_val_info = train_test_split(X_train, y_train_info, random_state=42, test_size=0.4, shuffle=True)

    # extract only OHE labels
    y_train = y_train_info[:, 5:]
    y_val = y_val_info[:, 5:]
    y_test = y_test_info[:, 5:]

    # ensure correct float type to conversion to tensor
    y_train = np.asarray(y_train).astype(np.float64)
    y_val = np.asarray(y_val).astype(np.float64)
    y_test = np.asarray(y_test).astype(np.float64)

    # delete old data to spare memory
    del(X)
    del(y_info)

    # subsample X_train for training purposes
    print("\nSampling training data..." , flush=True)
    if X_train.shape[0] > 400000:
        X_train = X_train[100000:400000]
        y_train = y_train[100000:400000]
        y_train_info = y_test_info[100000:400000]

    # print some information
    print("X train shape:", X_train.shape, flush=True)
    print("y_train shape:", y_train.shape, flush=True)
    print("X train region composition:", np.unique(y_train_info[:,0], return_counts=True), flush=True)
    
    print("\nX val shape:", X_val.shape, flush=True)
    print("y_val shape:", y_val.shape, flush=True)
    print("X val region composition:", np.unique(y_val_info[:,0], return_counts=True), flush=True)
      
    print("\nX_test shape", X_test.shape, flush=True)
    print("y test shape:", y_test.shape, flush=True)
    print("X test region composition:", np.unique(y_test_info[:,0], return_counts=True), flush=True)

    print(flush=True)
    print("y_train with OHE:\n", y_train, flush=True)
    print(flush=True)
    print("\ny_val with OHE:\n", y_val, flush=True)
    print(flush=True)
    print("\ny_test with OHE:\n", y_test, flush=True)

    # select a subset of the validation set
    if X_val.shape[0] > 150000:
        X_val_subset = X_val[100000:150000]
        y_val_subset = y_val[100000:150000]
    else:
        X_val_subset = X_val[:X_val.shape[0]]
        y_val_subset = y_val[:y_val.shape[0]]

    ###--- START TRAINING THE MODEL ---###

    # define model hyperparameters #################################################################################################################################
    n_epoch = 1000
    batch_size = 1024
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.08)
    loss="categorical_crossentropy"
    n_filters=20
    kernel_size=2
    dropout=0.2
    patience = 100
    min_delta=0.01
    train_metric = "accuracy"
    val_metric = "val_accuracy"

    # retrieve number of classes
    n_classes = np.unique(y_train).shape[0]

    # WaveNet 1dCNN classifier with causal dilation for two classes
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.InputLayer(input_shape=[X_train.shape[1], X_train.shape[2]]))
    for rate in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512) * 1:
        model.add(tf.keras.layers.Conv1D(filters=n_filters, kernel_size=kernel_size, activation="relu", 
                                    padding="causal", dilation_rate=rate))
    model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.MaxPooling1D(pool_size=2))
    model.add(tf.keras.layers.Conv1D(filters=50, kernel_size=2))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=2))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(50, activation='relu'))
    model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(20, activation='relu'))
    model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dense(5, activation='relu'))
    model.add(tf.keras.layers.BatchNormalization())
    model.add(tf.keras.layers.Dropout(dropout))
    model.add(tf.keras.layers.Dense(n_classes, activation='softmax'))
    model.compile(loss=loss, optimizer=optimizer, metrics=[train_metric])
    # print model summary
    print(flush=True)
    print(model.summary(), flush=True)

    # setting callbacks
    single_model_path = os.path.join(model_path, model_name)

    my_callbacks = [keras.callbacks.EarlyStopping(patience=patience,  monitor=val_metric, min_delta=min_delta),
                    keras.callbacks.ModelCheckpoint(single_model_path,
                    monitor=val_metric, verbose=1, 
                    save_best_only=True, mode='max')]

    # train the model
    history = model.fit(X_train, y_train, epochs=n_epoch, batch_size=batch_size, verbose=1, callbacks=my_callbacks,
                    validation_data=(X_val_subset, y_val_subset))

    # save the last model and the standard scaler instance
    dump(sc, os.path.join(model_path, model_name.split(".")[:-1][0]+"standard_scaler_final.joblib"))
    model.save(os.path.join(model_path, model_name.split(".")[:-1][0]+"_final.h5"))

    # load the best saved model
    model_best = keras.models.load_model(single_model_path)

    # plot history
    f, ax = plt.subplots(1,2, figsize=(15,5))
    ax[0].plot(history.history["loss"], label="Train. Loss")
    ax[0].plot(history.history["val_loss"], label="Val. Loss")
    ax[0].set_xlabel("Epochs")
    ax[0].set_ylabel("Loss")
    ax[0].set_title("Loss in Training and Validation datasets")
    ax[0].set_yscale("log")
    ax[0].legend()
    ax[1].plot(history.history[train_metric], label=f"Train. {train_metric}")
    ax[1].plot(history.history[val_metric], label=f"Val. {val_metric}")
    ax[1].set_xlabel("Epochs")
    ax[1].set_ylabel("Accuracy")
    ax[1].set_ylim(-0.1 ,1.1)
    ax[1].set_title("Accuracy in Training and Validation datasets")
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(model_path, f"{model_name}_Loss_Accuracy_BEST.tiff"))
    plt.show()



    # predict classes on training set and evaluate model perfomances
    labels = ["True Neg", "False Pos", "False Neg", "True Pos"]
    categories = ["CCcontext", "CTcontext"]

    y_train_hat = model_best.predict(X_train)
    print("\nX_train predicted classes count:", np.unique([np.argmax(i) for i in y_train_hat], return_counts=True), flush=True)
    print("\nAccuracy on training set is:", accuracy_score([np.argmax(i) for i in y_train], [np.argmax(i) for i in y_train_hat] ), flush=True)
    print("\nConfusion Matrix on training set:\n", confusion_matrix([np.argmax(i) for i in y_train], [np.argmax(i) for i in y_train_hat] ), flush=True)
    print("\nClassification Report on training set:\n", classification_report([np.argmax(i) for i in y_train], [np.argmax(i) for i in y_train_hat], zero_division=1), flush=True)

    make_confusion_matrix(confusion_matrix([np.argmax(i) for i in y_train], [np.argmax(i) for i in y_train_hat]), 
                        group_names=labels,
                        categories=categories,
                        title = "Confusion Matrix on Training Set",
                        figsize=(6,5),
                        path=os.path.join(model_path, f"{model_name}_Conf_Matrix_BEST_TRAIN.tiff"))



    ###--- Evaluate on validation set ---###
    y_val_hat = model_best.predict(X_val)
    print("\nX_val predicted classes count:", np.unique([np.argmax(i) for i in y_val_hat], return_counts=True), flush=True)
    print("\nAccuracy on validation set is:", accuracy_score([np.argmax(i) for i in y_val], [np.argmax(i) for i in y_val_hat] ), flush=True)
    print("\nConfusion Matrix on validation set:\n", confusion_matrix([np.argmax(i) for i in y_val], [np.argmax(i) for i in y_val_hat] ), flush=True)
    print("\nClassification Report on validation set:\n", classification_report([np.argmax(i) for i in y_val], [np.argmax(i) for i in y_val_hat], zero_division=1), flush=True)

    make_confusion_matrix(confusion_matrix([np.argmax(i) for i in y_val], [np.argmax(i) for i in y_val_hat]), 
                        group_names=labels,
                        categories=categories,
                        title = "Confusion Matrix on validation Set",
                        figsize=(6,5),
                        path=os.path.join(model_path, f"{model_name}_Conf_Matrix_BEST_VAL.tiff"))



    ###--- Evaluate on Test set ---###
    y_test_hat = model_best.predict(X_test)
    print("\nX_test predicted classes count:", np.unique([np.argmax(i) for i in y_test_hat], return_counts=True), flush=True)
    print("\nAccuracy on test set is:", accuracy_score([np.argmax(i) for i in y_test], [np.argmax(i) for i in y_test_hat] ), flush=True)
    print("\nConfusion Matrix on test set:\n", confusion_matrix([np.argmax(i) for i in y_test], [np.argmax(i) for i in y_test_hat] ), flush=True)
    print("\nClassification Report on test set:\n", classification_report([np.argmax(i) for i in y_test], [np.argmax(i) for i in y_test_hat], zero_division=1), flush=True)

    make_confusion_matrix(confusion_matrix([np.argmax(i) for i in y_test], [np.argmax(i) for i in y_test_hat]), 
                        group_names=labels,
                        categories=categories,
                        title = "Confusion Matrix on Test Set",
                        figsize=(6,5),
                        path=os.path.join(model_path, f"{model_name}_Conf_Matrix_BEST_TEST.tiff"))

    print("Training and Testing finished.", flush=True)



    # ### Evaluate on CT context reads
    # standardize and transform CT reads to 3D tensor

    X_CT_std = sc.transform(X_CT)
    X_CT_std = np.stack([X_CT_std[:,:7], X_CT_std[:,7:]], axis=2)
    X_CT_std



    # predict CT context reads as CC or TT context related read
    y_CT_hat_proba = model_best.predict(X_CT_std)
    y_CT_hat_proba



    # consolidate prediction results (Cproba and Tproba) of CT reads sites into a pandas dataframe
    df_CT = pd.DataFrame(y_info_CT[:,:-1])
    df_CT.columns = ["region", "position", "read_name", "strand"]
    df_CT["Cproba"] = y_CT_hat_proba[:,0]
    df_CT["Tproba"] = y_CT_hat_proba[:,1]
    print("\nResulting CT dataframe:\n", df_CT, flush=True)

    # saving to disk final df_CT
    final_df_path = os.path.join(model_path, 'df_CT_predicted.tsv')
    print(f"Saving final predicted dataframe to {final_df_path}." ,flush=True)
    df_CT.to_csv(final_df_path, sep="\t", header=True, index=None)

    # aggregating results on genomic space and saving aggregated results to disk
    aggregate_results_genomic_dimension(df_CT_predicted_filepath=final_df_path, 
                                        bam_filepath=bam_filepath, 
                                        Tproba_thr=Tproba_thr,
                                        ref_filepath=ref_filepath,
                                        model_type="CNN",
                                        organism=organism)
                                        
    print(f"Total Elapsed Time: {datetime.now() - start_time}", flush=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Train a CNN Wavenet model from scratch using CC and TT contexts reads currents and dwell times.")
    parser.add_argument("-mn", 
                        required=True,
                        type=str,
                        help="-model_name: \t a <str> indicating the model name to use when saving final and best models.")
    parser.add_argument("-ffp",
                        required=True,
                        type=str,
                        help="-features_folderpath: \t a <str> indicationg the full path for the folder produced by currents_dwells_retriever script where the currents and dwells time of each read are stored per position.")
    parser.add_argument("-bam",
                        required=True,
                        type=str,
                        help="-bam_file_fullpath: \t a <str> with the full_path for the input bam file")

    parser.add_argument("-Tproba_thr",
                        required=False,
                        type=float,
                        default=0.8,
                        help="-T_proba_threshold_value: \t a <float> indicating a probability threshold to filter out CT reads with lower probability given by model.")

    parser.add_argument("-ref_filepath",
                        required=False,
                        default=None,
                        type=str,
                        help="-reference_filepath: \t a <str> indicating the full-path for reference used for reads alignment.")
    
    parser.add_argument("-O",
                        "--organism",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-organism: \t a <str> indicating the organims analyzed in order to produce a p-value of the APOBEC1 signature resemblance after aggregation of genome space of the sites.")
    
    
    args = parser.parse_args()
    model_name = args.mn
    features_folderpath = args.ffp
    bam_filepath = args.bam
    Tproba_thr = args.Tproba_thr
    ref_filepath = args.ref_filepath
    organism = args.organism
    
    # print some starting info related to the C_to_U_classifier version, to the used program and to the input arguments
    print(f"[{datetime.now()}] C_to_U_classifier version: {__version__}", flush=True)
    print(f"[{datetime.now()}] Train CNN-Wavenet inner model from scratch on currents features. Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    train_model_from_scratch(model_name=model_name, 
                            features_folderpath=features_folderpath, 
                            bam_filepath=bam_filepath,
                            Tproba_thr=Tproba_thr,
                            ref_filepath=ref_filepath,
                            organism=organism)