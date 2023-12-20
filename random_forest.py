"""
Implement the random forest method used in the me-class paper
for the comparison
"""

import sys, os, argparse, math
import numpy as np
import pickle
import pandas as pd
import scipy.stats
import scipy.interpolate
import sklearn.ensemble
from sklearn import metrics
import joblib
from matplotlib import pyplot as plt
import seaborn as sns

# Down sample the data by averaging the value in size r window
def down_sample(x, r):
    # x: input 2D numpy array
    # r: ratio of down sampling
    ns = x.shape[1] // r       # number of samples
    #print(f"ns = {ns}")
    # Remove remainder columns more than ns * r
    x = x[:, 0:int(ns*r)]
    # Average every r elements
    x = x.reshape(x.shape[0], -1, r).mean(axis=2)
    return x


# Down sample the data by getting the max value in size r window
def down_sample_max(x, r):
    # x: input 2D numpy array
    # r: ratio of down sampling
    ns = x.shape[1] // r  # number of samples
    # Remove remainder columns more than ns * r
    x = x[:, 0:int(ns * r)]
    # Average every r elements
    x = x.reshape(x.shape[0], -1, r).max(axis=2)
    return x


# Down sample the data by getting the min value in size r window
def down_sample_min(x, r):
    # x: input 2D numpy array
    # r: ratio of down sampling
    ns = x.shape[1] // r  # number of samples
    # Remove remainder columns more than ns * r
    x = x[:, 0:int(ns * r)]
    # Average every r elements
    x = x.reshape(x.shape[0], -1, r).min(axis=2)
    return x


def data_preprocess(x, y, r=1):
    np.nan_to_num(x, copy=False)
    # take average of every ? bases
    if r > 1:
        x_avg = down_sample(x, r)
        #x_max = down_sample_max(x, r)
        #x_min = down_sample_min(x, r)
        #x = np.concatenate([x_avg, x_max, x_min], axis=1)
        x = x_avg
        #x = np.expand_dims(x, 2)

    # Quantize the x values
    #x = (x*10.0 + 0.5).astype(np.int32)
    #x = x.astype(np.float32) / 10.
    x_q = (np.abs(x) * 10 + 0.5).astype(np.int32)
    x_q = x_q.astype(np.float32) / 10.
    x = np.sign(x) * x_q
    # covert negative LFC to cat 0 and positive to cat 1
    y = (
        np.right_shift((np.sign(y).astype(np.int32) + 1), 1)
            .flatten()
    )
    return x, y


# me-class: Read .label file
def load_multilabel_expr(filename):
    fh = open(filename,'r')
    Ex = list()
    L = list()
    E = list()
    Y = list()
    I = list()
    D = list()
    header = ""
    for line in fh:
        if line.startswith("#"):
            ll = line.lstrip("#").rstrip().split()
            ii = ll.index("GENE_ID")
            ei = ll.index("EXPR")
            xi = ll.index("NUM_EXONS")
            pli = ll.index("POS_LOW")
            phi = ll.index("POS_HIGH")
            header = "#PROB_DW\tPROB_UP\t"+line.lstrip("#")
            continue
        ll = line.strip().split()
        D.append(ll)
        id_name = ll[ii]
        expr = float(ll[ei])
        E.append(int(ll[xi]))
        L.append(math.log(int(ll[phi])-int(ll[pli]),10))
        Ex.append(expr)
        I.append(id_name)
        if expr > 0:
            Y.append(1)
        else:
            Y.append(-1)
    return Y, L, E, I, D, header


def load_vals(filename):
    fh = open(filename, 'r')
    C = list()
    for line in fh:
        if line.startswith("#"):
            continue
        ll = [float(x) for x in line.strip().split()]
        C.append(ll)
    return C

def load_meclass_sample(sample_name, meclass_basedir="../me-class/example_dataset/"):
    sys.stderr.write("Loading sample "+sample_name+"\n")
    label_file = meclass_basedir+sample_name+".label"
    Y, L, E, I, D, header = load_multilabel_expr(label_file)
    data_file = meclass_basedir+sample_name+".tss.meth.dat"
    X = load_vals(data_file)
    return X, Y


if __name__ == "__main__":
    print("Random Forest Classifier")

    parser = argparse.ArgumentParser(
        description="Random forest model for analyzing methylation and LFC data from human epigenmoics project"
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="The input .pickle file that stores the data",
        default="./data_linear.pickle",
    )
    parser.add_argument(
        "--down-sampling",
        dest="down_sampling",
        type=int,
        default=1,
        help="Ratio of down sampling the methylation data",
    )
    # load prepared data
    args = parser.parse_args()
    pickle_file = args.input_file
    with open(pickle_file, "rb") as fp:
        data = pickle.load(fp)
    # covert NaN values to 0
    x_train = data["train_dataset"].astype(np.float32)
    x_valid = data["valid_dataset"].astype(np.float32)
    x_test = data["test_dataset"].astype(np.float32)
    y_train = data["train_labels"]
    y_valid = data["valid_labels"]
    y_test = data["test_labels"]

    x_train, y_train = data_preprocess(x_train, y_train, args.down_sampling)
    x_valid, y_valid = data_preprocess(x_valid, y_valid, args.down_sampling)
    x_test, y_test = data_preprocess(x_test, y_test, args.down_sampling)

    print(x_train.shape, y_train.shape)
    print(x_train[0:5, 0:5])
    print(np.min(x_train), np.max(x_train))
    print(np.min(y_train), np.max(y_train))

    # Read example data from me-class
    # We will use E095_E096 as testing, and the rest as training
    # in order to see how well it performs in a random forest classifier
    """sample_names = ["E096_E097", "E071_E079", "E094_E095", "E095_E096"]
    test_sample = "E094_E095"
    x_train = list()
    y_train = list()
    x_valid = list()
    y_valid = list()
    for s in sample_names:
        X, Y = load_meclass_sample(s, meclass_basedir="./out/")
        if s == test_sample:
            x_valid += X
            y_valid += Y
        else:
            x_train += X
            y_train += Y
    x_train = np.asarray(x_train)
    y_train = np.asarray(y_train)
    x_valid = np.asarray(x_valid)
    y_valid = np.asarray(y_valid)
"""
    rf = sklearn.ensemble.RandomForestClassifier(n_estimators=101, n_jobs=-1)
    #rf = sklearn.tree.DecisionTreeClassifier()
    #rf = joblib.load('rf_model.joblib')
    # Train the random forest classifier on the training data
    # and test it on the validation data
    rf.fit(x_train, y_train)

    # save the RF model to a file
    #joblib.dump(rf, 'rf_model.joblib')

    y_train_predict = rf.predict(x_train)
    print(metrics.confusion_matrix(y_train, y_train_predict))
    print("training accuracy = ", metrics.accuracy_score(y_train, y_train_predict))

    y_valid_predict = rf.predict(x_valid)
    confusion = metrics.confusion_matrix(y_valid, y_valid_predict)
    print(confusion)
    print("validation accuracy = ", metrics.accuracy_score(y_valid, y_valid_predict))

    # Test model on the test dataset
    y_test_predict = rf.predict(x_test)
    print(metrics.confusion_matrix(y_test, y_test_predict))
    print("test accuracy = ", metrics.accuracy_score(y_test, y_test_predict))

    import matplotlib.pyplot as plt

    plt.plot(rf.feature_importances_)
    plt.ylabel('feature importance')
    plt.show()
