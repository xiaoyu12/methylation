import numpy as np
import pandas as pd
import pickle
import glob
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt

# Main fuction
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Prepare methylation and LFC data for Ehux 1516 and 217 strains"
    )
    parser.add_argument(
        '-i',
        "--input",
        dest="in_file",
        default="data_EH1516_EH217_filtered.pickle",
        help="Input pickle file of methylation differences and LFC",
    )
    parser.add_argument(
        '-o',
        "--output",
        dest="out_file",
        default="data_prepared.pickle",
        help="Output pickle file",
    )

    args = parser.parse_args()

    # load the input pickle file
    with open(args.in_file, "rb") as fp:
        save = pickle.load(fp)
        meth_diff = save["meth_diff"]
        expr_lfc = save["expr_lfc"]

    print(meth_diff.shape)
    print(expr_lfc.shape)
    print(expr_lfc.head())

    # Use a seed for reproducibility
    np.random.seed(42)
    # shuffle the dataset randomly
    permut = np.random.permutation(meth_diff.shape[0])
    print("permutations:", permut[0:10])
    meth_diff = meth_diff.iloc[permut]
    expr_lfc = expr_lfc.iloc[permut]
    print(expr_lfc.head())
    print(np.sign(expr_lfc.head()))
    print(meth_diff.head())

    # split data into training, validation and testing sets
    l = meth_diff.shape[0]
    lv = int(l / 5)  # number of validation data
    lt = int(l / 5)  # number of testing data
    ln = l - lv - lt  # number of training data

    start, end = 0, ln
    print(f"start: {start}, end: {end}")
    mdf_train = np.asarray(meth_diff.iloc[start:end])
    lfc_train = np.asarray(expr_lfc.iloc[start:end])
    print(mdf_train.shape)
    print(mdf_train[0:5, 0:5])
    start, end = ln, ln + lv
    print(f"start: {start}, end: {end}")
    mdf_valid = np.asarray(meth_diff.iloc[start:end])
    lfc_valid = np.asarray(expr_lfc.iloc[start:end])
    start, end = ln + lv, ln + lv + lt
    print(f"start: {start}, end: {end}")
    mdf_test = np.asarray(meth_diff.iloc[start:end])
    lfc_test = np.asarray(expr_lfc.iloc[start:end])

    print(mdf_train.shape, lfc_train.shape)
    print(mdf_train[0:5, 0:5])
    print(np.min(mdf_train), np.max(mdf_train))
    print(np.min(lfc_train), np.max(lfc_train))

    pickle_file = args.out_file
    try:
        f = open(pickle_file, "wb")
        save = {
            "train_dataset": mdf_train,
            "train_labels": lfc_train,
            "valid_dataset": mdf_valid,
            "valid_labels": lfc_valid,
            "test_dataset": mdf_test,
            "test_labels": lfc_test,
        }
        pickle.dump(save, f)
        # pickle.dump(save, f)
        f.close()
    except Exception as e:
        print("Unable to save data to", pickle_file, ":", e)
        raise