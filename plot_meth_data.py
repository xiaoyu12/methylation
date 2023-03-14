import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)

import pickle
import numpy as np
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="The input .pickle file that stores the data",
    )

    # load prepared data
    args = parser.parse_args()
    pickle_file = args.input_file
    with open(pickle_file, "rb") as fp:
        save = pickle.load(fp)

    meth_diff = save["meth_diff"]
    expr_lfc = save["expr_lfc"]

    print(meth_diff.shape)
    print(expr_lfc.shape)
    print(expr_lfc.head())
    """x_train = data["train_dataset"].astype(np.float32)
    x_valid = data["valid_dataset"].astype(np.float32)
    x_test = data["test_dataset"].astype(np.float32)
    y_train = data["train_labels"]
    y_valid = data["valid_labels"]
    y_test = data["test_labels"]

    print(x_train.shape, y_train.shape)

    # choose 10 random data from x_train and plot them
    np.random.seed(42)
    sp = np.random.randint(x_train.shape[0], size=10)

    fig, axes = plt.subplots(5, 2, figsize=(25, 12))
    for i in range(5):
        for j in range(2):
            axes[i, j].plot(x_train[sp[i*2+j]], )

    plt.show()"""

# plot the methylation difference data for first 10 genes
    fig, axes = plt.subplots(5, 2, figsize=(25, 18))
    for i in range(5):
        for j in range(2):
            axes[i, j].plot(meth_diff.iloc[i * 2 + j], )
            axes[i, j].set_title(str(meth_diff.index[i*2+j]))

    plt.show()

    plt.plot(meth_diff.loc[431797], )
    plt.show()