import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse, sys, os
sys.path.append(os.getcwd())
import pycstm

def plot_scatter(data, filename):
    fig = plt.gcf()
    fig.set_size_inches(16.0, 16.0)
    plt.clf()
    plt.scatter(data[:, 0], data[:, 1], s=20, marker="o", edgecolors="none", color="blue")
    plt.savefig('./data/fig/'+filename+'.png')

if __name__ == "__main__":
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()
    doc_vectors = np.asarray(cstm.get_doc_vectors(), dtype=np.float32)
    # print(doc_vectors)
    for i in range(ndim_d-1):
        plot_scatter(doc_vectors[:, i:], "doc_scatter_{}_{}".format(i, i+1))