import pylab
import numpy as np
import seaborn as sns
import argparse, sys, os
sys.path.append(os.getcwd())
import pycstm

# def plot_kde(data, filename):
#     fig = pylab.gcf()
#     fig.set_size_inches(16.0, 16.0)
#     pylab.clf()
#     ax = sns.kdeplot(data[:, 0], data[:, 1], shade=True, cmap="Greens", n_levels=30, clip=[[-4, 4]] * 2)
#     kde = ax.get_figure()
#     kde.savefig('./data/fig/'+filename+'.png')

def plot_scatter(data, filename):
    fig = pylab.gcf()
    fig.set_size_inches(16.0, 16.0)
    pylab.clf()
    pylab.scatter(data[:, 0], data[:, 1], s=20, marker="o", edgecolors="none", color="blue")
    pylab.savefig('./data/fig/'+filename+'.png')

if __name__ == "__main__":
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()
    doc_vectors = np.asarray(cstm.get_doc_vectors(), dtype=np.float32)
    # print(doc_vectors)
    for i in range(ndim_d-1):
        plot_scatter(doc_vectors[:, i:], "doc_scatter_{}_{}".format(i, i+1))