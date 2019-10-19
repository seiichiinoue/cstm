import pylab
import numpy as np
import seaborn as sns
import argparse, sys, os
sys.path.append(os.getcwd())
import pycstm

def plot_words(words, ndim_vector):
    with sns.axes_style("white", {"font.family": ["MS Gothic"]}):
        for i in range(ndim_vector-1):
            fig = pylab.gcf()
            fig.set_size_inches(45.0, 45.0)
            pylab.clf()
            for meta in words:
                word_id, word, count, vector = meta
                pylab.text(vector[i], vector[i+1], word, fontsize=8)
            pylab.xlim(-4, 4)
            pylab.ylim(-4, 4)
            pylab.savefig("./data/fig/"+"words_{}_{}.png".format(i, i+1))

if __name__ == "__main__":
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()
    common_words = cstm.get_high_freq_words(1000)
    plot_words(common_words, ndim_d)