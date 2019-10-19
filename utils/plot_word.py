import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse, sys, os
sys.path.append(os.getcwd())
import pycstm
# we downloaded ipa font from https://ipafont.ipa.go.jp/node17#jp
plt.rcParams['font.family'] = 'IPAGothic'
sns.set(font=["IPAGothic"], font_scale=2)

def plot_words(words, ndim_vector):
    with sns.axes_style("white", {"font.family": ["IPAGothic"]}):
        for i in range(ndim_vector-1):
            fig = plt.gcf()
            fig.set_size_inches(45.0, 45.0)
            plt.clf()
            for meta in words:
                word_id, word, count, vector = meta
                plt.text(vector[i], vector[i+1], word, fontsize=8)
            plt.xlim(-4, 4)
            plt.ylim(-4, 4)
            plt.savefig("./data/fig/"+"words_{}_{}.png".format(i, i+1))

if __name__ == "__main__":
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()
    common_words = cstm.get_high_freq_words(1000)
    plot_words(common_words, ndim_d)