import sys, os, argparse
import numpy as np
sys.path.append(os.getcwd())
import pycstm

def levenshtein(s1, s2):
    n, m = len(s1), len(s2)
    # initialization
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i
    for i in range(m + 1):
        dp[0][i] = i
    # dp
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1,         # insertion
                           dp[i][j - 1] + 1,         # deletion
                           dp[i - 1][j - 1] + cost)  # replacement
    return dp[n][m]

def find_sim_words(tar):
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()

    words = cstm.get_words_similar_to_word(tar, 30)
    for meta in words:
        word_id, word, count, vector, cosine = meta
        vector = np.asarray(vector, dtype=np.float32)
        # word = word.encode(sys.stdout.encoding)
        print("word id: {} word: {} count: {} sim: {}".format(word_id, word, count, cosine))

def find_sim_docs(tar):
    tar = tar
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()

    filenames = cstm.get_doc_filenames()
    tmp = [0] * len(filenames)
    for i in range(len(filenames)):
        tmp[i] = levenshtein(filenames[i], tar)
    ret = filenames[np.argmin(tmp)]
    if tar == ret:
        print("target file: {}".format(tar))
    else:
        tar = ret
        print("ERROR: target file does not found in filename list")
        print("alternatively, we selected target file: {}".format(tar))
    if not tar.endswith('.txt'):
        tar += '.txt'

    docs = cstm.get_docs_similar_to_file(tar, 20)
    for meta in docs:
        doc_id, doc, vector, cosine = meta
        vector = np.asarray(vector, dtype=np.float32)
        # word = word.encode(sys.stdout.encoding)
        print("word id: {} word: {} sim: {}".format(doc_id, doc, cosine))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='this script for extract similar documents.')
    parser.add_argument('-word', help='target word')
    # parser.add_argument('-doc', help='target document file name')
    args = parser.parse_args()
    find_sim_words(args.word)
    # find_sim_docs(args.doc)
