import argparse, sys, os
sys.path.append(os.path.split(os.getcwd())[0])
import pycstm

if __name__ == "__main__":
    cstm = pycstm.cstm("./model/cstm.model")
    ndim_d = cstm.get_ndim_d()

    doc_vectors = cstm.get_doc_vectors()