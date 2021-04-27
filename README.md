# Continuous Space Topic Model

## Environment

- C++ 14+
- clang++ 9.0
- boost 1.71.0
- glog 0.4.0
- gflag 2.2.2
- boost-python3
- python3

## Usage

- process text data with mecab-python3

```bash
$ python3 utils/process.py --tar_path data/raw/ --save_path data/train/
```

- training CSTM with MH algorithm

```bash
$ make
$ ./cstm -ndim_d=20 -ignore_word_count=4 -epoch=100 -num_threads=1 -data_path=./data/train/ -validation_data_path=./data/validation/ -model_path=./model/cstm.model
```


- load CSTM model and plot vector

```bash
$ make install
$ python3 utils/plot_doc.py
$ python3 utils/plot_word.py
```

- caluculation cosine similarity between words/docs

```bash
$ python3 utils/similarity.py -word TARGET_WORD
```

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)
- [musyoku.github.io](http://musyoku.github.io/)
