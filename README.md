# Continuous Space Topic Model

## Description

the implementation of [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf) with C++

ABSTRACT:
>This paper proposes a continous space text model based on Gaussian processes. Introducing latent coordinates of words over which the Gaussian process is defined, we can encode word correlations directly and lead to a model that performs better than mixture models. Our model would serve as a foundation of more complex text models and also as a statistical visualization of texts.

trained latent word vector(corpus: NIPS, dim: 20)

![cstm_plot](https://seiichiinoue.github.io/img/cstm_result.png)

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
$ ./cstm -n_dim=20 -ignore_word_count=0 -epoch=100 -data_path=./data/train/ -model_path=./model/cstm.model
```


- load CSTM model and plot vector

```bash
$ make install
$ python3 utils/plot_doc.py
$ python3 utils/plot_word.py
```

- caluculation cosine similarity between words/docs

```bash
$ python3 utils/similarity.py
```

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)
- [musyoku.github.io](http://musyoku.github.io/)
