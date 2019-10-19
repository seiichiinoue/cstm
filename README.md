# Continuous Space Topic Model

## Todo

- python wrapper with boost-python
- script for plot vector

## Description

the implementation of [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf) with C++

ABSTRACT:
>This paper proposes a continous space text model based on Gaussian processes. Introducing latent coordinates of words over which the Gaussian process is defined, we can encode word correlations directly and lead to a model that performs better than mixture models. Our model would serve as a foundation of more complex text models and also as a statistical visualization of texts.

trained latent word vector(corpus: NIPS, dim: 20)

![cstm_plot](https://seiichiinoue.github.io/img/cstm_result.png)

## Environment

- C++ 17+
- clang++ 9.0
- boost 1.71.0
- glog 0.4.0
- gflag 2.2.2
- boost-python3
- python3

## Usage

- prepare text data with mecab-python3

```bash
$ python utils/text.py --lang_type japanese --tar_path data/raw/ --save_path data/train/
```

- training CSTM with MCMC

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

## Experimental Result

experimental setting as follows:

- dimention of latent coordinates = 20
- vocabulary size =  15349
- num of documents = 2 (too small)
- num of words = 292584

### word vectors of 1-2 dimentions

![](https://seiichiinoue.github.io/img/word_plot_1_2.png)

### word vectors of 9-10 dimentions

![](https://seiichiinoue.github.io/img/word_plot_9_10.png)

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)
- [musyoku.github.io](http://musyoku.github.io/)
