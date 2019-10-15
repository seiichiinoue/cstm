# Continuous Space Topic Model

## Todo

- model saving, vector saving
- script for plot vector

## Description

the implementation of [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)

ABSTRACT:
>This paper proposes a continous space text model based on Gaussian processes. Introducing latent coordinates of words over which the Gaussian process is defined, we can encode word correlations directly and lead to a model that performs better than mixture models. Our model would serve as a foundation of more complex text models and also as a statistical visualization of texts.

trained latent word vector(corpus: NIPS, dim: 20)

![cstm_plot](https://seiichiinoue.github.io/img/cstm_result.png)

## Environment

- C++ 11
- clang++ 10.0
- boost 1.71.0
- glog 0.4.0
- gflag 2.2.2

## Usage

- prepare text data with mecab-python3

```bash
$ python utils/text.py --lang_type japanese --tar_path data/kokoro.txt --wakati_path data/kokoro-wakati.txt
```

- training CSTM

```bash
$ make
$ ./cstm -n_dim=20 -ignore_word_count=0 -epoch=100 -data_path=[string: your data] -model_path=[string: saveplace]
```

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)
- [musyoku.github.io](http://musyoku.github.io/)
