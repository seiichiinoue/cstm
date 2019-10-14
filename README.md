# Continuous Space Topic Model

## Todo

- model saving, vector saving
- script for plot vector

## Description

the implementation of [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)

ABSTRACT:
>This paper proposes a continous space text model based on Gaussian processes. Introducing latent coordinates of words over which the Gaussian process is defined, we can encode word correlations directly and lead to a model that performs better than mixture models. Our model would serve as a foundation of more complex text models and also as a statistical visualization of texts.

## Environment

- C++ 11
- clang++ 10.0
- boost 1.71.0
- glog 0.4.0
- gflag 2.2.2

## Usage

- prepare text data with mecab-python3

```bash
$ python utils/text.py --tar_path data/kokoro.txt --wakati_path data/kokoro-wakati.txt
```

- training CSTM

```bash
$ make
$ ./cstm -NDIM=20 -IGNORE=0 -EPOCH=100
```

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)
- [musyoku.github.io](http://musyoku.github.io/)
