# Continuous Space Topic Model

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

```bash
$ make
$ ./cstm -NDIM=20 -IGNORE=10 -EPOCH=100
```

## Reference

- [Modeling Text through Gaussian Processes](http://chasen.org/~daiti-m/paper/nl213cstm.pdf)

- [musyoku.github.io](http://musyoku.github.io/)
